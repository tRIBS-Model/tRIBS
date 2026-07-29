# Testing Needed: In-Process Reach Partitioning (v6.0.0)

Status: **not yet implemented** — this file documents what the black-box suite
should cover for the in-process METIS partitioning added in v6.0.0, once the
suite is revived. The Big Spring serial/parallel equivalence test is the
natural home for most of it.

## What changed

Parallel tRIBS no longer needs the external MeshBuilder + gpmetis + perl
workflow to produce a `.reach` graphfile. METIS is vendored into the codebase
(`src/metis_builds/`, compiled into the binary) and the reach graph is
partitioned in-process (`src/tPartition/`, called from `tGraph::partition`).

New keyword semantics (breaking change from v5):

- `GRAPHOPTION`: partitioning method — `0` = SF (surface flow edges only),
  `1` = SSF (flow + subsurface flux edges), `2` = SSFH (SSF + headwater
  balancing constraint). Values outside 0–2 are a fatal error.
  (Old meaning was 0 = default split, 1 = reach file, 2 = inlet/outlet file.)
- `GRAPHFILE`: now **optional**. If set, the file must exist and is read
  (with validation). If blank/absent, the partition is generated every run
  and written to `<OUTFILENAME>_<SF|SSF|SSFH>_<np>nodes.reach`. The output
  path is never searched for an existing file.
- `PARALLELMODE`: gains value `2` = partition-only (tRIBSpar only): build
  the mesh and stream network, write the `.reach` graphfile, print partition
  statistics, and exit without simulating.

## 1. The golden invariant (highest priority)

**Partitioning must never change the hydrology.** Any run below should match
the serial reference within the tolerance already used by the Big Spring
serial/parallel test. Matrix per benchmark watershed, blank `GRAPHFILE`:

| Run | Tests |
|---|---|
| serial `tRIBS` | reference answer |
| `tRIBSpar -np 1` | parallel code path, no actual partitioning |
| `tRIBSpar -np 2`, `-np 4` | partition-count independence |
| `-np 3` with `GRAPHOPTION` 0, 1, and 2 | method independence |

Compare: outlet hydrograph (`_Outlet.qout` — preferred over the `.mrf` `Srf`
column), `.mrf` water balance, and integrated spatial output at the final
timestep. The GRAPHOPTION 0/1/2 runs must also match *each other*; any
divergence points directly at the partitioner.

Run on every available benchmark — mesh topology (reach count, confluence
structure, headwaters) is what stresses the graph construction. A watershed
with fewer reaches than `np` is a useful stress case.

## 2. Partition sanity

- The integrated spatial output includes processor rank: mapped, partitions
  should look like contiguous sub-basins, not salt-and-pepper.
- **Edge-cut identity** (cheap, high-value check): the edge-cut reported on
  the `Partitioned ...` line must equal the *channel* crossings for
  `GRAPHOPTION 0` and the *subsurface* crossings for `GRAPHOPTION 1` and `2`.
  Channel connections are a subset of subsurface neighbor pairs (reaches
  joined end-to-end share a junction node, so they are mesh neighbors too),
  so the SSF/SSFH graph is exactly the flux-pair set and the SF graph is the
  flow edges alone. If SF's edge-cut ever matches the subsurface number, SF
  is wrongly being fed flux edges. Channel crossings must always be ≤
  subsurface crossings.
- Console output: `Partitioned N reaches into P partitions (SF); E of T
  channel connections cross processor boundaries` and `Partition p: W nodes`.
  Node counts should be reasonably even — but judge against the floor the
  statistics report prints (set by the largest indivisible reach), not
  against a fixed percentage: with few reaches per partition, perfect
  balance is often unreachable. One partition with many times the nodes of
  another is still a failure.
- The partition messages must appear **once**, not once per MPI rank. All
  ranks run METIS; only the master reports.
- If a legacy MeshBuilder-generated `.reach` file exists for a benchmark,
  run with `GRAPHFILE` pointing at it: it should pass validation and
  reproduce the same hydrology (confirms reach numbering consistency
  between tRIBS and the retired MeshBuilder).

## 3. Feature mechanics and error paths

Each should exit with a clear message — never crash, hang, or run with a
bad partition:

1. Blank `GRAPHFILE` → "No parallel partitioning graph file (GRAPHFILE)
   provided..." message; `<basename>_<method>_<np>nodes.reach` written next
   to model outputs with `<partition> <reachID>` pairs covering every reach.
2. `GRAPHFILE` set to that written file → "Partitioning read from reach
   file..." and results identical to the generated run.
3. `GRAPHFILE` set to a nonexistent path → fatal "not found" error
   (no silent generation at the typo'd path).
4. A `.reach` file generated for `-np 3` run with `-np 2` → fatal
   partition-range error.
5. A hand-truncated `.reach` file → fatal reach-count error.
6. A duplicated line in the `.reach` file → fatal duplicate-assignment error.
7. `GRAPHOPTION` outside 0–2 → fatal error during input preprocessing
   (fails in seconds, before mesh building).
8. `PARALLELMODE 2` under `mpirun -np N` → `.reach` file written, partition
   statistics printed (per-partition table, load balance, edge cuts), exits
   after Part 3b without simulating; the `.reach` file must be byte-identical
   to the one a full `PARALLELMODE 1` run with the same config generates.
9. `PARALLELMODE 2` in the serial `tRIBS` binary → fatal error pointing at
   `mpirun -np <N> tRIBSpar`.
10. `PARALLELMODE 2` with `GRAPHFILE` set → validates and reports statistics
    for the existing file instead of generating a new one.

## 4. Determinism and performance

- Same parallel config run twice → the two generated `.reach` files must be
  byte-identical (all MPI ranks run METIS independently and rely on
  deterministic agreement).
- Wall time serial vs `-np 2/4`: parallel should not be slower than serial;
  a slowdown suggests a pathological partition.
- Cross-platform note: a `.reach` generated on macOS need not be
  byte-identical to one generated on Linux (compiler-dependent tie-breaks
  inside METIS). Each is valid; do not assert cross-platform file equality.

## Known gaps / future hardening (not test failures)

- Every MPI rank runs METIS independently and trusts identical results;
  master-computes-then-broadcast would be more robust.
- CI (`.github/workflows/`) compiles the vendored METIS on Linux and both
  macOS architectures but never *runs* a parallel partition — item 1 above
  is the runtime coverage that matters most.
- Restart/hotstart interacting with generated partitions is untested.
