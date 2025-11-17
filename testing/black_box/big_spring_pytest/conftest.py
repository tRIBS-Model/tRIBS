import json
import sys
import os
import pytest
import getpass
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed

extra_path = f"{os.environ['HOME']}/Documents/Repos/Forked/pytRIBS"
if extra_path not in sys.path:
    sys.path.append(extra_path)

from pytRIBS.classes import Model as model
from pytRIBS.classes import Results as results


@pytest.fixture(scope='module')
def setup_data():
    # Setup conditional for refining or adding tests because running big_springs for a month simulation can take between
    # 9-20 minutes depending on optimization. For operational testing run_flag should always be set to True
    run_flag = True

    # read in config.json
    with open('../config.json', 'r') as openfile:
        json_obj = json.load(openfile)

    bs_path = json_obj['bs_path']
    binary_par_path = json_obj['bin_parallel']
    binary_ser_path = json_obj['bin_serial']
    input_par_file = f'{bs_path}/src/in_files/big_spring_par.in'
    input_ser_file = f'{bs_path}/src/in_files/big_spring.in'

    if run_flag:

        # Get the current username
        username = getpass.getuser()

        # Get the current date
        current_date = datetime.now()

        print(f'\ntRIBS Integrated Test: Big Spring Benchmark\nUser: {username} \nDate: {current_date} \n')

        # Use ThreadPoolExecutor to run tasks in parallel
        with ThreadPoolExecutor(max_workers=2) as executor:
            futures = {
                executor.submit(run_parallel_task, binary_par_path, input_par_file): 'parallel',
                executor.submit(run_serial_task, binary_ser_path, input_ser_file): 'serial'
            }

            results_dict = {}
            for future in as_completed(futures):
                task_name = futures[future]
                try:
                    mrf_result = future.result()
                    results_dict[task_name] = mrf_result
                    print(f'{task_name} simulation completed')
                except Exception as e:
                    print(f'{task_name} simulation generated an exception: {e}')

        # Access the results
        parallel_result = results_dict.get('parallel')
        serial_result = results_dict.get('serial')

    else:
        parallel_result = get_existing_results(input_par_file)
        serial_result = get_existing_results(input_ser_file)

    r_par = results(input_par_file)

    yield parallel_result, serial_result, r_par


def run_parallel_task(binary_path, input_file):
    # setup pytribs Model class, read input file, run model, instantiate Results class, return mrf file.
    m = model()
    m.read_input_file(input_file)
    m.run(binary_path, input_file, mpi_command='mpirun -n 3', verbose=False)
    r = results(input_file)
    r.get_mrf_results()
    return r.mrf['mrf']


def run_serial_task(binary_path, input_file):
    # setup pytribs Model class, read input file, run model, instantiate Results class, return mrf file.
    m = model()
    m.read_input_file(input_file)
    m.run(binary_path, input_file, verbose=False)
    r = results(input_file)
    r.get_mrf_results()
    return r.mrf['mrf']


def get_existing_results(input_file):
    r = results(input_file)
    r.get_mrf_results()
    return r.mrf['mrf']
