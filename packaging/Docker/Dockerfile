FROM alpine:3.24 AS builder

RUN apk update && \
    apk add --no-cache \
        git \
        cmake \
        build-base \
        openmpi-dev \
        zlib-dev \
    && rm -rf /var/cache/apk/*

WORKDIR /build

# Build a released tag rather than a moving branch, so an image tag always
# rebuilds to the same source. Override for a one-off build with
# --build-arg TRIBS_VERSION=<tag>.
ARG TRIBS_VERSION=v6.0.0

RUN git clone --depth 1 --branch ${TRIBS_VERSION} https://github.com/tRIBS-Model/tRIBS.git .

RUN cmake -B cmake-build-parallel -S . -Dparallel=ON -DCMAKE_BUILD_TYPE=Release && \
    cmake --build cmake-build-parallel --target all

RUN cmake -B cmake-build-serial -S . -Dparallel=OFF -DCMAKE_BUILD_TYPE=Release && \
    cmake --build cmake-build-serial --target all


# Create the Final Image
FROM alpine:3.24

RUN apk update && \
    apk add --no-cache \
        bash \
        libstdc++ \
        openmpi \
        openssh-client \
        zlib \
    && rm -rf /var/cache/apk/*

# Create a dedicated group first.
RUN addgroup -S tribsgroup

# Create the user and add them to the new group.
RUN adduser -S tRIBSuser -G tribsgroup

# Create directories and set ownership to the user and their new group.
RUN mkdir -p /tribs/bin /tribs/shared && \
    chown -R tRIBSuser:tribsgroup /tribs

# Copy the compiled binaries to the correct location
COPY --from=builder /build/cmake-build-serial/tRIBS /tribs/bin/tRIBS
COPY --from=builder /build/cmake-build-parallel/tRIBSpar /tribs/bin/tRIBSpar

WORKDIR /tribs
USER tRIBSuser
CMD ["/bin/bash"]