FROM alpine:3.19 AS builder

RUN apk update && \
    apk add --no-cache \
        git \
        cmake \
        build-base \
        openmpi-dev \
        zlib-dev \
    && rm -rf /var/cache/apk/*

WORKDIR /build

RUN git clone --depth 1 --branch main https://github.com/tribshms/tRIBS.git .

RUN cmake -B cmake-build-parallel -S . -Dparallel=ON -DCMAKE_BUILD_TYPE=Release && \
    cmake --build cmake-build-parallel --target all

RUN cmake -B cmake-build-serial -S . -Dparallel=OFF -DCMAKE_BUILD_TYPE=Release && \
    cmake --build cmake-build-serial --target all


# Create the Final Image
FROM alpine:3.19

RUN apk update && \
    apk add --no-cache \
        bash \
        openmpi \
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