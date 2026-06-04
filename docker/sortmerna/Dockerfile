FROM --platform=linux/amd64 ubuntu:22.04 AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    cmake \
    git \
    libgflags-dev \
    && rm -rf /var/lib/apt/lists/*

# Build zlib
RUN git clone --depth 1 --branch v1.3.1 https://github.com/madler/zlib.git /zlib
RUN cmake -S /zlib -B /zlib/build \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_C_FLAGS="-march=x86-64" \
        -DCMAKE_INSTALL_PREFIX=/opt/zlib \
    && cmake --build /zlib/build --parallel $(nproc) \
    && cmake --install /zlib/build

# Build RocksDB v7.10.2 (required by sortmerna v4.3.7)
RUN git clone --depth 1 --branch v7.10.2 https://github.com/facebook/rocksdb.git /rocksdb
RUN cmake -S /rocksdb -B /rocksdb/build \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_CXX_FLAGS="-march=x86-64" \
        -DCMAKE_C_FLAGS="-march=x86-64" \
        -DCMAKE_INSTALL_PREFIX=/opt/rocksdb \
        -DZLIB_ROOT=/opt/zlib \
        -DROCKSDB_BUILD_SHARED=ON \
        -DWITH_TESTS=OFF \
        -DWITH_BENCHMARK_TOOLS=OFF \
        -DWITH_TOOLS=OFF \
        -DWITH_SNAPPY=OFF \
        -DWITH_LZ4=OFF \
        -DWITH_ZSTD=OFF \
        -DWITH_BZ2=OFF \
        -DWITH_ZLIB=ON \
    && cmake --build /rocksdb/build --parallel $(nproc) \
    && cmake --install /rocksdb/build

# Build SortMeRNA
RUN git clone --depth 1 --branch v4.3.7 https://github.com/biocore/sortmerna.git /sortmerna \
    && git clone --depth 1 https://github.com/cameron314/concurrentqueue.git /sortmerna/3rdparty/concurrentqueue
RUN cmake -S /sortmerna -B /sortmerna/build \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_CXX_FLAGS="-march=x86-64 -I/sortmerna/3rdparty/concurrentqueue" \
        -DCMAKE_INSTALL_PREFIX=/opt/sortmerna \
        -DRocksDB_DIR=/opt/rocksdb/lib/cmake/rocksdb \
        -DZLIB_ROOT=/opt/zlib \
    && cmake --build /sortmerna/build --parallel $(nproc) \
    && cmake --install /sortmerna/build


FROM --platform=linux/amd64 ubuntu:22.04

RUN apt-get update && apt-get install -y --no-install-recommends \
    libgflags2.2 \
    zlib1g \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /opt/sortmerna/bin/sortmerna /usr/local/bin/sortmerna
COPY --from=builder /opt/rocksdb/lib/librocksdb.so* /usr/local/lib/
COPY --from=builder /opt/zlib/lib/libz.so* /usr/local/lib/

RUN ldconfig

RUN sortmerna --version
