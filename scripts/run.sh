#!/usr/bin/env bash
set -euo pipefail
BUILD_DIR=build
if [! -d "$BUILD_DIR"]; then
    mkdir -p $BUILD_DIR
    pushd $BUILD_DIR
    cmake ..
    make -j
    popd
fi
./${BUILD_DIR}/src/kaon_sullivan_gen examples/kaon_default.yaml