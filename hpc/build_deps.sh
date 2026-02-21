#!/bin/bash
# ==============================================================================
# Build PETSc + CGNS from bundled source tarballs (offline)
# ==============================================================================
#
# This script compiles PETSc and CGNS from the source archives in deps/,
# requiring NO internet access.  It is intended for HPC environments where
# outbound network is unavailable.
#
# Prerequisites:
#   - GCC (gcc, gfortran)
#   - MPI (mpicc, mpifort, mpicxx)  -- system or module-loaded
#   - cmake >= 3.8 (for CGNS; if unavailable, see manual CGNS build below)
#   - make
#
# Usage:
#   export MPI_CC=mpicc          # or absolute path
#   export MPI_FC=mpifort        # or absolute path
#   export MPI_CXX=mpicxx        # or absolute path
#   bash hpc/build_deps.sh [install_prefix]
#
# Default install_prefix: ./deps/install
#
# After completion, set these before compiling ADflow:
#   export PETSC_DIR=<install_prefix>/petsc
#   export PETSC_ARCH=""
#   export CGNS_HOME=<install_prefix>/cgns
#
# ==============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
DEPS_DIR="$REPO_ROOT/deps"

# Install prefix
PREFIX="${1:-$DEPS_DIR/install}"
mkdir -p "$PREFIX"
PREFIX="$(cd "$PREFIX" && pwd)"

# Compiler defaults
MPI_CC="${MPI_CC:-mpicc}"
MPI_FC="${MPI_FC:-mpifort}"
MPI_CXX="${MPI_CXX:-mpicxx}"

NPROC=$(nproc 2>/dev/null || echo 4)

echo "======================================================================"
echo "  ADflow Offline Dependency Build"
echo "======================================================================"
echo "  DEPS_DIR:   $DEPS_DIR"
echo "  PREFIX:     $PREFIX"
echo "  MPI_CC:     $MPI_CC"
echo "  MPI_FC:     $MPI_FC"
echo "  MPI_CXX:    $MPI_CXX"
echo "  NPROC:      $NPROC"
echo "======================================================================"

# ==============================================================================
# 1. Build PETSc
# ==============================================================================
echo ""
echo "[1/2] Building PETSc 3.15.5 ..."
echo "----------------------------------------------------------------------"

PETSC_TAR="$DEPS_DIR/petsc-3.15.5.tar.gz"
SUPERLU_TAR="$DEPS_DIR/superlu_dist-6.4.0.tar.gz"
PETSC_BUILD="$DEPS_DIR/build/petsc-3.15.5"
PETSC_INSTALL="$PREFIX/petsc"

if [ -f "$PETSC_INSTALL/lib/petsc/conf/variables" ]; then
    echo "  PETSc already installed at $PETSC_INSTALL, skipping."
else
    if [ ! -f "$PETSC_TAR" ]; then
        echo "  ERROR: $PETSC_TAR not found"
        exit 1
    fi

    # Extract
    mkdir -p "$DEPS_DIR/build"
    if [ ! -d "$PETSC_BUILD" ]; then
        echo "  Extracting PETSc ..."
        tar -xzf "$PETSC_TAR" -C "$DEPS_DIR/build"
    fi

    # Configure
    cd "$PETSC_BUILD"

    PETSC_CONFIGURE_ARGS=(
        --prefix="$PETSC_INSTALL"
        --PETSC_ARCH=gcc-opt
        --with-cc="$MPI_CC"
        --with-fc="$MPI_FC"
        --with-cxx="$MPI_CXX"
        --with-debugging=0
        --with-shared-libraries=1
        --with-fortran-bindings=1
        --with-scalar-type=real
        COPTFLAGS="-O2"
        CXXOPTFLAGS="-O2"
        FOPTFLAGS="-O2"
    )

    # Add SuperLU_DIST if available
    if [ -f "$SUPERLU_TAR" ]; then
        echo "  Including SuperLU_DIST from local tarball"
        PETSC_CONFIGURE_ARGS+=(--download-superlu_dist="$SUPERLU_TAR")
    fi

    echo "  Configuring PETSc ..."
    python3 configure "${PETSC_CONFIGURE_ARGS[@]}"

    echo "  Compiling PETSc (nproc=$NPROC) ..."
    make PETSC_DIR="$PETSC_BUILD" PETSC_ARCH=gcc-opt all -j "$NPROC"

    echo "  Installing PETSc ..."
    make PETSC_DIR="$PETSC_BUILD" PETSC_ARCH=gcc-opt install

    echo "  PETSc installed: $PETSC_INSTALL"
fi

# ==============================================================================
# 2. Build CGNS
# ==============================================================================
echo ""
echo "[2/2] Building CGNS 4.2.0 ..."
echo "----------------------------------------------------------------------"

CGNS_TAR="$DEPS_DIR/CGNS-4.2.0.tar.gz"
CGNS_BUILD="$DEPS_DIR/build/CGNS-4.2.0"
CGNS_INSTALL="$PREFIX/cgns"

if [ -f "$CGNS_INSTALL/include/cgnslib.h" ] && [ -f "$CGNS_INSTALL/lib/libcgns.a" ]; then
    echo "  CGNS already installed at $CGNS_INSTALL, skipping."
else
    if [ ! -f "$CGNS_TAR" ]; then
        echo "  ERROR: $CGNS_TAR not found"
        exit 1
    fi

    # Extract
    mkdir -p "$DEPS_DIR/build"
    if [ ! -d "$CGNS_BUILD" ]; then
        echo "  Extracting CGNS ..."
        tar -xzf "$CGNS_TAR" -C "$DEPS_DIR/build"
    fi

    # Try cmake first, fall back to manual build
    if command -v cmake &>/dev/null; then
        echo "  Building CGNS with cmake ..."
        mkdir -p "$CGNS_BUILD/build_dir"
        cd "$CGNS_BUILD/build_dir"
        cmake .. \
            -DCMAKE_INSTALL_PREFIX="$CGNS_INSTALL" \
            -DCGNS_ENABLE_FORTRAN=ON \
            -DCGNS_BUILD_SHARED=OFF \
            -DCMAKE_C_COMPILER=gcc \
            -DCMAKE_Fortran_COMPILER=gfortran
        make -j "$NPROC"
        make install
    else
        echo "  cmake not found, building CGNS manually ..."

        cd "$CGNS_BUILD/src"
        # Configure with basic options
        ./configure \
            --prefix="$CGNS_INSTALL" \
            --enable-fortran \
            --disable-shared \
            --without-hdf5 \
            CC=gcc FC=gfortran
        make -j "$NPROC"
        make install
    fi

    echo "  CGNS installed: $CGNS_INSTALL"
fi

# ==============================================================================
# Summary
# ==============================================================================
echo ""
echo "======================================================================"
echo "  Dependencies built successfully!"
echo "======================================================================"
echo ""
echo "  Add to your environment before compiling ADflow:"
echo ""
echo "    export PETSC_DIR=$PETSC_INSTALL"
echo "    export PETSC_ARCH=\"\""
echo "    export CGNS_HOME=$CGNS_INSTALL"
echo ""
echo "  Then run:"
echo "    bash hpc/deploy.sh"
echo ""
echo "======================================================================"
