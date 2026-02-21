#!/bin/bash
# ==============================================================================
# ADflow GCC Environment for Paracloud HPC (BSCC)
# ==============================================================================
#
# Usage:
#   在 Slurm 脚本或交互式 SSH 中 source 此文件:
#     source /public3/home/m6s001303/software-m6s001303/adflow-gcc/env_gcc.sh
#
#   或在本地开发时:
#     source env_gcc.sh
#
# Background:
#   Intel ifort 编译的 ADflow 存在 Fortran module 数据重复问题,
#   导致 f2py 修改的湍流模型系数无法传递到计算内核.
#   GCC gfortran 编译版本已验证系数修改正确生效 (2026-02-21).
#
# ==============================================================================

# --- Step 1: Module system ---
source /public3/soft/modules/module.sh
module purge
module load gcc/12.2 mpi/openmpi/4.1.5-gcc12.2 miniforge/24.11

# --- Step 2: Conda environment ---
source activate adflow-gcc

# --- Step 3: Library paths ---
# PETSc + system OpenMPI
export LD_LIBRARY_PATH=/public3/soft/openmpi/4.1.5-gcc12.2/lib:/public3/home/m6s001303/software-m6s001303/adflow-gcc/petsc-install/lib:$LD_LIBRARY_PATH

# Force system OpenMPI (avoid conda's OpenMPI 5.x conflict)
export LD_PRELOAD=/public3/soft/openmpi/4.1.5-gcc12.2/lib/libmpi.so

# --- Step 4: ADflow Python path ---
# pip install 不可用 (HPC 无外网), 通过 PYTHONPATH + symlink 替代
export PYTHONPATH=/public3/home/m6s001303/software-m6s001303/adflow-gcc:$PYTHONPATH

# --- Step 5: Unbuffered output ---
export PYTHONUNBUFFERED=1

# ==============================================================================
# Key directories:
#   ADflow source:  /public3/home/m6s001303/software-m6s001303/adflow-gcc/
#   PETSc install:  /public3/home/m6s001303/software-m6s001303/adflow-gcc/petsc-install/
#   CGNS install:   /public3/home/m6s001303/software-m6s001303/adflow-gcc/cgns-install/
#   M6 mesh (ADF):  /public3/home/m6s001303/M6/ONERA_M6.cgns
#
# Limitation:
#   CGNS compiled without HDF5 - only ADF format meshes work.
#   NACA0012.cgns is HDF5 format and cannot be read.
# ==============================================================================
