# ADflow — Turbulence Model Coefficient Calibration Fork

Based on [MDO Lab ADflow](https://github.com/mdolab/adflow), modified to support **runtime modification of SA and SST turbulence closure coefficients** for Bayesian calibration and uncertainty quantification.

## How It Works

Original ADflow declares SA and SST closure coefficients as Fortran `parameter` (compile-time constants). This fork:

1. **Removes `parameter`** from SA and SST coefficients in `paramTurb.F90`, making them mutable module variables with default initial values
2. **Adds setter subroutines** (`setSAConstants`, `setSSTConstants`, etc.) in Fortran
3. **Replaces hardcoded constants** (e.g., `0.09` for beta*) with the module variables in `SST.F90` and `turbUtils.F90`
4. **Provides a Python ctypes API** (`adflow_turb_ctypes.py`) that directly writes to the Fortran symbol addresses in the loaded `libadflow.so`

## Python API

The Python interface uses ctypes to directly access Fortran module symbols, bypassing f2py (see [Known Issues](#known-issues) for why).

```python
from adflow import ADFLOW
from adflow_turb_ctypes import (
    set_sa_constants, set_sa_defaults, get_sa_constants,
    set_sst_constants, set_sst_defaults, get_sst_constants,
)

solver = ADFLOW(options=aeroOptions)

# --- SA: set 9 calibration parameters (cw1 auto-recomputed) ---
set_sa_constants(
    cb1=0.1355, cb2=0.622, sigma=2.0/3.0, kappa=0.41,
    cv1=7.1, cw2=0.3, cw3=2.0, ct3=1.2, ct4=0.5,
)

# Read back current values
print(get_sa_constants())
# {'cb1': 0.1355, 'cb2': 0.622, 'sigma': 0.6667, 'kappa': 0.41,
#  'cv1': 7.1, 'cw1': 3.2391, 'cw2': 0.3, 'cw3': 2.0, 'ct3': 1.2, 'ct4': 0.5}

# --- SST: set 9 independent parameters ---
set_sst_constants(
    sstk=0.41, a1=0.31, betas=0.09,
    sigk1=0.85, sigw1=0.5, beta1=0.075,
    sigk2=1.0, sigw2=0.856, beta2=0.0828,
)

# Reset to defaults
set_sa_defaults()
set_sst_defaults()
```

### Important Usage Notes

- **Create ONE solver instance** and switch coefficients between `AeroProblem` runs. Creating multiple `ADFLOW` instances causes a PETSc crash.
- **Re-apply coefficients before each solve**: ADflow's internal `referenceState` may be called on AeroProblem switch. Set coefficients right before `solver(ap)`.
- Copy `adflow_turb_ctypes.py` to your working directory or add this repo's root to `PYTHONPATH`.

### Coefficient Reference

**SA Model** — 9 calibration parameters:

| Parameter | Variable | Default | Description |
|-----------|----------|---------|-------------|
| cb1 | `rsaCb1` | 0.1355 | Production coefficient |
| cb2 | `rsaCb2` | 0.622 | Diffusion coefficient |
| sigma | `rsaCb3` | 2/3 | Diffusion ratio (**stored as rsaCb3, not c_b3**) |
| kappa | `rsaK` | 0.41 | von Karman constant |
| cv1 | `rsaCv1` | 7.1 | Wall damping coefficient |
| cw2 | `rsaCw2` | 0.3 | Destruction coefficient |
| cw3 | `rsaCw3` | 2.0 | Destruction coefficient |
| ct3 | `rsaCt3` | 1.2 | Transition coefficient |
| ct4 | `rsaCt4` | 0.5 | Transition coefficient |

Derived: `cw1 = cb1/kappa^2 + (1+cb2)/sigma` (auto-recomputed by setter)

Auxiliary (modifiable but not in standard calibration): `rsaCt1` (1.0), `rsaCt2` (2.0), `rsaCrot` (2.0)

**SST Model** — 9 independent parameters:

| Parameter | Variable | Default |
|-----------|----------|---------|
| kappa | `rSSTK` | 0.41 |
| a1 | `rSSTA1` | 0.31 |
| beta* | `rSSTBetas` | 0.09 |
| sigma_k1 | `rSSTSigk1` | 0.85 |
| sigma_w1 | `rSSTSigw1` | 0.5 |
| beta_1 | `rSSTBeta1` | 0.075 |
| sigma_k2 | `rSSTSigk2` | 1.0 |
| sigma_w2 | `rSSTSigw2` | 0.856 |
| beta_2 | `rSSTBeta2` | 0.0828 |

## Build & Deploy

This repository does **not** provide pre-compiled binaries. The workflow is: patch the ADflow source with `patch_adflow_turb.py`, then compile from source.

> **Note**: The `adflow/libadflow.so` in this repo is a build artifact and may be corrupted by git line-ending conversion. Always compile from source on the target platform.

### Option A: Docker (Recommended for Development)

#### 1. Pull the MDO Lab Docker image

```bash
docker pull mdolab/public:u22-gcc-ompi-stable
```

#### 2. Start container

```bash
docker run -it --name adflow-turb mdolab/public:u22-gcc-ompi-stable bash
```

#### 3. Clone and patch

```bash
git clone https://github.com/bominwang/adflow-turb-calibration.git /tmp/repo
python3 /tmp/repo/patch_adflow_turb.py
```

Expected output (4 files patched):

```
[1/4] Patching paramTurb.F90 (SA + SST mutable with initialisers)...
  Patched: .../paramTurb.F90
[2/4] Checking initializeFlow.F90 (no setter calls in referenceState)...
  .../initializeFlow.F90: no setter calls found, OK.
[3/4] Patching SST.F90 (f1 blending: 0.09 -> rSSTBetas)...
  Patched: .../SST.F90
[4/4] Patching turbUtils.F90 (eddy viscosity: 0.09 -> rSSTBetas)...
  Patched: .../turbUtils.F90
```

#### 4. Compile and install

```bash
cd /home/mdolabuser/repos/adflow
make clean && make
pip install -e .
```

Compilation warnings about `Type mismatch` are from upstream ADflow CGNS bindings and are harmless. Success is indicated by:

```
Testing if module libadflow can be imported...
Module libadflow was successfully imported
```

#### 5. Verify

```bash
cp /tmp/repo/adflow_turb_ctypes.py /home/mdolabuser/
python3 -c "
import sys; sys.path.insert(0, '/home/mdolabuser')
from adflow_turb_ctypes import set_sa_constants, get_sa_constants, set_sa_defaults

set_sa_defaults()
c = get_sa_constants()
assert abs(c['cb1'] - 0.1355) < 1e-6, f'default failed: {c}'

set_sa_constants(0.5, 0.622, 2/3, 0.41, 7.1, 0.3, 2.0, 1.2, 0.5)
c = get_sa_constants()
assert abs(c['cb1'] - 0.5) < 1e-6, f'set failed: {c}'

set_sa_defaults()
c = get_sa_constants()
assert abs(c['cb1'] - 0.1355) < 1e-6, f'reset failed: {c}'

print('ctypes interface: OK')
print('SA read/write/reset: OK')
"
```

### Option B: HPC Native Compilation

For HPC systems without Docker/Singularity. Requires GCC gfortran + OpenMPI.

See `hpc/deploy.sh` for an automated deployment script, or follow the manual steps below.

#### Prerequisites

- GCC gfortran (tested with 12.2)
- OpenMPI (tested with 4.1.5) or Intel MPI
- PETSc 3.15+ compiled with shared libraries
- CGNS 4.x compiled with Fortran support
- Python 3 with numpy and mpi4py

#### Manual steps

```bash
# 1. Clone
git clone https://github.com/bominwang/adflow-turb-calibration.git
cd adflow-turb-calibration

# 2. Set environment
export PETSC_DIR=/path/to/petsc-install
export PETSC_ARCH=""
export CGNS_HOME=/path/to/cgns-install

# 3. Patch
export ADFLOW_SRC=$(pwd)/src
python3 patch_adflow_turb.py

# 4. Create config/config.mk (adjust paths for your HPC)
cat > config/config.mk << 'EOF'
FF90 = mpifort
CC   = mpicc

FF77_FLAGS = -fPIC -fdefault-real-8 -fdefault-double-8 -g -O2 -fallow-argument-mismatch
FF90_FLAGS = $(FF77_FLAGS) -std=f2008
FFXX_OPT_FLAGS = -O2
C_FLAGS   = -fPIC -O2 -std=c99

AR       = ar
AR_FLAGS = -rvs

LINKER       = $(FF90)
LINKER_FLAGS =

CGNS_INCLUDE_FLAGS=-I$(CGNS_HOME)/include
CGNS_LINKER_FLAGS=-L$(CGNS_HOME)/lib -lcgns

include ${PETSC_DIR}/lib/petsc/conf/variables
PETSC_INCLUDE_FLAGS=${PETSC_CC_INCLUDES} -I$(PETSC_DIR)
PETSC_LINKER_FLAGS=${PETSC_LIB}

FF90_PRECISION_FLAGS = $(FF90_INTEGER_PRECISION_FLAG)$(FF90_REAL_PRECISION_FLAG)
CC_PRECISION_FLAGS   = $(CC_INTEGER_PRECISION_FLAG) $(CC_REAL_PRECISION_FLAG)

PYTHON = python3
PYTHON-CONFIG = python3-config
F2PY = f2py
EOF

# 5. Compile
make clean && make

# 6. Install (if network available)
pip install -e .

# If no network, use PYTHONPATH:
# export PYTHONPATH=$(pwd):$PYTHONPATH
```

#### HPC Common Issues

| Problem | Cause | Solution |
|---------|-------|----------|
| `cgns.mod: Unexpected EOF` | `.mod` from Intel, used with GCC | Recompile CGNS with GCC |
| `u_int64_t` undefined | `adStack.c` non-standard type | `#include <stdint.h>`, replace with `uint64_t` |
| `__intel_sse2_*` undefined | CGNS static lib has Intel objects | Recompile CGNS C sources with GCC |
| conda MPI conflict | conda OpenMPI 5.x overrides system 4.x | Use absolute compiler paths or `LD_PRELOAD` |
| `pip install` fails (no network) | Can't download setuptools | Use `PYTHONPATH` instead |

### Option C: Singularity (HPC with Container Support)

```bash
singularity pull adflow-turb.sif docker://mdolab/public:u22-gcc-ompi-stable
```

Or build a custom image with Dockerfile:

```dockerfile
FROM mdolab/public:u22-gcc-ompi-stable
COPY patch_adflow_turb.py adflow_turb_ctypes.py /tmp/
RUN python3 /tmp/patch_adflow_turb.py \
    && cd /home/mdolabuser/repos/adflow \
    && make clean && make \
    && pip install -e .
```

## Patched Files (4)

| File | Change |
|------|--------|
| `src/modules/paramTurb.F90` | SA + SST constants become mutable; 4 setter subroutines added |
| `src/initFlow/initializeFlow.F90` | Ensure `referenceState` does not call setters (would clobber user coefficients) |
| `src/turbulence/SST.F90` | Hardcoded `0.09` (beta*) in f1 blending replaced with `rSSTBetas` |
| `src/turbulence/turbUtils.F90` | Hardcoded `0.09` (beta*) in eddy viscosity replaced with `rSSTBetas` |

## Verification

### HPC A/B Test (ONERA M6, Paracloud BSCC, Job 36920821)

Three-case test with SA model on ONERA M6 wing (M=0.8395, alpha=3.06 deg). Single ADFLOW solver instance, coefficients switched between runs via ctypes API:

| Case | cb1 | CL | CD |
|------|-----|--------|---------|
| A (default) | 0.1355 | 0.2606478129 | 0.0187979813 |
| B (+48%) | 0.2 | 0.2625755177 | 0.0193669582 |
| C (-41%) | 0.08 | 0.2573790690 | 0.0178077544 |

All three cases produce **different** CL/CD. Trends are physically consistent:
- cb1 up → more turbulent production → thicker boundary layer → higher CL/CD
- cb1 down → less production → thinner boundary layer → lower CL/CD

### Docker NACA 0012 Cp Validation

SA and SST models each tested with 4 coefficient sets (1 default + 3 random). All cases converge and show coefficient-dependent Cp variations.

See `examples/NACA0012/` for scripts and results.

## Known Issues

### f2py Module Variable Bug (Critical)

**Do NOT use the f2py module interface (`pt = solver.adflow.paramturb`) for coefficient modification.** This is the reason this fork uses ctypes instead.

**Symptom**: Setting coefficients via f2py reads back correctly, but has **zero effect** on CFD computation. All CL/CD values are identical to machine precision regardless of coefficient changes.

**Root cause**: When f2py wraps a Fortran shared library, the module variable data exists in **two separate memory locations**: one accessed by the f2py Python wrapper, and one used by the Fortran computation subroutines (`sa_block`, `blockette`, etc.). Python writes to copy A; Fortran reads from copy B.

**Affected compilers**: Both Intel ifort/ifx AND GCC gfortran. This is NOT a compiler-specific issue.

**Verification**: Hardcoding `rsaCb1 = 0.5` in source code and recompiling produces different CL/CD (Job 36920761), confirming the Fortran module variables work correctly. The bug is exclusively in how f2py accesses them.

**Solution**: The ctypes API (`adflow_turb_ctypes.py`) uses `ctypes.c_double.in_dll()` to locate the actual Fortran symbol address (e.g., `__paramturb_MOD_rsacb1`) in the loaded `libadflow.so`. This writes to the same memory that computation code reads, and is verified to produce different CL/CD (Job 36920788, 36920821).

### CGNS HDF5 Format

If CGNS was compiled without HDF5, only ADF-format CGNS files can be read. Compile CGNS with `-DCGNS_ENABLE_HDF5=ON` for HDF5 support.

### PETSc Reinitialization

Creating multiple `ADFLOW` instances in one Python process causes a PETSc `MPI_ABORT`. Use a single solver instance and create new `AeroProblem` objects for different coefficient sets.

## Upstream

Based on [mdolab/adflow](https://github.com/mdolab/adflow) v2.12.1.

## License

Same as upstream: GNU Lesser General Public License (LGPL), version 2.1. See [LICENSE.md](LICENSE.md).
