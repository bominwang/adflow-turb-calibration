# ADflow — Turbulence Calibration Edition

Modified version of [MDO Lab ADflow](https://github.com/mdolab/adflow) with **runtime-modifiable SA and SST turbulence closure coefficients** for Bayesian calibration and uncertainty quantification.

## What's Changed

The original ADflow declares SA and SST closure coefficients as Fortran `parameter` (compile-time constants), making them immutable at runtime. This fork removes the `parameter` attribute and exposes all coefficients through the f2py Python interface.

### Modified Files (3 files)

| File | Change |
|------|--------|
| `src/modules/paramTurb.F90` | SA + SST constants made mutable; added 4 setter subroutines |
| `src/f2py/adflow.pyf` | Added `paramturb` module: 22 variables + 4 subroutines |
| `src/initFlow/initializeFlow.F90` | Added `setSADefaults()` and `setSSTDefaults()` initialization calls |

### Exposed Variables

**SA model** — 13 variables (9 calibration + 4 auxiliary):

| Variable | Parameter | Default | Note |
|----------|-----------|---------|------|
| `rsacb1` | c_b1 | 0.1355 | Production coefficient |
| `rsacb2` | c_b2 | 0.622 | Diffusion coefficient |
| `rsacb3` | **sigma** | 0.6667 | **Naming trap: stores sigma (=2/3), NOT c_b3** |
| `rsak` | kappa | 0.41 | von Karman constant |
| `rsacv1` | c_v1 | 7.1 | Wall damping |
| `rsacw1` | c_w1 | 3.2391 | Derived: cb1/k^2 + (1+cb2)/sigma |
| `rsacw2` | c_w2 | 0.3 | Destruction coefficient |
| `rsacw3` | c_w3 | 2.0 | Destruction coefficient |
| `rsact1` | c_t1 | 1.0 | Trip function |
| `rsact2` | c_t2 | 2.0 | Trip function |
| `rsact3` | c_t3 | 1.2 | Trip function |
| `rsact4` | c_t4 | 0.5 | Trip function |
| `rsacrot` | c_rot | 2.0 | Rotation correction |

**SST model** — 9 independent variables:

| Variable | Parameter | Default |
|----------|-----------|---------|
| `rsstk` | kappa | 0.41 |
| `rssta1` | a1 | 0.31 |
| `rsstbetas` | beta* | 0.09 |
| `rsstsigk1` | sigma_k1 | 0.85 |
| `rsstsigw1` | sigma_w1 | 0.5 |
| `rsstbeta1` | beta_1 | 0.075 |
| `rsstsigk2` | sigma_k2 | 1.0 |
| `rsstsigw2` | sigma_w2 | 0.856 |
| `rsstbeta2` | beta_2 | 0.0828 |

### Subroutines

| Subroutine | Description |
|------------|-------------|
| `setsadefaults()` | Reset all 13 SA variables to defaults |
| `setsaconstants(cb1, cb2, sigma, kappa, cv1, cw2, cw3, ct3, ct4)` | Set 9 SA parameters, auto-recompute c_w1 |
| `setsstdefaults()` | Reset all 9 SST variables to defaults |
| `setsstconstants(sstk, a1, betas, sigk1, sigw1, beta1, sigk2, sigw2, beta2)` | Set 9 SST parameters |

## Usage

```python
from adflow import ADFLOW

solver = ADFLOW(options=aeroOptions)
pt = solver.adflow.paramturb

# --- SA: set 9 calibration parameters (c_w1 auto-recomputed) ---
pt.setsaconstants(
    0.1355,    # cb1
    0.622,     # cb2
    2.0/3.0,   # sigma (stored as rsacb3)
    0.41,      # kappa
    7.1,       # cv1
    0.3,       # cw2
    2.0,       # cw3
    1.2,       # ct3
    0.5        # ct4
)

# --- SST: set 9 independent parameters ---
pt.setsstconstants(
    0.41,      # sstk
    0.31,      # a1
    0.09,      # betas
    0.85,      # sigk1
    0.5,       # sigw1
    0.075,     # beta1
    1.0,       # sigk2
    0.856,     # sigw2
    0.0828     # beta2
)

# --- Or set individually ---
pt.rsacw2 = 0.055       # SA
pt.rsstbeta1 = 0.06     # SST
```

## Pre-built Binary

A pre-compiled `libadflow.so` is included in `adflow/libadflow.so`, built with:
- Docker image: `mdolab/public:u22-gcc-ompi-stable`
- Ubuntu 22.04, GCC, OpenMPI

To use the pre-built binary, the runtime environment must match the build environment (same OS, MPI, and Python version as the MDO Lab Docker image).

## Building from Source

### Using the patch script

If you have the original ADflow source:

```bash
python patch_adflow_turb.py   # Patches 3 source files (idempotent)
cd /path/to/adflow
make clean && make
pip install -e .
```

### Using this repository directly

```bash
git clone https://github.com/bominwang/adflow-turb-calibration.git
cd adflow-turb-calibration
make clean && make
pip install -e .
```

## Validation

All 7 tests pass:

```
discover            : PASS   (22 variables + 4 subroutines)
sa_defaults         : PASS
sst_defaults        : PASS
sa_readwrite        : PASS   (13/13)
sst_readwrite       : PASS   (9/9)
sa_setter           : PASS   (9 params + c_w1 auto-recompute)
sst_setter          : PASS   (9 params)
```

## Upstream

Based on [mdolab/adflow](https://github.com/mdolab/adflow) v2.12.1.

## License

Same as upstream: GNU Lesser General Public License (LGPL), version 2.1. See [LICENSE.md](LICENSE.md).
