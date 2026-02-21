"""
Patch ADflow to expose SA and SST closure coefficients as runtime-modifiable variables.
=======================================================================================

Changes:
  1. paramTurb.F90:       Remove `parameter` from SA & SST constants, add default
                          initial values and setter subroutines
  2. initializeFlow.F90:  NO calls to setSADefaults/setSSTDefaults in referenceState
                          (defaults come from module-level initialisers; referenceState
                          is called on every AeroProblem switch and would clobber
                          user-supplied coefficients)
  3. SST.F90:             Replace hardcoded 0.09 (beta*) with rSSTBetas in f1 blending
  4. turbUtils.F90:       Replace hardcoded 0.09 (beta*) with rSSTBetas in f2/eddy-visc

NOTE: The f2py interface (.pyf patch) has been removed.  f2py module variable access
has a known memory duplication bug: Python writes go to f2py's copy of the module data,
but Fortran computation code reads from its own copy.  This makes coefficient
modifications appear to succeed (read-back matches) while having ZERO effect on the
actual CFD computation.  This bug affects BOTH Intel ifort AND GCC gfortran.

The recommended Python API is via ctypes (adflow_turb_ctypes.py), which directly
writes to the Fortran symbol addresses in the loaded shared library.

Usage (inside MDO Lab Docker or HPC):
    python3 patch_adflow_turb.py
    cd /path/to/adflow
    make clean && make
    pip install -e .

Python API after rebuild (via ctypes):
    from adflow_turb_ctypes import (
        set_sa_constants, set_sa_defaults, get_sa_constants,
        set_sst_constants, set_sst_defaults, get_sst_constants,
    )

    # SA coefficients
    set_sa_constants(cb1=0.1355, cb2=0.622, sigma=2/3, kappa=0.41,
                     cv1=7.1, cw2=0.3, cw3=2.0, ct3=1.2, ct4=0.5)

    # SST coefficients
    set_sst_constants(sstk=0.41, a1=0.31, betas=0.09,
                      sigk1=0.85, sigw1=0.5, beta1=0.075,
                      sigk2=1.0, sigw2=0.856, beta2=0.0828)
"""

import os
import re
import shutil

ADFLOW_SRC = os.environ.get(
    "ADFLOW_SRC", "/home/mdolabuser/repos/adflow/src"
)
PARAMTURB_F90 = os.path.join(ADFLOW_SRC, "modules", "paramTurb.F90")
INIT_FLOW_F90 = os.path.join(ADFLOW_SRC, "initFlow", "initializeFlow.F90")
SST_F90 = os.path.join(ADFLOW_SRC, "turbulence", "SST.F90")
TURBUTILS_F90 = os.path.join(ADFLOW_SRC, "turbulence", "turbUtils.F90")


# ============================================================
# 1. Patch paramTurb.F90
# ============================================================
def patch_paramturb():
    """Replace paramTurb.F90: SA & SST constants become mutable with
    default initial values, add setter subroutines."""

    backup = PARAMTURB_F90 + ".bak"
    if not os.path.exists(backup):
        shutil.copy2(PARAMTURB_F90, backup)
        print(f"  Backup: {backup}")

    new_content = r"""module paramTurb
!
!       Module that contains the constants for the turbulence models.
!
!       PATCHED: SA and SST constants are mutable (no 'parameter') to allow
!       runtime modification for Bayesian calibration of closure coefficients.
!       Default values are set via Fortran initialisers so they survive
!       repeated referenceState calls without being clobbered.
!       Other turbulence model constants (K-omega, K-tau, V2-f) remain as
!       compile-time parameters.
!
    use constants, only: realType, intType
    implicit none
    save
!
! ======================================================================
!       Spalart-Allmaras constants (MUTABLE).
!       Initialised here so defaults survive repeated referenceState calls.
! ======================================================================
!
!       9 calibration parameters (matching literature):
    real(kind=realType) :: rsaK   = 0.41_realType      ! kappa
    real(kind=realType) :: rsaCb1 = 0.1355_realType    ! production coeff
    real(kind=realType) :: rsaCb2 = 0.622_realType     ! diffusion coeff
    real(kind=realType) :: rsaCb3 = 0.66666666667_realType ! sigma (2/3)
    real(kind=realType) :: rsaCv1 = 7.1_realType       ! wall damping coeff
    real(kind=realType) :: rsaCw2 = 0.3_realType       ! destruction coeff
    real(kind=realType) :: rsaCw3 = 2.0_realType       ! destruction coeff
    real(kind=realType) :: rsaCt3 = 1.2_realType       ! transition coeff
    real(kind=realType) :: rsaCt4 = 0.5_realType       ! transition coeff
!
!       Derived constant (auto-recomputed by setter):
    real(kind=realType) :: rsaCw1 = 3.239067055837563_realType
!
!       Auxiliary (not in standard calibration set, but modifiable):
    real(kind=realType) :: rsaCt1  = 1.0_realType      ! transition coeff
    real(kind=realType) :: rsaCt2  = 2.0_realType      ! transition coeff
    real(kind=realType) :: rsaCrot = 2.0_realType      ! rotation correction
!
! ======================================================================
!       K-omega SST constants (MUTABLE).
!       Initialised here so defaults survive repeated referenceState calls.
! ======================================================================
!
    real(kind=realType) :: rSSTK     = 0.41_realType   ! kappa
    real(kind=realType) :: rSSTA1    = 0.31_realType   ! limiter constant a1
    real(kind=realType) :: rSSTBetas = 0.09_realType   ! beta*

    real(kind=realType) :: rSSTSigk1 = 0.85_realType   ! sigma_k1
    real(kind=realType) :: rSSTSigw1 = 0.5_realType    ! sigma_omega1
    real(kind=realType) :: rSSTBeta1 = 0.075_realType  ! beta_1

    real(kind=realType) :: rSSTSigk2 = 1.0_realType    ! sigma_k2
    real(kind=realType) :: rSSTSigw2 = 0.856_realType  ! sigma_omega2
    real(kind=realType) :: rSSTBeta2 = 0.0828_realType ! beta_2
!
! ======================================================================
!       K-omega constants (unchanged, compile-time parameters).
! ======================================================================
!
    real(kind=realType), parameter :: rkwK = 0.41_realType
    real(kind=realType), parameter :: rkwSigk1 = 0.5_realType
    real(kind=realType), parameter :: rkwSigw1 = 0.5_realType
    real(kind=realType), parameter :: rkwSigd1 = 0.5_realType
    real(kind=realType), parameter :: rkwBeta1 = 0.0750_realType
    real(kind=realType), parameter :: rkwBetas = 0.09_realType
!
! ======================================================================
!       K-tau constants (unchanged).
! ======================================================================
!
    real(kind=realType), parameter :: rktK = 0.41_realType
    real(kind=realType), parameter :: rktSigk1 = 0.5_realType
    real(kind=realType), parameter :: rktSigt1 = 0.5_realType
    real(kind=realType), parameter :: rktSigd1 = 0.5_realType
    real(kind=realType), parameter :: rktBeta1 = 0.0750_realType
    real(kind=realType), parameter :: rktBetas = 0.09_realType
!
! ======================================================================
!       V2-f constants (unchanged).
! ======================================================================
!
    real(kind=realType), parameter :: rvfC1 = 1.4_realType
    real(kind=realType), parameter :: rvfC2 = 0.3_realType
    real(kind=realType), parameter :: rvfBeta = 1.9_realType
    real(kind=realType), parameter :: rvfSigk1 = 1.0_realType
    real(kind=realType), parameter :: rvfSige1 = 0.7692307692_realType
    real(kind=realType), parameter :: rvfSigv1 = 1.00_realType
    real(kind=realType), parameter :: rvfCn = 70.0_realType

    real(kind=realType), parameter :: rvfN1Cmu = 0.190_realType
    real(kind=realType), parameter :: rvfN1A = 1.300_realType
    real(kind=realType), parameter :: rvfN1B = 0.250_realType
    real(kind=realType), parameter :: rvfN1Cl = 0.300_realType
    real(kind=realType), parameter :: rvfN6Cmu = 0.220_realType
    real(kind=realType), parameter :: rvfN6A = 1.400_realType
    real(kind=realType), parameter :: rvfN6B = 0.045_realType
    real(kind=realType), parameter :: rvfN6Cl = 0.230_realType

    real(kind=realType) :: rvfLimitK, rvfLimitE, rvfCl
    real(kind=realType) :: rvfCmu
!
! ======================================================================
!       Wall function fit variables (unchanged).
! ======================================================================
!
    integer(kind=intType) :: nFit

    real(kind=realType), dimension(:), allocatable :: ypT, reT
    real(kind=realType), dimension(:), allocatable :: up0, up1
    real(kind=realType), dimension(:), allocatable :: up2, up3

    real(kind=realType), dimension(:, :), allocatable :: tup0, tup1
    real(kind=realType), dimension(:, :), allocatable :: tup2, tup3
#ifndef USE_TAPENADE
    real(kind=realType), dimension(:), allocatable :: ypTb, reTb
    real(kind=realType), dimension(:), allocatable :: up0b, up1b
    real(kind=realType), dimension(:), allocatable :: up2b, up3b
#endif

    logical, dimension(:), allocatable :: tuLogFit

contains

! ======================================================================
!   SA setter subroutines
! ======================================================================

    subroutine setSADefaults()
        !
        ! Reset all SA closure coefficients to their standard default values.
        !
        implicit none

        rsaK   = 0.41_realType
        rsaCb1 = 0.1355_realType
        rsaCb2 = 0.622_realType
        rsaCb3 = 0.66666666667_realType   ! sigma = 2/3
        rsaCv1 = 7.1_realType
        rsaCw2 = 0.3_realType
        rsaCw3 = 2.0_realType
        rsaCt1 = 1.0_realType
        rsaCt2 = 2.0_realType
        rsaCt3 = 1.2_realType
        rsaCt4 = 0.5_realType
        rsaCrot = 2.0_realType

        ! Derived constant: c_w1 = c_b1/kappa^2 + (1+c_b2)/sigma
        rsaCw1 = rsaCb1 / (rsaK * rsaK) &
                 + (1.0_realType + rsaCb2) / rsaCb3

    end subroutine setSADefaults

    subroutine setSAConstants(cb1, cb2, sigma, kappa, cv1, cw2, cw3, ct3, ct4)
        !
        ! Set the 9 independent SA closure coefficients at runtime.
        ! Automatically recomputes the derived constant c_w1.
        !
        implicit none
        real(kind=realType), intent(in) :: cb1, cb2, sigma, kappa
        real(kind=realType), intent(in) :: cv1, cw2, cw3, ct3, ct4

        rsaCb1 = cb1
        rsaCb2 = cb2
        rsaCb3 = sigma    ! NOTE: rsaCb3 stores sigma, not c_b3
        rsaK   = kappa
        rsaCv1 = cv1
        rsaCw2 = cw2
        rsaCw3 = cw3
        rsaCt3 = ct3
        rsaCt4 = ct4

        ! CRITICAL: recompute derived constant c_w1
        rsaCw1 = rsaCb1 / (rsaK * rsaK) &
                 + (1.0_realType + rsaCb2) / rsaCb3

    end subroutine setSAConstants

! ======================================================================
!   SST setter subroutines
! ======================================================================

    subroutine setSSTDefaults()
        !
        ! Reset all SST closure coefficients to their standard default values.
        !
        implicit none

        rSSTK     = 0.41_realType
        rSSTA1    = 0.31_realType
        rSSTBetas = 0.09_realType

        rSSTSigk1 = 0.85_realType
        rSSTSigw1 = 0.5_realType
        rSSTBeta1 = 0.0750_realType

        rSSTSigk2 = 1.0_realType
        rSSTSigw2 = 0.856_realType
        rSSTBeta2 = 0.0828_realType

    end subroutine setSSTDefaults

    subroutine setSSTConstants(sstk, a1, betas, sigk1, sigw1, beta1, &
                               sigk2, sigw2, beta2)
        !
        ! Set the 9 SST closure coefficients at runtime.
        ! All 9 coefficients are independent (no derived constants).
        !
        implicit none
        real(kind=realType), intent(in) :: sstk, a1, betas
        real(kind=realType), intent(in) :: sigk1, sigw1, beta1
        real(kind=realType), intent(in) :: sigk2, sigw2, beta2

        rSSTK     = sstk
        rSSTA1    = a1
        rSSTBetas = betas

        rSSTSigk1 = sigk1
        rSSTSigw1 = sigw1
        rSSTBeta1 = beta1

        rSSTSigk2 = sigk2
        rSSTSigw2 = sigw2
        rSSTBeta2 = beta2

    end subroutine setSSTConstants

end module paramTurb
"""

    with open(PARAMTURB_F90, "w") as f:
        f.write(new_content)
    print(f"  Patched: {PARAMTURB_F90}")


# ============================================================
# 2. Ensure initializeFlow.F90 does NOT call setters
# ============================================================
def patch_initializeflow():
    """Ensure initializeFlow.F90 does NOT call setSADefaults/setSSTDefaults
    inside referenceState.  referenceState is invoked on every AeroProblem
    switch and would overwrite user-supplied coefficients.  Defaults are
    provided by module-level initialisers in paramTurb.F90 instead."""

    if not os.path.exists(INIT_FLOW_F90):
        print(f"  WARNING: {INIT_FLOW_F90} not found, skipping.")
        return

    backup = INIT_FLOW_F90 + ".bak"
    if not os.path.exists(backup):
        shutil.copy2(INIT_FLOW_F90, backup)
        print(f"  Backup: {backup}")

    with open(INIT_FLOW_F90, "r") as f:
        content = f.read()

    # Remove any existing calls to setSADefaults / setSSTDefaults
    # that a previous version of the patch may have inserted.
    if "call setSADefaults" in content or "call setSSTDefaults" in content:
        # Remove the call lines and any surrounding comment block we added
        content = content.replace(
            "        ! Initialize turbulence model closure coefficients to default values.\n"
            "        ! Must be called before any turbulence constants are used.\n"
            "        call setSADefaults()\n"
            "        call setSSTDefaults()\n\n",
            "        ! Turbulence closure coefficients are initialised by their module-level\n"
            "        ! default values in paramTurb.F90.  Do NOT call setSADefaults /\n"
            "        ! setSSTDefaults here, because referenceState is invoked every time\n"
            "        ! the AeroProblem changes, which would overwrite any user-supplied\n"
            "        ! custom coefficients set via setSAConstants / setSSTConstants.\n\n"
        )
        # Also handle the case where only the calls exist without our comment
        content = re.sub(
            r" *call setSADefaults\(\)\n", "", content
        )
        content = re.sub(
            r" *call setSSTDefaults\(\)\n", "", content
        )
        print(f"  Patched: {INIT_FLOW_F90} (removed setter calls from referenceState)")
    else:
        print(f"  {INIT_FLOW_F90}: no setter calls found, OK.")


# ============================================================
# 3. Patch SST.F90  (replace hardcoded 0.09 with rSSTBetas)
# ============================================================
def patch_sst():
    """In the f1SST subroutine, replace hardcoded 0.09_realType (beta*)
    with rSSTBetas from paramTurb, and add the necessary use statement."""

    if not os.path.exists(SST_F90):
        print(f"  WARNING: {SST_F90} not found, skipping.")
        return

    backup = SST_F90 + ".bak"
    if not os.path.exists(backup):
        shutil.copy2(SST_F90, backup)
        print(f"  Backup: {backup}")

    with open(SST_F90, "r") as f:
        content = f.read()

    changed = False

    # Add 'use paramTurb, only: rSSTBetas' to the f1SST subroutine
    # (it already has use constants, blockPointers, etc. but NOT paramTurb)
    old_imports = (
        "        use constants\n"
        "        use blockPointers\n"
        "        use inputTimeSpectral\n"
        "        use iteration\n"
        "        use turbMod\n"
        "        use utils, only: setPointers\n"
        "        use turbUtils, only: kwCDTerm"
    )
    new_imports = (
        "        use constants\n"
        "        use blockPointers\n"
        "        use inputTimeSpectral\n"
        "        use iteration\n"
        "        use turbMod\n"
        "        use paramTurb, only: rSSTBetas\n"
        "        use utils, only: setPointers\n"
        "        use turbUtils, only: kwCDTerm"
    )
    if old_imports in content and "use paramTurb, only: rSSTBetas" not in content:
        content = content.replace(old_imports, new_imports, 1)
        changed = True

    # Replace the hardcoded 0.09_realType with rSSTBetas in the f1 blending
    if "0.09_realType * w(i, j, k, itu2) * d2Wall(i, j, k))" in content:
        content = content.replace(
            "0.09_realType * w(i, j, k, itu2) * d2Wall(i, j, k))",
            "rSSTBetas * w(i, j, k, itu2) * d2Wall(i, j, k))"
        )
        changed = True

    if changed:
        with open(SST_F90, "w") as f:
            f.write(content)
        print(f"  Patched: {SST_F90}")
    else:
        print(f"  {SST_F90}: already patched or pattern not found.")


# ============================================================
# 4. Patch turbUtils.F90  (replace hardcoded 0.09 with rSSTBetas)
# ============================================================
def patch_turbutils():
    """In SSTEddyViscosity, replace hardcoded 0.09_realType (beta*) with
    rSSTBetas.  paramTurb is already imported in this subroutine."""

    if not os.path.exists(TURBUTILS_F90):
        print(f"  WARNING: {TURBUTILS_F90} not found, skipping.")
        return

    backup = TURBUTILS_F90 + ".bak"
    if not os.path.exists(backup):
        shutil.copy2(TURBUTILS_F90, backup)
        print(f"  Backup: {backup}")

    with open(TURBUTILS_F90, "r") as f:
        content = f.read()

    if "0.09_realType * w(i, j, k, itu2) * d2Wall(i, j, k))" in content:
        content = content.replace(
            "0.09_realType * w(i, j, k, itu2) * d2Wall(i, j, k))",
            "rSSTBetas * w(i, j, k, itu2) * d2Wall(i, j, k))"
        )
        with open(TURBUTILS_F90, "w") as f:
            f.write(content)
        print(f"  Patched: {TURBUTILS_F90}")
    else:
        print(f"  {TURBUTILS_F90}: already patched or pattern not found.")


# ============================================================
# Main
# ============================================================
def main():
    print("=" * 70)
    print("Patching ADflow for runtime SA + SST coefficient modification")
    print("=" * 70)

    print("\n[1/4] Patching paramTurb.F90 (SA + SST mutable with initialisers)...")
    patch_paramturb()

    print("\n[2/4] Checking initializeFlow.F90 (no setter calls in referenceState)...")
    patch_initializeflow()

    print("\n[3/4] Patching SST.F90 (f1 blending: 0.09 -> rSSTBetas)...")
    patch_sst()

    print("\n[4/4] Patching turbUtils.F90 (eddy viscosity: 0.09 -> rSSTBetas)...")
    patch_turbutils()

    print("\n" + "=" * 70)
    print("Patch complete!")
    print("=" * 70)
    print()
    print("Next steps:")
    print("  1. Rebuild ADflow:")
    print("     cd /path/to/adflow")
    print("     make clean && make")
    print("     pip install -e .")
    print()
    print("  2. Copy adflow_turb_ctypes.py to your working directory")
    print()
    print("Python API (via ctypes -- recommended):")
    print()
    print("  from adflow_turb_ctypes import (")
    print("      set_sa_constants, set_sa_defaults, get_sa_constants,")
    print("      set_sst_constants, set_sst_defaults, get_sst_constants,")
    print("  )")
    print()
    print("  # --- SA model (9 calibration params + auto cw1) ---")
    print("  set_sa_constants(")
    print("      cb1=0.1355, cb2=0.622, sigma=2./3.,")
    print("      kappa=0.41, cv1=7.1, cw2=0.3, cw3=2.0,")
    print("      ct3=1.2, ct4=0.5)")
    print()
    print("  # --- SST model (9 independent params) ---")
    print("  set_sst_constants(")
    print("      sstk=0.41, a1=0.31, betas=0.09,")
    print("      sigk1=0.85, sigw1=0.5, beta1=0.075,")
    print("      sigk2=1.0, sigw2=0.856, beta2=0.0828)")
    print()
    print("WARNING: Do NOT use the f2py module variable interface")
    print("  (pt = solver.adflow.paramturb; pt.rsacb1 = ...) for coefficient")
    print("  modification.  f2py has a memory duplication bug that makes")
    print("  modifications appear to succeed while having zero effect on")
    print("  the CFD computation.  Use the ctypes API above instead.")


if __name__ == "__main__":
    main()
