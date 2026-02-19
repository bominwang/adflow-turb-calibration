"""
Patch ADflow to expose SA and SST closure coefficients as runtime-modifiable variables.
=======================================================================================

Changes:
  1. paramTurb.F90:       Remove `parameter` from SA & SST constants, add default
                          initial values and setter subroutines
  2. adflow.pyf:          Add paramturb module to f2py interface (SA + SST vars & subs)
  3. initializeFlow.F90:  NO calls to setSADefaults/setSSTDefaults in referenceState
                          (defaults come from module-level initialisers; referenceState
                          is called on every AeroProblem switch and would clobber
                          user-supplied coefficients)
  4. SST.F90:             Replace hardcoded 0.09 (beta*) with rSSTBetas in f1 blending
  5. turbUtils.F90:       Replace hardcoded 0.09 (beta*) with rSSTBetas in f2/eddy-visc

Usage (inside MDO Lab Docker):
    python3 patch_adflow_turb.py
    cd /home/mdolabuser/repos/adflow
    make clean && make
    pip install -e .

Python API after rebuild:
    pt = solver.adflow.paramturb

    # SA coefficients
    pt.setsadefaults()
    pt.setsaconstants(cb1, cb2, sigma, kappa, cv1, cw2, cw3, ct3, ct4)

    # SST coefficients
    pt.setsstdefaults()
    pt.setsstconstants(sstk, a1, betas, sigk1, sigw1, beta1, sigk2, sigw2, beta2)
"""

import os
import re
import shutil

ADFLOW_SRC = os.environ.get(
    "ADFLOW_SRC", "/home/mdolabuser/repos/adflow/src"
)
PARAMTURB_F90 = os.path.join(ADFLOW_SRC, "modules", "paramTurb.F90")
ADFLOW_PYF = os.path.join(ADFLOW_SRC, "f2py", "adflow.pyf")
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
# 2. Patch adflow.pyf
# ============================================================
def patch_pyf():
    """Patch adflow.pyf: add paramturb module with SA + SST vars and subroutines."""

    backup = ADFLOW_PYF + ".bak"
    if not os.path.exists(backup):
        shutil.copy2(ADFLOW_PYF, backup)
        print(f"  Backup: {backup}")

    with open(ADFLOW_PYF, "r") as f:
        content = f.read()

    # Check if already patched
    if "module paramturb" in content.lower():
        print("  adflow.pyf already contains paramturb module, skipping.")
        return

    # f2py interface block for paramturb module
    paramturb_pyf = """
       module paramturb ! in :adflow:../../modules/paramTurb.F90
         use constants
         ! ----- SA constants (13 variables) -----
         real(kind=realtype) :: rsak
         real(kind=realtype) :: rsacb1
         real(kind=realtype) :: rsacb2
         real(kind=realtype) :: rsacb3
         real(kind=realtype) :: rsacv1
         real(kind=realtype) :: rsacw1
         real(kind=realtype) :: rsacw2
         real(kind=realtype) :: rsacw3
         real(kind=realtype) :: rsact1
         real(kind=realtype) :: rsact2
         real(kind=realtype) :: rsact3
         real(kind=realtype) :: rsact4
         real(kind=realtype) :: rsacrot

         ! ----- SST constants (9 variables) -----
         real(kind=realtype) :: rsstk
         real(kind=realtype) :: rssta1
         real(kind=realtype) :: rsstbetas
         real(kind=realtype) :: rsstsigk1
         real(kind=realtype) :: rsstsigw1
         real(kind=realtype) :: rsstbeta1
         real(kind=realtype) :: rsstsigk2
         real(kind=realtype) :: rsstsigw2
         real(kind=realtype) :: rsstbeta2

         ! ----- SA subroutines -----
         subroutine setsadefaults()
         end subroutine setsadefaults

         subroutine setsaconstants(cb1, cb2, sigma, kappa, cv1, cw2, cw3, ct3, ct4)
           real(kind=realtype), intent(in) :: cb1
           real(kind=realtype), intent(in) :: cb2
           real(kind=realtype), intent(in) :: sigma
           real(kind=realtype), intent(in) :: kappa
           real(kind=realtype), intent(in) :: cv1
           real(kind=realtype), intent(in) :: cw2
           real(kind=realtype), intent(in) :: cw3
           real(kind=realtype), intent(in) :: ct3
           real(kind=realtype), intent(in) :: ct4
         end subroutine setsaconstants

         ! ----- SST subroutines -----
         subroutine setsstdefaults()
         end subroutine setsstdefaults

         subroutine setsstconstants(sstk, a1, betas, sigk1, sigw1, beta1, sigk2, sigw2, beta2)
           real(kind=realtype), intent(in) :: sstk
           real(kind=realtype), intent(in) :: a1
           real(kind=realtype), intent(in) :: betas
           real(kind=realtype), intent(in) :: sigk1
           real(kind=realtype), intent(in) :: sigw1
           real(kind=realtype), intent(in) :: beta1
           real(kind=realtype), intent(in) :: sigk2
           real(kind=realtype), intent(in) :: sigw2
           real(kind=realtype), intent(in) :: beta2
         end subroutine setsstconstants

       end module paramturb
"""

    # Insert after "end module constants"
    marker = "       end module constants"
    if marker not in content:
        print("  ERROR: Could not find 'end module constants' in adflow.pyf")
        return

    content = content.replace(marker, marker + "\n" + paramturb_pyf)

    with open(ADFLOW_PYF, "w") as f:
        f.write(content)
    print(f"  Patched: {ADFLOW_PYF}")


# ============================================================
# 3. Patch initializeFlow.F90  (do NOT add setter calls)
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
# 4. Patch SST.F90  (replace hardcoded 0.09 with rSSTBetas)
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
# 5. Patch turbUtils.F90  (replace hardcoded 0.09 with rSSTBetas)
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

    print("\n[1/5] Patching paramTurb.F90 (SA + SST mutable with initialisers)...")
    patch_paramturb()

    print("\n[2/5] Patching adflow.pyf (f2py interface)...")
    patch_pyf()

    print("\n[3/5] Patching initializeFlow.F90 (no setter calls in referenceState)...")
    patch_initializeflow()

    print("\n[4/5] Patching SST.F90 (f1 blending: 0.09 -> rSSTBetas)...")
    patch_sst()

    print("\n[5/5] Patching turbUtils.F90 (eddy viscosity: 0.09 -> rSSTBetas)...")
    patch_turbutils()

    print("\n" + "=" * 70)
    print("Patch complete!")
    print("=" * 70)
    print()
    print("Next steps:")
    print("  1. Rebuild ADflow:")
    print("     cd /home/mdolabuser/repos/adflow")
    print("     make clean && make")
    print("     pip install -e .")
    print()
    print("Python API after rebuild:")
    print()
    print("  pt = solver.adflow.paramturb")
    print()
    print("  # --- SA model (9 calibration params + auto cw1) ---")
    print("  pt.setsaconstants(")
    print("      cb1=0.1355, cb2=0.622, sigma=2./3.,")
    print("      kappa=0.41, cv1=7.1, cw2=0.3, cw3=2.0,")
    print("      ct3=1.2, ct4=0.5)")
    print()
    print("  # --- SST model (9 independent params) ---")
    print("  pt.setsstconstants(")
    print("      sstk=0.41, a1=0.31, betas=0.09,")
    print("      sigk1=0.85, sigw1=0.5, beta1=0.075,")
    print("      sigk2=1.0, sigw2=0.856, beta2=0.0828)")
    print()
    print("  # --- Or set individually ---")
    print("  pt.rsacw2 = 0.055       # SA")
    print("  pt.rsstbeta1 = 0.06     # SST")


if __name__ == "__main__":
    main()
