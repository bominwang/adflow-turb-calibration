"""
Patch ADflow to expose SA and SST closure coefficients as runtime-modifiable variables.
=======================================================================================

Changes:
  1. paramTurb.F90: Remove `parameter` from SA & SST constants, add setter subroutines
  2. adflow.pyf:    Add paramturb module to f2py interface (SA + SST vars & subroutines)
  3. initializeFlow.F90: Add calls to setSADefaults() and setSSTDefaults()

Usage (inside MDO Lab Docker):
    python3 patch_adflow_turb.py
    cd /home/mdolabuser/repos/adflow
    make clean && make
    pip install -e .

Python API after rebuild:
    # SA coefficients
    pt = solver.adflow.paramturb
    pt.setsadefaults()
    pt.setsaconstants(cb1, cb2, sigma, kappa, cv1, cw2, cw3, ct3, ct4)

    # SST coefficients
    pt.setsstdefaults()
    pt.setsstconstants(sstk, a1, betas, sigk1, sigw1, beta1, sigk2, sigw2, beta2)
"""

import os
import shutil

ADFLOW_SRC = os.environ.get(
    "ADFLOW_SRC", "/home/mdolabuser/repos/adflow/src"
)
PARAMTURB_F90 = os.path.join(ADFLOW_SRC, "modules", "paramTurb.F90")
ADFLOW_PYF = os.path.join(ADFLOW_SRC, "f2py", "adflow.pyf")
INIT_FLOW_F90 = os.path.join(ADFLOW_SRC, "initFlow", "initializeFlow.F90")


# ============================================================
# 1. Patch paramTurb.F90
# ============================================================
def patch_paramturb():
    """Replace paramTurb.F90: SA & SST constants become mutable, add setter subroutines."""

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
!       Other turbulence model constants (K-omega, K-tau, V2-f) remain as
!       compile-time parameters.
!
    use constants, only: realType, intType
    implicit none
    save
!
! ======================================================================
!       Spalart-Allmaras constants (MUTABLE).
!       Default values set by setSADefaults().
! ======================================================================
!
!       9 calibration parameters (matching literature):
    real(kind=realType) :: rsaK       ! kappa (von Karman constant, default 0.41)
    real(kind=realType) :: rsaCb1     ! production coeff (default 0.1355)
    real(kind=realType) :: rsaCb2     ! diffusion coeff (default 0.622)
    real(kind=realType) :: rsaCb3     ! *** stores SIGMA (= 2/3), NOT c_b3 ***
    real(kind=realType) :: rsaCv1     ! wall damping coeff (default 7.1)
    real(kind=realType) :: rsaCw2     ! destruction coeff (default 0.3)
    real(kind=realType) :: rsaCw3     ! destruction coeff (default 2.0)
    real(kind=realType) :: rsaCt3     ! transition coeff (default 1.2)
    real(kind=realType) :: rsaCt4     ! transition coeff (default 0.5)
!
!       Derived constant (auto-recomputed by setter):
    real(kind=realType) :: rsaCw1     ! = rsaCb1/(rsaK^2) + (1+rsaCb2)/rsaCb3
!
!       Auxiliary (not in standard calibration set, but modifiable):
    real(kind=realType) :: rsaCt1     ! transition coeff (default 1.0)
    real(kind=realType) :: rsaCt2     ! transition coeff (default 2.0)
    real(kind=realType) :: rsaCrot    ! rotation correction (default 2.0)
!
! ======================================================================
!       K-omega SST constants (MUTABLE).
!       Default values set by setSSTDefaults().
! ======================================================================
!
    real(kind=realType) :: rSSTK      ! kappa (default 0.41)
    real(kind=realType) :: rSSTA1     ! limiter constant a1 (default 0.31)
    real(kind=realType) :: rSSTBetas  ! beta* (default 0.09)

    real(kind=realType) :: rSSTSigk1  ! sigma_k1 (default 0.85)
    real(kind=realType) :: rSSTSigw1  ! sigma_omega1 (default 0.5)
    real(kind=realType) :: rSSTBeta1  ! beta_1 (default 0.075)

    real(kind=realType) :: rSSTSigk2  ! sigma_k2 (default 1.0)
    real(kind=realType) :: rSSTSigw2  ! sigma_omega2 (default 0.856)
    real(kind=realType) :: rSSTBeta2  ! beta_2 (default 0.0828)
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
        ! Set all SA closure coefficients to their standard default values.
        ! Must be called once during solver initialization.
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
        ! Arguments (matching the standard SA calibration literature):
        !   cb1   - production coefficient      (default 0.1355)
        !   cb2   - diffusion coefficient        (default 0.622)
        !   sigma - diffusion ratio              (default 2/3, stored as rsaCb3)
        !   kappa - von Karman constant          (default 0.41)
        !   cv1   - wall damping coefficient     (default 7.1)
        !   cw2   - destruction coefficient      (default 0.3)
        !   cw3   - destruction coefficient      (default 2.0)
        !   ct3   - transition coefficient       (default 1.2)
        !   ct4   - transition coefficient       (default 0.5)
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
        ! Set all SST closure coefficients to their standard default values.
        ! Must be called once during solver initialization.
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
        ! Arguments:
        !   sstk  - von Karman constant kappa    (default 0.41)
        !   a1    - limiter constant              (default 0.31)
        !   betas - beta*                         (default 0.09)
        !   sigk1 - sigma_k1                     (default 0.85)
        !   sigw1 - sigma_omega1                 (default 0.5)
        !   beta1 - beta_1                       (default 0.075)
        !   sigk2 - sigma_k2                     (default 1.0)
        !   sigw2 - sigma_omega2                 (default 0.856)
        !   beta2 - beta_2                       (default 0.0828)
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
# 3. Patch initializeFlow.F90
# ============================================================
def patch_initializeflow():
    """Patch initializeFlow.F90: add calls to setSADefaults() and setSSTDefaults()."""

    if not os.path.exists(INIT_FLOW_F90):
        print(f"  WARNING: {INIT_FLOW_F90} not found.")
        print("  Manually add 'call setSADefaults()' and 'call setSSTDefaults()'.")
        return

    backup = INIT_FLOW_F90 + ".bak"
    if not os.path.exists(backup):
        shutil.copy2(INIT_FLOW_F90, backup)
        print(f"  Backup: {backup}")

    with open(INIT_FLOW_F90, "r") as f:
        content = f.read()

    # Check if already patched
    if "setSADefaults" in content and "setSSTDefaults" in content:
        print("  initializeFlow.F90 already patched, skipping.")
        return

    marker = "        ! Compute the dimensional viscosity from Sutherland's law"
    if marker not in content:
        print("  WARNING: Could not find insertion marker in initializeFlow.F90")
        print("  Manually add calls to setSADefaults() and setSSTDefaults().")
        return

    insert_code = """        ! Initialize turbulence model closure coefficients to default values.
        ! Must be called before any turbulence constants are used.
        call setSADefaults()
        call setSSTDefaults()

"""
    content = content.replace(marker, insert_code + marker)

    with open(INIT_FLOW_F90, "w") as f:
        f.write(content)
    print(f"  Patched: {INIT_FLOW_F90}")


# ============================================================
# Main
# ============================================================
def main():
    print("=" * 70)
    print("Patching ADflow for runtime SA + SST coefficient modification")
    print("=" * 70)

    print("\n[1/3] Patching paramTurb.F90 (SA + SST mutable)...")
    patch_paramturb()

    print("\n[2/3] Patching adflow.pyf (f2py interface)...")
    patch_pyf()

    print("\n[3/3] Patching initializeFlow.F90 (initialization calls)...")
    patch_initializeflow()

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
