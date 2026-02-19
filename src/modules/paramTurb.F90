module paramTurb
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
!       Initialised here so defaults survive repeated referenceState calls.
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
!       Default values set by setSSTDefaults().
! ======================================================================
!
!       Initialised here so defaults survive repeated referenceState calls.
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
