module turbMod
    !
    !       This local module contains variables used when the turbulence
    !       equations are solved.
    !
    use precision
       use complexify 
    implicit none
    save

    ! secondOrd:  whether or not a second order discretization for
    !             the advective terms must be used.
    ! sig1, sig2: Sigma coefficients in the diffusion terms of the
    !             different turbulence models.

    logical :: secondOrd
    complex(kind=realType) :: sig1, sig2

    ! dvt:     Pointer, which points to an unused part of dw. It is
    !          used for temporary storage of residual.
    ! vort:    Pointer, which points to an unused part of dw. It is
    !          used for temporary storage of the magnitude of
    !          vorticity squared.
    ! prod:    Pointer, which points to an unused part of dw. It is
    !          used for temporary storage of the unscaled production
    !          term.
    ! f1:      F1 blending function in the SST model.
    ! kwCD:   Cross diffusion term in the k-omega type models.
    ! ktCD:   Cross diffusion term in the k-tau model
    ! sct:     Time scale in the v2-f model.
    ! scl2:    Length scale in the v2-f model.
    ! strain2: Square of the strain.

    complex(kind=realType), dimension(:, :, :, :), pointer :: dvt
    complex(kind=realType), dimension(:, :, :), pointer :: vort
    complex(kind=realType), dimension(:, :, :), pointer :: prod
    complex(kind=realType), dimension(:, :, :), pointer :: f1
    complex(kind=realType), dimension(:, :, :), pointer :: kwCD
    complex(kind=realType), dimension(:, :, :), pointer :: ktCD
    complex(kind=realType), dimension(:, :, :), pointer :: sct
    complex(kind=realType), dimension(:, :, :), pointer :: scl2
    complex(kind=realType), dimension(:, :, :), pointer :: strain2

#ifndef USE_TAPENADE
    complex(kind=realType), dimension(:, :, :, :), pointer :: dvtd
    complex(kind=realType), dimension(:, :, :), pointer :: vortd
    complex(kind=realType), dimension(:, :, :), pointer :: prodd
    complex(kind=realType), dimension(:, :, :), pointer :: f1d
    complex(kind=realType), dimension(:, :, :), pointer :: kwCDd
    complex(kind=realType), dimension(:, :, :), pointer :: ktCDd
    complex(kind=realType), dimension(:, :, :), pointer :: sctd
    complex(kind=realType), dimension(:, :, :), pointer :: scl2d
    complex(kind=realType), dimension(:, :, :), pointer :: strain2d
#endif
end module turbMod
