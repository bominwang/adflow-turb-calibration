
module BCPointers

! Thiss module contains data structures used to apply BCs.

    use constants, only: intType, realType
       use complexify 
    implicit none
    save

    complex(kind=realType), dimension(:, :, :), pointer :: ww0, ww1, ww2, ww3
    complex(kind=realType), dimension(:, :), pointer :: pp0, pp1, pp2, pp3
    complex(kind=realType), dimension(:, :), pointer :: rlv0, rlv1, rlv2, rlv3
    complex(kind=realType), dimension(:, :), pointer :: rev0, rev1, rev2, rev3
    complex(kind=realType), dimension(:, :), pointer :: gamma0, gamma1, gamma2, gamma3
    complex(kind=realType), dimension(:, :, :), pointer :: ssi, ssj, ssk
    complex(kind=realType), dimension(:, :, :), pointer :: ss, xx
    complex(kind=realType), dimension(:, :), pointer :: dd2wall, sFace
    integer(kind=intType), dimension(:, :), pointer :: gcp

    integer(kind=intType) :: iStart, iEnd, iSize
    integer(kind=intType) :: jStart, jEnd, jSize

#ifndef USE_TAPENADE
    complex(kind=realType), dimension(:, :, :), pointer :: ww0d, ww1d, ww2d, ww3d
    complex(kind=realType), dimension(:, :), pointer :: pp0d, pp1d, pp2d, pp3d
    complex(kind=realType), dimension(:, :), pointer :: rlv0d, rlv1d, rlv2d, rlv3d
    complex(kind=realType), dimension(:, :), pointer :: rev0d, rev1d, rev2d, rev3d
    complex(kind=realType), dimension(:, :), pointer :: gamma0d, gamma1d, gamma2d, gamma3d
    complex(kind=realType), dimension(:, :, :), pointer :: ssid, ssjd, sskd, xxd
    complex(kind=realType), dimension(:, :, :), pointer :: ssd
    complex(kind=realType), dimension(:, :), pointer :: dd2walld, sFaced
#endif

end module BCPointers

#ifndef USE_TAPENADE
module BCPointers_d
    use BCPointers
       use complexify 
 end module BCPointers_d

module BCPointers_b
    use BCPointers
       use complexify 
 end module BCPointers_b

module BCPointers_fast_b
    use BCPointers
       use complexify 
 end module BCPointers_fast_b
#endif
