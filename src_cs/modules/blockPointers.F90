module blockPointers
    !
    !       This module contains the pointers for all variables inside a
    !       block. The pointers are set via the subroutine setPointers,
    !       which can be found in the utils directory. In this way the
    !       code becomes much more readable. The relation to the original
    !       multiblock grid is not copied, because it does not affect the
    !       computation.
    !       See the module block for the meaning of the variables.
    !       Note that the dimensions are not pointers, but integers.
    !       Consequently changing dimensions of a block must be done only
    !       with the variables of floDoms.
    !
    use constants, only: intType, realType, porType
    use block, only: fringeType, BCDataType, viscSubFaceType, flowDoms, nDom
       use complexify 
#ifndef USE_TAPENADE
    use block, only: flowDomsd
#endif
    implicit none
    !
    !       Additional info, such that it is known to which block the data
    !       inside this module belongs.
    !
    ! sectionID:   the section to which this block belongs.
    ! nbkLocal :   local block number.
    ! nbkGlobal:   global block number in the original cgns grid.
    ! mgLevel:     the multigrid level.
    ! spectralSol: the spectral solution index of this block.

    integer(kind=intType) :: sectionID
    integer(kind=intType) :: nbkLocal, nbkGlobal, mgLevel
    integer(kind=intType) :: spectralSol
    !
    !       Variables, which are either copied or the pointer is set to
    !       the correct variable in the block. See the module block for
    !       meaning of the variables.
    !
    integer(kind=intType) :: nx, ny, nz, il, jl, kl
    integer(kind=intType) :: ie, je, ke, ib, jb, kb
    integer(kind=intType) :: maxDim, imaxDim, jmaxDim

    logical :: rightHanded

    integer(kind=intType) :: iBegOr, iEndOr, jBegOr, jEndOr
    integer(kind=intType) :: kBegOr, kEndOr

    integer(kind=intType) :: nSubface, n1to1, nBocos, nViscBocos

    integer(kind=intType), dimension(:), pointer :: BCType
    integer(kind=intType), dimension(:), pointer :: BCFaceID

    integer(kind=intType), dimension(:), pointer :: cgnsSubface

    integer(kind=intType), dimension(:), pointer :: inBeg, inEnd
    integer(kind=intType), dimension(:), pointer :: jnBeg, jnEnd
    integer(kind=intType), dimension(:), pointer :: knBeg, knEnd

    integer(kind=intType), dimension(:), pointer :: dinBeg, dinEnd
    integer(kind=intType), dimension(:), pointer :: djnBeg, djnEnd
    integer(kind=intType), dimension(:), pointer :: dknBeg, dknEnd

    integer(kind=intType), dimension(:), pointer :: icBeg, icEnd
    integer(kind=intType), dimension(:), pointer :: jcBeg, jcEnd
    integer(kind=intType), dimension(:), pointer :: kcBeg, kcEnd

    integer(kind=intType), dimension(:), pointer :: neighBlock
    integer(kind=intType), dimension(:), pointer :: neighProc
    integer(kind=intType), dimension(:), pointer :: l1, l2, l3
    integer(kind=intType), dimension(:), pointer :: groupNum

    integer(kind=intType), dimension(:, :, :), pointer :: iblank
    integer(kind=intType), dimension(:, :, :), pointer :: status
    integer(kind=intType), dimension(:, :, :), pointer :: forcedRecv
    type(fringeType), dimension(:), pointer :: fringes
    integer(kind=intType), pointer :: nDonors
    integer(kind=intType), dimension(:, :, :, :), pointer :: fringePtr
    integer(kind=intType), dimension(:, :, :, :), pointer :: gInd
    integer(kind=intType), dimension(:, :), pointer :: orphans
    integer(kind=intType) :: nOrphans

    integer(kind=intType), dimension(:), pointer :: neighBlockOver
    integer(kind=intType), dimension(:), pointer :: neighProcOver

    type(BCDataType), dimension(:), pointer :: BCData
    type(viscSubfaceType), dimension(:), pointer :: viscSubface

    integer(kind=intType), dimension(:, :), pointer :: viscIMinPointer
    integer(kind=intType), dimension(:, :), pointer :: viscIMaxPointer
    integer(kind=intType), dimension(:, :), pointer :: viscJMinPointer
    integer(kind=intType), dimension(:, :), pointer :: viscJMaxPointer
    integer(kind=intType), dimension(:, :), pointer :: viscKMinPointer
    integer(kind=intType), dimension(:, :), pointer :: viscKMaxPointer

    complex(kind=realType), dimension(:, :, :, :), pointer :: x
    complex(kind=realType), dimension(:, :, :, :, :), pointer :: xOld
    complex(kind=realType), dimension(:, :, :, :), pointer :: sI, sJ, sK
    complex(kind=realType), dimension(:, :, :), pointer :: vol
    complex(kind=realType), dimension(:, :, :), pointer :: volref
    complex(kind=realType), dimension(:, :, :, :), pointer :: volOld
    complex(kind=realType), dimension(:, :, :), pointer :: skew
    complex(kind=realType), dimension(:, :, :, :), pointer :: dadidata

    integer(kind=porType), dimension(:, :, :), pointer :: porI, porJ, porK

    integer(kind=intType), dimension(:, :, :), pointer :: indFamilyI
    integer(kind=intType), dimension(:, :, :), pointer :: indFamilyJ
    integer(kind=intType), dimension(:, :, :), pointer :: indFamilyK

    integer(kind=intType), dimension(:, :, :), pointer :: factFamilyI
    integer(kind=intType), dimension(:, :, :), pointer :: factFamilyJ
    integer(kind=intType), dimension(:, :, :), pointer :: factFamilyK

    complex(kind=realType), dimension(:, :, :, :, :), pointer :: rotMatrixI
    complex(kind=realType), dimension(:, :, :, :, :), pointer :: rotMatrixJ
    complex(kind=realType), dimension(:, :, :, :, :), pointer :: rotMatrixK

    logical :: blockIsMoving, addGridVelocities

    complex(kind=realType), dimension(:, :, :), pointer :: sFaceI, sFaceJ, sfaceK
    complex(kind=realType), dimension(:, :, :, :), pointer :: w
    complex(kind=realType), dimension(:, :, :, :, :), pointer :: wOld

    complex(kind=realType), dimension(:, :, :), pointer :: p, gamma, aa
    complex(kind=realType), dimension(:, :, :), pointer :: shockSensor
    complex(kind=realType), dimension(:, :, :), pointer :: rlv, rev
    complex(kind=realType), dimension(:, :, :, :), pointer :: s
    complex(kind=realType), dimension(:, :, :), pointer :: p1
    complex(kind=realType), dimension(:, :, :, :), pointer :: dw, fw
    complex(kind=realType), dimension(:, :, :, :), pointer :: scratch
    complex(kind=realType), dimension(:, :, :, :, :), pointer :: dwOldRK
    complex(kind=realType), dimension(:, :, :, :), pointer :: w1, wr
    complex(kind=realType), dimension(:, :, :), pointer :: ux, uy, uz
    complex(kind=realType), dimension(:, :, :), pointer :: vx, vy, vz
    complex(kind=realType), dimension(:, :, :), pointer :: wx, wy, wz
    complex(kind=realType), dimension(:, :, :), pointer :: qx, qy, qz

    integer(kind=intType), dimension(:, :), pointer :: mgIFine
    integer(kind=intType), dimension(:, :), pointer :: mgJFine
    integer(kind=intType), dimension(:, :), pointer :: mgKFine

    complex(kind=realType), dimension(:), pointer :: mgIWeight
    complex(kind=realType), dimension(:), pointer :: mgJWeight
    complex(kind=realType), dimension(:), pointer :: mgKWeight

    integer(kind=intType), dimension(:, :), pointer :: mgICoarse
    integer(kind=intType), dimension(:, :), pointer :: mgJCoarse
    integer(kind=intType), dimension(:, :), pointer :: mgKCoarse

    complex(kind=realType), dimension(:, :, :, :), pointer :: wn
    complex(kind=realType), dimension(:, :, :), pointer :: pn
    complex(kind=realType), dimension(:, :, :), pointer :: dtl
    complex(kind=realType), dimension(:, :, :), pointer :: radI, radJ, radK

    complex(kind=realType), dimension(:, :, :), pointer :: d2Wall
    complex(kind=realType), dimension(:, :, :), pointer :: intermittency
    complex(kind=realType), dimension(:, :, :), pointer :: filterDES  ! eran-des
    complex(kind=realType), dimension(:, :, :, :), pointer :: bmti1
    complex(kind=realType), dimension(:, :, :, :), pointer :: bmti2
    complex(kind=realType), dimension(:, :, :, :), pointer :: bmtj1
    complex(kind=realType), dimension(:, :, :, :), pointer :: bmtj2
    complex(kind=realType), dimension(:, :, :, :), pointer :: bmtk1
    complex(kind=realType), dimension(:, :, :, :), pointer :: bmtk2
    complex(kind=realType), dimension(:, :, :), pointer :: bvti1, bvti2
    complex(kind=realType), dimension(:, :, :), pointer :: bvtj1, bvtj2
    complex(kind=realType), dimension(:, :, :), pointer :: bvtk1, bvtk2

    integer(kind=intType), dimension(:, :, :), pointer :: globalNode
    integer(kind=intType), dimension(:, :, :), pointer :: globalCell
    complex(kind=realType), dimension(:, :, :, :), pointer :: xSeed
    integer(kind=intType), dimension(:, :, :), pointer :: wallInd

    complex(kind=realType), dimension(:, :, :, :), pointer :: w_offTimeInstance
    complex(kind=realType), dimension(:, :, :), pointer :: vol_offTimeInstance

    ! Added by HDN
    complex(kind=realType), dimension(:, :, :, :), pointer :: xALE
    complex(kind=realType), dimension(:, :, :, :), pointer :: sVeloIALE, sVeloJALE, sVeloKALE
    complex(kind=realType), dimension(:, :, :, :, :), pointer :: sIALE, sJALE, sKALE
    complex(kind=realType), dimension(:, :, :, :), pointer :: sFaceIALE, sFaceJALE, sFaceKALE
    complex(kind=realType), dimension(:, :, :, :, :), pointer :: dwALE, fwALE

#ifndef USE_TAPENADE
    TYPE(VISCSUBFACETYPE), DIMENSION(:), POINTER :: viscsubfaced

    complex(kind=realtype), DIMENSION(:, :, :, :), POINTER :: xd
    complex(kind=realtype), DIMENSION(:, :, :, :), POINTER :: sid, sjd, skd

    complex(kind=realType), dimension(:, :, :), pointer :: vold

    complex(kind=realtype), DIMENSION(:, :, :, :, :), POINTER :: rotmatrixid
    complex(kind=realtype), DIMENSION(:, :, :, :, :), POINTER :: rotmatrixjd
    complex(kind=realtype), DIMENSION(:, :, :, :, :), POINTER :: rotmatrixkd

    complex(kind=realtype), DIMENSION(:, :, :), POINTER :: sfaceid
    complex(kind=realtype), DIMENSION(:, :, :), POINTER :: sfacejd
    complex(kind=realtype), DIMENSION(:, :, :), POINTER :: sfacekd

    complex(kind=realtype), DIMENSION(:, :, :, :), POINTER :: wd
    complex(kind=realtype), DIMENSION(:, :, :, :, :), POINTER :: woldd

    complex(kind=realtype), DIMENSION(:, :, :), POINTER :: uxd
    complex(kind=realtype), DIMENSION(:, :, :), POINTER :: uyd
    complex(kind=realtype), DIMENSION(:, :, :), POINTER :: uzd

    complex(kind=realtype), DIMENSION(:, :, :), POINTER :: vxd
    complex(kind=realtype), DIMENSION(:, :, :), POINTER :: vyd
    complex(kind=realtype), DIMENSION(:, :, :), POINTER :: vzd

    complex(kind=realtype), DIMENSION(:, :, :), POINTER :: wxd
    complex(kind=realtype), DIMENSION(:, :, :), POINTER :: wyd
    complex(kind=realtype), DIMENSION(:, :, :), POINTER :: wzd

    complex(kind=realtype), DIMENSION(:, :, :), POINTER :: qxd
    complex(kind=realtype), DIMENSION(:, :, :), POINTER :: qyd
    complex(kind=realtype), DIMENSION(:, :, :), POINTER :: qzd

    complex(kind=realtype), DIMENSION(:, :, :), POINTER :: pd, gammad, aad
    complex(kind=realtype), DIMENSION(:, :, :), POINTER :: rlvd, revd

    complex(kind=realtype), DIMENSION(:, :, :, :), POINTER :: sd

    complex(kind=realtype), DIMENSION(:, :, :, :), POINTER :: dwd, fwd
    complex(kind=realtype), DIMENSION(:, :, :, :), POINTER :: w1d, wrd
    complex(kind=realtype), DIMENSION(:, :, :, :), POINTER :: scratchd

    complex(kind=realtype), DIMENSION(:, :, :), POINTER :: dtld
    complex(kind=realtype), DIMENSION(:, :, :), POINTER :: radid
    complex(kind=realtype), DIMENSION(:, :, :), POINTER :: radjd
    complex(kind=realtype), DIMENSION(:, :, :), POINTER :: radkd

    complex(kind=realtype), DIMENSION(:, :, :, :), POINTER :: bmti1d
    complex(kind=realtype), DIMENSION(:, :, :, :), POINTER :: bmti2d
    complex(kind=realtype), DIMENSION(:, :, :, :), POINTER :: bmtj1d
    complex(kind=realtype), DIMENSION(:, :, :, :), POINTER :: bmtj2d
    complex(kind=realtype), DIMENSION(:, :, :, :), POINTER :: bmtk1d
    complex(kind=realtype), DIMENSION(:, :, :, :), POINTER :: bmtk2d

    complex(kind=realtype), DIMENSION(:, :, :), POINTER :: bvti1d, bvti2d
    complex(kind=realtype), DIMENSION(:, :, :), POINTER :: bvtj1d, bvtj2d
    complex(kind=realtype), DIMENSION(:, :, :), POINTER :: bvtk1d, bvtk2d

    complex(kind=realtype), DIMENSION(:, :, :), POINTER :: d2walld

    complex(kind=realType), dimension(:, :, :, :), pointer :: w_offTimeInstanced
    complex(kind=realType), dimension(:, :, :), pointer :: vol_offTimeInstanced

    type(BCDataType), dimension(:), pointer :: BCDatad

    complex(kind=realType), dimension(:, :, :, :, :), pointer :: PCMat
    complex(kind=realType), dimension(:, :, :, :), pointer :: PCvec1, PCvec2

    complex(kind=realType), dimension(:, :, :, :), pointer :: i_D_fact, j_D_fact, k_D_fact
    complex(kind=realType), dimension(:, :, :, :), pointer :: i_L_Fact, j_L_Fact, k_L_Fact
    complex(kind=realType), dimension(:, :, :, :), pointer :: i_U_Fact, j_U_Fact, k_U_Fact
    complex(kind=realType), dimension(:, :, :, :), pointer :: i_U2_Fact, j_U2_Fact, k_U2_Fact
    integer(kind=intType), dimension(:, :, :, :), pointer :: i_ipiv, j_ipiv, k_ipiv
#endif

end module blockPointers
