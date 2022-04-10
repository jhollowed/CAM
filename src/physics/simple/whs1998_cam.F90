module whs1998_cam

    !-----------------------------------------------------------------------
    !
    ! --JH--
    !
    ! Purpose: implement forcingd from Modified Held/Suarez IDEALIZED physics algorithm
    !          (modified with Williamson stratosphere):
    !
    !   Williamson, D. L., J. G. Olson and B. A. Boville, 1998: A comparison
    !   of semi--Lagrangian and Eulerian tropical climate simulations.
    !   Mon. Wea. Rev., vol 126, pp. 1001-1012.
    !
    !-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pcols, pver

  implicit none
  private
  save

  public :: whs1998_init, whs1998_tend

  real(r8), parameter :: efoldf  =  1._r8  ! efolding time for wind dissipation
  real(r8), parameter :: efolda  = 40._r8  ! efolding time for T dissipation
  real(r8), parameter :: efolds  =  4._r8  ! efolding time for T dissipation
  real(r8), parameter :: sigmab  =  0.7_r8 ! threshold sigma level
  real(r8), parameter :: t00     = 200._r8 ! minimum reference temperature
  real(r8), parameter :: kf      = 1._r8/(86400._r8*efoldf) ! 1./efolding_time for wind dissipation

  real(r8), parameter :: onemsig = 1._r8 - sigmab ! 1. - sigma_reference

  real(r8), parameter :: ka      = 1._r8/(86400._r8 * efolda) ! 1./efolding_time for temperature diss.
  real(r8), parameter :: ks      = 1._r8/(86400._r8 * efolds)

!======================================================================= 
contains
!======================================================================= 

  subroutine whs1998_init(pbuf2d)
    use physics_buffer,     only: physics_buffer_desc
    use cam_history,        only: addfld, add_default
    use physconst,          only: cappa, cpair
    use ref_pres,           only: pref_mid_norm, psurf_ref
    use whs1998,            only: whs1998_mod_init

    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    character(len=512)            :: errmsg
    integer                       :: errflg

    ! Set model constant values
    call whs1998_mod_init(pver, cappa, cpair, psurf_ref, pref_mid_norm, errmsg, errflg)

    ! This field is added by radiation when full physics is used
    call addfld('QRS', (/ 'lev' /), 'A', 'K/s', &
         'Temperature tendency associated with the relaxation toward the equilibrium temperature profile')
    call add_default('QRS', 1, ' ')
 end subroutine whs1998_init

  subroutine whs1998_tend(state, ptend, ztodt)
    use physconst,          only: cpairv
    use phys_grid,          only: get_rlat_all_p
    use physics_types,      only: physics_state, physics_ptend
    use physics_types,      only: physics_ptend_init
    use cam_abortutils,     only: endrun
    use cam_history,        only: outfld
    use whs1998,   only: whs1998_mod_run

    !
    ! Input arguments
    !
    type(physics_state), intent(inout) :: state
    real(r8),            intent(in)    :: ztodt            ! Two times model timestep (2 delta-t)
                                                           !
                                                           ! Output argument
                                                           !
    type(physics_ptend), intent(out)   :: ptend            ! Package tendencies
                                                           !
    !---------------------------Local workspace-----------------------------

    integer                            :: lchnk            ! chunk identifier
    integer                            :: ncol             ! number of atmospheric columns

    real(r8)                           :: clat(pcols)      ! latitudes(radians) for columns
    real(r8)                           :: pmid(pcols,pver) ! mid-point pressure
    integer                            :: i, k             ! Longitude, level indices

    character(len=512)            :: errmsg
    integer                       :: errflg

    !
    !-----------------------------------------------------------------------
    !

    lchnk = state%lchnk
    ncol  = state%ncol

    call get_rlat_all_p(lchnk, ncol, clat)
    do k = 1, pver
      do i = 1, ncol
        pmid(i,k) = state%pmid(i,k)
      end do
    end do

    ! initialize individual parameterization tendencies
    call physics_ptend_init(ptend, state%psetcols, 'whs1998', ls=.true., lu=.true., lv=.true.)

    call whs1998_mod_run(pver, ncol, clat, state%pmid, &
         state%u, state%v, state%t, ptend%u, ptend%v, ptend%s, errmsg, errflg)

    ! Note, we assume that there are no subcolumns in simple physics
    pmid(:ncol,:) = ptend%s(:ncol, :)/cpairv(:ncol,:,lchnk)
    if (pcols > ncol) then
      pmid(ncol+1:,:) = 0.0_r8
    end if
    call outfld('QRS', pmid, pcols, lchnk)

  end subroutine whs1998_tend

end module whs1998_cam
