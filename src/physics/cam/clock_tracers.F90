!===============================================================================
! Age of air "clock" test tracers
! provides dissipation rate for diagnostic constituents
!
! Joe Hollowed
! 04/10/2022
! This module being written as a stripped-down version of aoa_tracers, 
! which enables the advection of just 2 clock tracers rather than 4, and removes
! all surface fluxes
!===============================================================================

module clock_tracers

  use shr_kind_mod, only: r8 => shr_kind_r8
  use spmd_utils,   only: masterproc
  use ppgrid,       only: pcols, pver
  use constituents, only: pcnst, cnst_add, cnst_name, cnst_longname
  use cam_logfile,  only: iulog
  use ref_pres,     only: pref_mid_norm

  implicit none
  private
  save

  ! Public interfaces
  public :: clock_tracers_register         ! register constituents
  public :: clock_tracers_implements_cnst  ! true if named constituent is implemented by this package
  public :: clock_tracers_init_cnst        ! initialize constituent field
  public :: clock_tracers_init             ! initialize history fields, datasets
  public :: clock_tracers_timestep_init    ! place to perform per timestep initialization
  public :: clock_tracers_timestep_tend    ! calculate tendencies
  public :: clock_tracers_readnl           ! read namelist options

  ! Private module data

  integer, parameter :: ncnst=2  ! number of constituents implemented by this module

  ! constituent names
  character(len=8), parameter :: c_names(ncnst) = (/'CLOCK1', 'CLOCK2'/)

  ! constituent source/sink names
  character(len=8), parameter :: src_names(ncnst) = (/'CLOCK1SRC', 'CLOCK2SRC'/)

  integer :: ifirst ! global index of first constituent
  integer :: ixclock1 ! global index for CLOCK1 tracer
  integer :: ixclock2 ! global index for CLOCK2 tracer

  ! Data from namelist variables
  logical :: clock_tracers_flag  = .false.    ! true => turn on test tracer code, namelist variable
  logical :: clock_read_from_ic_file = .true. ! true => tracers initialized from IC file


!===============================================================================
contains
!===============================================================================


subroutine clock_tracers_readnl(nlfile)

    use namelist_utils,     only: find_group_name
    use units,              only: getunit, freeunit
    use mpishorthand
    use cam_abortutils,     only: endrun

    implicit none

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'clock_tracers_readnl'


    namelist /clock_tracers_nl/ clock_tracers_flag, clock_read_from_ic_file

    !-----------------------------------------------------------------------------

    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'clock_tracers_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, clock_tracers_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    call mpibcast(clock_tracers_flag, 1, mpilog,  0, mpicom)
    call mpibcast(clock_read_from_ic_file, 1, mpilog,  0, mpicom)
#endif

  endsubroutine clock_tracers_readnl


!================================================================================


  subroutine clock_tracers_register
    !-----------------------------------------------------------------------
    !
    ! Purpose: register advected constituents
    !
    !-----------------------------------------------------------------------
    use physconst,  only: cpair, mwdry
    !-----------------------------------------------------------------------

    if (.not. clock_tracers_flag) return

    call cnst_add(c_names(1), mwdry, cpair, 0._r8, ixclock1, readiv=clock_read_from_ic_file, &
                  longname='Age-of_air clock tracer 1')
    ifirst = ixclock1
    call cnst_add(c_names(2), mwdry, cpair, 0._r8, ixclock2, readiv=clock_read_from_ic_file, &
                  longname='Age-of_air clock tracer 2')

  end subroutine clock_tracers_register


!===============================================================================


  function clock_tracers_implements_cnst(name)
    !-----------------------------------------------------------------------
    !
    ! Purpose: return true if specified constituent is implemented by this package
    !
    !-----------------------------------------------------------------------

    character(len=*), intent(in) :: name   ! constituent name
    logical :: clock_tracers_implements_cnst        ! return value

    !---------------------------Local workspace-----------------------------
    integer :: m
    !-----------------------------------------------------------------------

    clock_tracers_implements_cnst = .false.

    if (.not. clock_tracers_flag) return

    do m = 1, ncnst
       if (name == c_names(m)) then
          clock_tracers_implements_cnst = .true.
          return
       end if
    end do

  end function clock_tracers_implements_cnst


!===============================================================================


  subroutine clock_tracers_init_cnst(name, latvals, lonvals, mask, q)

    !-----------------------------------------------------------------------
    !
    ! Purpose: initialize test tracers mixing ratio fields
    !  This subroutine is called at the beginning of an initial run ONLY
    !
    !-----------------------------------------------------------------------

    character(len=*), intent(in)  :: name
    real(r8),         intent(in)  :: latvals(:) ! lat in degrees (ncol)
    real(r8),         intent(in)  :: lonvals(:) ! lon in degrees (ncol)
    logical,          intent(in)  :: mask(:)    ! Only initialize where .true.
    real(r8),         intent(out) :: q(:,:)     ! kg tracer/kg dry air (gcol, plev)

    integer :: m
    !-----------------------------------------------------------------------

    if (.not. clock_tracers_flag) return

    do m = 1, ncnst
       if (name ==  c_names(m))  then
          ! pass global constituent index
          call init_cnst_3d(ifirst+m-1, latvals, lonvals, mask, q)
       endif
    end do

  end subroutine clock_tracers_init_cnst


!===============================================================================


  subroutine clock_tracers_init

    !-----------------------------------------------------------------------
    !
    ! Purpose: initialize age of air constituents
    !          (declare history variables)
    !-----------------------------------------------------------------------

    use cam_history,    only: addfld, add_default

    integer :: m, mm, k
    !-----------------------------------------------------------------------

    if (.not. clock_tracers_flag) return

    ! Set names of tendencies and declare them as history variables

    do m = 1, ncnst
       mm = ifirst+m-1
       call addfld(cnst_name(mm), (/ 'lev' /), 'A', 'kg/kg', cnst_longname(mm))
       call addfld(src_names(m),  (/ 'lev' /), 'A', 'kg/kg/s', trim(cnst_name(mm))//' source/sink')

       call add_default (cnst_name(mm), 1, ' ')
       call add_default (src_names(m),  1, ' ')
    end do

  end subroutine clock_tracers_init


!===============================================================================


  subroutine clock_tracers_timestep_tend(state, ptend, cflx, landfrac, dt)

    use physics_types, only: physics_state, physics_ptend, physics_ptend_init
    use cam_history,   only: outfld
    use time_manager,  only: get_nstep
    
    !--JH--
    use ref_pres,     only: pref_mid_norm
    use time_manager, only: get_curr_time

    ! Arguments
    type(physics_state), intent(inout) :: state           ! state variables --JH--
    type(physics_ptend), intent(out)   :: ptend              ! package tendencies
    real(r8),            intent(inout) :: cflx(pcols,pcnst)  ! Surface constituent flux (kg/m^2/s)
    real(r8),            intent(in)    :: landfrac(pcols)    ! Land fraction
    real(r8),            intent(in)    :: dt                 ! timestep

    !----------------- Local workspace-------------------------------

    integer :: i, k
    integer :: lchnk                          ! chunk identifier
    integer :: ncol                           ! no. of column in chunk
    integer :: nstep                          ! current timestep number
    logical  :: lq(pcnst)
    
    !--JH--
    integer  :: day,sec           ! date variables
    real(r8) :: t                 ! tracer boundary condition
    real(r8) :: clock1_scaling    ! scale CLOCK1 from nstep to time
    real(r8) :: clock2_tau        ! relaxation timescale for CLOCK2
    
    !------------------------------------------------------------------

    if (.not. clock_tracers_flag) then
       call physics_ptend_init(ptend,state%psetcols,'none') !Initialize an empty ptend for use with physics_update
       return
    end if

    lq(:)      = .FALSE.
    lq(ixclock1) = .TRUE.
    lq(ixclock2) = .TRUE.
    call physics_ptend_init(ptend,state%psetcols, 'clock_tracers', lq=lq)

    nstep = get_nstep()
    lchnk = state%lchnk
    ncol  = state%ncol
    
    !--JH--
    ! get current time in days
    call get_curr_time(day,sec)
    t=day + sec/86400.0              ! current time in days
    !compute CLOCK1 time scaling
    clock1_scaling = 1/86400.0
    clock2_tau = (1.0/10.0) * 86400.0  ! 1/10 day in seconds

    do k = 1, pver
       do i = 1, ncol

          ! =======================================
          !
          ! UPDATE TENDENCIES/TRACER CONCENTRATION
          !
          ! =======================================

          ! CLOCK1
          ! --JH--: This tracer will be used as a clock tracer with a source of
          ! 1 everywhere above ~700hPa
          if (pref_mid_norm(k) <= 0.7) then
              ptend%q(i,k,ixclock1) = 1.0_r8 * clock1_scaling
          else
              ptend%q(i,k,ixclock1) = 0.0_r8
              state%q(i,k,ixclock1) = 0.0_r8
          end if

          ! CLOCK2
          ! --JH--: This tracer will be used as a clock tracer which assumes the
          ! value of the model time eveywhere below ~700hPa
          if (pref_mid_norm(k) >= 0.7) then
              ptend%q(i,k,ixclock2) = 0.0_r8
              state%q(i,k,ixclock2) = t
              !ptend%q(i,k,ixclock2) = (t-state%q(i,k,ixclock2)) / clock2_tau 
                    ! ^ relaxation rather than eplicitly set; not used
          else
              ptend%q(i,k,ixclock2) = 0.0_r8
          end if

       end do
    end do

    ! record tendencies on history files
    call outfld (src_names(1), ptend%q(:,:,ixclock1), pcols, lchnk)
    call outfld (src_names(2), ptend%q(:,:,ixclock2), pcols, lchnk)

    ! Set tracer surface fluxes to zero
    do i = 1, ncol
       ! CLOCK1
       cflx(i,ixclock1) = 1.e-6_r8
       ! CLOCK2
       cflx(i,ixclock2) = 1.e-6_r8
    end do

  end subroutine clock_tracers_timestep_tend


!===========================================================================


  subroutine init_cnst_3d(m, latvals, lonvals, mask, q)

    integer,  intent(in)  :: m          ! global constituent index
    real(r8), intent(in)  :: latvals(:) ! lat in degrees (ncol)
    real(r8), intent(in)  :: lonvals(:) ! lon in degrees (ncol)
    logical,  intent(in)  :: mask(:)    ! Only initialize where .true.
    real(r8), intent(out) :: q(:,:)     ! kg tracer/kg dry air (gcol,plev)

    integer :: j, k, gsize
    !-----------------------------------------------------------------------

    if (masterproc) then
      write(iulog,*) 'AGE-OF-AIR CONSTITUENTS: INITIALIZING ',cnst_name(m),m
    end if

    if (m == ixclock1) then

       q(:,:) = 0.0_r8

    else if (m == ixclock2) then

       q(:,:) = 0.0_r8

    end if

  end subroutine init_cnst_3d


!=====================================================================


end module clock_tracers
