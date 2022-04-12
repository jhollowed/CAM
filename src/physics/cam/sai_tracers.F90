!===============================================================================
! Stratospheric aerosol injection test tracers
! provides dissipation rate for diagnostic constituents
!
! Joe Hollowed
! 04/10/2022
! This module being written based on the structure of aoa_tracers.
! Enables the advection of 5 tracer species:
! - SAI_SO2:    SO2 species coincident with plume
! - SAI_ASH:    ash species coincident with plume
! - SAI_THETA:  uniform potential temperature dynamic tracer
! - SAI_PV:     uniform potential vorticity dynamic tracer
! - SAI_CLOCK:  sai tracer pulse coincident with plume
!===============================================================================

module sai_tracers

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
  public :: sai_tracers_register         ! register constituents
  public :: sai_tracers_implements_cnst  ! true if named constituent is implemented by this package
  public :: sai_tracers_init_cnst        ! initialize constituent field
  public :: sai_tracers_init             ! initialize history fields, datasets
  public :: sai_tracers_timestep_init    ! place to perform per timestep initialization
  public :: sai_tracers_timestep_tend    ! calculate tendencies
  public :: sai_tracers_readnl           ! read namelist options

  ! Private module data

  integer, parameter :: ncnst=5  ! number of constituents implemented by this module

  ! constituent names
  character(len=8), parameter :: c_names(ncnst) = (/'SAI_SO2', 'SAI_ASH',&
                                                    'SAI_THETA', 'SAI_PV', 'SAI_CLOCK'/)

  integer :: ifirst ! global index of first constituent
  integer :: ixsai1 ! global index for SAI_SO2 tracer
  integer :: ixsai2 ! global index for SAI_ASH tracer
  integer :: ixsai3 ! global index for SAI_THETA tracer
  integer :: ixsai4 ! global index for SAI_PV tracer
  integer :: ixsai5 ! global index for SAI_CLOCK tracer

  ! Data from namelist variables
  logical  :: sai_tracers_flag  = .false.      ! true => turn on test tracer code, namelist variable
  logical  :: sai_read_from_ic_file  = .false. ! this is not updated, and is always false (only analytic)
  real(r8) :: sai_duration                     ! sai duration, in hours
  integer  :: sai_decay_type                   ! sai decay type, either exponential (1) or flat (2)
  real(r8) :: sai_so2_efold                    ! e-folding timescale for SO2
  real(r8) :: sai_ash_efold                    ! e-folding timescale for ash


!===============================================================================
contains
!===============================================================================


subroutine sai_tracers_readnl(nlfile)

    use namelist_utils,     only: find_group_name
    use units,              only: getunit, freeunit
    use mpishorthand
    use cam_abortutils,     only: endrun

    implicit none

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'sai_tracers_readnl'


    namelist /sai_tracers_nl/ sai_tracers_flag, sai_duration, sai_decay_type, &
                              sai_so2_efold, sai_ash_efold

    !-----------------------------------------------------------------------------

    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'sai_tracers_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, sai_tracers_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    call mpibcast(sai_tracers_flag, 1, mpilog,  0, mpicom)
    call mpibcast(sai_duration, 1, mpilog,  0, mpicom)
    call mpibcast(sai_decay_type, 1, mpilog,  0, mpicom)
    call mpibcast(sai_so2_efold, 1, mpilog,  0, mpicom)
    call mpibcast(sai_ash_efold, 1, mpilog,  0, mpicom)
#endif

  endsubroutine sai_tracers_readnl


!================================================================================


  subroutine sai_tracers_register
    !-----------------------------------------------------------------------
    !
    ! Purpose: register advected constituents
    !
    !-----------------------------------------------------------------------
    use physconst,  only: cpair, mwdry
    !-----------------------------------------------------------------------

    if (.not. sai_tracers_flag) return

    call cnst_add(c_names(1), mwdry, cpair, 0._r8, ixsai1, readiv=sai_read_from_ic_file, &
                  longname='SAI SO2 tracer')
    ifirst = ixsai1
    call cnst_add(c_names(2), mwdry, cpair, 0._r8, ixsai2, readiv=sai_read_from_ic_file, &
                  longname='SAI ash tracer')
    call cnst_add(c_names(3), mwdry, cpair, 0._r8, ixsai3, readiv=sai_read_from_ic_file, &
                  longname='SAI potential temperature tracer')
    call cnst_add(c_names(4), mwdry, cpair, 0._r8, ixsai4, readiv=sai_read_from_ic_file, &
                  longname='SAI potential vorticity tracer')
    call cnst_add(c_names(5), mwdry, cpair, 0._r8, ixsai5, readiv=sai_read_from_ic_file, &
                  longname='SAI clock pulse tracer')

  end subroutine sai_tracers_register


!===============================================================================


  function sai_tracers_implements_cnst(name)
    !-----------------------------------------------------------------------
    !
    ! Purpose: return true if specified constituent is implemented by this package
    !
    !-----------------------------------------------------------------------

    character(len=*), intent(in) :: name   ! constituent name
    logical :: sai_tracers_implements_cnst        ! return value

    !---------------------------Local workspace-----------------------------
    integer :: m
    !-----------------------------------------------------------------------

    sai_tracers_implements_cnst = .false.

    if (.not. sai_tracers_flag) return

    do m = 1, ncnst
       if (name == c_names(m)) then
          sai_tracers_implements_cnst = .true.
          return
       end if
    end do

  end function sai_tracers_implements_cnst


!===============================================================================


  subroutine sai_tracers_init_cnst(name, latvals, lonvals, mask, q)

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

    if (.not. sai_tracers_flag) return

    do m = 1, ncnst
       if (name ==  c_names(m))  then
          ! pass global constituent index
          call init_cnst_3d(ifirst+m-1, latvals, lonvals, mask, q)
       endif
    end do

  end subroutine sai_tracers_init_cnst


!===============================================================================


  subroutine sai_tracers_init

    !-----------------------------------------------------------------------
    !
    ! Purpose: initialize age of air constituents
    !          (declare history variables)
    !-----------------------------------------------------------------------

    use cam_history,    only: addfld, add_default

    integer :: m, mm, k
    !-----------------------------------------------------------------------

    if (.not. sai_tracers_flag) return

    ! Set names of tendencies and declare them as history variables

    do m = 1, ncnst
       mm = ifirst+m-1
       call addfld(cnst_name(mm), (/ 'lev' /), 'A', 'kg/kg', cnst_longname(mm))

       call add_default (cnst_name(mm), 1, ' ')
    end do

  end subroutine sai_tracers_init


!===============================================================================


  subroutine sai_tracers_timestep_tend(state, ptend, cflx, landfrac, dt)

    use physics_types, only: physics_state, physics_ptend, physics_ptend_init
    use cam_history,   only: outfld
    use time_manager,  only: get_nstep
    
    !--JH--
    use ref_pres,     only: pref_mid_norm
    use time_manager, only: get_curr_time
    use physconst,    only: pi, rearth, cpair, rair 

    ! Arguments
    type(physics_state), intent(inout) :: state              ! state variables --JH--
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
    real(r8), parameter :: deg2rad = pi/180._r8
    real(r8), parameter :: rad2deg = 180._r8/pi
    integer  :: day,sec
    real(r8) :: td                 
    real(r8) :: ts                 
    real(r8) :: lat                 
    real(r8) :: lon                 
    real(r8) :: zz                 
    real(r8) :: plat                 
    real(r8) :: plon                 
    real(r8) :: pz                 
    real(r8) :: pzm                 
    real(r8) :: pxm                 
    real(r8) :: ptau                 
    real(r8) :: pdur             
    real(r8) :: pkSO2                 
    real(r8) :: pkAsh                 
    real(r8) :: pASO2
    real(r8) :: pAAsh
    real(r8) :: P0             
    
    !------------------------------------------------------------------

    if (.not. sai_tracers_flag) then
       call physics_ptend_init(ptend,state%psetcols,'none') !Initialize an empty ptend for use with physics_update
       return
    end if

    lq(:)      = .FALSE.
    lq(ixsai1) = .TRUE.
    lq(ixsai2) = .TRUE.
    lq(ixsai3) = .TRUE.
    lq(ixsai4) = .TRUE.
    lq(ixsai5) = .TRUE.
    call physics_ptend_init(ptend,state%psetcols, 'sai_tracers', lq=lq)

    nstep = get_nstep()
    lchnk = state%lchnk
    ncol  = state%ncol
    
    !--JH--
    call get_curr_time(day,sec)
    td = day + sec/86400.0              ! current time in days
    ts = (day*24*60*60) + sec          ! current time in s
    
    !-------- compute SAI_CLOCK time scaling --------
    clock_scaling = 1/86400.0
    
    ! -------- SAI plume parameters --------
    plat   = 15.15 * deg2rad       ! meridional center of plume in rad
    plon   = 120.35 * deg2rad      ! zonal center of plume in rad
    pz     = 25000.0               ! vertical center of plume in meters
    pzm    = 15000.0               ! width of plume in vertical in meters
    pxm    = 200000.0              ! width of plume in horizonal in meters
    pASO2  = 3000.0                ! SO2 amplitude normalization
    pAAsh  = 3000.0                ! Ash amplitude normalization
    pkSO2  =  1.0/sai_so2_efold    ! timescale for SO2 removal
    pkAsh  = 1.0/sai_ash_efold     ! timescale for ash removal
    pdur   = sai_duration * 3600   ! injection duration in seconds
    ptau   = -LOG(0.05)/pdur       ! exponential time decay constant
    P0  = 100000                   ! reference pressure

    do k = 1, pver
       do i = 1, ncol

          ! =======================================
          !
          ! UPDATE TENDENCIES/TRACER CONCENTRATION
          !
          ! =======================================
          
          lat = state%lat(i)    ! latitude of column in rad
          lon = state%lon(i)    ! longitude of column in rad
          zz = state%zm(i, k)   ! geopotential height above surface at column midpoints in meters


          ! =============================                        
          ! ========== SAI_SO2 ==========
          ptend%q(i,k,ixsai1) = -pkSO2 * state%q(i, k, ixaoa1) + &
                                EXP(-(1.0_r8/2.0_r8) * &
                                     (((lat-plat)/(pxm/(4.0_r8*rearth)))**2.0_r8 + &
                                      ((lon-plon)/(pxm/(4.0_r8*rearth)))**2.0_r8)) * &
                                EXP(-(1.0_r8/2.0_r8) * &
                                     ((zz-pz)/(pzm/4.0_r8))**2.0_r8)
          ! ---- add time dependency, based on chosen decay type
          if (sai_decay_type == 1) then
              ! ---- exponential decay
              ptend%q(i,k,ixsai1) = ptend%q(i,k,ixsai1) * EXP(-ptau * t)
          else if (sai_decay_type == 2) then
              ! ---- constant amplitude, zero after duration elapsed
              if(ts > tdur):
                  ptend%q(i,k,ixsai1) = 0
              end if
          end if
          ! ---- normalize plume
          ptend%q(i,k,ixsai1) = ptend%q(i,k,ixsai1) * pASO2
          

          ! =============================                        
          ! ========== SAI_ASH ==========
          ptend%q(i,k,ixsai2) = -pkAsh * state%q(i, k, ixaoa2) + &
                                EXP(-(1.0_r8/2.0_r8) * &
                                     (((lat-plat)/(pxm/(4.0_r8*rearth)))**2.0_r8 + &
                                      ((lon-plon)/(pxm/(4.0_r8*rearth)))**2.0_r8)) * &
                                EXP(-(1.0_r8/2.0_r8) * &
                                     ((zz-pz)/(pzm/4.0_r8))**2.0_r8)
          ! ---- add time dependency, based on chosen decay type
          if (sai_decay_type == 1) then
              ! ---- exponential decay
              ptend%q(i,k,ixsai2) = ptend%q(i,k,ixsai2) * EXP(-ptau * t)
          else if (sai_decay_type == 2) then
              ! ---- constant amplitude, zero after duration elapsed
              if(ts > tdur):
                  ptend%q(i,k,ixsai1) = 0
              end if
          end if
          ! ---- normalize plume
          ptend%q(i,k,ixsai2) = ptend%q(i,k,ixsai2) * pAAsh
          

          ! ===============================                        
          ! ========== SAI_THETA ==========
          ! Potential Temperature
          ! initialize within the first minute of the injection
          if (t < 60.0_r8) then
              state%q(i,k,ixsai3) = state%t(i,k) * (P0 / state%pmid(i, k))**(rair/cpair)
          end if
          ptend%q(i,k,ixsai3) = 0.0_r8
                                      

          ! ============================                        
          ! ========== SAI_PV ==========
          ! Potential vorticity
          ptend%q(i,k,ixsai4) = 0.0_r8

          
          ! --=============================                        
          ! ========== SAI_CLOCK ==========
          ! clock tracer with a source of 1 everywhere above ~700hPa
          if (pref_mid_norm(k) <= 0.7) then
              ptend%q(i,k,ixsai5) = 1.0_r8 * sai1_scaling
          else
              ptend%q(i,k,ixsai5) = 0.0_r8
              state%q(i,k,ixsai5) = 0.0_r8
          end if


       end do
    end do

    ! Set tracer surface fluxes to zero
    do i = 1, ncol
       cflx(i,ixsai1) = 0.0_r8
       cflx(i,ixsai2) = 0.0_r8
       cflx(i,ixsai3) = 0.0_r8
       cflx(i,ixsai4) = 0.0_r8
       cflx(i,ixsai5) = 0.0_r8
    end do

  end subroutine sai_tracers_timestep_tend


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

    if (m == ixsai1) then

       q(:,:) = 0.0_r8

    else if (m == ixsai2) then

       q(:,:) = 0.0_r8

    end if

  end subroutine init_cnst_3d


!=====================================================================


end module sai_tracers
