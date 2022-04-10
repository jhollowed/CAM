module ic_whs1998

  !-----------------------------------------------------------------------
  !
  ! Purpose: Set Held-Suarez initial conditions based on input coordinates
  !
  !
  !-----------------------------------------------------------------------
  use cam_logfile,         only: iulog
  use shr_kind_mod,        only: r8 => shr_kind_r8
  use cam_abortutils,      only: endrun
  use spmd_utils,          only: masterproc
  use shr_sys_mod,         only: shr_sys_flush

  implicit none
  private

  ! Public interface
  public :: whs1998_set_ic

!==============================================================================
CONTAINS
!==============================================================================

  subroutine whs1998_set_ic(latvals, lonvals, U, V, T, PS, PHIS,           &
       Q, m_cnst, mask, verbose)
    use const_init,    only: cnst_init_default
    use constituents,  only: cnst_name

    !--JH--
    use physconst,   only: gravit, pi, rair

    !-----------------------------------------------------------------------
    !
    ! Purpose: Set Held-Suarez initial values for dynamics state variables
    !
    !-----------------------------------------------------------------------

    ! Dummy arguments
    
    ! --JH--: While comments here specify that these should be in degrees, the
    ! ic's in SE and FV3 are set in dyn_comp.F90 via a call to "analytic_ic_set_ic", an 
    ! interface defined in dynamics/tests/inic_analytic.F90 which aliases 
    ! "dyn_set_inic_cblock". 
    ! The comments in that function also state that it expects lat,lon in
    ! degrees, but still, dyn_comp for FV3 and SE pass them as radians.
    ! So here they also show up as radians, at least for FV3 and SE. My
    ! implementation will make this assumption, and so this will not work for
    ! any dycores for which the dyn_comp implementation actually *does* respect
    ! the expected input units for analytic_ic_set_ic...
    
    real(r8),           intent(in)    :: latvals(:) ! lat in degrees (ncol)
    real(r8),           intent(in)    :: lonvals(:) ! lon in degrees (ncol)

    real(r8), optional, intent(inout) :: U(:,:)     ! zonal velocity
    real(r8), optional, intent(inout) :: V(:,:)     ! meridional velocity
    real(r8), optional, intent(inout) :: T(:,:)     ! temperature
    real(r8), optional, intent(inout) :: PS(:)      ! surface pressure
    real(r8), optional, intent(out)   :: PHIS(:)    ! surface geopotential
    real(r8), optional, intent(inout) :: Q(:,:,:)   ! tracer (ncol, lev, m)
    integer,  optional, intent(in)    :: m_cnst(:)  ! tracer indices (reqd. if Q)
    logical,  optional, intent(in)    :: mask(:)    ! Only init where .true.
    logical,  optional, intent(in)    :: verbose    ! For internal use

    ! Local variables
    logical, allocatable              :: mask_use(:)
    logical                           :: verbose_use
    integer                           :: i, k, m
    integer                           :: ncol
    integer                           :: nlev
    integer                           :: ncnst
    character(len=*), parameter       :: subname = 'HS94_SET_IC'
    
    !--JH--
    ! adding topography modification from Gerber+Polvani 2008 (GB08)
    real(r8), parameter :: deg2rad = pi / 180.0_r8
    real(r8)            :: phi0     = 25.0_r8 * deg2rad   ! GP08 parameter phi_0 in radians
    real(r8)            :: phi1     = 65.0_r8 * deg2rad   ! GP08 parameter phi_1 in radians
    real(r8)            :: h0       = 3000.0_r8           ! GP08 parameter h0, in m
    real(r8)            :: mm       = 2.0_r8              ! GP08 parameter m; wavenumber of topo
    real(r8)            :: sinlat    
    real(r8)            :: sinlatsq
    real(r8)            :: surface_pressure(size(latvals))
    real(r8)            :: surface_height(size(latvals))
    real(r8), parameter :: p00 = 1.e5_r8
    real(r8), parameter :: T0 = 250._r8


    allocate(mask_use(size(latvals)))
    if (present(mask)) then
      if (size(mask_use) /= size(mask)) then
        call endrun('cnst_init_default: input, mask, is wrong size')
      end if
      mask_use = mask
    else
      mask_use = .true.
    end if

    if (present(verbose)) then
      verbose_use = verbose
    else
      verbose_use = .true.
    end if

    ncol = size(latvals, 1)
    nlev = -1
    if (present(U)) then
      nlev = size(U, 2)
      do k = 1, nlev
        where(mask_use)
          U(:,k) = 0.0_r8
        end where
      end do
      if(masterproc .and. verbose_use) then
        write(iulog,*) '          U initialized by "',subname,'"'
      end if
    end if

    if (present(V)) then
      nlev = size(V, 2)
      do k = 1, nlev
        where(mask_use)
          V(:,k) = 0.0_r8
        end where
      end do
      if(masterproc .and. verbose_use) then
        write(iulog,*) '          V initialized by "',subname,'"'
      end if
    end if
    
    ! ===================== GP08 HS init =========================
    !--JH--
    ! insert GB08 topography in the NH, SH
    ! assuming latvals, lonvals in radians, see comment above

    if (present(T)) then
      nlev = size(T, 2)
      do k = 1, nlev
        where(mask_use)
          !--JH--
          !T(:,k) = 250._r8
          T(:,k) = T0
        end where
      end do
      if(masterproc .and. verbose_use) then
        write(iulog,*) '          T initialized by "',subname,'"'
      end if
    end if

    ! --------------- init topography --------------- 
    do i = 1, ncol
        if (mask_use(i)) then
            if(abs(latvals(i)) > phi0 .and. abs(latvals(i)) < phi1) then
                sinlat = sin( (abs(latvals(i)) - phi0)/(phi1 - phi0) * pi )
                sinlatsq = sinlat * sinlat
                surface_height(i) = h0 * sinlatsq * cos(mm * lonvals(i))
                surface_pressure(i) = p00 * exp( -gravit*surface_height(i)/(rair*T0))
            else
                surface_height(i) = 0.0_r8
                surface_pressure(i) = p00
            end if
        end if
    end do 
    
   ! --------------- init PHIS --------------- 
    if (present(PHIS)) then  
      where(mask_use)
        ! surface geopotential: topo height times gravity 
        PHIS(:) = gravit*surface_height(:)
      end where
      
      if(masterproc .and. verbose_use) then
        write(iulog,*) '          PHIS initialized by "',subname,'"'
      end if
    end if
    
   ! --------------- init PS --------------- 
    if (present(PS)) then
      where(mask_use)
        PS(:) = surface_pressure(:)
      end where
      
      if(masterproc .and. verbose_use) then
        write(iulog,*) '          PS initialized by "',subname,'"'
      end if
    end if
    
    ! ===================== End GP08 HS init =========================

    if (present(Q)) then
      nlev = size(Q, 2)
      ncnst = size(m_cnst, 1)
      do m = 1, ncnst
        if (m_cnst(m) == 1) then
          ! No water vapor in Held-Suarez
          do k = 1, nlev
            where(mask_use)
              Q(:,k,m_cnst(m)) = 0.0_r8
            end where
          end do
          if(masterproc .and. verbose_use) then
            write(iulog,*) '          ', trim(cnst_name(m_cnst(m))), ' initialized by "',subname,'"'
          end if
        else
          call cnst_init_default(m_cnst(m), latvals, lonvals, Q(:,:,m_cnst(m)),&
               mask=mask_use, verbose=verbose_use, notfound=.false.)
        end if
      end do
    end if
    

    deallocate(mask_use)

  end subroutine whs1998_set_ic

end module ic_whs1998
