module whs1998
    
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

  use ccpp_kinds, only: kind_phys

  !--JH--
  use physconst,   only: gravit, rair
  use cam_logfile, only: iulog
  use ref_pres,           only: ptop_ref  !--JH--

  implicit none
  private
  save

  public :: whs1998_mod_init
  public :: whs1998_mod_run

  !!
  !! Forcing parameters
  !!
  real(kind_phys), parameter :: efoldf  =  1._kind_phys  ! efolding time for wind dissipation
  real(kind_phys), parameter :: efolda  = 40._kind_phys  ! efolding time for T dissipation
  real(kind_phys), parameter :: efolds  =  4._kind_phys  ! efolding time for T dissipation
  real(kind_phys), parameter :: sigmab  =  0.7_kind_phys ! threshold sigma level
  real(kind_phys), parameter :: t00     = 200._kind_phys ! minimum reference temperature
  real(kind_phys), parameter :: kf      = 1._kind_phys/(86400._kind_phys*efoldf) ! 1./efolding_time for wind dissipation

  real(kind_phys), parameter :: onemsig = 1._kind_phys - sigmab ! 1. - sigma_reference

  real(kind_phys), parameter :: ka      = 1._kind_phys/(86400._kind_phys * efolda) ! 1./efolding_time for temperature diss.
  real(kind_phys), parameter :: ks      = 1._kind_phys/(86400._kind_phys * efolds)

  !!
  !! Model constants, reset in init call
  !!
  real(kind_phys)              :: cappa = 2.0_kind_phys / 7.0_kind_phys  ! R/Cp
  real(kind_phys)              :: cpair = 1004.0_kind_phys        ! specific heat of dry air (J/K/kg)
  real(kind_phys)              :: psurf_ref = 0.0_kind_phys       ! Surface pressure
  ! pref_mid_norm are layer midpoints normalized by surface pressure ('eta' coordinate)
  real(kind_phys), allocatable :: pref_mid_norm(:)



!======================================================================= 
contains
!======================================================================= 

  subroutine whs1998_mod_init(pver_in, cappa_in, cpair_in, &
                                   psurf_ref_in, pref_mid_norm_in, errmsg, errflg)
    !! Dummy arguments
    integer,           intent(in) :: pver_in
    real(kind_phys),   intent(in) :: cappa_in
    real(kind_phys),   intent(in) :: cpair_in
    real(kind_phys),   intent(in) :: psurf_ref_in
    real(kind_phys),   intent(in) :: pref_mid_norm_in(:)
    character(len=512),intent(out):: errmsg
    integer,           intent(out):: errflg

    integer               :: pver                     ! Num vertical levels

    errmsg = ' '
    errflg = 0

    pver = pver_in
    allocate(pref_mid_norm(pver))
    cappa         = cappa_in
    cpair         = cpair_in
    psurf_ref     = psurf_ref_in
    pref_mid_norm = pref_mid_norm_in   ! Layer midpoints normalized by surface pressure 
                                       ! ('eta' coordinate)

  end subroutine whs1998_mod_init

  subroutine whs1998_mod_run(pver, ncol, clat, pmid, &
       u, v, t, du, dv, s, errmsg, errflg)

    !
    ! Input arguments
    !
    integer,  intent(in)  :: pver      ! Num vertical levels
    integer,  intent(in)  :: ncol      ! Num active columns
    real(kind_phys), intent(in)  :: clat(:)   ! latitudes(radians) for columns
    real(kind_phys), intent(in)  :: pmid(:,:) ! mid-point pressure
    real(kind_phys), intent(in)  :: u(:,:)    ! Zonal wind (m/s)
    real(kind_phys), intent(in)  :: v(:,:)    ! Meridional wind (m/s)
    real(kind_phys), intent(in)  :: t(:,:)    ! Temperature (K)
    !
    ! Output arguments
    !
    real(kind_phys),   intent(out) :: du(:,:)   ! Zonal wind tend
    real(kind_phys),   intent(out) :: dv(:,:)   ! Meridional wind tend
    real(kind_phys),   intent(out) :: s(:,:)    ! Heating rate
    character(len=512),intent(out):: errmsg
    integer,           intent(out):: errflg
    !
    !---------------------------Local workspace-----------------------------
    !
    integer  :: i, k          ! Longitude, level indices

    real(kind_phys) :: kv            ! 1./efolding_time (normalized) for wind
    real(kind_phys) :: kt            ! 1./efolding_time for temperature diss.
    real(kind_phys) :: trefa         ! "radiative equilibrium" T
    real(kind_phys) :: trefc         ! used in calc of "radiative equilibrium" T
    real(kind_phys) :: cossq(ncol)   ! coslat**2
    real(kind_phys) :: cossqsq(ncol) ! coslat**4
    real(kind_phys) :: sinsq(ncol)   ! sinlat**2
    real(kind_phys) :: coslat(ncol)  ! cosine(latitude)

    ! --JH--
    
    ! add forcing parameters for WHS modification
    real(kind_phys)   :: efoldaa          ! efolding time for T dissipation
    real(kind_phys)   :: kaa              ! 1./efolding_time for temperature diss.
    real(kind_phys)   :: pi               ! pi
    real(kind_phys)   :: p0strat          ! threshold pressure in Pa
    real(kind_phys)   :: phi0             ! threshold latitude
    real(kind_phys)   :: dphi0            ! del-latitude
    real(kind_phys)   :: a0               ! coefficent
    real(kind_phys)   :: aeq              ! 100 mb in sigma coordinate
    real(kind_phys)   :: apole            ! 2 mb in sigma coordinate
    real(kind_phys)   :: lapsew           ! lapse rate
    real(kind_phys)   :: constw           ! constant
    real(kind_phys)   :: lapsec           ! lapse rate
    real(kind_phys)   :: constc           ! constant
    real(kind_phys)   :: acoslat          ! abs(acos(coslat))
    real(kind_phys)   :: pih                 ! 0.5*pi

    
    !
    !-----------------------------------------------------------------------
    !
    
    ! --JH-- 
    ! initialize parameter values for WHS mod 
    efoldaa = 40._kind_phys
    kaa = 1._kind_phys/(86400._kind_phys*efoldaa)
    pi = 4._kind_phys*atan(1._kind_phys)
    pih = 2._kind_phys*atan(1._kind_phys)
    phi0   = 60._kind_phys*pi/180._kind_phys
    dphi0  = 15._kind_phys*pi/180._kind_phys
    a0     = 2.65_kind_phys/dphi0
    aeq    = 10000._kind_phys / psurf_ref 
    !apole  = 200._kind_phys / psurf_ref            ! original value from Williamson+ modificaiton
    apole  = (ptop_ref * 0.9_kind_phys) / psurf_ref ! place apole slightly above model top
    lapsew = -3.345e-03_kind_phys
    constw = rair*lapsew/gravit
    lapsec =  2.00e-03_kind_phys
    constc = rair*lapsec/gravit


    errmsg = ' '
    errflg = 0

    do i = 1, ncol
      coslat (i) = cos(clat(i))
      sinsq  (i) = sin(clat(i))*sin(clat(i))
      cossq  (i) = coslat(i)*coslat(i)
      cossqsq(i) = cossq (i)*cossq (i)
    end do

    ! --JH--
    ! Compute idealized radiative heating rates (as dry static energy)

    ! based on source modofications at
    ! /glade/u/home/cjablono/cam3_5_41/models/atm/cam/src/physics/cam/tphysidl.F90 (henceforth TPF)
    ! usages of ptend there are replaced with the local var s
    ! usages of tmp (for implicit time stepping?) were removed

    do k = 1, pver
        do i = 1, ncol
          if (pref_mid_norm(k) > sigmab) then
            ! --- apply when pressure greater than sigmab, as in HS+1994
            kt = ka + (ks - ka)*cossqsq(i)*(pref_mid_norm(k) - sigmab)/onemsig
          else
            ! --- apply when pressure below sigmab, as in HS+1994
            kt =  ka
          endif
          
          acoslat = abs(acos(coslat(i)))
          trefc   = 315._kind_phys - (60._kind_phys * sinsq(i))

          ! --- apply standard HS temperature forcing below aeq Pa
          trefa=(trefc - 10._kind_phys*cossq(i)*log((pmid(i,k)/psurf_ref)))*(pmid(i,k)/psurf_ref)**cappa
          trefa   = max(t00,trefa)

          ! --- apply first term of WHS+1998 Eq1 Apdx A for p < aeq
          if (pmid(i,k) < (aeq * psurf_ref)) then    ! scale aeq to Pa
             trefa = t00*((pmid(i,k)/10000._kind_phys))**constc
          endif
          
          ! --- apply second term of WHS+1998 Eq1 Apdx A for p < p0strat
          p0strat = aeq - (aeq - apole)*0.5_kind_phys*(1._kind_phys + tanh(a0*(acoslat - phi0)))
          p0strat = p0strat * psurf_ref   ! scale to Pa
          if (pmid(i,k) < p0strat) then
             trefa = trefa + t00*( ((pmid(i,k)/p0strat))**constw - 1._kind_phys )
          endif
          
          ! update heating 
          s(i,k)  = (trefa - t(i,k))*kt*cpair
        end do
    end do
    
    !
    ! Add diffusion near the surface for the wind fields
    !
    du(:,:) = 0._kind_phys
    dv(:,:) = 0._kind_phys

    ! this included in TPF, not sure why
    !kf = 1._r8/(86400._r8*efoldf)

    do k = 1, pver
      if (pref_mid_norm(k) > sigmab) then
        kv  = kf*(pref_mid_norm(k) - sigmab)/onemsig
        do i = 1, ncol
          du(i,k) = -kv*u(i,k)
          dv(i,k) = -kv*v(i,k)
        end do
      end if
    end do

  end subroutine whs1998_mod_run

end module whs1998
