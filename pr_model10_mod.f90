!--------------------------------------------------------------------------------
module model_mod
!--------------------------------------------------------------------------------
use ODE_mod
use carbonate_chemistry_mod
use spline_mod
use ext_driver_mod
use resize_mod

implicit none
integer,public  , parameter :: nstate = 5, naux=9  ! size of state vector + aux variables  
double precision, parameter :: pi=3.141592654, eps=5d-9


double precision,allocatable,dimension(:) :: mut, muc,     &
                                             Pt,Pc,        & !
                                             Rt,Rc,        &
                                             kla1t, kla1c, &     ! rate of O2 production
                                             kla2t, kla2c, &
                                             kM1t,  kM1c,  &
                                             kM2t,  kM2c,  &
                                             kM3t,  kM3c,  &
                                             fP1t,  fP1c,  &
                                             fP3t,  fP3c,  &
                                             tat,   tac,   &
                                             tauRt, tauRc, &
                                             tauPt, tauPc, &
                                             flert, flerc, &
                                             PQdt,  PQdc,  &  ! photosynthetic quotient day
                                             PQnt,  PQnc,  &  !                         night
                                             NCdt,  NCdc      ! Nitrogen-Carbon ratio


! carbonate chemistry parameters
double precision :: S =33.,                     &  ! Salinity
                    TK=273.15+22.,              &  ! Temp in kelvin
                    TA =2300.,                  &  ! alkinity in uM
                    KM=500.,                    &  ! Michealis-Menten term for inorg C uptake
                    xO2 (2) = (/0.2094,0.0000/),&  ! partial pressure of O2
                    xCO2(2) = (/400d-6,1.00d0/),&  ! partial pressure of CO2
                    K1f, K2f 
contains

function indices( bool_array, n ) result (inds)
logical,intent(in) :: bool_array(:)
integer,intent(in) :: n
integer            :: inds(n)

integer :: i,ii

inds = 0

ii=1
do i = 1, size(bool_array)
   if( bool_array(i) )then
      inds(ii)=i
      ii = ii+1
      if( ii>n ) return
   endif
enddo
return 
end function indices

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function sim_PR_de( theta, theta_typ, time_steps ) result( y )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double precision, intent(in) :: theta(:)

double precision, intent(in) :: time_steps(:)
character(len=*),    intent(in) :: theta_typ(:)
double precision             :: y(nstate+naux,size(time_steps))

double precision :: y0(nstate), y1(nstate)

integer          :: n, nstep
double precision :: dt, t0, t1, tt  ! time step info

double precision :: X_0       ! initial biomass
double precision :: O2_0      ! initial O2 concentration
double precision :: DIC_0     ! initial O2 concentration

double precision :: A(nstate,nstate)

double precision :: O2,DIC ,H, pH, CO3, HCO3, CO2,& ! PR=rate of photosynthesis or respiration
                    PM,P,P_CO2,P_HCO3,P_CO3,R,PQ,KM1,KM2,fccm, f_photo, X, mu, kLa(2)
integer :: i,nc, m,msub
integer,allocatable :: ind_tmp(:)
character(len=10) :: str_tmp,comp
logical :: mask(size(theta))

nstep = size(time_steps)

y = 0.d0


! unpack the parameters
! initial parameters
!X_0     = theta( 1)
O2_0    = theta( 1)   ! oxygen
DIC_0   = theta( 2)
TA      = theta( 3)

mask = theta_typ=='O20'  ; ind_tmp = indices(mask,count(mask) )  ; O2_0 = theta(ind_tmp(1))

! load the spline points
mask = theta_typ=='R'    ; Rc    = theta( indices(mask,count(mask) ) )
mask = theta_typ=='P'    ; Pc    = theta( indices(mask,count(mask) ) )
mask = theta_typ=='kla1' ; kla1c = theta( indices(mask,count(mask) ) )
mask = theta_typ=='kla2' ; kla2c = theta( indices(mask,count(mask) ) )
mask = theta_typ=='km1'  ; km1c = theta( indices(mask,count(mask) ) )
mask = theta_typ=='km2'  ; km2c = theta( indices(mask,count(mask) ) )
mask = theta_typ=='km3'  ; km3c = theta( indices(mask,count(mask) ) )
mask = theta_typ=='ta'   ; tac  = theta( indices(mask,count(mask) ) )
mask = theta_typ=='tauP' ; tauPc = theta( indices(mask,count(mask) ) )
mask = theta_typ=='tauR' ; tauRc = theta( indices(mask,count(mask) ) )
mask = theta_typ=='PQd'  ; PQdc  = theta( indices(mask,count(mask) ) )
mask = theta_typ=='PQn'  ; PQnc  = theta( indices(mask,count(mask) ) )
mask = theta_typ=='NCd'  ; NCdc  = theta( indices(mask,count(mask) ) )

! transform the units
Pc = Pc*24.
Rc = Rc*24.
kla1c   = log(2.)*60.*24./kla1c
kla2c   = log(2.)*60.*24./kla2c
tac     = tac*24.
tauPc   = log(2.)*60.*24./tauPc
tauRc   = log(2.)*60.*24./tauRc

! load the intial conditaions
y0 = 0.d0
y0(1) = O2_0
y0(2) = DIC_0
y0(3) = TA
!y0(4) = 0.d0
y0(4) = 1.d0

call CC_solve_DIC_Talk( DIC_0*1d-6, TA*1.d-6, TK, S, CO2, H,  pH=pH, HCO3=HCO3,CO3=CO3, &
                const=10 )
CO2 =   CO2*1d6
HCO3 = HCO3*1d6
CO3 =   CO3*1d6
y =0.d0

t0 = time_steps(1)
call calc_XPR(t0,y0, mu,PM,P,P_CO2,P_HCO3,P_CO3,R,PQ, pH,CO2,HCO3,CO3)
y(        :nstate,1) = y0
y(nstate+1:      ,1) =(/ pH,CO2,HCO3,CO3, PM,P_CO2,P_HCO3,P_CO3,R /)


do n = 2, nstep
   t0 = time_steps(n-1)+eps
   t1 = time_steps(n  )-eps
   dt = t1-t0

   ! step the state vector
   call RK4 ( y0, t0, dt, dy_PR,  y1 )
   !y1 =  magnus4_nl_real( y0, t0, dt, dA_PR )
   y(        :nstate,n) = y1

   ! calculate the auxilary varaiable and store them
   call calc_XPR(t1,y1, mu,PM,P,P_CO2,P_HCO3,P_CO3,R,PQ, pH,CO2,HCO3,CO3)
   y(nstate+1:      ,n) =(/ pH,CO2,HCO3,CO3, PM,P_CO2,P_HCO3,P_CO3,R /)

   y0 = y1     ! Don't delete!
enddo


return
end function sim_PR_de

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine calc_XPR( t,y, mu, PM,P,P_CO2,P_HCO3,P_CO3,R,PQ, pH,CO2, HCO3, CO3 )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double precision,intent(in) :: t, y(:) 
double precision,intent(out):: mu,PM,P,P_HCO3,P_CO2, P_CO3, R,PQ , pH,CO2, HCO3, CO3

double precision :: t_dawn, t_dusk, X, O2,DIC,TA, I ,S,TK, KM,fler,t_resp,t_photo,&
                    fphoto,  fP1,fP3 , H, KM1,KM2,KM3


t_dawn = 24*(t-floor(t))-8.
t_dusk = 24*(t-floor(t))-20.

if( t_dusk <  0. ) t_dusk=12.
if( t_dawn <  0. ) t_dawn=12.
if( t_dusk >12.0 ) t_dusk=12.
if( t_dawn >12.0 ) t_dawn=12.


O2 = y(1)
DIC= y(2)
TA = y(3)

I  = ext_light(t)              ! irradiance at top of column (uEin)
TK = ext_temp(t)+273.15d0
S  = ext_sal( t )

K1f = ext_K1f( t )
K2f = ext_K2f( t )

call CC_solve_DIC_Talk( DIC*1d-6, TA*1.d-6, TK, S, CO2, H,  pH=pH,HCO3=HCO3,CO3=CO3, &
               K1_f=K1f, K2_f=K2f, const=10 )
CO2  =  CO2*1d6
HCO3 = HCO3*1d6
CO3  =  CO3*1d6

fler   = 0.d0 ; if(size(flert)>0) fler   = spline_hc( t,flert, flerc )
t_resp = 1d6  ; if(size(tauRt)>0) t_resp = spline_hc( t,tauRt, tauRc )
t_photo= 1d6  ; if(size(tauPt)>0) t_photo= spline_hc( t,tauPt, tauPc )

mu =  0.0
!if( nmu1-nmu0 >-1 )then
!   if( mu_t0 <=t .and. t<=mu_t0+mu_dt ) mu = spline_hc( t, mut, muc )
!endif

if( t_dawn < 10.0 )then
   R = -spline_hc( t,Rt,Rc)!*(1.0+fler*(1.-exp(-t_dawn*t_resp)))
else
   R = -spline_hc( t,Rt,Rc)!*(1.0+fler*exp(-t_dusk*t_resp))
endif

P  = 0.d0
PQ = 1.d0
P_HCO3 = 0.d0
P_CO2  = 0.d0
P_CO3  = 0.d0
if ( I>0.1) then
   KM1=0.
   KM2=0.
   KM3=0.
   PM   = spline_hc( t, Pt,Pc)
   fP1  = 0.
   fP3  = 0.
   if(size(fP1t)>0) fP1  = spline_hc( t, fP1t, fP1c )
   if(size(fP3t)>0) fP3  = spline_hc( t, fP3t, fP3c )
   if(size(Km1t)>0) Km1  = spline_hc( t, kM1t ,kM1c )
   if(size(Km2t)>0) Km2  = spline_hc( t, kM2t ,kM2c )
   if(size(Km3t)>0) Km3  = spline_hc( t, kM3t ,kM3c )

   fphoto= (1.-exp(-t_dawn*t_photo))
   P_CO2  = fP1*PM* CO2/( CO2+KM1) ! *fphoto
   P_HCO3 =     PM*HCO3/(HCO3+KM2) !*fphoto 
   P_CO3  = fP3*PM* CO3/( CO2+KM3) 
   P=P_HCO3 + P_CO2 + P_CO3

   if(size(PQdt)>0) PQ  = spline_hc( t, PQdt ,PQdc  )
else
   if(size(PQnt)>0) PQ  = spline_hc( t, PQnt ,PQnc  )
endif

return
end subroutine calc_XPR

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function dy_PR( t, y ) result( dy )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double precision,intent(in) :: t,  y(:)
double precision            ::    dy(size(y))

double precision :: DIC,TA, X,O2, PR ,H, pH, HCO3, CO2, CO3,& ! PR=rate of photosynthesis or respiration
                     I,spg(2), O2_H(2), CO2_H(2), kLa(2),  PM,P,P_HCO3,P_CO2,P_CO3,R,PQ, dTA, NC

double precision :: t_dawn, t_dusk , mu=0.d0

t_dawn = t-floor((t)/24.)*24-8.
t_dusk = t-floor((t)/24.)*24-20.

O2 = y(1)
DIC= y(2)
TA = y(3)

I       = ext_light(t)              ! irradiance at top of column (uEin)
spg     = ext_gas  (t)              ! sparging rate


! calculate the Henry's law level of dissolved gases for the input streams
S  = ext_sal( t )
TK = ext_temp(t)+273.15d0
O2_H  = xO2  * K0_O2 ( TK, S ) * rho_sw( TK, S ) * 1000.
CO2_H = xCO2 * K0_CO2( TK, S ) * rho_sw( TK, S ) * 1000.

kLa = 0.d0
kLa(1) = spline_hc( t, kla1t,kla1c )
kLa(1) = lint_1D  ( t, kla1t,kla1c ,3 )
if( size(kla2t)>0) kLa(2) = lint_1D( t, kla2t,kla2c,3 )

! calculate the growth, P and R rates
call calc_XPR( t, y, mu,PM,P,P_CO2,P_HCO3,P_CO3,R,PQ, pH, CO2,HCO3,CO3) 

dTA = 0.
if(size(tat)>0) dTA = lint_1D( t, tat, tac,1 )
NC  = 0.
if(size(NCdt)>0) NC = spline_hc( t, NCdt, NCdc )
if( I < 0.1 ) NC= 0

dy =  0.
dy(1) =     (P+R)/PQ         +  dot_product( kLa       * spg, ( O2_H- O2) )   
dy(2) = -   (P+R) +  0*dTA   +  dot_product( kLa*0.893 * spg, (CO2_H-CO2) ) 
dy(3) =  NC*(P+R) +2.0*dTA
dy(4) =    +(P+R)       ! C fixed in Cell
return
end function dy_PR

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function dA_PR( t, y ) result( dA )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double precision,intent(in) :: t,  y(:)
double precision            ::    dA(size(y),size(y))


double precision :: DIC,TA, X,O2, PR ,H, pH, HCO3, CO2, CO3,& ! PR=rate of photosynthesis or respiration
                     I,spg(2), O2_H(2), CO2_H(2), kLa(2),  PM,P,P_HCO3,P_CO2,P_CO3,R,PQ, dTA, NC

double precision :: t_dawn, t_dusk , mu=0.d0

integer :: ny

ny = size(y)


O2 = y(1)
DIC= y(2)
TA = y(3)

I       = ext_light(t)              ! irradiance at top of column (uEin)
spg     = ext_gas  (t)              ! sparging rate


! calculate the Henry's law level of dissolved gases for the input streams
S  = ext_sal( t )
TK = ext_temp(t)+273.15d0
O2_H  = xO2  * K0_O2 ( TK, S ) * rho_sw( TK, S ) * 1000.
CO2_H = xCO2 * K0_CO2( TK, S ) * rho_sw( TK, S ) * 1000.

kla=0.d0
!kLa(1) = spline_hc( t, kla1t,kla1c )
kLa(1) = lint_1D( t, kla1t,kla1c ,2 )
if( size(kla2t)>0 ) kLa(2) = lint_1D( t, kla2t,kla2c,1 )


! calculate the growth, P and R rates
call calc_XPR( t, y, mu,PM,P,P_CO2,P_HCO3,P_CO3,R,PQ, pH, CO2,HCO3,CO3)

dTA = 0.
if( size(tat)>0 ) dTA = spline_hc( t, tat, tac )
NC  = 0.
if( size(NCdt)>0 ) NC = spline_hc( t, NCdt, NCdc )


dA = 0
dA(1,1 ) =           - dot_product( kLa       , spg              )
dA(1,2 ) =  (P+R)/DIC
dA(1,ny) =           + dot_product( kLa       * spg, ( O2_H    ) )


dA(2,2 ) = -(P+R)/DIC - dot_product( kLa*0.893 , spg)*CO2/DIC
dA(2,ny) =            + dot_product( kLa*0.893 * spg, (CO2_H    ) )

if(I>0.1) dA(3,ny) =  NC*(P+R)


return
end function dA_PR

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine initialise()
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

call set_size( mut,   0 ) ; call set_size( muc,   0 )
call set_size( Pt,    0 ) ; call set_size( Pc,    0 )
call set_size( Rt,    0 ) ; call set_size( Rc,    0 )
call set_size( kla1t, 0 ) ; call set_size( kla1c, 0 )
call set_size( kla2t, 0 ) ; call set_size( kla2c, 0 )
call set_size( km1t,  0 ) ; call set_size( km1c,  0 )
call set_size( km2t,  0 ) ; call set_size( km2c,  0 )
call set_size( km3t,  0 ) ; call set_size( km3c,  0 )
call set_size( fP1t , 0 ) ; call set_size( fP1c, 0 )
call set_size( fP3t , 0 ) ; call set_size( fP3c, 0 )
call set_size( tat,   0 ) ; call set_size( tac,   0 )
call set_size( tauPt, 0 ) ; call set_size( tauPc, 0 )
call set_size( tauRt, 0 ) ; call set_size( tauRc, 0 )
call set_size( flert, 0 ) ; call set_size( flerc, 0 )
call set_size( PQdt , 0 ) ; call set_size( PQdc , 0 )
call set_size( PQnt , 0 ) ; call set_size( PQnc , 0 )

return
end subroutine initialise

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine set_spline_var( varname, Cp )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! sets the spline control points and saves the variables for the spline setup
character*(*)   , intent(in)   :: varname
double precision,intent(in) :: Cp(:)      ! control points of spline

integer :: nc

select case( trim(varname))
case("mu")   ; mut    = Cp
case("P")    ; Pt     = Cp 
case("R")    ; Rt     = Cp   
case("kla1") ; kla1t  = Cp    
case("kla2") ; kla2t  = Cp   
case("km1" ) ; km1t   = Cp   
case("km2" ) ; km2t   = Cp   
case("km3" ) ; km3t   = Cp
case("fP1" ) ; fP1t   = Cp   
case("fP3" ) ; fP3t   = Cp   
case("dta" ) ; tat    = Cp
case("tauR") ; tauRt  = Cp
case("tauP") ; tauPt  = Cp   
case("fler") ; flert  = Cp   
case("PQd")  ; PQdt   = Cp   
case("PQn")  ; PQnt   = Cp
case("NCd")  ; NCdt   = Cp   
end select

return
end subroutine set_spline_var

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine set_size( arr, n )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double precision, allocatable,intent(inout) :: arr(:)
integer,intent(in) :: n
if( allocated(arr) )then
   if( size(arr) .ne. n )then
      deallocate(arr)
      allocate(arr(n))
   endif
else
   allocate(arr(n))
endif
end subroutine set_size

end module model_mod

