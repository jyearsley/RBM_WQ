module Block_Energy
!
!   Energy budget variables
!
!   Incoming short wave radiation, kcal/m**2/sec
!
    real(sp), dimension(:), allocatable::q_ns
!
!   Incoming atmospheric radiation, kcal/m**2/sec
!
    real(sp), dimension(:), allocatable::q_na
!
!   Air temperature at surface, deg. C
!
    real(sp), dimension(:), allocatable::dbt
!  
!   Wind speed, m/sec
!
    real(sp), dimension(:), allocatable::wind
!
!   Vapor pressure of air at surface, mb
!
    real(sp), dimension(:), allocatable::ea
!
!   Air pressure at surface, mb
!
    real(sp), dimension(:), allocatable::press 

!
    real(sp), dimension (:), allocatable::mu,alphamu,beta,gmma,smooth_param

!   Some important constants
!
      real(sp)             :: lvp,rb,rho
      real(sp),parameter   :: evap_coeff=1.5e-9,pf=0.640,pi=3.14159,rfac=304.8       
!
end module Block_Energy  
