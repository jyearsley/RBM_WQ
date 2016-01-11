MODULE KINETICS
use 
IMPLICIT NONE
CONTAINS
!**********************************************************************
!                                                                     * 
! Subroutine for Water Quality Constituents                           *
!                                                                     *
!**********************************************************************
SUBROUTINE WQ(No_Const,WQ)
!**********************************************************************
! Select cases from WQ                                                *
!**********************************************************************
do nq=1,No_Const
  select case (CONST_Table(nq))
!**********************************************************************
! Water Temperature using energy budget method                        *
!**********************************************************************
    case (1)
    use Block_Energy
    implicit none
    integer::i,ncell,nd
    real::A,B,e0,q_surf,q_conv,q_evap,q_ws,td,T_surf
    real, dimension(2):: q_fit, T_fit
!
    td=nd
    T_fit(1)=T_surf-1.0
    T_fit(2)=T_surf+1.0
    do i=1,2
      e0=2.1718E8*EXP(-4157.0/(T_fit(i)+239.09))
      rb=pf*(dbt(ncell)-T_fit(i))
      lvp=597.0-0.57*T_fit(i)
      q_evap=1000.*lvp*(evap_a+evap_b*wind(ncell)+evap_c*(wind(cell)**2)
      if(q_evap.lt.0.0) q_evap=0.0
      q_conv=rb*q_evap
      q_evap=q_evap*(e0-ea(ncell))
      q_ws=6.693E-2+1.471E-3*T_fit(i)
      q_fit(i)=(1.-r_short)*q_ns(ncell)+(1.-r_atmo)*q_na(ncell)-q_ws-q_evap+q_conv
    end do
!
!     q=AT+B
!
!     Linear fit over the range of 2.0 deg C.
!     These results can be used to estimate the "equilibrium" 
!     temperature and linear rate constant.
!
    A=(q_fit(1)-q_fit(2))/(T_fit(1)-T_fit(2))
    q_surf=0.5*(q_fit(1)+q_fit(2))
    B=(q_surf/A)-(T_fit(1)+T_fit(2))/2.
!
!     ******************************************************
!               Return to Subroutine SYSTMM
!     ******************************************************!**********************************************************************
! Dissolved Oxygen                                                    *
!**********************************************************************
    case (2)
!**********************************************************************

!**********************************************************************
! ELEVATION CORRECTION FOR GASES                                      *
    real             :: z_fctr
    real, intent(in) :: z      ! Elevation of surface above MSL, meters
!
    z_fctr = (1.0-0.0001148*z)
!
!
! DISSOLVED OXYGEN SATURATION
!
!     real            :: do_sat  
    real,intent(in) :: t 
!    
!
! D I S S O L V E D   O X Y G E N
!
        k2 = (reaer_a+reaer_b*wind+reaer_c*wind*wind)
  do_sat = z_fctr*EXP(7.7117-1.31403*(LOG(t+45.93)))

!
! Rate of change
!
  CONST_rate(2)     = (-k_bod*BOD+k2*(do_sat- DO0)/z)*delta_t
!
!**********************************************************************
! B I O C H E M I C A L   O 2   D E M A N D                           *
!**********************************************************************
! 
!**********************************************************************
    case (3)                          
!**********************************************************************
!
    real            ::cbod
    real,intent(in) :: k_bod,bod,t
! 
    CONST_rate(3) = -k_bod*(Q10_bod**(t-20.))*CONST(3,nr,nncell)*delta_t
!
!
!
!************************************************************************
!**                   D E C A Y   C O N S T A N T S                     *
!************************************************************************
!
!**********************************************************************
!  C O L I F O R M                                                    *
!**********************************************************************
!
!**********************************************************************
  CASE (4)
!**********************************************************************
    real,intent(in)      :: k_coli,Q10_coli,t,coli0
    real                 :: coli
!
    coli = -k_coli*(Q10_coli**(t-20.0))*coli0
!
!**********************************************************************
!  S U S P E N D E D   S O L I D S                                    *
!**********************************************************************
  CASE (5)
!**********************************************************************
!
    real,intent(in) :: W_Settl,Delta_T,Depth,tss_0
    real            :: tss
!
    tss = (1.-w_settl*delta_t/z)*tss_0!

  END SELECT
!
END MODULE KINETICS
