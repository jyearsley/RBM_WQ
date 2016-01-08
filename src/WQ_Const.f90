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
      q_evap=1000.*lvp*evap_coeff*wind(ncell)
      if(q_evap.lt.0.0) q_evap=0.0
      q_conv=rb*q_evap
      q_evap=q_evap*(e0-ea(ncell))
      q_ws=6.693E-2+1.471E-3*T_fit(i)
      q_fit(i)=q_ns(ncell)+q_na(ncell)-q_ws-q_evap+q_conv
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
!               Return to Subroutine RIVMOD
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
!     real            :: do_x  
    real,intent(in) :: T 
!
    do_x = EXP(7.7117-1.31403*(LOG(T+45.93)))
!    
!
! D I S S O L V E D   O X Y G E N
!
!        O2EX = (0.5+0.05*WIND*WIND)/86400.0
!          DOSS(KT,I) = 0.0
!          DO K=KT,KB(I)
!            APR       = (O2AG*AGR(K,I)-O2AR*ARR(K,I))*ALGAE(K,I)
!            DOSS(K,I) = APR-O2NH4*NH4D(K,I)-O2OM*(LPOMD(K,I)
!     .                  +SEDD(K,I))-SODD(K,I)*A3(K,I)-O2OM*DOMD(K,I)
!     .                  -CBODD(K,I)*CBOD(K,I)*RBOD
!          END DO
  DO_SAT = DO_SATURATION(T)
!
! Rate of change
!
  CONST_rate(2)     = (-K_bod*BOD+K2*(DO_SAT- DO0)/Z)*delta_t
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
    real,intent(in) :: K_bod,bod,T
! 
    CONST_rate(3) = -K_bod*(1.047**(T-20.))*CONST(3,nr,nncell)*delta_t
!
!
!                   R A T E   M U L T I P L I E R S                  **
!
!      ENTRY RATE_MULTIPLIERS
!        DO I=IU,ID
!            LAM1       = FR(T1(K,I),NH4T1,NH4T2,NH4K1,NH4K2)
!            NH4RM(K,I) = LAM1/(1.0+LAM1-NH4K1)
!            LAM1       = FR(T1(K,I),NO3T1,NO3T2,NO3K1,NO3K2)
!            NO3RM(K,I) = LAM1/(1.0+LAM1-NO3K1)
!            LAM1       = FR(T1(K,I),OMT1,OMT2,OMK1,OMK2)
!            OMRM(K,I)  = LAM1/(1.0+LAM1-OMK1)
!            LAM1       = FR(T1(K,I),AT1,AT2,AK1,AK2)
!            LAM2       = FF(T1(K,I),AT3,AT4,AK3,AK4)
!            ARMR(K,I)  = LAM1/(1.0+LAM1-AK1)
!            ARMF(K,I)  = LAM2/(1.0+LAM2-AK4)
!          END DO
!        END DO
!      RETURN

!************************************************************************
!**                   D E C A Y   C O N S T A N T S                     *
!************************************************************************
!
!      ENTRY DECAY_CONSTANTS
!        DO I=IU,ID
!          DO K=KT,KB(I)
!            A1(K,I)    = (1.0+SIGN(1.0,DO(K,I)-O2LIM))*0.5
!            A2(K,I)    = (1.0+SIGN(1.0,O2LIM-DO(K,I)))*0.5
!            A3(K,I)    = (1.0+SIGN(1.0,DO(K,I)-1.E-10))*0.5
!            DOMD(K,I)  = OMRM(K,I)*(LDOMDK*LDOM(K,I)+RDOMDK
!     .                   *RDOM(K,I))*A3(K,I)
!            NH4D(K,I)  = NH4DK*NH4RM(K,I)*NH4(K,I)*A1(K,I)
!        CTION SUSP_SOLIDS
!**********************************************************************
!  C O L I F O R M                                                    *
!**********************************************************************
!
!**********************************************************************
  CASE (4)
!**********************************************************************
    real,intent(in)      :: K_coli,Q10_coli,T,coli0
    real                 :: coli
!
    Coli = -K_coli*(Q10_coli**(T-20.0))*coli0
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
    tss = (1-W_Settl*Delta_T/Depth)*tss_0!

  END SELECT
!
END MODULE KINETICS


!************************************************************************
!**                        L A B I L E   D O M                         **
!************************************************************************

!      ENTRY LABILE_DOM
!        D    NO3D(K,I)  = NO3DK*NO3RM(K,I)*NO3(K,I)*A2(K,I)
!            CBODD(K,I) = KBOD*TBOD**(T1(K,I)-20.0)*A3(K,I)
!            LPOMD(K,I) = LPOMDK*OMRM(K,I)*LPOM(K,I)*A3(K,I)
!            SEDD(K,I)  = SDK*OMRM(K,I)*SED(K,I)*A3(K,I)
!          END DO
!        END DO
!        DO I=IU,ID
!          FPSS(KT,I)   = PARTP*SS(KT,I)/(PARTP*(SS(KT,I)+FE(KT,I))+1.0)
!          FPFE(KT,I)   = PARTP*FE(KT,I)/(PARTP*(SS(KT,I)+FE(KT,I))+1.0)
!          SETOUT(KT,I) = (SSS*FPSS(KT,I)+FES*FPFE(KT,I))/HKT2(I)
!     .                   *A1(KT,I)
!          DO K=KT+1,KB(I)
!            FPSS(K,I)   = PARTP*SS(K,I)/(PARTP*(SS(K,I)+FE(K,I))+1.0)
!            FPFE(K,I)   = PARTP*FE(K,I)/(PARTP*(SS(K,I)+FE(K,I))+1.0)
!            SETIN(K,I)  = SETOUT(K-1,I)
!            SETOUT(K,I) = (SSS*FPSS(K,I)+FES*FPFE(K,I))/H(K)*A1(K,I)
!          END DO
!        END DO
!        DO I=IU,ID
!          SODD(KT,I) = SOD(I)/BHKT2(I)*OMRM(KT,I)*(B(KTI(I),I)
!     .                 -B(KT+1,I))
!          DO K=KT+1,KB(I)-1
!            SODD(K,I) = SOD(I)/BH(K,I)*OMRM(K,I)*(B(K,I)-B(K+1,I))
!          END DO
!          SODD(KB(I),I) = SOD(I)/BH(KB(I),I)*OMRM(KB(I),I)*B(KB(I),I)
!        END DO
!      RETURN
!

!************************************************************************
!**                        L A B I L E   D O M                         **
!************************************************************************

!      ENTRY LABILE_DOM
!        DO I=IU,ID
!          DO K=KT,KB(I)
!            DECAY       = OMRM(K,I)*A3(K,I)*(LDOMDK+LRDDK)*LDOM(K,I)
!            APR         = (AER(K,I)+(1.0-APOM)*AMR(K,I))*ALGAE(K,I)
!            LDOMSS(K,I) = APR-DECAY
!          END DO
!        END DO
!      RETURN

!************************************************************************
!**                     R E F R A C T O R Y   D O M                    **
!************************************************************************

!      ENTRY REFRACTORY_DOM
!        DO I=IU,ID
!          DO K=KT,KB(I)
!            RDOMSS(K,I) = OMRM(K,I)*(LRDDK*LDOM(K,I)-RDOMDK
!     .                    *RDOM(K,I))*A3(K,I)
!          END DO
!        END DO
!      RETURN

!************************************************************************
!**                      P H Y T O P L A N K T O N                     **
!************************************************************************

!      ENTRY PHYTOPLANKTON
!        LTCOEF = (1.0-BETA)*SRO*4.186E6/ASAT
!        DO I=IU,ID
!
!********* Limiting factor


!          GAMMA = EXH2O+EXSS*SS(KT,I)+EXOM*(ALGAE(KT,I)+LPOM(KT,I))
!          LAM1  = LTCOEF
!          LAM2  = LTCOEF*EXP(-GAMMA*DEPTHM(KT))
!          LLIM  = 2.718282*(EXP(-LAM2)-EXP(-LAM1))/(GAMMA*HKT2(I))
!          FDPO4 = 1.0-FPSS(KT,I)-FPFE(KT,I)
!          PLIM  = FDPO4*PO4(KT,I)/(FDPO4*PO4(KT,I)+AHSP)
!          NLIM  = (NH4(KT,I)+NO3(KT,I))/(NH4(KT,I)+NO3(KT,I)+AHSN)
!          LIMIT = MIN(PLIM,NLIM,LLIM)
!          IF (LIMIT.EQ.PLIM) THEN
!            WRITE (LFAC(KT,I),'(F4.3)') PLIM
!            LF         = ' P '
!            LFPR(KT,I) = LF//LFAC(KT,I)
!          ELSE IF (LIMIT.EQ.NLIM) THEN
!            WRITE (LFAC(KT,I),'(F4.3)') NLIM
!            LF         = ' N '
!            LFPR(KT,I) = LF//LFAC(KT,I)
!          ELSE IF (LIMIT.EQ.LLIM) THEN
!            WRITE (LFAC(KT,I),'(F4.3)') LLIM
!            LF         = ' L '
!            LFPR(KT,I) = LF//LFAC(KT,I)
!          END IF

!********* Sources/sinks!
!          ARR(KT,I)   = ARMR(KT,I)*ARMF(KT,I)*AR*A3(KT,I)
!          AMR(KT,I)   = (ARMR(KT,I)+1.0-ARMF(KT,I))*AM
!          AGR(KT,I)   = ARMR(KT,I)*ARMF(KT,I)*AG*LIMIT
!          AGR(KT,I)   = MIN(AGR(KT,I),PO4(KT,I)/(BIOP*DLT*ALGAE(KT,I)
!     .                  +NONZERO),(NH4(KT,I)+NO3(KT,I))/(BION*DLT
!     .                  *ALGAE(KT,I)+NONZERO))
!          AER(KT,I)   = MIN((1.0-LLIM)*AE,AGR(KT,I))
!          GROWTH      = (AGR(KT,I)-ARR(KT,I)-AER(KT,I)-AMR(KT,I))
!     .                  *ALGAE(KT,I)
!          NETS        = -AS*ALGAE(KT,I)/HKT2(I)
!          ALGSS(KT,I) = GROWTH+NETS
!          DO K=KT+1,KB(I)


!*********** Limiting factor


!            GAMMA = EXH2O+EXSS*SS(K,I)+EXOM*(ALGAE(K,I)+LPOM(K,I))
!            LAM1  = LTCOEF*EXP(-GAMMA*DEPTHM(K))
!            LAM2  = LTCOEF*EXP(-GAMMA*DEPTHM(K+1))
!            LLIM  = 2.718282*(EXP(-LAM2)-EXP(-LAM1))/(GAMMA*H(K))
!            FDPO4 = 1.0-FPSS(K,I)-FPFE(K,I)
!            PLIM  = FDPO4*PO4(K,I)/(FDPO4*PO4(K,I)+AHSP)
!            NLIM  = (NH4(K,I)+NO3(K,I))/(NH4(K,I)+NO3(K,I)+AHSN)
!            LIMIT = MIN(PLIM,NLIM,LLIM)
!            IF (LIMIT.EQ.PLIM) THEN
!              WRITE (LFAC(K,I),'(F4.3)') PLIM
!              LF        = ' P '
!              LFPR(K,I) = LF//LFAC(K,I)
!            ELSE IF (LIMIT.EQ.NLIM) THEN
!              WRITE (LFAC(K,I),'(F4.3)') NLIM
!              LF        = ' N '
!              LFPR(K,I) = LF//LFAC(K,I)
!            ELSE IF (LIMIT.EQ.LLIM) THEN
!              WRITE (LFAC(K,I),'(F4.3)') LLIM
!              LF        = ' L '
!              LFPR(K,I) = LF//LFAC(K,I)
!            END IF

!*********** Sources/sinks

!            ARR(K,I)   = ARMR(K,I)*ARMF(K,I)*AR*A3(K,I)
!            AMR(K,I)   = (ARMR(K,I)+1.0-ARMF(K,I))*AM
!            AGR(K,I)   = ARMR(K,I)*ARMF(K,I)*AG*LIMIT
!            AGR(K,I)   = MIN(AGR(K,I),PO4(K,I)/(BIOP*DLT*ALGAE(K,I)
!     .                   +NONZERO),(NH4(K,I)+NO3(K,I))/(BION*DLT
!     .                   *ALGAE(K,I)+NONZERO))
!            AER(K,I)   = MIN((1.0-LLIM)*AE,AGR(K,I))
!            GROWTH     = (AGR(K,I)-ARR(K,I)-AER(K,I)-AMR(K,I))
!     .                   *ALGAE(K,I)
!            NETS       = AS*(ALGAE(K-1,I)-ALGAE(K,I))/H(K)
!            ALGSS(K,I) = GROWTH+NETS
!          END DO
!        END DO
!      RETURN

!************************************************************************
!**                         L A B I L E   P O M                        **
!************************************************************************

!      ENTRY LABILE_POM
!        DO I=IU,ID
!          APR          = APOM*AMR(KT,I)*ALGAE(KT,I)
!          NETS         = -POMS*LPOM(KT,I)/HKT2(I)
!          LPOMSS(KT,I) = APR-LPOMD(KT,I)+NETS
!          DO K=KT+1,KB(I)
!            APR         = APOM*AMR(K,I)*ALGAE(K,I)
!            NETS        = POMS*(LPOM(K-1,I)-LPOM(K,I))/H(K)
!            LPOMSS(K,I) = APR-LPOMD(K,I)+NETS
!          END DO
!        END DO
!      RETURN

!************************************************************************
!**                         P H O S P H O R U S                        **
!************************************************************************

!      ENTRY PHOSPHORUS
!        DO I=IU,ID
!          DO K=KT,KB(I)
!            APR        = (ARR(K,I)-AGR(K,I))*ALGAE(K,I)
!            PO4SS(K,I) = BIOP*(APR+LPOMD(K,I)+DOMD(K,I)+SEDD(K,I))
!     .                   +PO4R*SODD(K,I)*A2(K,I)+SETIN(K,I)
!     .                   *PO4(K-1,I)-SETOUT(K,I)*PO4(K,I)
!          END DO
!        END DO
!      RETURN

!************************************************************************
!**                          A M M O N I U M                           **
!************************************************************************

!      ENTRY AMMONIUM
!        DO I=IU,ID
!          DO K=KT,KB(I)
!            APR        = (ARR(K,I)-AGR(K,I)*NH4(K,I)/(NH4(K,I)
!     .                   +NO3(K,I)+NONZERO))*ALGAE(K,I)
!            NH4SS(K,I) = BION*(APR+LPOMD(K,I)+DOMD(K,I)+SEDD(K,I))
!     .                   +NH4R*SODD(K,I)*A2(K,I)-NH4D(K,I)
!          END DO
!        END DO
!      RETURN

!************************************************************************
!**                            N I T R A T E                           **
!************************************************************************

!      ENTRY NITRATE
!        DO I=IU,ID
!          DO K=KT,KB(I)
!            ALGC       = BION*(1.0-NH4(K,I)/(NH4(K,I)+NO3(K,I)+NONZERO))
!     .                   *AGR(K,I)*ALGAE(K,I)
!            NO3SS(K,I) = NH4D(K,I)-NO3D(K,I)-ALGC 
!          END DO
!        END DO
!      RETURN


!************************************************************************
!**                           S E D I M E N T                          **
!************************************************************************

!      ENTRY SEDIMENT
!        DO I=IU,ID
!          SETTLE    = (AS*ALGAE(KT,I)+POMS*LPOM(KT,I))*DLT
!     .                /HKT2(I)*(1.0-B(KT+1,I)/BKT(I))
!          SED(KT,I) = MAX(SED(KT,I)+SETTLE-SEDD(KT,I)*DLT,0.0)
!          DO K=KT+1,KB(I)-1
!            SETTLE   = (AS*ALGAE(K,I)+POMS*LPOM(K,I))*DLT/H(K)
!     .                 *(1.0-B(K+1,I)/B(K,I))
!            SED(K,I) = MAX(SED(K,I)+SETTLE-SEDD(K,I)*DLT,0.0)
!          END DO
!          SETTLE       = (AS*ALGAE(KB(I),I)+POMS*LPOM(KB(I),I))
!     .                   *DLT/H(KB(I))
!          SED(KB(I),I) = MAX(SED(KB(I),I)+SETTLE-SEDD(KB(I),I)
!     .                   *DLT,0.0)
!        END DO
!      RETURN

!************************************************************************
!**                   I N O R G A N I C   C A R B O N                  **
!************************************************************************
!      ENTRY INORGANIC_CARBON
!        CO2EX = (0.5+0.05*WIND*WIND)/86400.0
!        DO I=IU,ID
!          TICSS(KT,I) = 0.0
!          DO K=KT,KB(I)
!            APR        = (ARR(K,I)-AGR(K,I))*ALGAE(K,I)
!            TICSS(K,I) = BIOC*(APR+DOMD(K,I)+LPOMD(K,I)+SEDD(K,I))
!     .                   +CO2R*SODD(K,I)*A3(K,I)
!          END DO
!          IF (.NOT.ICE(I)) TICSS(KT,I) = TICSS(KT,I)+CO2EX*(0.286
!     .                                   *EXP(-0.0314*(T2(KT,I))*PALT)
!     .                                   -CO2(KT,I))/HKT2(I)
!        END DO
!      RETURN

!************************************************************************
!**                             P H   C O 2                            **
!************************************************************************

!      ENTRY PH_CO2
!        DO I=IU,ID
!          DO K=KT,KB(I)
!            CARB = TIC(K,I)/1.2E4
!            ALK  = ALKAL(K,I)/5.0E4
!            T1K  = T1(K,I)+273.15

!*********** Activity equilibrium constants

!            KW = 10.0**(35.3944-5242.39/T1K-0.00835*T1K-11.8261
!     .           *LOG10(T1K))
!            K1 = 10.0**(14.8435-3404.71/T1K-0.032786*T1K)
!            K2 = 10.0**(6.498-2902.39/T1K-0.02379*T1K)

!*********** Ionic strength

!            IF (FRESH_WATER) THEN
!              S2 = 2.5E-05*TDS(K,I)
!            ELSE
!              S2 = 0.00147+0.019885*TDS(K,I)+0.000038*TDS(K,I)**2
!            END IF

!*********** Debye-Huckel terms and activity coefficients

!            SQRS2  = SQRT(S2)
!            HCO3T  = 10.0**(-0.5085*SQRS2/(1.0+1.3124*SQRS2)+4.745694E-03
!     .               +4.160762E-02*S2-9.284843E-03*S2**2)
!            CO3T   = 10.0**(-2.0340*SQRS2/(1.0+1.4765*SQRS2)+1.205665E-02
!     .               +9.715745E-02*S2-2.067746E-02*S2**2)
!            H2CO3T = 0.0755*S2
!            KW     = KW/HCO3T
!            K1     = K1*H2CO3T/HCO3T
!            K2     = K2*HCO3T/CO3T

!*********** pH evaluation

!            PH_LEFT  = -14.0
!            PH_RIGHT = 0.0
!            HION     = 10.0**PH_LEFT
!            BICARB   = CARB*K1*HION/(K1*HION+K1*K2+HION**2)
!            F_LEFT   = BICARB*(HION+2.0*K2)/HION+KW/HION-ALK-HION
!            IF (F_LEFT.LT.0) THEN
!             PH_START = PH_LEFT
!             PH_DIFF  = PH_RIGHT-PH_LEFT
!            ELSE
!             PH_START = PH_RIGHT
!             PH_DIFF  = PH_LEFT-PH_RIGHT
!            ENDIF
!            J = 0
!            DO WHILE (J.LT.50)
!              J       = J+1
!              PH_DIFF = PH_DIFF*0.5
!              PH_MID  = PH_START+PH_DIFF
!              HION    = 10.0**PH_MID
!              BICARB  = CARB*K1*HION/(K1*HION+K1*K2+HION**2)
!              FMID    = BICARB*(HION+2.0*K2)/HION+KW/HION-ALK-HION
!              IF (FMID.LT.0) PH_START = PH_MID
!              IF (ABS(PH_DIFF).LT.0.01.OR.FMID.EQ.0.) J = 51!              IF (ABS(PH_DIFF).LT.0.01.OR.FMID.EQ.0.) J = 51
!            ENDDO

!*********** pH, carbon dioxide, bicarbonate, and carbonate concentrations

!            PH(K,I)   = -PH_MID
!            CO2(K,I)  = TIC(K,I)/(1.0+K1/HION+K1*K2/HION**2)
!            HCO3(K,I) = TIC(K,I)/(1.0+HION/K1+K2/HION)
!            CO3(K,I)  = TIC(K,I)/((HION*HION)/(K1*K2)+HION/K2+1.0)
!          END DO
!        END DO
!      RETURN

