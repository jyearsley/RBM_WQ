Module Read_Input
!
USE Block_Energy
use Block_File_Names
USE Block_Hydro
USE Block_Network
USE Block_WQ
! 
IMPLICIT NONE
!
contains
SUBROUTINE READ_PARAMETERS
!
! 
integer(I4B), private             :: nc,nq,nr
integer(I4B),dimension(:),private :: nwq
real(SP),dimension(:,:),private     :: dummy_C_init,dummy_Source
!

do nr=1,NREACH
!
! Read Mohseni paramters for headwaters, if water temperature is to 
! be simulated
!
  if (do_TEMP) then
    read(40,*) (MU(nr),ALPHAMU(nr),BETA(nr),GMMA(nr),SMOOTH_PARAM(nr)
  end if
  
  do nc=1,NO_CELLS(nr)
    read(50,*) (nwq(nq),dummy_C_init(nq),nq=1,NQ_total)
    do nq = 1,NQ_total
      CONST_init(nq)=dummy_C_init(nwq(nq))
    end do
    read(60,*) (nwq(nq),dummy_Source(nq),nq=1,NQ_total)
    do nq = 1,NQ_total
      CONST_src(nq)=dummy_Source(nwq(nq))
    end do
!
! Read kinetics parameters here, but postpone until non-conservative 
! constituents like DO, BOD, N, and P are included
!
! READ KINETICS parameters next
!
!
  end do
END SUBROUTINE READ_PARAMETERS
!
! Read the forcing files for hydrology and hydraulics and for the
! heat budget, if temperature is being simulated
!  
SUBROUTINE READ_FORCING
!
integer(I4B),private :: nc,ncell,nnd,no_flow,no_heat,nr,nrec_flow,nrec_heat

! DO I REALLY NEED Q_avg???????????????????????????????????????
real(SP):: Q_avg
!
no_flow=0
no_heat=0
!
do nr=1,NREACH
  do nc=1,NO_CELLS(nr)-1 
    no_flow=no_flow+1
    no_heat=no_heat+1
!
    nrec_flow=FLOW_CELLS*(ndays-1)+no_flow
    nrec_heat=HEAT_CELLS*(ndays-1)+no_heat
!
    read(35,'(2i5,2f10.1,2f6.1,f7.1,f6.2)' &
           ,rec=nrec_flow) nnd,ncell &
           ,Q_In(no_heat),Q_Out(no_heat),Q_Diff(no_heat) &  
           ,DEPTH(no_heat),WIDTH(no_heat),U(no_heat)

!
    if(u(no_heat).lt.0.01) u(no_heat)=0.01
    if(ncell.ne.no_heat) write(*,*) 'Flow file error',ncell,no_heat 
!
!

    delta_n=NDELTA(ncell)
! 
    Q_avg=0.5*(Q_IN(no_heat)+Q_OUT(no_heat))
    Q_DIFF(no_heat)=Q_DIFF(no_heat)/delta_n
    dt(no_heat)=DX(no_heat)/U(no_heat)

!
!  Added check to see if travel time of parcel exceeds the
!  computational interval.  If so, it writes to file fort.45.
!  One should check to see if there are NaN's at that node.
!  (UW_JRY_2011/03/15)
!
    if(dt(no_heat).gt.dt_comp) write(45,*) &
           'Travel time=',dt(no_heat) &
            , '> dt_comp at node -',no_heat
    if (do_TEMP) then
      read(36,'(i5,2f6.1,2f7.4,f6.3,f7.1,f5.1)'                        &
               ,rec=nrec_heat) ncell                                   &
               ,DBT(no_heat),EA(no_heat)                               &
               ,Q_NS(no_heat),Q_NA(no_heat),RHO                        & 
               ,PRESS(no_heat),WIND(no_heat)
    end if
  end do
!
!      Tributary flow is Q_out from the next to the last cell
!
!
!       Read the meteorology for the last cell, but not the flow
!
  if (do_TEMP)
    no_heat=no_heat+1 
    Q_IN(no_heat)=Q_OUT(no_heat-1)
    Q_OUT(no_heat)=Q_IN(no_heat)
    Q_trib(nr)=Q_OUT(no_heat)    
    nrec_heat=heat_cells*(ndays-1)+no_heat
    read(36,'(i5,2f6.1,2f7.4,f6.3,f7.1,f5.1)'                          &
             ,rec=nrec_heat) ncell                                     &
             ,DBT(no_heat),EA(no_heat)                                 &   
             ,Q_NS(no_heat),Q_NA(no_heat),RHO                          &
             ,PRES(no_heat),WIND(no_heat)
  end if
!
!  The flow and hydraulics for the last cell has to be 
!  modified so they do not
!  take the values of the segment to which it is tributary
!
  Q_in(ncell)=Q_out(ncell-1)
!  Q_out(ncell)=Q_in(ncell-1)
  Q_diff(no_heat)=0.0
  u(no_heat)=u(no_heat-1)
  depth(no_heat)=depth(no_heat-1)
  width(no_heat)=width(no_heat-1)
  dt(no_heat)=0.5*dx(ncell)/u(no_heat)
end do
END SUBROUTINE Read_Forcing
!
! End of module
!
END MODULE Read_Input
