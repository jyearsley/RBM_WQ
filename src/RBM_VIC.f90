!
      PROGRAM RBM10_VIC
!
!     Dynamic river basin model for simulating water quality in
!     branching river systems with freely-flowing river segments. 
!
!     This version uses Reverse Particle Tracking in the Lagrangian
!     mode and Lagrangian interpolation in the Eulerian mode.
!
!     This version of the model software has the following limited features:
!     
!     1. Hydrologic and meteorologic forcings are from a direct access file
!        prepared by the routing program.  Hydraulic parameters are estimated
!        by the routing program using Leopold relationships describing depth
!        and speed as a function of flow.
!
!     2. Output is written to Unit 20 (specified in the command line for executing
!        RBM) and includes the simulated temperatures for the first computational cell
!        in each segment in the same order as that of the *.network file and at the 
!        same time step.  
!
!     Topology and routing is set up to be consistent with output
!     from the Variable Infiltration Capacity (VIC) model developed by the
!     Land Surface Hydrology Group at the University of Washington.
!     Model details are described in
!    
!     Yearsley, J. (2012), A grid-based approach for simulating stream temperature,
!     Water Resour. Res., 48, W03506, doi:10.1029/2011WR011515
!
!     For additional information contact:
!
!     John Yearsley
!     Land Surface Hydrology Group
!     Dept. of Civil and Environmental Engineering
!     Box 352700
!     University of Washington
!     Seattle, Washington
!     98195-2700
!
use Block_File_Names
!
!
implicit none
!
!

integer(i2b) iargc
integer(i2b) numarg

!
! Command line input
!
numarg = iargc ( )
if (numarg .lt. 2) then
  write (*,*) 'Too few arguments were given'
  write (*,*) ' '
  write (*,*) 'First:  Location and prefix of input files'
  write (*,*) '        (networkfile and parameterfile)'
  write (*,*) 'Second: Location and prefix of output files'
  write (*,*) ' '
  write (*,*) 'eg: $ <program-name> ./input/Salmon_0.50 ./output/Salmon_Test'
  write (*,*) ' '
  stop
end if
 call getarg ( 1, inPrefix )
 call getarg ( 2, outPrefix )
e
!
!     Open file with hydrologic data
!
open(unit=35,FILE=TRIM(flow_file) ,FORM='FORMATTED',ACCESS='DIRECT' ,RECL=60,STATUS='old')
!
! If we're simulating water temperature, open the file with the heat budget
!
if (do_TEMP) then
  read(90,'(A)') heat_file
!
!     Open file with meteorologic data
!     
  open(unit=36,FILE=TRIM(heat_file) ,FORM='FORMATTED',ACCESS='DIRECT' ,RECL=50,STATUS='old')
end if
!
!     Call systems programs to get started
!
!     SUBROUTINE BEGIN reads control file, sets up topology and
!     important properties of reaches
!
write(*,*) 'Calling BEGIN'
!
!     SUBROUTINE BEGIN reads in river system information from the NETWORK file
!
CALL BEGIN
!
!     SUBROUTINE SYSTMM performs the simulations
!
CALL SYSTMM(temp_file,param_file) ! (WUR_WF_MvV_2011/01/05)
!
!     Close files after simulation is complete
!
write(*,*) ' Closing files after simulation'

CLOSE(35)
CLOSE(36)
CLOSE(90)
STOP
END PROGRAM RBM10_VIC
