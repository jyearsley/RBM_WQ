Subroutine BEGIN(param_file,spatial_file)
!
use Block_Energy
use Block_File_Names
use Block_Hydro
use Block_Network
!
implicit none
!    
    character(ch)                  :: end_date,start_date     
    character(ch)                  :: lat
    character(len=10)              :: long
    character(len=200)             :: param_file,source_file,spatial_file
    character(len=5),dimension(20) :: WQ_const
!
    integer(i2b)                   :: Julian
    integer(i2b)                   :: head_name,trib_cell
    integer(i2b)                   :: jul_start,main_stem,nyear1,nyear2,nc,ncell,nseg
    integer(i2b)                   :: ns_max_test,nndlta,node,nr,cum_sgmnt
    integer(i2b)                   :: ns,no_sources
    integer(i2b)                   :: nq,nnq
!
    logical(lg)                    :: first_cell
!
    real(sp)                       :: rmile0,rmile1,xwpd
!
!
! Identify input/output files
!
! Identify the network file
!
  Net_File       = TRIM(inPrefix)//'_Network'
  write(*,*) 'Network file:     '      ,Net_File
!
! Identify other necessary files
!
  Param_File     = TRIM(inPrefix)//'_Parameters'
  write(*,*) 'Parameter file:   '      ,Param_File
!
  Source_File    = TRIM(inPrefix)//'_Source'
  write(*,*) 'Source file:            ',Source_File
!
  Spatial_File   = TRIM(outPrefix)//'_Spat'
  write(*,*) 'Spatial file:           ',Spatial_File
!
  CONST_file     = TRIM(outPrefix)//'_CONST'
  write(*,*) 'Constituent output file: ',CONST_File
!
!
OPEN(UNIT=90,FILE=TRIM(Net_File),STATUS='OLD')
!
!     Read header information from control file
!
read(90,*)
read(90,'(A)') flow_file
!
! Define constituents to be modeled
! 
NQ_total=0
!
read(90,'(20a5)',ADVANCE='NO',EOR=100,END=100) WQ_const
!
100 continue
!
do nq=1,20
  do nnq=1,20
    if (TRIM(WQ_const(nq)) .eq. WQ_token(nnq)) then
!
! Test here to see if this will include a simulation of water temperature
! Note: The constituent number for water temperature (nnq) = 1
!
      if (nnq .eq. 1) do_TEMP = .TRUE.
      write(*,*) 'WQ_test ',nnq,WQ_const(nq)
      NQ_total=NQ_total+1
!
! Establish constituent table
!
      CONST_Table(nq)=nnq
!      
    end if
  end do
end do
!
! Open the hydrologic forcing files
 

! Open the heat forcing files if temperature is to be simulated
!
! Fetch the name of the directory with constituent input files

!
!   Mohseni parameters, if used
!
!
!
!     Card Group I
!
read(90,*) start_date,end_date
read(start_date,'(i4,2i2)') start_year,start_month,start_day
read(end_date,  '(i4,2i2)') end_year,end_month,end_day
nyear1=start_year
nyear2=end_year
write(*,'(2(2x,i4,2i2))')  &
 start_year,start_month,start_day,end_year,end_month,end_day
!
!     Establish the Julian day for which simulations begin
!
jul_start = Julian(start_year,start_month,start_day)
!
read(90,*) nreach,flow_cells,heat_cells,no_rows,no_cols
!
! Allocate dynamic arrays
!
 allocate(ndelta(heat_cells))
 allocate(mu(nreach))
 allocate(alphamu(nreach))
 allocate(beta(nreach))
 allocate(gmma(nreach))
 allocate (smooth_param(nreach))
 allocate(dx(heat_cells))
 allocate(no_celm(nreach))
 no_celm=0
 allocate(no_cells(nreach))
 no_cells=0
 allocate(no_tribs(heat_cells))
 no_tribs=0
 allocate(trib(heat_cells,10))
 trib=0
 allocate(head_cell(nreach))
 allocate(segment_cell(nreach,ns_max))
 allocate(x_dist(nreach,0:ns_max))

!
! Read the source files using the constituent table (CONST_Table)
!
   read(90,'(A)') source_file ! (WUR_WF_MvV_2011/05/23)
   print *,'source file: ', source_file ! (WUR_WF_MvV_2011/05/23)
   open(40,file=TRIM(source_file),status='old')
!
! Initialize headwaters values for those constituents that
! are being simulated
!
do nr = 1,nreach
  read(40,*) (CONST_init(nr,nq)=1,NQ_Total)
end do
!
end if
!
!     Start reading the reach date and initialize the reach index, NR
!     and the cell index, NCELL
!
ncell=0
!
ns_max_test=-1
!
!     Card Group IIb. Reach characteristics
!s
do nr=1,nreach
!
!     Initialize NSEG, the total number of segments in this reach
!
  nseg=0
  write(*,*) ' Starting to read reach ',nr
!
!     Read the number of cells in this reach, the headwater #,
!     the number of the cell where it enters the next higher order stream,
!     the headwater number of the next higher order stream it enters, and
!     the river mile of the headwaters.
!
  read(90,'(i5,11x,i4,10x,i5,15x,i5,15x,f10.0,i5)') no_cells(nr) &
      ,head_name,trib_cell,main_stem,rmile0
!
!     If this is reach that is tributary to cell TRIB_CELL, give it the
!     pointer TRIB(TRIB_CELL) the index of this reach for further use.
!     Also keep track of the total number of tributaries for this cell
!
  if (trib_cell.gt.0) then
    no_tribs(trib_cell)=no_tribs(trib_cell)+1
    trib(trib_cell,no_tribs(trib_cell))=nr
  end if
*****************************************************************************
!   DON"T FORGET TO REMOVE THIS
!     Reading Mohseni parameters for each headwaters (UW_JRY_2011/06/18)
!
!  read(90,*) alphaMu(nr),beta(nr) &
!            ,gmma(nr),mu(nr),smooth_param(nr)
*****************************************************************************
!
!     Reading Reach Element information
!
  first_cell=.true.
  do nc=1,no_cells(nr)
    ncell=ncell+1
!
!   Read the data for point sources
!
    if (source) then
!
!  Place holder for point source input
!
    end if 
!
!     The headwaters index for each cell in this reach is given
!     in the order the cells are read
!
!     Card Type 3. Cell indexing #, Node # Row # Column Lat Long RM
!
!     Variable ndelta read in here.  At present, number of elements
!     is entered manually into the network file (UW_JRY_2011/03/15)
!
    read(90,'(5x,i5,5x,i5,8x,i5,6x,a8,6x,a10,7x,f10.0,i5)')  &
              node,nrow(ncell),ncol(ncell),lat,long,rmile1,ndelta(ncell)
!
!    Set the number of segments of the default, if not specified
!
    if (ndelta(ncell).lt.1) ndelta(ncell)=n_default
    if(first_cell) then
      first_cell=.false.
      head_cell(nr)=ncell
      x_dist(nr,0)=5280.*rmile0
    end if
!
! Added variable ndelta (UW_JRY_2011/03/15)
!
    dx(ncell)=5280.*(rmile0-rmile1)/ndelta(ncell)
    rmile0=rmile1
    nndlta=0
200 continue
    nndlta=nndlta+1
    nseg=nseg+1
    segment_cell(nr,nseg)=ncell
    x_dist(nr,nseg)=x_dist(nr,nseg-1)-dx(ncell)
!
!   Write Segment List for mapping to temperature output (UW_JRY_2008/11/19)
!
    open(22,file=TRIM(spatial_file),status='unknown') ! (changed by WUR_WF_MvV_2011/01/05)
    write(22,'(4i6,1x,a8,1x,a10,i5)') nr,ncell,nrow(ncell),ncol(ncell),lat,long,nndlta
!
! 
!
!  Added variable ndelta  (UW_JRY_2011/03/15)
!
    if(nndlta.lt.ndelta(ncell)) go to 200  
    no_celm(nr)=nseg
    segment_cell(nr,nseg)=ncell
    x_dist(nr,nseg)=5280.*rmile1
!
! End of segment loop
!
  end do
if(ns_max_test.lt.nseg) ns_max_test=nseg
!
! End of reach loop
!
end do
if(ns_max_test.gt.ns_max) then
  write(*,*) 'RBM is terminating because'
  write(*,*) 'NS_MAX exceeded. Change NS_MAX in Block_Network to: ',ns_max_test
  stop
end if
!
! NOTE!!: This is a hardwired limitation on the number of daily time steps
!
nwpd=1
!
xwpd=nwpd
dt_comp=86400./xwpd
!
!     ******************************************************
!                         Return to RMAIN
!     ******************************************************
!
900 continue
!
!
end subroutine BEGIN
