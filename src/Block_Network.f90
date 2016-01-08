Module Block_Network
!
! Module with stream topology variables
!
! integer(i2b) arrays
!
    integer(i2b), dimension(:), allocatable   :: no_celm,no_cells,ndelta,no_tribs
    integer(i2b), dimension(:), allocatable   :: head_cell
    integer(i2b), dimension(:), allocatable   :: ncol,nrow,no_rows,no_cols
    integer(i2b), dimension(:,:), allocatable :: segment_cell,trib
!
! integer(i2b) variables 
!
    integer(i2b)                              :: flow_cells,heat_cells
    integer(i2b)                              :: ndays,nreach,ntrb,nwpd
    integer(i2b),parameter                    :: ns_max=200
    integer(i2b)                              :: start_year,start_month,start_day
    integer(i2b)                              :: end_year,end_month,end_day
    integer(i2b)                              :: n_default=2
!
! Real(sp) variables
!
    real(sp)                                   :: delta_n,dt_comp
end module Block_Network
