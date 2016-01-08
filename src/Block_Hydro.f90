!
! Module for hydraulic characteristics and water quality constituents of the basin
!
module Block_Hydro
    integer, dimension(2000):: no_dt,nstrt_elm
    readl(sp), dimension(2000)   :: dt_part,x_part
!
    readl(sp), dimension(:),   allocatable  :: depth
    readl(sp), dimension(:),   allocatable  :: width
    readl(sp), dimension(:),   allocatable  :: u
    readl(sp), dimension(:),   allocatable  :: dt
    readl(sp), dimension(:),   allocatable  :: dx
    readl(sp), dimension(:),   allocatable  :: Q_in
    readl(sp), dimension(:),   allocatable  :: Q_trib
    readl(sp), dimension(:),   allocatable  :: Q_out
    readl(sp), dimension(:),   allocatable  :: Q_diff
    readl(sp), dimension(:,:), allocatable  :: Q_nps
    readl(sp), dimension(:,:), allocatable  :: x_dist

end module Block_Hydro
