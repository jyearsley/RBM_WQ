module block_wq_params
!
! Module with parameters for constituent kinetics
!
! Energy budget parameters: Bowen ratio, reflectivity (shortwave and atmospheric)
! and evaporation coefficients
!
     real(sp)       :: evap_a,evap_b,evap_c,rb,r_short,r_atmo 
!
! Dissolved oxygen parameters: Q10,deaeration coefficients
!
     real(sp)       :: Q10_do,reaer_a,reaer_b,reaer_c 
!
! Biochemical oxygen demand parameters: Q10, decay rate
!
     real(sp)       :: Q10_bod,k_bod
!
! Coliform bacteria parameters: Q10, decay rate
!
     real(sp)       :: Q10_colif,k_colif
!
! Total suspended solids parameters: settling velocity
!
     real(sp)       :: w_settl
!
end module block_wq_params
