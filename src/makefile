#
# Makefile for the grid-based semi-Lagrangian water temperature model, RBM10_VIC
#
# Start of the makefile
#
# Defining variables
#
objects = RBM_VIC.o BEGIN.o SYSTMM.o PARTICLE_TRACK.o \
          ENERGY.o JULIAN.o TNTRP.o READ_INPUT.o \
          Block_Energy.o Block_Hydro.o Block_Network.o Block_WQ.o\
	  WRITE.o	  
f90comp = gfortran
# Makefile
RBM_VIC: $(objects)
	$(f90comp) -o RBM_VIC $(objects)
#bgin.mod: Begin.o Begin.f90
#	$(f90comp) -c BEGIN.f90
Block_WQ.o: Block_WQ.f90
	$(f90comp) -c Block_WQ.f90
block_wq.mod: Block_WQ.o Block_WQ.f90
	$(f90comp) -c Block_WQ.f90
Block_Energy.o: Block_Energy.f90
	$(f90comp) -c Block_Energy.f90
block_energy.mod: Block_Energy.f90
	$(f90comp) -c Block_Energy.f90
Block_Hydro.o: Block_Hydro.f90
	$(f90comp) -c Block_Hydro.f90
block_hydro.mod: Block_Hydro.o Block_Hydro.f90
	$(f90comp) -c Block_Hydro.f90
Block_Network.o: Block_Network.f90
	$(f90comp) -c Block_Network.f90
block_network.mod: Block_Network.o Block_Network.f90
	$(f90comp) -c Block_Network.f90
BEGIN.o: block_energy.mod block_network.mod block_hydro.mod Begin.f90
	$(f90comp) -c BEGIN.f90
READ_INPUT.o: block_energy.mod block_hydro.mod block_network.mod READ_INPUT.f90
	$(f90comp) -c READ_INPUT.f90 
systm.mod: SYSTMM.o SYSTMM.f90
	$(f90comp) -c SYSTMM.f90
SYSTMM.o: block_network.mod block_energy.mod block_hydro.mod SYSTMM.f90
	$(f90comp) -c SYSTMM.f90
ENERGY.o: block_energy.mod ENERGY.f90
	$(f90comp) -c ENERGY.f90
PARTICLE_TRACK.o: block_hydro.mod block_network.mod PARTICLE_TRACK.f90
	$(f90comp) -c PARTICLE_TRACK.f90
WRITE.o: WRITE.f90
	$(f90comp) -c WRITE.f90
JULIAN.o: JULIAN.f90
	$(f90comp) -c JULIAN.f90
TNTRP.o: TNTRP.f90
	$(f90comp) -c TNTRP.f90

RBM_VIC.o: bgin.mod systm.mod RBM_VIC.f90
	$(f90comp) -c RBM_VIC.f90

# Cleaning everything
clean:
#	rm bgin.mod block_energy.mod block_hydro.mod block_network.mod\
	rm block_energy.mod block_hydro.mod block_network.mod\
           systm.mod RBM_VIC
	rm $(objects)
