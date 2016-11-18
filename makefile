#-----------------------------------------------------------------------
#
#  Makefile for Solver
#  Friedrich Ulrich, 2016
#
#-----------------------------------------------------------------------
# Linker that is usedf
LL	= caf
# linker flags
# for single image compilation
#LFLAGS= -fcoarray=single -ffree-line-length-none 
# compiler that is used
CC	= caf
# compiler flags
# for single image compilation
#CFLAGS= -fcoarray=single -ffree-line-length-none -g 
# for multi-image compilation
# CFLAGS= 
# All Source Files
SRCS	= main.f90 inputfile.f90 VARIABLEN.f90 INI_MOD.f90 energy_tot_mod.f90 solver_mod.f90 vel_feld_mod.f90 druck_feld_mod.f90 temp_feld_mod.f90 xmom_mod.f90 ymom_mod.f90 energy_mod.f90 check_conti_mod.f90 conti_mod.f90 zmom_mod.f90 sync_mod.f90 derivatives_x_mod.f90 sutherland_mod.f90 wang_johnson_mod.f90 derivatives_y_mod.f90 derivatives_z_mod.f90 spannungstensor_mod.f90 print_all_mod.f90 boundry_condition_mod.f90 neumann_mod.f90 energie-feld_mod.f90 lagrange_intpol.f nin.f test_mod.f90 derivatives_y2_mod.f90 char_bc.f
# wall_clock.f90 INI_MOD.f90 print_all_mod.f90 wang_johnson_mod.f90 explizit_mod.f90 sync_mod.f90 implicit_mod.f90 navier_mod.f90 output_mod.f90 wang_johnson_2_mod.f90 nav_sim_mod.f90
# all Ojects
OBJECTS	= main.o inputfile.o VARIABLEN.o INI_MOD.o energy_tot_mod.o solver_mod.o vel_feld_mod.o druck_feld_mod.o temp_feld_mod.o xmom_mod.o ymom_mod.o energy_mod.o check_conti_mod.o conti_mod.o zmom_mod.o sync_mod.o derivatives_x_mod.o sutherland_mod.o wang_johnson_mod.o derivatives_y_mod.o derivatives_z_mod.o spannungstensor_mod.o print_all_mod.o boundry_condition_mod.o neumann_mod.o energie-feld_mod.o lagrange_intpol.o nin.o test_mod.o derivatives_y2_mod.o char_bc.o
#wall_clock.o  INI_MOD.o print_all_mod.o explizit_mod.o sync_mod.o implicit_mod.o output_mod.o wang_johnson_mod.o navier_mod.o nav_sim_mod.o wang_johnson_2_mod.o
MODS = inputfile.mod variablen.mod ini_mod.mod energy_tot_mod.mod solver_mod.mod vel_feld_mod.mod druck_feld_mod.mod temp_feld_mod.mod xmom_mod.mod ymom_mod.mod energy_mod.mod check_conti_mod.mod conti_mod.mod zmom_mod.mod sync_mod.mod derivatives_x_mod.mod sutherland_mod.mod wang_johnson_mod.mod derivatives_y_mod.mod derivatives_z_mod.mod spannungstensor_mod.mod print_all_mod.mod boundry_condition_mod.mod neumann_mod.mod energie-feld_mod.mod lagrange_intpol.mod nin.mod test_mod.mod derivatives_y2_mod.mod char_bc.mod
#wall_clock.mod  ini_mod.mod print_all_mod.mod explizit_mod.mod sync_mod.mod implicit_mod.mod output_mod.mod wang_johnson_mod.mod navier_mod.mod wang_johnson_2_mod.mod nav_sim_mod.mod
# Name of the executable file
PROGRAM = eulerrun.x
# file suffixes
.SUFFIXES:
.SUFFIXES: .f .f90 .o  .mod
# compile but not linking
.f.o:	
	$(CC) $(CFLAGS) -c $< 
	
.f90.o:	
	$(CC) $(CFLAGS) -c $< 



all:	$(PROGRAM)
$(PROGRAM):	$(OBJECTS)
	$(LL)  $(LFLAGS) $(OBJECTS) -o $(PROGRAM) 
# dependencies
main.o: inputfile.o VARIABLEN.o INI_MOD.o energy_tot_mod.o solver_mod.o vel_feld_mod.o druck_feld_mod.o temp_feld_mod.o print_all_mod.o boundry_condition_mod.o energie-feld_mod.o test_mod.o
#wall_clock.o  explizit_mod.o implicit_mod.o output_mod.o INI_MOD.o navier_mod.o nav_sim_mod.o
sync_mod.o: VARIABLEN.o  
wang_johnson_mod.o: VARIABLEN.o 
#explizit_mod.o: VARIABLEN.o sync_mod.o
#implicit_mod.o: wang_johnson_mod.o VARIABLEN.o sync_mod.o
#output_mod.o: VARIABLEN.o print_all_mod.o
print_all_mod.o: VARIABLEN.o
INI_MOD.o: VARIABLEN.o lagrange_intpol.o derivatives_x_mod.o derivatives_y_mod.o nin.o
lagrange_intpol.o: VARIABLEN.o
energy_tot_mod.o: VARIABLEN.o
solver_mod.o: VARIABLEN.o xmom_mod.o ymom_mod.o energy_mod.o check_conti_mod.o conti_mod.o zmom_mod.o sync_mod.o sutherland_mod.o spannungstensor_mod.o
sutherland_mod.o: VARIABLEN.o
druck_feld_mod.o: VARIABLEN.o
vel_feld_mod.o: VARIABLEN.o
temp_feld_mod.o: VARIABLEN.o
energie-feld_mod.o: VARIABLEN.o
nin.o: VARIABLEN.o
test_mod.o: VARIABLEN.o derivatives_x_mod.o derivatives_y_mod.o derivatives_z_mod.o 
xmom_mod.o: VARIABLEN.o derivatives_x_mod.o derivatives_y_mod.o derivatives_z_mod.o neumann_mod.o
ymom_mod.o: VARIABLEN.o derivatives_x_mod.o derivatives_y_mod.o derivatives_z_mod.o neumann_mod.o
zmom_mod.o: VARIABLEN.o derivatives_x_mod.o derivatives_y_mod.o derivatives_z_mod.o neumann_mod.o
energy_mod.o: VARIABLEN.o derivatives_x_mod.o derivatives_y_mod.o derivatives_z_mod.o neumann_mod.o
check_conti_mod.o: VARIABLEN.o neumann_mod.o
conti_mod.o: VARIABLEN.o derivatives_x_mod.o derivatives_y_mod.o derivatives_z_mod.o neumann_mod.o
spannungstensor_mod.o: VARIABLEN.o derivatives_x_mod.o derivatives_y_mod.o derivatives_z_mod.o neumann_mod.o
derivatives_x_mod.o: VARIABLEN.o sync_mod.o wang_johnson_mod.o
derivatives_y_mod.o: VARIABLEN.o sync_mod.o wang_johnson_mod.o derivatives_y2_mod.o
derivatives_y2_mod.o: VARIABLEN.o sync_mod.o wang_johnson_mod.o
derivatives_z_mod.o: VARIABLEN.o sync_mod.o derivatives_y_mod.o derivatives_y2_mod.o
boundry_condition_mod.o: VARIABLEN.o sync_mod.o char_bc.o
neumann_mod.o: VARIABLEN.o sync_mod.o
#wang_johnson_2_mod.o: VARIABLEN.o
#nav_sim_mod.o: VARIABLEN.o implicit_mod.o  navier_mod.o
# cleaing all unneeded files
clean:
	-/bin/rm -f $(OBJECTS) $(MODS)
