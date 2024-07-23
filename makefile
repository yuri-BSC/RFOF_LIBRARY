# Variables
COMPILER = gfortran -c
LOAD = gfortran
FLAGS = -g -Wall
EXEC = rftest
OBJ_DIR = ./

OBJFILES = \
$(OBJ_DIR)/itm_types.o \
$(OBJ_DIR)/dum_magnetic_field.o \
$(OBJ_DIR)/RFOF_magnetic_field.o \
$(OBJ_DIR)/RFOF_constants.o \
$(OBJ_DIR)/RFOF_parameters.o \
$(OBJ_DIR)/RFOF_markers.o \
$(OBJ_DIR)/RFOF_local_magnetic_field.o \
$(OBJ_DIR)/Local.o \
$(OBJ_DIR)/Fullwave_mod.o \
$(OBJ_DIR)/GlobalParam_mod.o \
$(OBJ_DIR)/coherentwave_mod.o \
$(OBJ_DIR)/euitm_schemas.o \
$(OBJ_DIR)/euitm_waves_interface.o \
$(OBJ_DIR)/RFOF_waves.o \
$(OBJ_DIR)/RFOF_numerics.o \
$(OBJ_DIR)/RFOF_resonance_memory.o \
$(OBJ_DIR)/RFOF_resonance_condition.o \
$(OBJ_DIR)/RFOF_random_numbers.o \
$(OBJ_DIR)/RFOF_diagnostics.o \
$(OBJ_DIR)/RFOF_types.o \
$(OBJ_DIR)/RFOF_kick.o \
$(OBJ_DIR)/RFOF_mpi_serial.o \
$(OBJ_DIR)/RFOF_Efield_update.o \
$(OBJ_DIR)/RFOF_main.o \
$(OBJ_DIR)/RFOF_wiener_sample_paths.o \
$(OBJ_DIR)/dummy_orbit.o \
$(OBJ_DIR)/rftest.o \
$(OBJ_DIR)/steinbrecher_template.o \
$(OBJ_DIR)/magnetic_field_interface_to_RFOF.o \
#$(OBJ_DIR)/RFOF_Stochastics.o \
#$(OBJ_DIR)/RFOF_BesselProcess0.o \
#$(OBJ_DIR)/RFOF_Bessel0_5_FET_TestMain.o \
#$(OBJ_DIR)/RFOF_BesselProcess0_3.o \
#$(OBJ_DIR)/RFOF_BesselProcess0_5.o \
#$(OBJ_DIR)/RFOF_random_numbers_SG2.o \

#Default target
all: $(EXEC)

#Compile the program
$(EXEC) : $(OBJFILES)
	$(LOAD) $(FLAGS) -o ./$(EXEC) $(OBJFILES)

#Dependencies

$(OBJ_DIR)/itm_types.o : $(OBJ_DIR)/itm_types.f90
	$(COMPILER) -o $(OBJ_DIR)/itm_types.o \
	$(OBJ_DIR)/itm_types.f90

	
$(OBJ_DIR)/RFOF_wiener_sample_paths.o : $(OBJ_DIR)/RFOF_wiener_sample_paths.F90
	$(COMPILER) -o $(OBJ_DIR)/RFOF_wiener_sample_paths.o \
	$(OBJ_DIR)/RFOF_wiener_sample_paths.F90

$(OBJ_DIR)/rftest.o : $(OBJ_DIR)/rftest.F90
	$(COMPILER) -o $(OBJ_DIR)/rftest.o \
	$(OBJ_DIR)/rftest.F90

$(OBJ_DIR)/dummy_orbit.o : $(OBJ_DIR)/dummy_orbit.F90
	$(COMPILER) -o $(OBJ_DIR)/dummy_orbit.o \
	$(OBJ_DIR)/dummy_orbit.F90

$(OBJ_DIR)/RFOF_numerics.o : $(OBJ_DIR)/RFOF_numerics.F90
	$(COMPILER) -o $(OBJ_DIR)/RFOF_numerics.o \
	$(OBJ_DIR)/RFOF_numerics.F90   

$(OBJ_DIR)/RFOF_types.o : $(OBJ_DIR)/RFOF_types.F90
	$(COMPILER) -o $(OBJ_DIR)/RFOF_types.o \
	$(OBJ_DIR)/RFOF_types.F90    

$(OBJ_DIR)/RFOF_parameters.o : $(OBJ_DIR)/RFOF_parameters.F90
	$(COMPILER) -o $(OBJ_DIR)/RFOF_parameters.o \
	$(OBJ_DIR)/RFOF_parameters.F90    

$(OBJ_DIR)/RFOF_waves.o : $(OBJ_DIR)/RFOF_waves.F90
	$(COMPILER) -o $(OBJ_DIR)/RFOF_waves.o \
	$(OBJ_DIR)/RFOF_waves.F90    

$(OBJ_DIR)/RFOF_magnetic_field.o : $(OBJ_DIR)/RFOF_magnetic_field.F90
	$(COMPILER) -o $(OBJ_DIR)/RFOF_magnetic_field.o \
	$(OBJ_DIR)/RFOF_magnetic_field.F90    

$(OBJ_DIR)/RFOF_markers.o : $(OBJ_DIR)/RFOF_markers.F90
	$(COMPILER) -o $(OBJ_DIR)/RFOF_markers.o \
	$(OBJ_DIR)/RFOF_markers.F90

$(OBJ_DIR)/RFOF_resonance_memory.o : $(OBJ_DIR)/RFOF_resonance_memory.F90
	$(COMPILER) -o $(OBJ_DIR)/RFOF_resonance_memory.o \
	$(OBJ_DIR)/RFOF_resonance_memory.F90

$(OBJ_DIR)/RFOF_resonance_condition.o : $(OBJ_DIR)/RFOF_resonance_condition.F90
	$(COMPILER) -o $(OBJ_DIR)/RFOF_resonance_condition.o \
	$(OBJ_DIR)/RFOF_resonance_condition.F90

$(OBJ_DIR)/RFOF_kick.o : $(OBJ_DIR)/RFOF_kick.F90
	$(COMPILER) -o $(OBJ_DIR)/RFOF_kick.o \
	$(OBJ_DIR)/RFOF_kick.F90

$(OBJ_DIR)/RFOF_local_magnetic_field.o : $(OBJ_DIR)/RFOF_local_magnetic_field.F90
	$(COMPILER) -o $(OBJ_DIR)/RFOF_local_magnetic_field.o \
	$(OBJ_DIR)/RFOF_local_magnetic_field.F90

$(OBJ_DIR)/RFOF_main.o : $(OBJ_DIR)/RFOF_main.F90
	$(COMPILER) -o $(OBJ_DIR)/RFOF_main.o \
	$(OBJ_DIR)/RFOF_main.F90

$(OBJ_DIR)/RFOF_constants.o : $(OBJ_DIR)/RFOF_constants.F90
	$(COMPILER) -o $(OBJ_DIR)/RFOF_constants.o \
	$(OBJ_DIR)/RFOF_constants.F90

$(OBJ_DIR)/dum_magnetic_field.o : $(OBJ_DIR)/dum_magnetic_field.F90
	$(COMPILER) -o $(OBJ_DIR)/dum_magnetic_field.o \
	$(OBJ_DIR)/dum_magnetic_field.F90

$(OBJ_DIR)/RFOF_diagnostics.o : $(OBJ_DIR)/RFOF_diagnostics.F90
	$(COMPILER) -o $(OBJ_DIR)/RFOF_diagnostics.o \
	$(OBJ_DIR)/RFOF_diagnostics.F90

$(OBJ_DIR)/euitm_waves_interface.o : $(OBJ_DIR)/euitm_waves_interface.f90
	$(COMPILER) -o $(OBJ_DIR)/euitm_waves_interface.o \
	$(OBJ_DIR)/euitm_waves_interface.f90

$(OBJ_DIR)/RFOF_Efield_update.o : $(OBJ_DIR)/RFOF_Efield_update.F90
	$(COMPILER) -o $(OBJ_DIR)/RFOF_Efield_update.o \
	$(OBJ_DIR)/RFOF_Efield_update.F90

$(OBJ_DIR)/steinbrecher_template.o : $(OBJ_DIR)/steinbrecher_template.F90
	$(COMPILER) -o $(OBJ_DIR)/steinbrecher_template.o \
	$(OBJ_DIR)/steinbrecher_template.F90

$(OBJ_DIR)/RFOF_random_numbers.o : $(OBJ_DIR)/RFOF_random_numbers.F90
	$(COMPILER) -o $(OBJ_DIR)/RFOF_random_numbers.o \
	$(OBJ_DIR)/RFOF_random_numbers.F90

$(OBJ_DIR)/RFOF_mpi_serial.o : $(OBJ_DIR)/RFOF_mpi_serial.F90
	$(COMPILER) -o $(OBJ_DIR)/RFOF_mpi_serial.o \
	$(OBJ_DIR)/RFOF_mpi_serial.F90

$(OBJ_DIR)/GlobalParam_mod.o : $(OBJ_DIR)/GlobalParam.f90
	$(COMPILER) -o $(OBJ_DIR)/GlobalParam_mod.o \
	$(OBJ_DIR)/GlobalParam.f90

$(OBJ_DIR)/coherentwave_mod.o : $(OBJ_DIR)/coherentwave.f90
	$(COMPILER) -o $(OBJ_DIR)/coherentwave_mod.o \
	$(OBJ_DIR)/coherentwave.f90

$(OBJ_DIR)/euitm_schemas.o : $(OBJ_DIR)/euitm_schemas.f90
	$(COMPILER) -o $(OBJ_DIR)/euitm_schemas.o \
	$(OBJ_DIR)/euitm_schemas.f90

$(OBJ_DIR)/Fullwave_mod.o : $(OBJ_DIR)/Fullwave.f90
	$(COMPILER) -o $(OBJ_DIR)/Fullwave_mod.o \
	$(OBJ_DIR)/Fullwave.f90

$(OBJ_DIR)/Local.o : $(OBJ_DIR)/Local.f90
	$(COMPILER) -o $(OBJ_DIR)/Local.o \
	$(OBJ_DIR)/Local.f90

$(OBJ_DIR)/magnetic_field_interface_to_RFOF.o : $(OBJ_DIR)/magnetic_field_interface_to_RFOF.F90
	$(COMPILER) -o $(OBJ_DIR)/magnetic_field_interface_to_RFOF.o \
	$(OBJ_DIR)/magnetic_field_interface_to_RFOF.F90

#$(OBJ_DIR)/type_waves.o : $(OBJ_DIR)/type_waves.F90
#	$(COMPILER) -o $(OBJ_DIR)/type_waves.o \
	$(OBJ_DIR)/type_waves.F90

$(OBJ_DIR)/RFOF_BesselProcess0.o : $(OBJ_DIR)/RFOF_BesselProcess0.F90
	$(COMPILER) -o $(OBJ_DIR)/RFOF_BesselProcess0.o \
	$(OBJ_DIR)/RFOF_BesselProcess0.F90

$(OBJ_DIR)/RFOF_Bessel0_5_FET_TestMain.o : $(OBJ_DIR)/RFOF_Bessel0_5_FET_TestMain.F90
	$(COMPILER) -o $(OBJ_DIR)/RFOF_Bessel0_5_FET_TestMain.o \
	$(OBJ_DIR)/RFOF_Bessel0_5_FET_TestMain.F90

$(OBJ_DIR)/RFOF_BesselProcess0_3.o : $(OBJ_DIR)/RFOF_BesselProcess0_3.F90
	$(COMPILER) -o $(OBJ_DIR)/RFOF_BesselProcess0_3.o \
	$(OBJ_DIR)/RFOF_BesselProcess0_3.F90

$(OBJ_DIR)/RFOF_BesselProcess0_5.o : $(OBJ_DIR)/RFOF_BesselProcess0_5.F90
	$(COMPILER) -o $(OBJ_DIR)/RFOF_BesselProcess0_5.o \
	$(OBJ_DIR)/RFOF_BesselProcess0_5.F90

$(OBJ_DIR)/RFOF_random_numbers_SG2.o : $(OBJ_DIR)/RFOF_random_numbers_SG2.F90
	$(COMPILER) -o $(OBJ_DIR)/RFOF_random_numbers_SG2.o \
	$(OBJ_DIR)/RFOF_random_numbers_SG2.F90

$(OBJ_DIR)/RFOF_Stochastics.o : $(OBJ_DIR)/RFOF_Stochastics.F90
	$(COMPILER) -o $(OBJ_DIR)/RFOF_Stochastics.o \
	$(OBJ_DIR)/RFOF_Stochastics.F90

make clean:
	rm $(OBJ_DIR)/*.o
	rm $(OBJ_DIR)/*.mod
#rm $(EXEC)
	
