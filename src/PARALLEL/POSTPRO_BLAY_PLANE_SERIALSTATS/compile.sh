rm *.o *.mod *.exe
#mpif90 -C -fbacktrace -g -c module_postpro_subvolume.f90 main_postpro_blay_serial_stats.f90
mpif90 -g -c module_postpro_subvolume.f90 main_postpro_blay_serial_stats.f90
mpif90 module_postpro_subvolume.o main_postpro_blay_serial_stats.o -o postpro_blay_props_serial_stats.exe
