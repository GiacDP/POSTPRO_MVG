rm *.o *.mod *.exe
mpif90 -g -c module_postpro_subvolume.f90 main_addmom_subvolume_serial_stats.f90
mpif90 module_postpro_subvolume.o main_addmom_subvolume_serial_stats.o -o addmom_parallel_subvolume_serial_stats.exe
