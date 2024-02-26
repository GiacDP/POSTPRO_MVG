rm *.o *.mod *.exe
mpif90 \-g \-c module_postpro_subvolume.f90 main_vel_tau_2dstat.f90
mpif90 \-g module_postpro_subvolume.o main_vel_tau_2dstat.o -o postpro_vel_tau_2dstat.exe
