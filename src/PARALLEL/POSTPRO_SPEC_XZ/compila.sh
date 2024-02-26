rm *.mod *.o wallspec.exe
mpif90 -fdefault-real-8 -fdefault-double-8 -c modfftw.f90  -I/usr/local/include -O2
mpif90 -fdefault-real-8 -fdefault-double-8 -c welch.f90    -I/usr/local/include -O2
mpif90 -fdefault-real-8 -fdefault-double-8 -c main_wallspec.f90 -O2 #-fconvert=big-endian
mpif90 -fdefault-real-8 -fdefault-double-8 -o wallspec.exe *.o -O2 -L/usr/local/lib/ -lfftw3 #-fconvert=big-endian
