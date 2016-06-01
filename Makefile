
#FFLAGS += -I/opt/sw/x86_64/glibc-2.12/ivybridge-ep/hdf5/1.8.12/intel-14.0.2/include

# set flags
F90 = mpiifort
#F90 = mpif90
FFLAGS += -O3 -g -fpp -DMPI
#FFLAGS += -O0 -g -fbacktrace -C -DMPI
FINCLUDE += -I/opt/sw/x86_64/glibc-2.12/ivybridge-ep/hdf5/1.8.12/intel-14.0.2/include
#FINCLUDE += -I/opt/sw/x86_64/glibc-2.12/ivybridge-ep/hdf5/1.8.14/gnu-4.4.7/include
LD = $(F90)
LDFLAGS += -lhdf5_fortran -lhdf5hl_fortran -lmkl_rt  #-limf
LDINCLUDE += -L/opt/sw/x86_64/glibc-2.12/ivybridge-ep/hdf5/1.8.12/intel-14.0.2/lib -L/cm/shared/apps/intel/composer_xe_2015.2.164/mkl/lib/intel64
#LDINCLUDE += -L/opt/sw/x86_64/glibc-2.12/ivybridge-ep/hdf5/1.8.14/gnu-4.4.7/lib -L/cm/shared/apps/intel/composer_xe_2015.2.164/mkl/lib/intel64

# get source files
F90_MAIN_SOURCES := main.f90 parameters_module.f90 index.f90 lapack_module.f90 
F90_SOURCES := vertex_chann_sym.f90 hdf5_module.f90 
MAIN_OBJECTS := $(patsubst %.f90,%.o,$(F90_MAIN_SOURCES))
OBJECTS := $(patsubst %.f90,%.o,$(F90_SOURCES))


all: main vertex_chann_sym

main: $(MAIN_OBJECTS) 
	$(LD) $(MAIN_OBJECTS) -o $@ $(FFLAGS) $(LDINCLUDE) $(LDFLAGS)

vertex_chann_sym: $(OBJECTS) 
	$(LD) $(OBJECTS) -o $@ $(FFLAGS) $(LDINCLUDE) $(LDFLAGS)

main.o: parameters_module.o index.o lapack_module.o
vertex_chann_sym.o: hdf5_module.o

%.o: %.f90
	$(F90) $(FFLAGS) $(FINCLUDE) -c $< -o $@

%.o: %.f
	$(F90) $(FFLAGS) $(FINCLUDE) -c $< -o $@

clean:
	rm -f *.o *.mod $(PROG)
