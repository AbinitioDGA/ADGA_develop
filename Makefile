F90 = mpiifort
FFLAGS += -O3 -g -fpp -DMPI -traceback #-check all
FINCLUDE += -I/opt/sw/x86_64/glibc-2.12/ivybridge-ep/hdf5/1.8.12/intel-14.0.2/include
LD = $(F90)
LDFLAGS += -lhdf5_fortran -lhdf5hl_fortran -lmkl_rt  #-limf
LDINCLUDE += -L/opt/sw/x86_64/glibc-2.12/ivybridge-ep/hdf5/1.8.12/intel-14.0.2/lib -L/cm/shared/apps/intel/composer_xe_2015.2.164/mkl/lib/intel64

# get source files
F90_MAIN_SOURCES := main.f90 parameters_module.f90 vq_module.f90 kq_tools.f90 lapack_module.f90 one_particle_quant_module.f90 \
		    susc_module.f90 eom_module.f90 hdf5_module.f90 aux.f90 mpi_org.f90
F90_VERTEX_SOURCES := vertex_chann_sym.f90 hdf5_module.f90 parameters_module.f90 aux.f90
#F90_VQ_SOURCES := vq_check.f90 vq_module.f90 parameters_module.f90 kq_tools.f90 hdf5_module.f90
MAIN_OBJECTS := $(patsubst %.f90,%.o,$(F90_MAIN_SOURCES))
VERTEX_OBJECTS := $(patsubst %.f90,%.o,$(F90_VERTEX_SOURCES))
VQ_OBJECTS := $(patsubst %.f90,%.o,$(F90_VQ_SOURCES))

# define all objects
ALL_OBJECTS := $(MAIN_OBJECTS)
ALL_OBJECTS += $(VERTEX_OBJECTS)
ALL_OBJECTS += $(VQ_OBJECTS)

# remove duplicates
OBJECTS := $(sort $(ALL_OBJECTS))


.PHONY: all
all: main vertex_chann_sym
#vq_check

main: $(MAIN_OBJECTS)
	$(LD) $^ -o $@ $(FFLAGS) $(LDINCLUDE) $(LDFLAGS)

vertex_chann_sym: $(VERTEX_OBJECTS)
	$(LD) $^ -o $@ $(FFLAGS) $(LDINCLUDE) $(LDFLAGS)

#vq_check: $(VQ_OBJECTS)
#	$(LD) $(VQ_OBJECTS) -o $@ $(FFLAGS) $(LDINCLUDE) $(LDFLAGS)

main.o: parameters_module.o aux.o kq_tools.o lapack_module.o one_particle_quant_module.o susc_module.o \
				eom_module.o vq_module.o hdf5_module.o mpi_org.o
vertex_chann_sym.o: hdf5_module.o parameters_module.f90 aux.o
#vq_check.o: parameters_module.o vq_module.o kq_tools.o
eom_module.o: parameters_module.o one_particle_quant_module.o
susc_module.o: parameters_module.o hdf5_module.o
one_particle_quant_module.o: aux.o lapack_module.o parameters_module.o
kq_tools.o: parameters_module.o
vq_module.o: parameters_module.o hdf5_module.o aux.o
hdf5_module.o: parameters_module.o aux.o
mpi_org.o: parameters_module.o 

$(OBJECTS): %.o: %.f90
	$(LD) $(FFLAGS) $(FINCLUDE) -c $< -o $@

.PHONY: clean
clean:
	rm -f *.o *.mod main vertex_chann_sym
