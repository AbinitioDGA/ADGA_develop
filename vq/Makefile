VPROG=readV
INSTALL_DIR=../bin
MAINDIR=../src

include ../make_config 

FINCLUDE += -I$(MAINDIR)
OBJECTS = $(MAINDIR)/parameters_module.o $(MAINDIR)/aux.o $(MAINDIR)/hdf5_module.o readV_UJ_nonloc.o

all: $(INSTALL_DIR)/$(VPROG)

$(INSTALL_DIR)/$(VPROG): $(OBJECTS)
	$(LD) $(OBJECTS) -o $@ $(FFLAGS) $(LDINCLUDE) $(LDFLAGS)

readV_UJ_nonloc.o: $(MAINDIR)/hdf5_module.o
$(MAINDIR)/hdf5_module.o: $(MAINDIR)/parameters_module.o $(MAINDIR)/aux.o

%.o: %.f90
	$(F90) $(FFLAGS) $(FINCLUDE) -c $< -o $@

%.o: %.f
	$(F90) $(FFLAGS) $(FINCLUDE) -c $< -o $@

.PHONY: clean pristine

clean:
	rm -f *.o *.mod

pristine:
	rm -f *.o *.mod $(INSTALL_DIR)/$(VPROG)

