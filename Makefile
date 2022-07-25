# See the adapter documentation for getting the adapter dependencies:
# https://precice.org/adapter-calculix-get-calculix.html
# Set the following variables before building:
# Path to original CalculiX source (e.g. $(HOME)/ccx_2.xx/src )
CCX_VERSION			= 2.19
CCX             = $(HOME)/CalculiX/ccx_$(CCX_VERSION)/src

### Change these if you built SPOOLES, ARPACK, or yaml-cpp from source ###
# SPOOLES include flags (e.g. -I$(HOME)/SPOOLES.2.2 )
SPOOLES_INCLUDE   = -I/usr/include/spooles/
# SPOOLES library flags (e.g. $(HOME)/SPOOLES.2.2/spooles.a)
SPOOLES_LIBS      = -lspooles
#
# ARPACK include flags (e.g. -I$(HOME)/ARPACK)
ARPACK_INCLUDE    =
# ARPACK library flags (e.g. $(HOME)/ARPACK/libarpack_INTEL.a)
ARPACK_LIBS       = -larpack -llapack -lblas
#
# yaml-cpp include flags (e.g. -I$(HOME)/yaml-cpp/include)
YAML_INCLUDE      = -I/usr/include/
# yaml-cpp library flags (e.g. -L$(HOME)/yaml-cpp/build -lyaml-cpp)
YAML_LIBS         = -lyaml-cpp

# Get the CFLAGS and LIBS from pkg-config (preCICE version >= 1.4.0).
# If pkg-config cannot find the libprecice.pc meta-file, you may need to set the
# path where this is stored into PKG_CONFIG_PATH when building the adapter.
PKGCONF_CFLAGS  = $(shell pkg-config --cflags libprecice)
PKGCONF_LIBS    = $(shell pkg-config --libs libprecice)

# Specify where to store the generated .o files
OBJDIR = bin

# Includes and libs
INCLUDES = \
	-I./ \
	-I./adapter \
	-I$(CCX) \
	$(SPOOLES_INCLUDE) \
	$(PKGCONF_CFLAGS) \
	$(ARPACK_INCLUDE) \
	$(YAML_INCLUDE)

LIBS = \
	$(SPOOLES_LIBS) \
	$(PKGCONF_LIBS) \
	-lstdc++ \
	$(YAML_LIBS) \
	$(ARPACK_LIBS) \
	-lpthread -lm -lc

# Compilers and flags
#CFLAGS = -g -Wall -std=c++11 -O0 -fopenmp $(INCLUDES) -DARCH="Linux" -DSPOOLES -DARPACK -DMATRIXSTORAGE
#FFLAGS = -g -Wall -O0 -fopenmp $(INCLUDES)

CFLAGS = -Wall -O3 -fopenmp $(INCLUDES) -DARCH="Linux" -DSPOOLES -DARPACK -DMATRIXSTORAGE -DUSE_MT

# OS-specific options
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	CC = /usr/local/bin/gcc
else
	CC = mpicc
endif

FFLAGS = -Wall -O3 -fopenmp $(INCLUDES) ${ADDITIONAL_FFLAGS}
# Note for GCC 10 or newer: add -fallow-argument-mismatch in the above flags
FC = mpifort
# FC = mpif90
# FC = gfortran

# Include a list of all the source files
include $(CCX)/Makefile.inc
SCCXMAIN = ccx_$(CCX_VERSION).c

# Append additional sources
SCCXC += nonlingeo_precice.c dyna_precice.c CCXHelpers.c PreciceInterface.c
SCCXF += getflux.f getkdeltatemp.f



# Source files in this folder and in the adapter directory
$(OBJDIR)/%.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@
$(OBJDIR)/%.o : %.f
	$(FC) $(FFLAGS) -c $< -o $@
$(OBJDIR)/%.o : adapter/%.c
	$(CC) $(CFLAGS) -c $< -o $@
$(OBJDIR)/%.o : adapter/%.cpp
	g++ -std=c++11 $(YAML_INCLUDE) -c $< -o $@ $(LIBS)
	#$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@ $(LIBS)

# Source files in the $(CCX) folder
$(OBJDIR)/%.o : $(CCX)/%.c
	$(CC) $(CFLAGS) -c $< -o $@
$(OBJDIR)/%.o : $(CCX)/%.f
	$(FC) $(FFLAGS) -c $< -o $@

# Generate list of object files from the source files, prepend $(OBJDIR)
OCCXF = $(SCCXF:%.f=$(OBJDIR)/%.o)
OCCXC = $(SCCXC:%.c=$(OBJDIR)/%.o)
OCCXMAIN = $(SCCXMAIN:%.c=$(OBJDIR)/%.o)
OCCXC += $(OBJDIR)/ConfigReader.o $(OBJDIR)/2D3DCoupling.o



$(OBJDIR)/ccx_preCICE: $(OBJDIR) $(OCCXMAIN) $(OBJDIR)/ccx_$(CCX_VERSION).a
	$(FC) -fopenmp -Wall -O3 -o $@ $(OCCXMAIN) $(OBJDIR)/ccx_$(CCX_VERSION).a $(LIBS)

$(OBJDIR)/ccx_$(CCX_VERSION).a: $(OCCXF) $(OCCXC)
	ar vr $@ $?

$(OBJDIR):
	mkdir -p $(OBJDIR)

clean:
	rm -f $(OBJDIR)/*.o $(OBJDIR)/ccx_$(CCX_VERSION).a $(OBJDIR)/ccx_preCICE
