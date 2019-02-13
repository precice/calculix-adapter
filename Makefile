# Specify the locations of: the original CCX source, SPOOLES and ARPACK
CCX             = $(HOME)/Software/Calculix/Calculix/2.15/ccx_2.15/CalculiX/ccx_2.15/src
SPOOLES         = $(HOME)/Software/Calculix/SPOOLES.2.2
ARPACK          = $(HOME)/Software/Calculix/ARPACK
PRECICE_ROOT    = $(HOME)/Software/preCICE
YAML            = $(HOME)/Software/yaml-cpp

# Specify where to store the generated .o files
OBJDIR 		= bin

# Includes and libs
INCLUDES = \
	-I./ \
	-I./adapter \
	-I$(CCX) \
	-I$(SPOOLES) \
	-I$(PRECICE_ROOT)/src \
	-I$(YAML)/include

LIBS = \
	$(SPOOLES)/spooles.a \
    	-L$(PRECICE_ROOT)/build/last -lprecice \
    	-lboost_log \
    	-lboost_log_setup \
    	-lboost_thread \
    	-lboost_system \
    	-lboost_filesystem \
	-lboost_program_options \
    	-lpython2.7 \
    	-lstdc++ \
	-L/home/daviske/Software/petsc/x86_64/lib -lpetsc \
    	-lmpi_cxx \
    	-L$(YAML)/build -lyaml-cpp \
        -lxml2

# OS-specific options
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	LIBS += $(ARPACK)/libarpack_MAC.a
else
	LIBS += $(ARPACK)/libarpack_INTEL.a
    	LIBS += -lpthread -lm -lc
endif

# Compilers and flags
#CFLAGS = -g -Wall -std=c++11 -O0 -fopenmp $(INCLUDES) -DARCH="Linux" -DSPOOLES -DARPACK -DMATRIXSTORAGE
#FFLAGS = -g -Wall -O0 -fopenmp $(INCLUDES)

CFLAGS = -Wall -O3 -fopenmp $(INCLUDES) -DARCH="Linux" -DSPOOLES -DARPACK -DMATRIXSTORAGE

# OS-specific options
ifeq ($(UNAME_S),Darwin)
        CC = /usr/local/bin/gcc
else
        CC = mpicc
endif

FFLAGS = -Wall -O3 -fopenmp $(INCLUDES)
#FC = mpifort
#FC = mpif90
#FC = gfortran
F77 = f77

# Include a list of all the source files
include $(CCX)/Makefile.inc
SCCXMAIN = ccx_2.15.c

# Append additional sources
SCCXC += nonlingeo_precice.c CCXHelpers.c PreciceInterface.c
SCCXF += getflux.f getkdeltatemp.f



# Source files in this folder and in the adapter directory
$(OBJDIR)/%.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@
$(OBJDIR)/%.o : %.f
	$(FC) $(FFLAGS) -c $< -o $@
$(OBJDIR)/%.o : adapter/%.c
	$(CC) $(CFLAGS) -c $< -o $@
$(OBJDIR)/%.o : adapter/%.cpp
	g++ -std=c++11 -I$(YAML)/include -c $< -o $@ $(LIBS)
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
OCCXC += $(OBJDIR)/ConfigReader.o



$(OBJDIR)/ccx_preCICE: $(OBJDIR) $(OCCXMAIN) $(OBJDIR)/ccx_2.15.a
	$(FC) -fopenmp -Wall -O3 -o $@ $(OCCXMAIN) $(OBJDIR)/ccx_2.15.a $(LIBS)

$(OBJDIR)/ccx_2.15.a: $(OCCXF) $(OCCXC)
	ar vr $@ $?

$(OBJDIR):
	mkdir -p $(OBJDIR)

clean:
	rm -f $(OBJDIR)/*.o $(OBJDIR)/ccx_2.15.a $(OBJDIR)/ccx_preCICE
