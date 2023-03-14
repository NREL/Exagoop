AMREX_HOME ?= ../../amrex
MPM_HOME = /projects/fcicfppm/ndeak/Exagoop/

EBASE = mpm

DEBUG	= FALSE

DIM	= 3

COMP    = gnu

USE_MPI   = TRUE 
USE_OMP   = FALSE
USE_CUDA  = FALSE
USE_EB    = TRUE

TINY_PROFILE = FALSE
USE_PARTICLES = TRUE
USE_SUNDIALS = TRUE

ifeq ($(USE_SUNDIALS),TRUE)

DEFINES += -DAMREX_USE_SUNDIALS

SUNDIALS_ROOT ?= ../ThirdParty/sundials/$(COMP)

REL_SUNDIALS_LIB_DIR := $(MPM_HOME)/ThirdParty/sundials/$(COMP)/lib
SUNDIALS_LIB_DIR := $(abspath $(lastword $(REL_SUNDIALS_LIB_DIR)))
INCLUDE_LOCATIONS += $(SUNDIALS_LIB_DIR)/../include
LIBRARIES += -L$(SUNDIALS_LIB_DIR) -lsundials_arkode
LIBRARIES += -L$(SUNDIALS_LIB_DIR) -lsundials_nvecmanyvector
LIBRARIES += -Wl,-rpath=${SUNDIALS_LIB_DIR}


# USE_CVODE_LIBS ?= TRUE
# USE_ARKODE_LIBS ?= TRUE

# INCLUDE_LOCATIONS += $(SUNDIALS_ROOT)/include
# LIBRARY_LOCATIONS += $(SUNDIALS_LIB_DIR)

# LIBRARIES += -L$(SUNDIALS_LIB_DIR) -lsundials_cvode
# LIBRARIES += -L$(SUNDIALS_LIB_DIR) -lsundials_arkode
# LIBRARIES += -L$(SUNDIALS_LIB_DIR) -lsundials_nvecmanyvector
# LIBRARIES += -Wl,-rpath=$(SUNDIALS_LIB_DIR)
endif

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Particle/Make.package
include $(AMREX_HOME)/Src/EB/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package

ifeq ($(USE_SUNDIALS),TRUE)
include $(AMREX_HOME)/Src/Extern/SUNDIALS/Make.package
endif

include $(MPM_HOME)/Source/Make.package

INCLUDE_LOCATIONS += $(MPM_HOME)/Source
VPATH_LOCATIONS   += $(MPM_HOME)/Source

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
