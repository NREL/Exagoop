AMREX_HOME ?= ../../amrex
MPM_HOME = ..

EBASE = mpm

DEBUG	= FALSE

DIM	= 3

COMP    = gnu

USE_MPI   = TRUE 
USE_OMP   = FALSE
USE_CUDA  = FALSE
USE_EB    = TRUE
DO_POST   = FALSE

TINY_PROFILE = FALSE
USE_PARTICLES = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Particle/Make.package
include $(AMREX_HOME)/Src/EB/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package

include $(MPM_HOME)/Source/Make.package

INCLUDE_LOCATIONS += $(MPM_HOME)/Source
VPATH_LOCATIONS   += $(MPM_HOME)/Source

ifeq ($(DO_POST), TRUE)
  HYDE_HOME ?= /Users/mday2/src/Hyde
  INCLUDE_LOCATIONS += $(HYDE_HOME)/tgsl-share/library $(HYDE_HOME)/JUNK/partio/src/lib
  LIBRARY_LOCATIONS += $(HYDE_HOME)/JUNK/partio/build/src/lib $(HOME)/.brew/lib
  LIBRARIES += -lpartio -ltbb
endif

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
