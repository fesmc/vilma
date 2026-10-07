.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# === Directory layout ========================================================
objdir  = obj
srcdir  = src
testdir = tests
bindir  = bin

# === Command-line switches ===================================================
debug  ?= 0
openmp ?= 1
# Optional VILMA1 solid-earth backend (solver="v1"). OFF by default: VILMA1 is a
# hand-installed, precompiled library absent on most machines, so it must never
# become a build dependency of VILMA2. With vilma_v1=0 (the default) not one
# VILMA1 symbol is referenced, no VILMA1 include path is needed, and a namelist
# asking for solver="v1" aborts with an actionable message from the
# compiled-in stub (src/vilma_v1.f90, gated on -DVILMA_V1). See doc/vilma-v1-backend.md.
#   make vilma vilma_v1=1 VILMA_V1_ROOT=/path/to/vilma1
vilma_v1  ?= 0

defaulttarget: usage

# === Compiler configuration ==================================================
# Inserted by `configme` (or the legacy config/legacy/config.py). The fragment
# defines the compiler (FC), flag sets (FFLAGS_BASE, FFLAGS_OPENMP, MODFLAGS,
# CPPFLAGS_PP, DFLAGS_*), the netCDF wiring (INC_NC, LIB_NC), and any
# host-specific link extras (LFLAGS_EXTRA).

<COMPILER_CONFIGURATION>

# === Shared dependency wiring ================================================
# Pulled in after the compiler fragment so it can reference the variables above.
# Defines the fesm-utils / FFTW / SHTns paths and the final FFLAGS_VILMA / LFLAGS_VILMA
# / CPPFLAGS_VILMA flag sets, honouring the openmp= switch.
include config/common.mk

# === Git version string ======================================================
GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always --tags 2>/dev/null || echo "unknown")
CPPFLAGS_VILMA := $(CPPFLAGS_VILMA) -DVERSION=\"$(GIT_VERSION)\"

# === Debug / profile selection ===============================================
DFLAGS = $(DFLAGS_NODEBUG)
ifeq ($(debug),1)
	DFLAGS = $(DFLAGS_DEBUG)
endif
ifeq ($(debug),2)
	DFLAGS = $(DFLAGS_PROFILE)
endif

# === Assemble the flags actually used by the rules ===========================
FFLAGS   = $(FFLAGS_VILMA)
LFLAGS   = $(LFLAGS_VILMA)
CPPFLAGS = $(CPPFLAGS_VILMA)
ifeq ($(openmp),1)
	FFLAGS := $(FFLAGS) $(FFLAGS_OPENMP)
endif

# === Source/rule lists and build targets =====================================
include config/Makefile_vilma.mk
