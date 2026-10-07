# Shared dependency wiring for VILMA2.
#
# Loaded by config/Makefile *after* the compiler fragment, so it may reference
# FFLAGS_BASE, MODFLAGS, CPPFLAGS_PP, INC_NC and LIB_NC defined there.
#
# All numerical dependencies (FFTW, SHTns) and the fesm-utils helper library are
# provided by the fesm-utils package, expected at the repo root as a symlink:
#
#     ln -s ../fesm-utils fesm-utils
#
# vilma_remap (conservative lon-lat -> Gauss remapping for the standalone driver) uses
# the `coords` module, which lives on the fesm-utils `coords-dev` branch. Point the
# symlink at a checkout on that branch and build its utils library:
#     (in the fesm-utils checkout)  configme config && make fesmutils-static
#
# SHTns (spherical-harmonic transforms) is built into fesm-utils with:
#     cd fesm-utils && make shtns openmp=0   # (openmp=1 for the omp variant)
# It links FFTW, which fesm-utils also builds. The OpenMP dependency variants
# are used by default (openmp=1); `make openmp=0` swaps in the serial variants.
# The swap is driven by the openmp= switch in the OpenMP section below.

# --- fesm-utils helper library (ncio, nml, mapping_scrip, ...) ---------------
# FESMUTILS_LIB is the on-disk libfesmutils.a (what `fesmutils-static` produces).
# Objects that `use` a fesm-utils module (coords) depend on it so a rebuilt
# fesm-utils (changed module interface) forces them to recompile against the new
# .mod, rather than silently relinking a stale object (ABI mismatch -> crash).
FESMUTILSROOT = fesm-utils
INC_FESMUTILS = -I$(FESMUTILSROOT)/include-serial
LIB_FESMUTILS = -L$(FESMUTILSROOT)/include-serial -lfesmutils
FESMUTILS_LIB = $(FESMUTILSROOT)/include-serial/libfesmutils.a

# --- FFTW --------------------------------------------------------------------
FFTWROOT = fesm-utils/fftw/fftw-serial
INC_FFTW = -I$(FFTWROOT)/include
LIB_FFTW = -L$(FFTWROOT)/lib -lfftw3 -lm

# --- SHTns (provides shtns.f03 Fortran 2003 interface + libshtns.a) ----------
SHTNSROOT = fesm-utils/SHTns/shtns-serial
INC_SHTNS = -I$(SHTNSROOT)/include
LIB_SHTNS = -L$(SHTNSROOT)/lib -lshtns

# --- (LIS removed) -----------------------------------------------------------
# The per-degree solve is now a dependency-free pivoted banded LU (vilma_band); LIS
# is no longer linked. Keeping INC_LIS / LIB_LIS empty so the flag lists below
# (and any external references) stay valid.
INC_LIS =
LIB_LIS =

# --- OpenMP build (make openmp=1) --------------------------------------------
# Two things happen for an OpenMP build:
#   1. The serial dependency builds above are swapped for their OpenMP variants
#      (done here): include-omp for fesm-utils, fftw-omp (-lfftw3_omp -lfftw3),
#      and shtns-omp (libshtns_omp.a, i.e. -lshtns_omp).
#   2. The compiler's OpenMP flag (-fopenmp / FFLAGS_OPENMP) is appended to
#      FFLAGS by config/Makefile, which also threads the per-degree loop in
#      vilma_response (begin_step / commit_step) over independent per-degree systems,
#      each solved by the re-entrant banded LU (vilma_band).
# (There is no LIS variant to reconcile: the iterative solver was removed in
# favour of the direct banded LU.)
ifeq ($(openmp),1)
	INC_FESMUTILS = -I$(FESMUTILSROOT)/include-omp
	LIB_FESMUTILS = -L$(FESMUTILSROOT)/include-omp -lfesmutils
	FESMUTILS_LIB = $(FESMUTILSROOT)/include-omp/libfesmutils.a

	FFTWROOT = fesm-utils/fftw/fftw-omp
	INC_FFTW = -I$(FFTWROOT)/include
	LIB_FFTW = -L$(FFTWROOT)/lib -lfftw3_omp -lfftw3 -lm

	SHTNSROOT = fesm-utils/SHTns/shtns-omp
	INC_SHTNS = -I$(SHTNSROOT)/include
	LIB_SHTNS = -L$(SHTNSROOT)/lib -lshtns_omp
endif

# --- VILMA1 backend (make vilma_v1=1 VILMA_V1_ROOT=...) -------------------------------
# OPTIONAL and OFF by default. VILMA1 (Martinec/Klemann; the CLIMBER-X i_geo=2
# backend) is a hand-installed, precompiled library: a `vega_pism.a` archive plus
# a directory of `.mod` files. It is absent on most machines, so it must never
# become a dependency of VILMA2.
#
# vilma_v1=0 (the default): CPPFLAGS_VILMA_V1 / INC_VILMA_V1 / LIB_VILMA_V1 are all EMPTY, so
#   the compile line carries no -DVILMA_V1 and no VILMA_V1_ROOT include, and the link line
#   no archive. src/vilma_v1.f90 then compiles to a pure-Fortran stub referencing
#   no VILMA1 symbol, which aborts with an actionable message if solver="v1" is
#   selected at runtime. The build is identical to a tree without this switch.
# vilma_v1=1: -DVILMA_V1 activates the real wrapper; VILMA_V1_ROOT must point at an install
#   containing include/*.mod and lib/vega_pism.a.
#
# Mirrors the vilma= toggle in CLIMBER-X's config/common.mk.
VILMA_V1_ROOT ?= vilma1
CPPFLAGS_VILMA_V1 =
INC_VILMA_V1 =
LIB_VILMA_V1 =
ifeq ($(vilma_v1),1)
	CPPFLAGS_VILMA_V1 = -DVILMA_V1
	INC_VILMA_V1      = -I$(VILMA_V1_ROOT)/include
	LIB_VILMA_V1      = $(VILMA_V1_ROOT)/lib/vega_pism.a
endif

# --- Final flag sets ---------------------------------------------------------
# MODFLAGS (-I/-J objdir) and FFLAGS_BASE come from the compiler fragment.
# INC_SHTNS is what lets `include 'shtns.f03'` in src/vilma_sht.f90 be found.
CPPFLAGS_VILMA = $(CPPFLAGS_PP) $(CPPFLAGS_VILMA_V1)
FFLAGS_VILMA   = $(FFLAGS_BASE) $(MODFLAGS) $(INC_NC) $(INC_FESMUTILS) $(INC_FFTW) $(INC_SHTNS) $(INC_LIS) $(INC_VILMA_V1)

# Static archives resolve left-to-right, so a library must precede the libraries
# it depends on: SHTns before FFTW (SHTns calls FFTW), fesm-utils before netCDF.
# LIB_VILMA_V1 is empty unless vilma_v1=1; VILMA1 calls netCDF, so it precedes LIB_NC.
LFLAGS_VILMA   = $(LIB_FESMUTILS) $(LIB_SHTNS) $(LIB_FFTW) $(LIB_LIS) $(LIB_VILMA_V1) $(LIB_NC) $(LFLAGS_EXTRA)
