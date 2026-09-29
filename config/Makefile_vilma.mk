# Source lists, compilation rules and targets for VILMA.
# Included by config/Makefile after the flag sets are assembled.

# --- Library object list (in module-dependency order) ------------------------
obj_vilma = \
	$(objdir)/vilma_precision.o \
	$(objdir)/vilma_constants.o \
	$(objdir)/vilma_params.o \
	$(objdir)/vilma_control.o \
	$(objdir)/vilma_sht.o \
	$(objdir)/vilma_tensor_sh.o \
	$(objdir)/vilma_field.o \
	$(objdir)/vilma_earth_structure.o \
	$(objdir)/vilma_radial_integrals.o \
	$(objdir)/vilma_band.o \
	$(objdir)/vilma_radial_fe.o \
	$(objdir)/vilma_viscoelastic.o \
	$(objdir)/vilma_response.o \
	$(objdir)/vilma_sle.o \
	$(objdir)/vilma_timestep.o \
	$(objdir)/vilma_rotation.o \
	$(objdir)/vilma_remap.o \
	$(objdir)/vilma_v1.o \
	$(objdir)/vilma_coupling.o \
	$(objdir)/vilma_io.o \
	$(objdir)/vilma_drive.o \
	$(objdir)/vilma.o

# --- Inter-module dependencies (so `make -j` stays correct) ------------------
$(objdir)/vilma_constants.o:        $(objdir)/vilma_precision.o
$(objdir)/vilma_params.o:           $(objdir)/vilma_precision.o $(objdir)/vilma_constants.o
$(objdir)/vilma_control.o:          $(objdir)/vilma_precision.o $(objdir)/vilma_constants.o $(objdir)/vilma_params.o
$(objdir)/vilma_sht.o:              $(objdir)/vilma_precision.o
$(objdir)/vilma_tensor_sh.o:        $(objdir)/vilma_precision.o $(objdir)/vilma_sht.o \
                                    $(objdir)/vilma_constants.o
$(objdir)/vilma_field.o:            $(objdir)/vilma_precision.o $(objdir)/vilma_sht.o
$(objdir)/vilma_earth_structure.o:  $(objdir)/vilma_precision.o $(objdir)/vilma_constants.o \
                                    $(objdir)/vilma_params.o $(objdir)/vilma_sht.o
$(objdir)/vilma_radial_integrals.o: $(objdir)/vilma_precision.o
$(objdir)/vilma_band.o:             $(objdir)/vilma_precision.o
$(objdir)/vilma_radial_fe.o:        $(objdir)/vilma_constants.o $(objdir)/vilma_earth_structure.o \
                                    $(objdir)/vilma_radial_integrals.o $(objdir)/vilma_band.o
$(objdir)/vilma_viscoelastic.o:     $(objdir)/vilma_radial_fe.o $(objdir)/vilma_earth_structure.o
$(objdir)/vilma_response.o:         $(objdir)/vilma_radial_fe.o $(objdir)/vilma_earth_structure.o \
                                    $(objdir)/vilma_sht.o $(objdir)/vilma_tensor_sh.o \
                                    $(objdir)/vilma_constants.o $(objdir)/vilma_viscoelastic.o
$(objdir)/vilma_sle.o:              $(objdir)/vilma_sht.o $(objdir)/vilma_constants.o \
                                    $(objdir)/vilma_response.o $(objdir)/vilma_rotation.o
$(objdir)/vilma_timestep.o:         $(objdir)/vilma_response.o $(objdir)/vilma_sle.o \
                                    $(objdir)/vilma_sht.o $(objdir)/vilma_viscoelastic.o \
                                    $(objdir)/vilma_precision.o $(objdir)/vilma_rotation.o
$(objdir)/vilma_rotation.o:         $(objdir)/vilma_sht.o $(objdir)/vilma_constants.o \
                                    $(objdir)/vilma_earth_structure.o $(objdir)/vilma_radial_fe.o \
                                    $(objdir)/vilma_viscoelastic.o
# Optional VILMA-v1 backend. Compiled ALWAYS; with vilma_v1=0 (the default) -DVILMA_V1 is
# absent and this is a pure-Fortran stub that references no VILMA-v1 symbol.
$(objdir)/vilma_v1.o:               $(objdir)/vilma_precision.o $(objdir)/vilma_constants.o \
                                    $(objdir)/vilma_params.o $(objdir)/vilma_sht.o
$(objdir)/vilma_coupling.o:         $(objdir)/vilma_v1.o \
                                    $(objdir)/vilma_response.o $(objdir)/vilma_sle.o \
                                    $(objdir)/vilma_rotation.o $(objdir)/vilma_earth_structure.o \
                                    $(objdir)/vilma_sht.o $(objdir)/vilma_params.o $(objdir)/vilma_remap.o \
                                    $(objdir)/vilma_timestep.o $(objdir)/vilma_viscoelastic.o
$(objdir)/vilma_remap.o:            $(objdir)/vilma_precision.o $(objdir)/vilma_sht.o
$(objdir)/vilma_io.o:               $(objdir)/vilma_coupling.o $(objdir)/vilma_response.o \
                                    $(objdir)/vilma_viscoelastic.o $(objdir)/vilma_sht.o \
                                    $(objdir)/vilma_constants.o
$(objdir)/vilma_drive.o:            $(objdir)/vilma_params.o $(objdir)/vilma_control.o \
                                    $(objdir)/vilma_sht.o \
                                    $(objdir)/vilma_coupling.o $(objdir)/vilma_io.o \
                                    $(objdir)/vilma_remap.o \
                                    $(objdir)/vilma_constants.o $(objdir)/vilma_precision.o
# The umbrella module re-exports every component, so it compiles last.
$(objdir)/vilma.o:                  $(objdir)/vilma_coupling.o $(objdir)/vilma_io.o \
                                    $(objdir)/vilma_drive.o $(objdir)/vilma_params.o \
                                    $(objdir)/vilma_control.o

# --- Pattern rule ------------------------------------------------------------
# Every object depends on FESMUTILS_LIB (libfesmutils.a): fesm-utils modules
# (coords/ncio/nml) are used across VILMA, and a rebuilt fesm-utils with a
# changed module interface must force a recompile against the new .mod rather than
# silently relinking a stale object (ABI mismatch -> segfault, the i_geo=3 coords
# crash). Blanket dep keeps this correct without tracking which objects use which
# fesm-utils module. FESMUTILS_LIB (config/common.mk) is required by every build.
$(objdir)/%.o: $(srcdir)/%.f90 $(FESMUTILS_LIB) | $(objdir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) -c $< -o $@

$(objdir):
	mkdir -p $(objdir)

$(bindir):
	mkdir -p $(bindir)

# --- Library -----------------------------------------------------------------
vilma-static: $(obj_vilma)
	ar rcs $(objdir)/libvilma.a $(obj_vilma)
	@echo ""
	@echo "    $(objdir)/libvilma.a is ready."
	@echo ""

# --- Standalone driver -------------------------------------------------------
vilma: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(srcdir)/vilma_main.f90 \
		-o $(bindir)/vilma.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/vilma.x is ready."

# --- Offline lon-lat -> Gauss remapper ---------------------------------------
vilma_remap: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(srcdir)/vilma_remap_main.f90 \
		-o $(bindir)/vilma_remap.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/vilma_remap.x is ready."

# --- Offline reference (bed + ice) generation onto the Gauss grid ------------
vilma_mkref: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(srcdir)/vilma_mkref_main.f90 \
		-o $(bindir)/vilma_mkref.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/vilma_mkref.x is ready."

# --- Tests -------------------------------------------------------------------
test_params: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_params.f90 \
		-o $(bindir)/test_params.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_params.x is ready."

test_drive: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_drive.f90 \
		-o $(bindir)/test_drive.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_drive.x is ready."

test_band: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_band.f90 \
		-o $(bindir)/test_band.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_band.x is ready."

test_sht: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_sht.f90 \
		-o $(bindir)/test_sht.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_sht.x is ready."

test_earth: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_earth.f90 \
		-o $(bindir)/test_earth.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_earth.x is ready."

test_mesh: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_mesh.f90 \
		-o $(bindir)/test_mesh.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_mesh.x is ready."

test_integrals: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_integrals.f90 \
		-o $(bindir)/test_integrals.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_integrals.x is ready."

test_assembly: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_assembly.f90 \
		-o $(bindir)/test_assembly.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_assembly.x is ready."

test_love: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_love.f90 \
		-o $(bindir)/test_love.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_love.x is ready."

# Reference dumper for the Julia port (writes NetCDF oracle data; not in TESTS).
dump_reference: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/dump_reference.f90 \
		-o $(bindir)/dump_reference.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/dump_reference.x is ready."

test_relax: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_relax.f90 \
		-o $(bindir)/test_relax.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_relax.x is ready."

test_tidal: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_tidal.f90 \
		-o $(bindir)/test_tidal.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_tidal.x is ready."

test_rotation: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_rotation.f90 \
		-o $(bindir)/test_rotation.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_rotation.x is ready."

test_rotation_sle: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_rotation_sle.f90 \
		-o $(bindir)/test_rotation_sle.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_rotation_sle.x is ready."

test_etd1: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_etd1.f90 \
		-o $(bindir)/test_etd1.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_etd1.x is ready."

# Coupling-order characterization (§3c): measures the convergence order of the
# strain<->memory coupling. Standalone like test_etd1 -- a dt-sweep diagnostic,
# intentionally NOT in TESTS / `make check`. Build + run directly.
test_couple_order: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_couple_order.f90 \
		-o $(bindir)/test_couple_order.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_couple_order.x is ready."

# SLE<->memory coupling-order characterization (§3c 3b): drives a fast-evolving load
# through the full sea-level driver and measures the order restored by co-converging
# the ocean load σ and the end-of-step memory τ. Standalone dt-sweep diagnostic, like
# test_couple_order -- intentionally NOT in TESTS / `make check`. Build + run directly.
test_sle_couple_order: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_sle_couple_order.f90 \
		-o $(bindir)/test_sle_couple_order.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_sle_couple_order.x is ready."

# Adaptive-Δt controller (§3c): field step-doubling estimate order + the adaptive
# stepper converging to a fine reference with far fewer steps. Standalone diagnostic,
# NOT in `make check`.
test_timestep: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_timestep.f90 \
		-o $(bindir)/test_timestep.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_timestep.x is ready."

test_response: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_response.f90 \
		-o $(bindir)/test_response.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_response.x is ready."

test_sle: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_sle.f90 \
		-o $(bindir)/test_sle.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_sle.x is ready."

test_flotation: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_flotation.f90 \
		-o $(bindir)/test_flotation.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_flotation.x is ready."

test_ve_response: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_ve_response.f90 \
		-o $(bindir)/test_ve_response.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_ve_response.x is ready."

test_tensor_sh: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_tensor_sh.f90 \
		-o $(bindir)/test_tensor_sh.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_tensor_sh.x is ready."

test_response_3d: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_response_3d.f90 \
		-o $(bindir)/test_response_3d.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_response_3d.x is ready."

test_visc_lid: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_visc_lid.f90 \
		-o $(bindir)/test_visc_lid.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_visc_lid.x is ready."

test_visc_load: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_visc_load.f90 \
		-o $(bindir)/test_visc_load.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_visc_load.x is ready."

test_rotinv: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_rotinv.f90 \
		-o $(bindir)/test_rotinv.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_rotinv.x is ready."

test_benchmark_lvz: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_benchmark_lvz.f90 \
		-o $(bindir)/test_benchmark_lvz.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_benchmark_lvz.x is ready."

test_sle_ve: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_sle_ve.f90 \
		-o $(bindir)/test_sle_ve.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_sle_ve.x is ready."

test_benchmark_love: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_benchmark_love.f90 \
		-o $(bindir)/test_benchmark_love.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_benchmark_love.x is ready."

test_coupling: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_coupling.f90 \
		-o $(bindir)/test_coupling.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_coupling.x is ready."

test_restart: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_restart.f90 \
		-o $(bindir)/test_restart.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_restart.x is ready."

test_couple_remap: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_couple_remap.f90 \
		-o $(bindir)/test_couple_remap.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_couple_remap.x is ready."

test_spinup: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_spinup.f90 \
		-o $(bindir)/test_spinup.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_spinup.x is ready."

test_benchmark_disc: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_benchmark_disc.f90 \
		-o $(bindir)/test_benchmark_disc.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_benchmark_disc.x is ready."

test_benchmark_martinec: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_benchmark_martinec.f90 \
		-o $(bindir)/test_benchmark_martinec.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_benchmark_martinec.x is ready."

test_field: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_field.f90 \
		-o $(bindir)/test_field.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_field.x is ready."

test_flotation_load: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_flotation_load.f90 \
		-o $(bindir)/test_flotation_load.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_flotation_load.x is ready."

test_marine_reference: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_marine_reference.f90 \
		-o $(bindir)/test_marine_reference.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_marine_reference.x is ready."

# Standalone SLE benchmark (Martinec-2018 case E2): ~750 steps at lmax=128, runs
# in minutes -- intentionally NOT in TESTS / `make check`. Build with `make
# openmp=1 test_benchmark_sle` and run $(bindir)/test_benchmark_sle.x directly.
test_benchmark_sle: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_benchmark_sle.f90 \
		-o $(bindir)/test_benchmark_sle.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_benchmark_sle.x is ready."

test_sle_eustatic: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_sle_eustatic.f90 \
		-o $(bindir)/test_sle_eustatic.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_sle_eustatic.x is ready."

test_sle_subgrid: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_sle_subgrid.f90 \
		-o $(bindir)/test_sle_subgrid.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_sle_subgrid.x is ready."

test_remap: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_remap.f90 \
		-o $(bindir)/test_remap.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_remap.x is ready."

test_toroidal: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_toroidal.f90 \
		-o $(bindir)/test_toroidal.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/test_toroidal.x is ready."

diag_tensor_grid: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/diag_tensor_grid.f90 \
		-o $(bindir)/diag_tensor_grid.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/diag_tensor_grid.x is ready."

diag_visc3d_paths: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/diag_visc3d_paths.f90 \
		-o $(bindir)/diag_visc3d_paths.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/diag_visc3d_paths.x is ready."

# 3-D lateral-viscosity benchmark (protocol output, no pass/fail; not in check).
bench_visc3d: vilma-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/bench_visc3d.f90 \
		-o $(bindir)/bench_visc3d.x $(objdir)/libvilma.a $(LFLAGS)
	@echo "    $(bindir)/bench_visc3d.x is ready."

TESTS = test_params test_drive test_band test_sht test_earth test_mesh test_integrals test_assembly test_love test_relax test_tidal test_rotation test_rotation_sle test_response test_sle test_flotation test_flotation_load test_marine_reference test_etd1 test_ve_response test_tensor_sh test_response_3d test_toroidal test_sle_ve test_benchmark_love test_coupling test_couple_remap test_spinup test_restart test_benchmark_disc test_benchmark_martinec test_field test_sle_subgrid test_visc_load test_visc_lid test_rotinv test_remap

check: $(TESTS)
	@echo ""
	@echo "=== Running VILMA test suite ==="
	@for t in $(TESTS); do \
		echo "--- $$t ---"; \
		$(bindir)/$$t.x || exit 1; \
	done
	@echo ""
	@echo "=== All tests passed ==="

# Slow quantitative benchmarks: full-resolution, multi-minute runs kept OUT of
# `make check`. Build with `make openmp=1 check-slow` (the per-step solve is ~100x
# faster threaded). The SLE benchmark sweeps all four Martinec migrating-coast cases;
# test_rotinv re-runs the off-pole rotational-invariance check at full resolution
# (lmax 128, vs lmax 16 in `make check`).
SLOW = test_benchmark_sle test_rotinv test_benchmark_lvz
# STACK: the 3-D tensor path (vilma_tensor_sh) holds several whole-grid automatic
# arrays per thread. At lmax >= 128 that overruns the default 8 MB stack and the
# run dies with SIGSEGV at an address just below the stack top -- a failure that
# looks like a code bug and is not one. Raise the limit here so the target works
# out of the box; an external user should not have to know this. (`make check`
# needs none of it: test_rotinv runs at lmax 16 there.)
BIGSTACK = ulimit -s unlimited 2>/dev/null || ulimit -s 262144 2>/dev/null || true

check-slow: $(SLOW)
	@echo ""
	@echo "=== Running VILMA slow benchmarks ==="
	@for c in C2 D3 E2 F1; do \
		echo "--- test_benchmark_sle $$c ---"; \
		( $(BIGSTACK); OMP_STACKSIZE=$${OMP_STACKSIZE:-256M} \
		  $(bindir)/test_benchmark_sle.x $$c ) || exit 1; \
	done
	@echo "--- test_rotinv (lmax 128) ---"
	@( $(BIGSTACK); OMP_STACKSIZE=$${OMP_STACKSIZE:-256M} \
	   $(bindir)/test_rotinv.x 128 ) || exit 1
	@echo "--- test_benchmark_lvz (Weerdesteijn 2023 LVZ, lmax 512) ---"
	@( $(BIGSTACK); OMP_STACKSIZE=$${OMP_STACKSIZE:-256M} \
	   $(bindir)/test_benchmark_lvz.x ) || exit 1
	@echo ""
	@echo "=== All slow benchmarks passed ==="

# --- Housekeeping ------------------------------------------------------------
.PHONY: usage check check-slow clean showconfig vilma vilma-static

usage:
	@echo ""
	@echo "    * VILMA build *"
	@echo ""
	@echo " make vilma-static : build libvilma.a"
	@echo " make vilma        : build the standalone driver (bin/vilma.x)"
	@echo " make check            : build + run the test suite"
	@echo " make openmp=1 check-slow : build + run the slow full-res benchmarks"
	@echo " make test_sht         : build the SHT round-trip test"
	@echo " make clean            : remove objects and binaries"
	@echo " make showconfig       : show the active build configuration"
	@echo ""
	@echo "   switches:  debug=0|1|2   openmp=0|1   vilma_v1=0|1"
	@echo ""
	@echo "   vilma_v1=1 additionally links the optional VILMA-v1 backend (solver=\"v1\"):"
	@echo "     make vilma vilma_v1=1 VILMA_V1_ROOT=/path/to/vilma1"
	@echo "   It is OFF by default and is not a dependency; see doc/vilma-v1-backend.md."
	@echo ""

showconfig:
	@echo "----------------------"
	@echo "VILMA build configuration"
	@echo "----------------------"
	@echo "compiler  : $(FC)"
	@echo "host      : $(shell hostname)"
	@echo "openmp    : $(openmp)"
	@echo "debug     : $(debug)"
	@echo "vilma_v1  : $(vilma_v1)   (VILMA_V1_ROOT=$(VILMA_V1_ROOT))"
	@echo "FFLAGS    : $(FFLAGS)"
	@echo "LFLAGS    : $(LFLAGS)"

clean:
	rm -f $(objdir)/*.o $(objdir)/*.mod $(objdir)/*.a
	rm -f $(bindir)/*.x
	rm -rf $(bindir)/*.dSYM
