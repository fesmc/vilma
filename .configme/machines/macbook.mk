# Machine configuration: macbook (intended compiler: gfortran).
#
# PROJECT-TIER fragment: VILMA2's own copy, which takes precedence over
# configme's shipped machines/macbook.mk (orchestrator > user > shipped). It
# reproduces the shipped fragment and additionally pins the production
# optimization flags, which configme's generic compilers/gfortran.mk cannot
# carry (it defaults to -O2 for every project and every architecture).
#
# netCDF is auto-detected by configme (nf-config/nc-config) — no NC_FROOT /
# NC_CROOT in .zshrc required. To pin it instead, assign INC_NC / LIB_NC here
# and configme will use those as an override.

# Disable the default -Wl,-zmuldefs: gfortran forwards it to Apple's ld (ld64),
# which rejects it ("ld: unknown options: -zmuldefs"). It is a GNU-ld/ELF flag.
LFLAGS_EXTRA =

# --- Production optimization (overrides the gfortran default -O2) ------------
# Rationale in doc/performance-assessment.md §1. The hot kernels
# (dissipative_rhs, advance_memory, the band solve) are tight double-precision
# loops that benefit directly.
#
# -mcpu=native is the aarch64 (Apple Silicon) idiom for "tune+ISA for this
# host" (the x86 -march=native equivalent); it enables NEON vectorization.
# -ffast-math matches the CLIMBER-X production build (the host turns it on and
# we cannot reliably configure the coupled build differently), so we keep it on
# for parity — NOTE it relaxes IEEE semantics, so re-run `make check` on the
# target machine to confirm the SLE mass-conservation tests (~1e-16) still pass
# after a toolchain or flag change.
#
# This is a per-machine override: other machines (levante, albedo, pik_hpc2024)
# still build at the shipped -O2 until they get an equivalent fragment.
DFLAGS_NODEBUG = -O3 -mcpu=native -funroll-loops -ffast-math
