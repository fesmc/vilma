module vilma_response
   !! Surface-load response operator: the abstraction the sea-level equation
   !! (vilma_sle) is built on. Given a spectral surface mass-density load σ_lm
   !! [kg m^-2], it returns the two fields the SLE needs,
   !!
   !!     u_lm  — radial displacement of the solid surface  [m]
   !!     n_lm  — geoid / sea-surface-equipotential height   [m]
   !!
   !! per spherical-harmonic coefficient. The SLE depends only on this interface,
   !! so the elastic and (later) viscoelastic earth responses are swappable.
   !!
   !! Geoid mapping. The per-degree solve returns U(a) and F(a) = φ₁(a), the
   !! surface coefficients of radial displacement and the perturbed gravitational
   !! potential. Martinec's φ₁ carries the load's own direct potential with the
   !! sign OPPOSITE to φ^L (φ₁ → −φ^L for a rigid sphere, k = −F/φ^L − 1; see
   !! vilma_radial_fe%loading_love). The geopotential perturbation is therefore −F,
   !! and Bruns' formula gives the geoid height
   !!
   !!     N(a) = −F(a)/g .
   !!
   !! This uses only U and F — NOT the horizontal Love number l, whose sign /
   !! normalization is still being calibrated — so the SLE is not blocked by
   !! that open item. Both U and the −F/g geoid are pinned by the validated
   !! rigid (U→0, 1+k→1) and fluid (U→−(2j+1)/3·φ^L/g, 1+k→0) limits.
   use vilma_precision,       only: wp
   use vilma_constants,       only: pi, grav_G
   use vilma_earth_structure, only: earth_gravity_at, earth_model, RHEOL_MAXWELL
   use vilma_radial_fe,       only: radial_operator_load_rhs, radial_operator_solve_vec, radial_operator_destroy, radial_operator_solve, radial_operator_assemble, radial_mesh_build, radial_mesh, radial_operator, &
                                 idx_u, idx_v, idx_f, ndof_of, &
                                 toroidal_operator, toroidal_operator_assemble, toroidal_operator_solve_vec, &
                                 toroidal_operator_destroy
   use vilma_viscoelastic,    only: NLAM, ve_strain_constants, dissipative_rhs, &
                                 advance_memory, strain_coeffs, scheme_is_implicit, &
                                 SCHEME_FE, SCHEME_TRAP, &
                                 NLAM_TOR, ve_strain_constants_tor, dissipative_rhs_tor, &
                                 advance_memory_tor, strain_coeffs_tor
   use vilma_sht,             only: sht_grid, sht_grid_lmidx, sht_grid_synthesis, sht_grid_analysis
   use vilma_tensor_sh,       only: tensor_sh, TLAM_SPH, tensor_sh_init, tensor_sh_thread_cfg, tensor_sh_synth, tensor_sh_analysis, tensor_sh_destroy
   use, intrinsic :: iso_c_binding, only: c_ptr
   implicit none
   private

   public :: response, RESP_NULL, RESP_ELASTIC, RESP_VE
   public :: response_init_null, response_init_elastic, response_init_ve
   public :: response_apply, response_horizontal, response_horizontal_toroidal, response_destroy
   public :: response_begin_step, response_commit_step
   public :: response_prepare_endpoint, response_advance_endpoint
   public :: response_endpoint_converged, response_finalize_step
   public :: response_set_dt, response_prime_sigma
   public :: response_save_state, response_restore_state, response_stash_coarse
   public :: response_coarse_fine_error, response_max_rate, response_memory_norm
   public :: response_enable_lateral_visc, response_enable_lateral_visc_from_nodes
   public :: response_radial_rate

   integer, parameter :: RESP_NULL = 0, RESP_ELASTIC = 1, RESP_VE = 2


   type :: response
      !! Surface-load response operator as a tagged union. `kind` selects the
      !! behaviour and the free functions response_* dispatch on it:
      !!   RESP_NULL    rigid, non-self-gravitating: u ≡ 0, N ≡ 0 (eustatic limit)
      !!   RESP_ELASTIC time-independent per-degree gains, precomputed once
      !!   RESP_VE      viscoelastic field driver (stateful Maxwell memory)
      !! The sea-level equation (vilma_sle) depends only on this interface, so the
      !! three earth responses are interchangeable. A default-initialised value
      !! (kind = RESP_NULL) is a valid null response with no construction needed.
      integer :: kind = RESP_NULL
      ! --- RESP_ELASTIC: per-degree elastic surface gains ---------------------
      real(wp), allocatable :: ugain(:)  !! (0:lmax) U(a) per unit σ_l  [m / (kg m^-2)]
      real(wp), allocatable :: ngain(:)  !! (0:lmax) N(a)=−F(a)/g per unit σ_l
      real(wp), allocatable :: vgain(:)  !! (0:lmax) V(a) per unit σ_l (horizontal)
      ! --- RESP_VE: viscoelastic field-driver state ---------------------------



      !! Viscoelastic field driver. Holds one per-degree saddle-point operator
      !! (assembled + factored once) shared across all orders m, plus an
      !! independent Maxwell memory-stress history per spectral coefficient
      !! (l,m) — each (l,m) load has its own time history, so memory cannot be
      !! collapsed across m. Because the operator and the M = μΔt/η factors are
      !! real, each complex (l,m) history is two real histories (re/im).
      !!
      !! The per-step response is affine: solving the unit load once per degree
      !! gives the elastic gains gu(l), gn(l) AND the nodal field used to update
      !! memory; begin_step solves the frozen memory forcing per (l,m) for the
      !! drift; apply combines them; commit_step advances the memory.
      integer  :: lmax = 0, nr = 0, ne = 0, ndof = 0, nlm = 0
      real(wp) :: g = 0.0_wp, a = 0.0_wp, dt = 0.0_wp, time = 0.0_wp
      !! begin_step skips the drift solve for a coefficient whose Maxwell memory is
      !! below skip_tol × (the largest memory over all coefficients): its drift is
      !! negligible, so it is set to zero rather than solved. Self-consistent (all
      !! memory is zero at t=0 ⇒ all skipped ⇒ exact elastic first step) and cheap
      !! to gate. skip_tol = 0 disables skipping (solve every coefficient).
      real(wp) :: skip_tol = 1.0e-4_wp
      real(wp), allocatable :: mnorm(:)                 !! (nk) max|memory| per slot
      type(radial_operator), allocatable :: ops(:)      !! (1:lmax) per-degree operator
      ! degree-independent element fields
      real(wp), allocatable :: r(:)                     !! node radii (nr)
      real(wp), allocatable :: mu(:), Mk(:)             !! (ne) shear, M=μΔt/η
      ! MkPerDt==0 marks an element with NO Maxwell memory (elastic or fluid): set in
      ! init by rheology, so it is the single predicate for skipping such elements in
      ! the memory advance and for leaving them untouched by a loaded 3D viscosity field.
      ! Rung 6 — laterally-varying viscosity (3D). When lat_visc is set the Maxwell
      ! memory advance goes pseudo-spectral: the lateral product M(θ,φ)·τ couples
      ! harmonics, so per (element, component, radial shape-coeff) the memory and
      ! the current strain are synthesised to the Gauss grid, advanced pointwise
      ! τ⁺=(1−M)τ−2μM·ε with the lateral M-field, and analysed back (advance_memory_3d).
      ! MkPerDt3 is the Δt-invariant rate μ/η_eff on the grid per element; set_dt
      ! rescales Mk3 = MkPerDt3·Δt exactly, as for the 1-D Mk. With a laterally
      ! UNIFORM field this reproduces the 1-D advance to SHT round-trip precision.
      logical  :: lat_visc = .false.                    !! 3D lateral viscosity active
      real(wp), allocatable :: Mk3(:,:,:)               !! (nphi,nlat,ne) M=μΔt/η_eff
      real(wp), allocatable :: MkPerDt3(:,:,:)          !! (nphi,nlat,ne) μ/η_eff
      type(tensor_sh) :: tsh                            !! dyadic tensor-SH transformer
      ! VILMA1-style 1-D/3-D layer split (mod_visc3d k1p/k2p). An element is treated as
      ! genuinely "3-D" — and pays the pseudo-spectral tensor-SH advance — only when its
      ! lateral log10(η) spread exceeds visc3d_tol; otherwise it collapses to a scalar
      ! effective rate (lateral mean) and advances on the cheap degree-diagonal spectral
      ! path, exactly like 1-D viscosity. With a laterally-uniform field NO element is 3-D,
      ! so the 3-D run costs the same as 1-D. e3d lists the 3-D elements; active1d(e) is
      ! .true. for the spectrally-advanced (1-D-effective + elastic/fluid) elements.
      real(wp) :: visc3d_tol = 1.0e-3_wp                !! lateral log10(η) spread → "3-D" (dex)
      logical  :: deg1_cm = .false.                    !! degree-1 in the CM frame (geocenter motion kept)
      integer  :: ne3d = 0                              !! # genuinely-3-D elements
      integer,  allocatable :: e3d(:)                   !! (ne3d) indices of the 3-D elements
      logical,  allocatable :: active1d(:)              !! (ne) advance spectrally (skip in 3-D path)
      ! Toroidal degree of freedom (design-toroidal.md). A laterally varying
      ! viscosity drives toroidal flow through the pointwise M(θ,φ)·τ product, and
      ! that flow feeds back into the spheroidal memory (Martinec 2000 after eq
      ! 110). It exists only when some element is genuinely 3-D: a radially
      ! symmetric or laterally uniform field forces no W at all. Then the memory
      ! arrays widen from NLAM to NLAM+NLAM_TOR channels (λ3, λ4 appended), each
      ! degree gets its tridiagonal W operator, and W rides along as pure drift —
      ! a surface load never forces it directly (eq 84 has no δW term), so it has
      ! no elastic gain and no xWn.
      integer  :: nlam = NLAM                           !! memory channels carried: NLAM or NLAM+NLAM_TOR
      ! Carry the toroidal field when a 3-D element appears (default). .false.
      ! reproduces the spheroidal-only model exactly — the control for measuring
      ! what the toroidal coupling does (design-toroidal.md V7), not a speed knob.
      logical  :: toroidal = .true.
      type(toroidal_operator), allocatable :: tops(:)   !! (1:lmax) per-degree W operator
      real(wp), allocatable :: nrmt(:,:)                !! (NLAM_TOR,1:lmax) Z³,Z⁴ norms
      real(wp), allocatable :: sat(:,:,:), sbt(:,:,:), sct(:,:,:)  !! (2,NLAM_TOR,1:lmax)
      complex(wp), allocatable :: dWa(:)                !! (nk) surface toroidal drift W(a)
      real(wp), allocatable :: dWn_re(:,:), dWn_im(:,:) !! (nr,nk) nodal W from τ_n
      real(wp), allocatable :: edWn_re(:,:), edWn_im(:,:) !! (nr,nk) nodal W from the trial τ_{n+1}
      complex(wp), allocatable :: dWa_prev(:)           !! (nk) prev-iterate W(a) (convergence)
      ! per-degree constants and unit-load response
      real(wp), allocatable :: Jr(:)                    !! (1:lmax) l(l+1)
      real(wp), allocatable :: nrmc(:,:)                !! (NLAM,1:lmax) Z:Z norms
      real(wp), allocatable :: sa(:,:,:), sb(:,:,:), sc(:,:,:)  !! (4,NLAM,1:lmax)
      real(wp), allocatable :: gu(:), gn(:), gv(:)      !! (0:lmax) elastic gains (gv=V(a))
      real(wp), allocatable :: xUn(:,:), xVn(:,:)       !! (nr,1:lmax) unit-load nodal U,V
      ! Per-(l,m) state is stored in DEGREE-GROUPED order k = 1..nk (all orders m of
      ! a degree l contiguous, l ascending; degree 0 carries no memory and is
      ! excluded). This makes begin_step/commit_step iterate k contiguously AND
      ! reuse the per-degree operator ops(l) across its orders (cache-hot), instead
      ! of the SHTns m-major lm order which switches operator on nearly every solve.
      integer :: nk = 0                                 !! # deforming coeffs (l>=1)
      integer, allocatable :: k2lm(:)                   !! (nk) slot k -> SHTns lm index
      integer, allocatable :: kdeg(:)                   !! (nk) degree l of slot k
      integer, allocatable :: kbeg(:)                   !! (1:lmax+1) first k of each degree
      ! per-(l,m) memory stress, split into real/imag (NLAM,ne,nk)
      real(wp), allocatable :: Are(:,:,:), Aim(:,:,:)
      real(wp), allocatable :: Bre(:,:,:), Bim(:,:,:)
      real(wp), allocatable :: Cre(:,:,:), Cim(:,:,:)
      ! frozen per-step drift (from the memory τ_n), set by begin_step
      complex(wp), allocatable :: dUa(:), dFa(:), dVa(:) !! (nk) surface drift (U,F,V)
      real(wp), allocatable :: dUn_re(:,:), dUn_im(:,:) !! (nr,nk) nodal drift U from τ_n (ε_n)
      real(wp), allocatable :: dVn_re(:,:), dVn_im(:,:) !! (nr,nk) nodal drift V from τ_n
      ! Memory time-integration scheme (see vilma_viscoelastic). FE (default) advances
      ! the memory once in commit_step from the report strain — the historical path,
      ! bit-identical. TRAP (2nd-order) is implicit in the end-of-step strain, so
      ! commit_step Picard-iterates the endpoint (frozen load: σ held at the converged
      ! value) — re-solving the drift against the trial τ_{n+1} each pass. §3c part 3a.
      integer  :: scheme = SCHEME_FE
      integer  :: max_couple_iter = 1        !! coupling-iteration cap for implicit schemes
      real(wp) :: couple_tol = 1.0e-6_wp     !! relative surface-drift change to stop iterating
      integer  :: couple_iters_last = 0      !! iterations taken last commit (diagnostic)
      ! PROFILE: phase wall-clock [s] + call counts, cumulative. The host (vilma_drive)
      ! zeroes these at the start of the transient so the reported per-step breakdown
      ! covers the transient only. t_drift = solve_drift (per-degree band LU);
      ! t_mem = the memory advance (vilma_advance / trapezoid_advance_all, which in the
      ! 3-D case is dominated by the dyadic SHT round-trip). The remainder of the
      ! solid_earth_update cost (SLE iteration, load/geoid SHTs) is t_upd - t_drift - t_mem.
      real(wp) :: t_drift = 0.0_wp, t_mem = 0.0_wp
      integer  :: n_drift = 0,      n_mem = 0
      ! Co-convergence state for the SLE driver's σ<->τ fixed point (§3c 3b). Each
      ! advance_endpoint does ONE trapezoid pass and refreshes the report drift to
      ! the new τ_{n+1}; the driver re-converges σ against it and repeats until done.
      logical  :: couple_done = .false.      !! co-convergence reached (drift settled)
      integer  :: couple_pass = 0            !! co-convergence passes taken this step
      ! Scratch for the implicit commit (allocated lazily on first TRAP commit):
      real(wp), allocatable :: Are0(:,:,:), Aim0(:,:,:)  !! (NLAM,ne,nk) τ_n snapshot
      real(wp), allocatable :: Bre0(:,:,:), Bim0(:,:,:)
      real(wp), allocatable :: Cre0(:,:,:), Cim0(:,:,:)
      real(wp), allocatable :: edUn_re(:,:), edUn_im(:,:) !! (nr,nk) endpoint nodal drift U (ε_{n+1})
      real(wp), allocatable :: edVn_re(:,:), edVn_im(:,:) !! (nr,nk) endpoint nodal drift V
      complex(wp), allocatable :: dUa_prev(:)             !! (nk) prev-iterate surface drift (convergence)
      ! Start-of-step load σ_n for the trapezoidal ε_n term. The rule is ½(ε_n+ε_{n+1})
      ! with ε_n = σ_n·xUn + drift(τ_n) and ε_{n+1} = σ_{n+1}·xUn + drift(τ_{n+1}); ε_n
      ! must use the load at t_n, not the current σ_{n+1}. For a held load σ_n=σ_{n+1}
      ! so this is invisible (the historical/3a path), but for a fast-evolving load
      ! (the SLE 3b driver) it sets the order. sigma_n carries the previous step's
      ! converged load; sigma_next stages σ_{n+1} until finalize commits it. Until the
      ! first step finalizes (sigma_primed=.false.) ε_n falls back to σ_{n+1}, which
      ! reproduces the historical first step exactly for a load present at t=0.
      complex(wp), allocatable :: sigma_n(:), sigma_next(:)  !! (nlm)
      logical :: sigma_primed = .false.
      ! Δt-invariant memory rate Mk/Δt = μ/η (set once in init): set_dt rescales the
      ! memory factor Mk = MkPerDt·Δt exactly, so the adaptive controller can change
      ! Δt with no operator re-factor and no drift from repeated rescaling.
      real(wp), allocatable :: MkPerDt(:)                 !! (ne)
      ! Controller state snapshots (§3c controller, lazily allocated): buffer A holds
      ! the entering state τ_n that a rejected/fine step restores to; buffer B holds the
      ! coarse τ_{n+1} for the step-doubling error estimate. Distinct from Are0 (which
      ! sle_solve overwrites internally each step). See save_state/stash_coarse.
      real(wp), allocatable :: Are_s(:,:,:), Aim_s(:,:,:), Bre_s(:,:,:), Bim_s(:,:,:)
      real(wp), allocatable :: Cre_s(:,:,:), Cim_s(:,:,:)            !! buffer A (τ_n)
      real(wp), allocatable :: Are_c(:,:,:), Aim_c(:,:,:), Bre_c(:,:,:), Bim_c(:,:,:)
      real(wp), allocatable :: Cre_c(:,:,:), Cim_c(:,:,:)            !! buffer B (τ_coarse)
      real(wp)    :: time_s = 0.0_wp                       !! saved time (buffer A)
      complex(wp), allocatable :: sigma_n_s(:)             !! saved σ_n (buffer A)
      logical     :: sigma_primed_s = .false.
   end type response

contains

   ! ===================================================================
   ! Public interface: free functions over `response`. The eight SLE-facing
   ! operators dispatch on self%kind; the constructors set it. RESP_NULL and
   ! RESP_ELASTIC share the stateless (no-op) step brackets.
   ! ===================================================================

   subroutine response_init_null(self)
      !! Rigid, non-self-gravitating response (u ≡ 0, N ≡ 0). Equivalent to a
      !! default-initialised value; provided for explicit construction.
      type(response), intent(out) :: self
      self%kind = RESP_NULL
   end subroutine response_init_null

   subroutine response_apply(self, sht, sigma_lm, u_lm, n_lm)
      type(response), intent(inout) :: self
      type(sht_grid), intent(in)    :: sht
      complex(wp),    intent(in)    :: sigma_lm(:)
      complex(wp),    intent(out)   :: u_lm(:)
      complex(wp),    intent(out)   :: n_lm(:)
      select case (self%kind)
      case (RESP_ELASTIC); call elastic_response_apply(self, sht, sigma_lm, u_lm, n_lm)
      case (RESP_VE);      call ve_response_apply(self, sht, sigma_lm, u_lm, n_lm)
      case default;        call null_response_apply(self, sht, sigma_lm, u_lm, n_lm)
      end select
   end subroutine response_apply

   subroutine response_horizontal(self, sht, sigma_lm, v_lm)
      type(response), intent(in)    :: self
      type(sht_grid), intent(in)    :: sht
      complex(wp),    intent(in)    :: sigma_lm(:)
      complex(wp),    intent(out)   :: v_lm(:)
      select case (self%kind)
      case (RESP_ELASTIC); call elastic_response_horizontal(self, sht, sigma_lm, v_lm)
      case (RESP_VE);      call ve_response_horizontal(self, sht, sigma_lm, v_lm)
      case default;        call response_horizontal_default(self, sht, sigma_lm, v_lm)
      end select
   end subroutine response_horizontal

   subroutine response_begin_step(self, sht)
      type(response), intent(inout) :: self
      type(sht_grid), intent(in)    :: sht
      select case (self%kind)
      case (RESP_VE);    call ve_response_begin(self, sht)
      case default;      call response_begin_default(self, sht)
      end select
   end subroutine response_begin_step

   subroutine response_commit_step(self, sht, sigma_lm)
      type(response), intent(inout) :: self
      type(sht_grid), intent(in)    :: sht
      complex(wp),    intent(in)    :: sigma_lm(:)
      select case (self%kind)
      case (RESP_VE);    call ve_response_commit(self, sht, sigma_lm)
      case default;      call response_commit_default(self, sht, sigma_lm)
      end select
   end subroutine response_commit_step

   subroutine response_prepare_endpoint(self, sht)
      type(response), intent(inout) :: self
      type(sht_grid), intent(in)    :: sht
      select case (self%kind)
      case (RESP_VE);    call ve_response_prepare_endpoint(self, sht)
      case default;      call response_prepare_default(self, sht)
      end select
   end subroutine response_prepare_endpoint

   subroutine response_advance_endpoint(self, sht, sigma_lm)
      type(response), intent(inout) :: self
      type(sht_grid), intent(in)    :: sht
      complex(wp),    intent(in)    :: sigma_lm(:)
      select case (self%kind)
      case (RESP_VE);    call ve_response_advance_endpoint(self, sht, sigma_lm)
      case default;      call response_advance_default(self, sht, sigma_lm)
      end select
   end subroutine response_advance_endpoint

   logical function response_endpoint_converged(self) result(done)
      type(response), intent(in) :: self
      select case (self%kind)
      case (RESP_VE);    done = ve_response_endpoint_converged(self)
      case default;      done = response_converged_default(self)
      end select
   end function response_endpoint_converged

   subroutine response_finalize_step(self, sht)
      type(response), intent(inout) :: self
      type(sht_grid), intent(in)    :: sht
      select case (self%kind)
      case (RESP_VE);    call ve_response_finalize_step(self, sht)
      case default;      call response_finalize_default(self, sht)
      end select
   end subroutine response_finalize_step

   subroutine response_destroy(self)
      !! Release whatever state the response holds (safe for any kind: every
      !! deallocation is allocated-guarded).
      type(response), intent(inout) :: self
      call elastic_response_destroy(self)
      call ve_response_destroy(self)
      self%kind = RESP_NULL
   end subroutine response_destroy


   subroutine response_begin_default(self, sht)
      !! No-op step bracket for stateless responses.
      type(response), intent(inout) :: self
      type(sht_grid),           intent(in)    :: sht
   end subroutine response_begin_default

   subroutine response_commit_default(self, sht, sigma_lm)
      type(response), intent(inout) :: self
      type(sht_grid),           intent(in)    :: sht
      complex(wp),              intent(in)    :: sigma_lm(:)
   end subroutine response_commit_default

   subroutine response_prepare_default(self, sht)
      !! No-op endpoint bracket for stateless / 1st-order responses.
      type(response), intent(inout) :: self
      type(sht_grid),           intent(in)    :: sht
   end subroutine response_prepare_default

   subroutine response_advance_default(self, sht, sigma_lm)
      !! No memory to advance (stateless response): nothing to do.
      type(response), intent(inout) :: self
      type(sht_grid),           intent(in)    :: sht
      complex(wp),              intent(in)    :: sigma_lm(:)
   end subroutine response_advance_default

   logical function response_converged_default(self) result(done)
      !! Stateless / 1st-order responses converge in a single pass.
      type(response), intent(in) :: self
      done = .true.
   end function response_converged_default

   subroutine response_finalize_default(self, sht)
      type(response), intent(inout) :: self
      type(sht_grid),           intent(in)    :: sht
   end subroutine response_finalize_default

   subroutine response_horizontal_default(self, sht, sigma_lm, v_lm)
      !! No horizontal displacement for a rigid / non-deforming response.
      type(response), intent(in)    :: self
      type(sht_grid),           intent(in)    :: sht
      complex(wp),              intent(in)    :: sigma_lm(:)  !! load [kg m^-2]
      complex(wp),              intent(out)   :: v_lm(:)      !! spheroidal V(a) [m]
      v_lm = (0.0_wp, 0.0_wp)
   end subroutine response_horizontal_default

   subroutine null_response_apply(self, sht, sigma_lm, u_lm, n_lm)
      type(response), intent(inout) :: self
      type(sht_grid),       intent(in)    :: sht
      complex(wp),          intent(in)    :: sigma_lm(:)
      complex(wp),          intent(out)   :: u_lm(:)
      complex(wp),          intent(out)   :: n_lm(:)
      u_lm = (0.0_wp, 0.0_wp)
      n_lm = (0.0_wp, 0.0_wp)
   end subroutine null_response_apply

   subroutine response_deg1_to_cm(op, x)
      !! Put a solved degree-1 state into the CENTRE-OF-MASS frame.
      !!
      !! Degree 1 has a rigid-translation null space, so the solve needs a gauge
      !! and the gauge IS the reference frame. radial_operator constrains
      !! w'd = 0 (Martinec 2000 eq 83): zero volume-integrated displacement, a
      !! centre-of-FIGURE-like frame in which the solid Earth does not translate.
      !! The CM frame instead holds the centre of mass of Earth+load fixed and
      !! lets the solid Earth translate — geocenter motion, which is a real part
      !! of the degree-1 sea-level fingerprint.
      !!
      !! Because the translation is a null direction, changing frame is a
      !! post-solve projection and cannot disturb the deformation: add c·n with
      !! n the null mode (stored by radial_operator_assemble).
      !!
      !! c is fixed by the CM condition on the POTENTIAL. Outside the Earth the
      !! degree-1 potential of the whole system vanishes in the CM frame, which
      !! in Martinec's variables is F(a) = 0 — equivalently 1 + k₁ = 0, k₁ = −1.
      !! So c = −F(a)/F_n(a), and F(a) comes out exactly zero afterwards.
      !!
      !! This is the same condition the code already ASSERTED for the geoid by
      !! setting ngain(1) = 0 while leaving the displacement in the w'd = 0
      !! frame. That mixture is what this removes: here the condition is imposed
      !! once, on the state, so the geoid and the displacement refer to one frame.
      type(radial_operator), intent(in)    :: op
      real(wp),              intent(inout) :: x(:)
      real(wp) :: fa_n, c
      integer  :: nr
      if (.not. allocated(op%nullmode)) return          ! not the bordered degree
      nr   = op%nr
      fa_n = op%nullmode(idx_f(nr))
      if (fa_n == 0.0_wp) return                        ! translation carries no potential: nothing to fix
      c = -x(idx_f(nr)) / fa_n
      x = x + c*op%nullmode
   end subroutine response_deg1_to_cm

   subroutine response_init_elastic(self, earth, lmax)
      !! Precompute the per-degree elastic surface gains for degrees 0..lmax.
      !!
      !!   l = 0 : incompressibility (Div u = 0) forbids degree-0 radial
      !!           deformation, so U(0)=0; the geoid feels only the load's own
      !!           monopole potential, N(0) = φ^L_0/g = 4πGa/g per unit σ.
      !!   l ≥ 1 : assemble the per-degree saddle-point operator, solve a unit
      !!           surface load, store U(a) and N(a) = −F(a)/g.
      type(response), intent(inout) :: self
      type(earth_model),       intent(in)    :: earth
      integer,                 intent(in)    :: lmax
      type(radial_mesh)     :: mesh
      type(radial_operator) :: op
      integer  :: l
      real(wp) :: ua, va, fa

      call elastic_response_destroy(self)
      self%kind = RESP_ELASTIC
      self%lmax = lmax
      self%a    = earth%r_earth
      self%g    = earth_gravity_at(earth, earth%r_earth)
      allocate(self%ugain(0:lmax), self%ngain(0:lmax), self%vgain(0:lmax))

      ! degree 0: no deformation, pure monopole geoid offset
      self%ugain(0) = 0.0_wp
      self%ngain(0) = 4.0_wp*pi*grav_G*self%a / self%g
      self%vgain(0) = 0.0_wp                    ! no horizontal at degree 0

      call radial_mesh_build(mesh, earth)
      do l = 1, lmax
         call radial_operator_assemble(op, earth, mesh, l)
         if (l == 1 .and. self%deg1_cm) then
            ! CM frame: solve for the whole state so the frame projection can be
            ! applied, then read the surface coefficients back off it.
            block
               real(wp), allocatable :: x1(:)
               allocate(x1(op%ndof))
               call radial_operator_solve_vec(op, radial_operator_load_rhs(op, 1.0_wp), x1)
               call response_deg1_to_cm(op, x1)
               ua = x1(idx_u(op%nr));  va = x1(idx_v(op%nr));  fa = x1(idx_f(op%nr))
               deallocate(x1)
            end block
         else
            call radial_operator_solve(op, 1.0_wp, ua, va, fa)  ! unit surface load coefficient
         end if
         self%ugain(l) = ua
         self%ngain(l) = -fa / self%g
         self%vgain(l) = va
         call radial_operator_destroy(op)
      end do

      ! degree-1 frame. Two conventions, `deg1_frame` in the namelist:
      !
      !   "cf": the per-degree solve fixes the displacement gauge
      !        (wᵀd=0, no volume-integrated translation, centre-of-figure-like)
      !        while the geoid is referenced to CM, where the degree-1 external
      !        potential vanishes ⇒ N₁≡0. Validated against the Spada-2011 disc
      !        n_disc, which matches once N₁ is dropped. The two halves are in
      !        DIFFERENT frames, so rsl carries no degree 1 at all.
      !   "cm" (default): response_deg1_to_cm has already put the state in the CM frame, so
      !        F(a)=0 came out of the projection and ngain(1) is zero as a RESULT,
      !        not an override — and the displacement carries geocenter motion.
      !
      ! VILMA1 runs in CM and reports the term in vega_deg1.dat; the F−V residual
      ! on the disc benchmark is 98–99.7 % degree 1 with "cf". See
      ! paper-fastearth3d-experiments notes/vilma-comparison.md §11.
      if (lmax >= 1 .and. .not. self%deg1_cm) self%ngain(1) = 0.0_wp
   end subroutine response_init_elastic

   subroutine elastic_response_apply(self, sht, sigma_lm, u_lm, n_lm)
      !! Spectral multiply: u_lm = ugain(l)·σ_lm, n_lm = ngain(l)·σ_lm. Degrees
      !! above the precomputed lmax are zeroed.
      type(response), intent(inout) :: self
      type(sht_grid),          intent(in)    :: sht
      complex(wp),             intent(in)    :: sigma_lm(:)
      complex(wp),             intent(out)   :: u_lm(:)
      complex(wp),             intent(out)   :: n_lm(:)
      integer :: l, m, lm, lcap

      u_lm = (0.0_wp, 0.0_wp)
      n_lm = (0.0_wp, 0.0_wp)
      lcap = min(self%lmax, sht%lmax)
      do m = 0, sht%mmax*sht%mres, sht%mres
         do l = m, lcap
            lm = sht_grid_lmidx(sht, l, m)
            u_lm(lm) = self%ugain(l) * sigma_lm(lm)
            n_lm(lm) = self%ngain(l) * sigma_lm(lm)
         end do
      end do
   end subroutine elastic_response_apply

   subroutine elastic_response_horizontal(self, sht, sigma_lm, v_lm)
      !! Spheroidal multiply: v_lm = vgain(l)·σ_lm (degree-1 left as solved, like
      !! ugain — the horizontal displacement is in the CE-like gauge, not the geoid
      !! CM frame). Synthesize ∇₁(Σ v_lm Y_lm) for (u_θ, u_φ).
      type(response), intent(in)    :: self
      type(sht_grid),          intent(in)    :: sht
      complex(wp),             intent(in)    :: sigma_lm(:)
      complex(wp),             intent(out)   :: v_lm(:)
      integer :: l, m, lm, lcap
      v_lm = (0.0_wp, 0.0_wp)
      lcap = min(self%lmax, sht%lmax)
      do m = 0, sht%mmax*sht%mres, sht%mres
         do l = m, lcap
            lm = sht_grid_lmidx(sht, l, m)
            v_lm(lm) = self%vgain(l) * sigma_lm(lm)
         end do
      end do
   end subroutine elastic_response_horizontal

   subroutine elastic_response_destroy(self)
      type(response), intent(inout) :: self
      if (allocated(self%ugain)) deallocate(self%ugain)
      if (allocated(self%ngain)) deallocate(self%ngain)
      if (allocated(self%vgain)) deallocate(self%vgain)
      self%lmax = 0
   end subroutine elastic_response_destroy

   ! --- viscoelastic field driver ---------------------------------------------

   subroutine response_init_ve(self, earth, sht, dt)
      !! Assemble the per-degree operators, precompute the unit-load response and
      !! Maxwell constants, and zero the per-(l,m) memory. Tied to the grid sht
      !! (sets lmax = sht%lmax and the coefficient layout).
      type(response), intent(inout) :: self
      type(earth_model),  intent(in)    :: earth
      type(sht_grid),     intent(in)    :: sht
      real(wp),           intent(in)    :: dt
      type(radial_mesh) :: mesh
      real(wp), allocatable :: x(:)
      real(wp) :: eta_e
      integer  :: l, m, e, lay, node, k

      call ve_response_destroy(self)
      self%kind = RESP_VE
      call radial_mesh_build(mesh, earth)
      self%lmax = sht%lmax;  self%nlm = sht%nlm
      self%nr = mesh%nr;  self%ne = mesh%ne;  self%ndof = ndof_of(mesh%nr)
      self%dt = dt;  self%time = 0.0_wp
      self%a  = earth%r_earth;  self%g = earth_gravity_at(earth, earth%r_earth)

      ! degree-independent element fields: node radii, shear, Maxwell factor
      allocate(self%r(self%nr));  self%r = mesh%r
      allocate(self%mu(self%ne), self%Mk(self%ne), self%MkPerDt(self%ne))
      allocate(self%active1d(self%ne));  self%active1d = .true.   ! all spectral until 3-D split
      do e = 1, self%ne
         lay = mesh%elem_layer(e)
         self%mu(e) = earth%layers(lay)%mu
         eta_e      = earth%layers(lay)%eta
         ! Only genuinely Maxwell layers carry memory. Classify by RHEOLOGY, not by
         ! eta>0: the elastic lithosphere is stored with eta=huge, which would give a
         ! tiny-but-NONZERO rate μ/huge≈3e-298 — close enough to be inert in the 1-D
         ! advance (Mk rounds away) but not exactly 0, which muddies the "no memory"
         ! test. Setting it (and the inviscid core) to exactly 0 makes MkPerDt==0 the
         ! universal "this element has no Maxwell memory" predicate, used to skip
         ! elastic/fluid elements in the advance and to leave them untouched by a
         ! loaded 3D viscosity field. (Observably identical: a μ/huge rate only ever
         ! produced ~1e-267 memory, below ULP in every force/uplift.)
         if (earth%layers(lay)%rheology == RHEOL_MAXWELL) then
            self%Mk(e)      = self%mu(e)*dt/eta_e     ! M=μΔt/η; set_dt rescales from MkPerDt
            self%MkPerDt(e) = self%mu(e)/eta_e        ! Δt-invariant rate μ/η
         else
            self%Mk(e)      = 0.0_wp                  ! elastic (η→∞) / fluid (μ=0): no memory
            self%MkPerDt(e) = 0.0_wp
         end if
      end do

      ! per-degree constants, operators, and unit-load response
      allocate(self%Jr(self%lmax), self%nrmc(NLAM,self%lmax))
      allocate(self%sa(4,NLAM,self%lmax), self%sb(4,NLAM,self%lmax), &
               self%sc(4,NLAM,self%lmax))
      allocate(self%gu(0:self%lmax), self%gn(0:self%lmax))
      allocate(self%xUn(self%nr,self%lmax), self%xVn(self%nr,self%lmax))
      allocate(self%gv(0:self%lmax))
      allocate(self%ops(self%lmax), x(self%ndof))

      ! degree 0: monopole geoid, no deformation, no memory (no operator).
      ! degree 1: geocenter motion — carried in the CM frame. The sparse KKT
      ! border in radial_operator removes the rigid-translation null space
      ! (wᵀd = 0, Blewitt 2003), so j=1 assembles and steps like any other
      ! degree; it joins the l-loop below.
      self%gu(0) = 0.0_wp
      self%gn(0) = 4.0_wp*pi*grav_G*self%a/self%g
      self%gv(0) = 0.0_wp                       ! no horizontal at degree 0

      do l = 1, self%lmax
         self%Jr(l) = real(l, wp)*real(l+1, wp)
         call ve_strain_constants(self%Jr(l), self%nrmc(:,l), &
                                  self%sa(:,:,l), self%sb(:,:,l), self%sc(:,:,l))
         call radial_operator_assemble(self%ops(l), earth, mesh, l)
         call radial_operator_solve_vec(self%ops(l), radial_operator_load_rhs(self%ops(l), 1.0_wp), x)
         if (l == 1 .and. self%deg1_cm) call response_deg1_to_cm(self%ops(l), x)
         self%gu(l) = x(idx_u(self%nr))
         self%gn(l) = -x(idx_f(self%nr))/self%g
         self%gv(l) = x(idx_v(self%nr))         ! surface horizontal V(a)
         do node = 1, self%nr
            self%xUn(node,l) = x(idx_u(node))
            self%xVn(node,l) = x(idx_v(node))
         end do
      end do

      ! degree-1 geoid frame (see response_init_elastic): the geoid is referenced
      ! to the CM frame ⇒ N₁≡0. Zero the degree-1 geoid gain here; the degree-1
      ! relaxation drift is likewise zeroed in begin_step. Displacement (gu(1),
      ! xUn/xVn) is left as solved (CE-like geocenter, h₁≈0).
      if (self%lmax >= 1 .and. .not. self%deg1_cm) self%gn(1) = 0.0_wp

      ! Degree-grouped coefficient map: slot k = 1..nk over (l>=1, m=0..min(l,mmax)),
      ! l ascending then m ascending. k2lm bridges back to the SHTns lm index for
      ! the load/uplift/geoid spectra; kdeg gives the degree (operator) per slot.
      self%nk = 0
      do l = 1, self%lmax
         do m = 0, min(l, sht%mmax*sht%mres), sht%mres
            self%nk = self%nk + 1
         end do
      end do
      allocate(self%k2lm(self%nk), self%kdeg(self%nk), self%kbeg(self%lmax+1))
      k = 0
      do l = 1, self%lmax
         self%kbeg(l) = k + 1                  ! first slot of degree l (contiguous)
         do m = 0, min(l, sht%mmax*sht%mres), sht%mres
            k = k + 1
            self%k2lm(k) = sht_grid_lmidx(sht, l, m)
            self%kdeg(k) = l
         end do
      end do
      self%kbeg(self%lmax+1) = k + 1           ! sentinel (one past the last slot)

      ! per-(l,m) memory + drift, all in degree-grouped k order (zeroed)
      self%nlam = NLAM                          ! spheroidal only until a 3-D element appears
      allocate(self%Are(self%nlam,self%ne,self%nk), self%Aim(self%nlam,self%ne,self%nk))
      allocate(self%Bre(self%nlam,self%ne,self%nk), self%Bim(self%nlam,self%ne,self%nk))
      allocate(self%Cre(self%nlam,self%ne,self%nk), self%Cim(self%nlam,self%ne,self%nk))
      ! Zeroed in PARALLEL, on the same schedule(static) partition over k that
      ! every consumer of these arrays uses (vilma_advance, trapezoid_advance_all,
      ! response_memory_norm, the save/restore buffers). This is the first touch,
      ! so it is what maps the pages onto NUMA domains: zeroing serially puts all
      ! ~6 x NLAM x ne x nk of it on the master thread's domain, and every threaded
      ! loop over it then reads across the socket interconnect. On a 2-socket EPYC
      ! 7763 that is the difference between local and remote bandwidth for the
      ! single largest array set in the model.
      !$omp parallel do default(shared) private(k) schedule(static)
      do k = 1, self%nk
         self%Are(:,:,k) = 0.0_wp;  self%Aim(:,:,k) = 0.0_wp
         self%Bre(:,:,k) = 0.0_wp;  self%Bim(:,:,k) = 0.0_wp
         self%Cre(:,:,k) = 0.0_wp;  self%Cim(:,:,k) = 0.0_wp
      end do
      allocate(self%dUa(self%nk), self%dFa(self%nk), self%dVa(self%nk))
      allocate(self%dUn_re(self%nr,self%nk), self%dUn_im(self%nr,self%nk))
      allocate(self%dVn_re(self%nr,self%nk), self%dVn_im(self%nr,self%nk))
      allocate(self%mnorm(self%nk))
      self%dUa = (0.0_wp,0.0_wp); self%dFa = (0.0_wp,0.0_wp); self%dVa = (0.0_wp,0.0_wp)
      self%dUn_re = 0.0_wp; self%dUn_im = 0.0_wp
      self%dVn_re = 0.0_wp; self%dVn_im = 0.0_wp
      self%mnorm = 0.0_wp                       ! zero memory ⇒ all slots skipped initially
   end subroutine response_init_ve

   subroutine ve_response_begin(self, sht)
      !! Freeze the per-(l,m) drift from the entering memory τ_n: solve the memory
      !! forcing −∫τ^V:δε with the load held at zero, storing surface + nodal drift.
      !! The nodal drift (the ε_n term) lands in self%dUn_*/dVn_*; the implicit
      !! commit re-uses the same solver against the trial τ_{n+1} (see solve_drift).
      type(response), intent(inout) :: self
      type(sht_grid),     intent(in)    :: sht
      if (self%nlam > NLAM) then
         call solve_drift(self, sht, self%dUn_re, self%dUn_im, self%dVn_re, self%dVn_im, &
                          self%dWn_re, self%dWn_im)
      else
         call solve_drift(self, sht, self%dUn_re, self%dUn_im, self%dVn_re, self%dVn_im)
      end if
   end subroutine ve_response_begin

   subroutine solve_drift(self, sht, Un_re, Un_im, Vn_re, Vn_im, Wn_re, Wn_im)
      !! Solve the per-(l,m) drift (load=0 memory forcing) from self's CURRENT memory
      !! arrays (self%Are…): surface drift → self%dUa/dFa/dVa, nodal drift → the four
      !! target arrays. begin_step passes self%dUn_* (drift from τ_n); the implicit
      !! commit passes self%edUn_* (drift from the trial τ_{n+1}). Targets are distinct
      !! components from everything read via self, so there is no argument aliasing.
      !! With the toroidal channels carried, Wn_re/Wn_im (required then) receive the
      !! nodal W from its own operator, and self%dWa its surface value.
      type(response), intent(inout) :: self
      type(sht_grid),     intent(in)    :: sht
      real(wp),           intent(out)   :: Un_re(:,:), Un_im(:,:), Vn_re(:,:), Vn_im(:,:)
      real(wp), optional, intent(out)   :: Wn_re(:,:), Wn_im(:,:)
      real(wp), allocatable :: fre(:), fim(:), xre(:), xim(:)
      real(wp), allocatable :: gre(:), gim(:), wre(:), wim(:)
      logical :: tor
      integer :: k, l, node, nb
      real(wp) :: thr
      integer(kind=8) :: pc0, pc1, prate           ! PROFILE: drift-solve wall-clock
      call system_clock(pc0, prate)

      ! Refresh the per-slot memory magnitude from the current memory (a cheap pass
      ! vs the solves) so it is always consistent — including after a restart, which
      ! reloads the memory arrays but not this derived cache. Explicit loop (no
      ! abs(slice) temporaries, which would each heap-allocate).
      ! Each (nlam, ne) slot block is contiguous, so it is scanned as one flat run
      ! of nb values (absmax6): a runtime-width channel loop inside the element
      ! loop does not vectorize, and cost the 1-D drift solve 5 % once nlam
      ! stopped being NLAM.
      nb = self%nlam*self%ne
      !$omp parallel do default(shared) private(k) schedule(static)
      do k = 1, self%nk
         self%mnorm(k) = absmax6(nb, self%Are(1,1,k), self%Bre(1,1,k), self%Cre(1,1,k), &
                                 self%Aim(1,1,k), self%Bim(1,1,k), self%Cim(1,1,k))
      end do
      !$omp end parallel do
      ! Skip the drift solve for coefficients with negligible memory (their drift is
      ! negligible too): zero their drift instead of solving. thr is relative to the
      ! largest memory present.
      thr = self%skip_tol * maxval(self%mnorm)
      tor = self%nlam > NLAM
      if (tor .and. .not. (present(Wn_re) .and. present(Wn_im))) &
         error stop 'solve_drift: the toroidal channels are carried, so W targets are required'

      ! Solve for the drift, PARALLEL OVER DEGREE l so each per-degree operator
      ! ops(l) is touched by a single thread. Safe because EVERY degree solves
      ! through the re-entrant banded LU (vilma_band) on threadprivate scratch,
      ! including degree 1, whose dense KKT border merely widens the band (see
      ! vilma_radial_fe). There is no LIS solver any more, so no degree has to be
      ! serialized for re-entrancy. The scratch vectors are per-thread; dynamic schedule balances the rising work
      ! per degree (l+1 orders). Inactive (ordinary serial loop) unless openmp=1.
      !$omp parallel default(shared) private(l, k, node, fre, fim, xre, xim, gre, gim, wre, wim)
      allocate(fre(self%ndof), fim(self%ndof), xre(self%ndof), xim(self%ndof))
      if (tor) allocate(gre(self%nr), gim(self%nr), wre(self%nr), wim(self%nr))
      !$omp do schedule(dynamic)
      do l = 1, self%lmax
         do k = self%kbeg(l), self%kbeg(l+1) - 1
            if (self%mnorm(k) <= thr) then       ! negligible memory ⇒ negligible drift
               self%dUa(k) = (0.0_wp,0.0_wp);  self%dFa(k) = (0.0_wp,0.0_wp)
               self%dVa(k) = (0.0_wp,0.0_wp)
               Un_re(:,k) = 0.0_wp;  Un_im(:,k) = 0.0_wp
               Vn_re(:,k) = 0.0_wp;  Vn_im(:,k) = 0.0_wp
               if (tor) then
                  self%dWa(k) = (0.0_wp,0.0_wp)
                  Wn_re(:,k) = 0.0_wp;  Wn_im(:,k) = 0.0_wp
               end if
               cycle
            end if
            fre = 0.0_wp;  fim = 0.0_wp
            call dissipative_rhs(self%ne, self%r, self%sa(:,:,l), self%sb(:,:,l), &
                 self%sc(:,:,l), self%nrmc(:,l), self%Are(:,:,k), self%Bre(:,:,k), &
                 self%Cre(:,:,k), fre)
            call dissipative_rhs(self%ne, self%r, self%sa(:,:,l), self%sb(:,:,l), &
                 self%sc(:,:,l), self%nrmc(:,l), self%Aim(:,:,k), self%Bim(:,:,k), &
                 self%Cim(:,:,k), fim)
            call radial_operator_solve_vec(self%ops(l), fre, xre)
            call radial_operator_solve_vec(self%ops(l), fim, xim)
            if (l == 1 .and. self%deg1_cm) then
               ! The relaxation drift is a degree-1 state like any other and needs
               ! the same frame projection; its F(a) then vanishes by itself.
               call response_deg1_to_cm(self%ops(l), xre)
               call response_deg1_to_cm(self%ops(l), xim)
            end if
            self%dUa(k) = cmplx(xre(idx_u(self%nr)), xim(idx_u(self%nr)), wp)
            self%dFa(k) = cmplx(xre(idx_f(self%nr)), xim(idx_f(self%nr)), wp)
            self%dVa(k) = cmplx(xre(idx_v(self%nr)), xim(idx_v(self%nr)), wp)
            if (l == 1 .and. .not. self%deg1_cm) self%dFa(k) = (0.0_wp, 0.0_wp)   ! N₁≡0 (see init)
            do node = 1, self%nr
               Un_re(node,k) = xre(idx_u(node))
               Un_im(node,k) = xim(idx_u(node))
               Vn_re(node,k) = xre(idx_v(node))
               Vn_im(node,k) = xim(idx_v(node))
            end do
            if (tor) then
               ! The toroidal drift: its own forcing, its own operator. No frame
               ! projection at l = 1 — the net-rotation border already fixes it.
               gre = 0.0_wp;  gim = 0.0_wp
               call dissipative_rhs_tor(self%ne, self%r, self%sat(:,:,l), self%sbt(:,:,l), &
                    self%sct(:,:,l), self%nrmt(:,l), self%Are(:,:,k), self%Bre(:,:,k), &
                    self%Cre(:,:,k), gre)
               call dissipative_rhs_tor(self%ne, self%r, self%sat(:,:,l), self%sbt(:,:,l), &
                    self%sct(:,:,l), self%nrmt(:,l), self%Aim(:,:,k), self%Bim(:,:,k), &
                    self%Cim(:,:,k), gim)
               call toroidal_operator_solve_vec(self%tops(l), gre, wre)
               call toroidal_operator_solve_vec(self%tops(l), gim, wim)
               Wn_re(:,k) = wre;  Wn_im(:,k) = wim
               self%dWa(k) = cmplx(wre(self%nr), wim(self%nr), wp)
            end if
         end do
      end do
      !$omp end do
      deallocate(fre, fim, xre, xim)
      if (tor) deallocate(gre, gim, wre, wim)
      !$omp end parallel
      call system_clock(pc1)
      self%t_drift = self%t_drift + real(pc1-pc0,wp)/prate;  self%n_drift = self%n_drift + 1
   end subroutine solve_drift

   pure real(wp) function absmax6(n, x1, x2, x3, x4, x5, x6) result(amax)
      !! max|x| over six contiguous blocks of n values each (0 for n = 0): the six
      !! memory arrays of one (l,m) slot. Exact, so the scan order cannot move a bit.
      !! One loop over all six, not six loops: at low thread counts the scan is
      !! latency-bound, and six concurrent streams draw more bandwidth than one
      !! at a time. Callers pass each block's first element, A(1,1,k), not the
      !! section A(:,:,k), which ifx copied into a temporary first.
      integer,  intent(in) :: n
      real(wp), intent(in) :: x1(n), x2(n), x3(n), x4(n), x5(n), x6(n)
      integer :: i
      amax = 0.0_wp
      do i = 1, n
         amax = max(amax, abs(x1(i)), abs(x2(i)), abs(x3(i)), &
                          abs(x4(i)), abs(x5(i)), abs(x6(i)))
      end do
   end function absmax6

   pure real(wp) function absdiffmax6(n, x1, x2, x3, x4, x5, x6, y1, y2, y3, y4, y5, y6) &
         result(dmax)
      !! max|x − y| over six pairs of contiguous blocks, as absmax6.
      integer,  intent(in) :: n
      real(wp), intent(in) :: x1(n), x2(n), x3(n), x4(n), x5(n), x6(n)
      real(wp), intent(in) :: y1(n), y2(n), y3(n), y4(n), y5(n), y6(n)
      integer :: i
      dmax = 0.0_wp
      do i = 1, n
         dmax = max(dmax, abs(x1(i) - y1(i)), abs(x2(i) - y2(i)), abs(x3(i) - y3(i)), &
                          abs(x4(i) - y4(i)), abs(x5(i) - y5(i)), abs(x6(i) - y6(i)))
      end do
   end function absdiffmax6

   subroutine ve_response_apply(self, sht, sigma_lm, u_lm, n_lm)
      !! Affine response at the frozen time: u = gu(l)·σ + drift_U,
      !! N = gn(l)·σ − drift_F/g.
      type(response), intent(inout) :: self
      type(sht_grid),     intent(in)    :: sht
      complex(wp),        intent(in)    :: sigma_lm(:)
      complex(wp),        intent(out)   :: u_lm(:)
      complex(wp),        intent(out)   :: n_lm(:)
      integer :: k, l, lm, lm0

      u_lm = (0.0_wp,0.0_wp);  n_lm = (0.0_wp,0.0_wp)
      ! degree 0: monopole geoid, no deformation, no memory
      lm0 = sht_grid_lmidx(sht, 0, 0)
      u_lm(lm0) = self%gu(0)*sigma_lm(lm0)
      n_lm(lm0) = self%gn(0)*sigma_lm(lm0)
      ! degrees l>=1, in degree-grouped k order (gn(1)=0 and dFa(k)=0 give N₁≡0)
      do k = 1, self%nk
         l  = self%kdeg(k)
         lm = self%k2lm(k)
         u_lm(lm) = self%gu(l)*sigma_lm(lm) + self%dUa(k)
         n_lm(lm) = self%gn(l)*sigma_lm(lm) - self%dFa(k)/self%g
      end do
   end subroutine ve_response_apply

   subroutine ve_response_horizontal(self, sht, sigma_lm, v_lm)
      !! Spheroidal V at the frozen time: v_lm = gv(l)·σ + drift_V. Uses the same
      !! frozen drift (dVa) as the last begin_step, so calling it after a converged
      !! step gives the horizontal consistent with apply()'s u/N. Degree 1 left as
      !! solved (CE-like gauge, like u — not the geoid CM frame).
      type(response), intent(in)    :: self
      type(sht_grid),     intent(in)    :: sht
      complex(wp),        intent(in)    :: sigma_lm(:)
      complex(wp),        intent(out)   :: v_lm(:)
      integer :: k, l, lm
      v_lm = (0.0_wp,0.0_wp)               ! degree 0 has no horizontal (gv(0)=0)
      do k = 1, self%nk
         l  = self%kdeg(k)
         lm = self%k2lm(k)
         v_lm(lm) = self%gv(l)*sigma_lm(lm) + self%dVa(k)
      end do
   end subroutine ve_response_horizontal

   subroutine response_horizontal_toroidal(self, sht, t_lm)
      !! The toroidal part of the surface horizontal displacement, as coefficients
      !! t_lm of u_h = e_r × ∇₁(Σ t_lm Y_lm) (vilma_sht%tor_synthesis): the surface W(a)
      !! of the last begin_step. Zero for every response that carries no toroidal
      !! field — null, elastic, and a VE response with no 3-D element — which
      !! is exact, not an approximation: nothing forces W there.
      type(response), intent(in)  :: self
      type(sht_grid), intent(in)  :: sht
      complex(wp),    intent(out) :: t_lm(:)
      integer :: k
      t_lm = (0.0_wp, 0.0_wp)
      if (self%kind /= RESP_VE) return
      if (self%nlam == NLAM) return
      do k = 1, self%nk
         t_lm(self%k2lm(k)) = self%dWa(k)
      end do
   end subroutine response_horizontal_toroidal

   subroutine ve_response_commit(self, sht, sigma_lm)
      !! Advance the memory with the converged load, frozen σ (held/slow-load step;
      !! §3c part 3a). Explicit (FE): total nodal strain = σ·(unit-load nodal) +
      !! drift(τ_n), one Maxwell update per (l,m). Implicit (TRAP): the endpoint is
      !! solved by Picard iteration — re-solve the drift against the trial τ_{n+1},
      !! form the endpoint strain, trapezoid-advance from τ_n, repeat to couple_tol.
      !! Advances time by Δt. For fast-evolving loads the SLE driver instead iterates
      !! prepare_endpoint/advance_endpoint/finalize_step so σ co-converges (3b).
      type(response), intent(inout) :: self
      type(sht_grid),     intent(in)    :: sht
      complex(wp),        intent(in)    :: sigma_lm(:)
      real(wp) :: cnorm, snorm
      integer  :: iter

      if (.not. scheme_is_implicit(self%scheme)) then
         call vilma_advance(self, sht, sigma_lm)       ! explicit: byte-for-byte historical
         self%couple_iters_last = 1
         self%time = self%time + self%dt
         return
      end if

      ! --- implicit (TRAP): iterate the endpoint to a consistent τ_{n+1} ------------
      call ensure_commit_scratch(self)
      call snapshot_taun(self)                      ! τ_n base for every trapezoid pass
      call keep_drift(self)
      do iter = 1, self%max_couple_iter
         ! Endpoint drift from the current τ_{n+1} estimate (self%Are…); also refreshes
         ! self%dUa (the surface drift, used as the convergence signal). Then reset to
         ! τ_n and trapezoid-advance with (ε_n, ε_{n+1}).
         call solve_endpoint_drift(self, sht)
         call trapezoid_advance_all(self, sht, sigma_lm)
         ! iter 1 re-solves drift from τ_n and so reproduces begin_step's drift; the
         ! fixed point only moves at iter 2 (never exit on the first pass).
         call drift_change(self, cnorm, snorm)
         call keep_drift(self)
         if (iter >= 2 .and. cnorm <= self%couple_tol*max(snorm, tiny(1.0_wp))) exit
      end do
      self%couple_iters_last = min(iter, self%max_couple_iter)
      self%sigma_n = sigma_lm;  self%sigma_primed = .true.   ! σ_n for the next step's ε_n
      self%time = self%time + self%dt
   end subroutine ve_response_commit

   subroutine vilma_advance(self, sht, sigma_lm)
      !! Explicit forward-Euler memory advance: one Maxwell update per (l,m) from the
      !! report strain σ·(unit-load nodal) + drift(τ_n). The historical path, shared by
      !! commit_step and advance_endpoint so both stay byte-identical for FE. With
      !! laterally-varying viscosity the genuinely-3-D elements (e3d) advance pseudo-
      !! spectrally (advance_memory_3d); the laterally-uniform + elastic/fluid elements
      !! (active1d) advance on this cheap spectral path, masked by `active`.
      type(response), intent(inout) :: self
      type(sht_grid),     intent(in)    :: sht
      complex(wp),        intent(in)    :: sigma_lm(:)
      real(wp), allocatable :: Ure(:), Uim(:), Vre(:), Vim(:)
      real(wp) :: sre, sim
      integer  :: k, l, lm, node
      integer(kind=8) :: pc0, pc1, prate           ! PROFILE: memory-advance wall-clock
      call system_clock(pc0, prate)

      ! 3-D elements first (no-op when ne3d == 0, i.e. 1-D or laterally-uniform field).
      if (self%lat_visc) call advance_memory_3d(self, sht, sigma_lm)

      !$omp parallel default(shared) private(k, l, lm, node, sre, sim, Ure, Uim, Vre, Vim)
      allocate(Ure(self%nr), Uim(self%nr), Vre(self%nr), Vim(self%nr))
      !$omp do schedule(static)
      do k = 1, self%nk                          ! degree-grouped order
         l  = self%kdeg(k)
         lm = self%k2lm(k)
         sre = real(sigma_lm(lm), wp);  sim = aimag(sigma_lm(lm))
         do node = 1, self%nr
            Ure(node) = sre*self%xUn(node,l) + self%dUn_re(node,k)
            Uim(node) = sim*self%xUn(node,l) + self%dUn_im(node,k)
            Vre(node) = sre*self%xVn(node,l) + self%dVn_re(node,k)
            Vim(node) = sim*self%xVn(node,l) + self%dVn_im(node,k)
         end do
         call advance_memory(self%ne, self%mu, self%Mk, Ure, Vre, self%Jr(l), &
              self%Are(:,:,k), self%Bre(:,:,k), self%Cre(:,:,k), active=self%active1d)
         call advance_memory(self%ne, self%mu, self%Mk, Uim, Vim, self%Jr(l), &
              self%Aim(:,:,k), self%Bim(:,:,k), self%Cim(:,:,k), active=self%active1d)
         if (self%nlam > NLAM) then              ! W is pure drift: no σ·xWn term
            call advance_memory_tor(self%ne, self%mu, self%Mk, self%dWn_re(:,k), self%Jr(l), &
                 self%Are(:,:,k), self%Bre(:,:,k), self%Cre(:,:,k), active=self%active1d)
            call advance_memory_tor(self%ne, self%mu, self%Mk, self%dWn_im(:,k), self%Jr(l), &
                 self%Aim(:,:,k), self%Bim(:,:,k), self%Cim(:,:,k), active=self%active1d)
         end if
      end do
      !$omp end do
      deallocate(Ure, Uim, Vre, Vim)
      !$omp end parallel
      call system_clock(pc1)
      self%t_mem = self%t_mem + real(pc1-pc0,wp)/prate;  self%n_mem = self%n_mem + 1
   end subroutine vilma_advance

   subroutine response_enable_lateral_visc(self, sht, pert_elem)
      !! Rung 6 — turn on laterally-varying viscosity. `pert_elem` is the log10
      !! viscosity perturbation per element on the Gauss grid, (nphi,nlat,ne):
      !! η_eff(θ,φ) = η_radial · 10^pert, so the Maxwell rate scales by 10^(−pert).
      !! Elastic/fluid elements (MkPerDt = 0) stay memory-free regardless — the
      !! lithosphere remains exactly elastic.
      !!
      !! 1-D/3-D split (VILMA1 mod_visc3d): an element is flagged genuinely 3-D only
      !! when its lateral log10(η) spread exceeds visc3d_tol; that subset (e3d) pays
      !! the pseudo-spectral tensor-SH advance. Every other Maxwell element collapses
      !! to a scalar effective rate (its lateral MEAN) and advances on the cheap
      !! degree-diagonal spectral path. A laterally-uniform field therefore flags NO
      !! element 3-D, so the cost equals the 1-D run (exactly, not just to SHT round-off).
      type(response), intent(inout) :: self
      type(sht_grid),     intent(in)    :: sht
      real(wp),           intent(in)    :: pert_elem(:,:,:)   !! (nphi,nlat,ne) log10 η perturbation
      integer :: e, ne_grid
      real(wp) :: spread, mean_pert
      if (size(pert_elem,1) /= sht%nphi .or. size(pert_elem,2) /= sht%nlat .or. &
          size(pert_elem,3) /= self%ne) &
         error stop 'enable_lateral_visc: pert_elem must be (nphi,nlat,ne)'
      call tensor_sh_init(self%tsh, sht)
      if (allocated(self%Mk3))      deallocate(self%Mk3)
      if (allocated(self%MkPerDt3)) deallocate(self%MkPerDt3)
      if (allocated(self%e3d))      deallocate(self%e3d)
      allocate(self%MkPerDt3(sht%nphi, sht%nlat, self%ne))
      allocate(self%Mk3(sht%nphi, sht%nlat, self%ne))
      ne_grid = sht%nphi*sht%nlat
      ! Build the per-grid rate from the ORIGINAL radial MkPerDt, THEN classify and (for
      ! 1-D-effective elements) overwrite the scalar MkPerDt with the lateral-mean rate.
      do e = 1, self%ne
         self%MkPerDt3(:,:,e) = self%MkPerDt(e) * 10.0_wp**(-pert_elem(:,:,e))
      end do
      self%active1d = .true.;  self%ne3d = 0
      do e = 1, self%ne
         if (self%MkPerDt(e) == 0.0_wp) cycle             ! elastic/fluid: spectral, memory-free
         spread = maxval(pert_elem(:,:,e)) - minval(pert_elem(:,:,e))
         if (spread > self%visc3d_tol) then
            self%active1d(e) = .false.;  self%ne3d = self%ne3d + 1   ! genuinely 3-D
         else
            mean_pert = sum(pert_elem(:,:,e))/real(ne_grid, wp)      ! collapse to scalar
            self%MkPerDt(e) = self%MkPerDt(e) * 10.0_wp**(-mean_pert)
         end if
      end do
      self%Mk  = self%MkPerDt * self%dt                   ! rescale 1-D scalar rates
      self%Mk3 = self%MkPerDt3 * self%dt
      allocate(self%e3d(self%ne3d))
      self%ne3d = 0
      do e = 1, self%ne
         if (.not. self%active1d(e)) then
            self%ne3d = self%ne3d + 1;  self%e3d(self%ne3d) = e
         end if
      end do
      self%lat_visc = .true.
      ! A genuinely 3-D element is what couples spheroidal and toroidal. Once the
      ! toroidal channels are on they stay on: a later laterally uniform field
      ! (the pre-spinup's mean field, say) no longer forces W, but whatever W and
      ! toroidal memory exist must still relax, not be dropped.
      if (self%toroidal .and. self%ne3d > 0 .and. self%nlam == NLAM) call enable_toroidal(self)
   end subroutine response_enable_lateral_visc

   subroutine enable_toroidal(self)
      !! Carry the toroidal degree of freedom from here on: per-degree W operators
      !! and constants, the W drift arrays, and every memory array and state buffer
      !! widened from NLAM to NLAM+NLAM_TOR channels with the new ones zeroed —
      !! the exact toroidal state of a run that has only ever been spheroidal.
      type(response), intent(inout) :: self
      integer :: l
      allocate(self%tops(self%lmax), self%nrmt(NLAM_TOR,self%lmax))
      allocate(self%sat(2,NLAM_TOR,self%lmax), self%sbt(2,NLAM_TOR,self%lmax), &
               self%sct(2,NLAM_TOR,self%lmax))
      do l = 1, self%lmax
         call toroidal_operator_assemble(self%tops(l), self%r, self%mu, l)
         call ve_strain_constants_tor(self%Jr(l), self%nrmt(:,l), &
                                      self%sat(:,:,l), self%sbt(:,:,l), self%sct(:,:,l))
      end do
      allocate(self%dWa(self%nk));  self%dWa = (0.0_wp, 0.0_wp)
      allocate(self%dWn_re(self%nr,self%nk), self%dWn_im(self%nr,self%nk))
      self%dWn_re = 0.0_wp;  self%dWn_im = 0.0_wp
      self%nlam = NLAM + NLAM_TOR
      call widen(self%Are);  call widen(self%Aim)
      call widen(self%Bre);  call widen(self%Bim)
      call widen(self%Cre);  call widen(self%Cim)
      if (allocated(self%Are0)) then
         call widen(self%Are0);  call widen(self%Aim0)
         call widen(self%Bre0);  call widen(self%Bim0)
         call widen(self%Cre0);  call widen(self%Cim0)
         call ensure_commit_scratch_tor(self)
      end if
      if (allocated(self%Are_s)) then
         call widen(self%Are_s);  call widen(self%Aim_s)
         call widen(self%Bre_s);  call widen(self%Bim_s)
         call widen(self%Cre_s);  call widen(self%Cim_s)
         call widen(self%Are_c);  call widen(self%Aim_c)
         call widen(self%Bre_c);  call widen(self%Bim_c)
         call widen(self%Cre_c);  call widen(self%Cim_c)
      end if
   contains
      subroutine widen(x)
         !! (NLAM,ne,nk) → (NLAM+NLAM_TOR,ne,nk), spheroidal channels kept, the
         !! new ones zero. Filled in parallel on the same schedule(static)
         !! partition over k as every loop that reads these arrays: this is their
         !! first touch, which places the pages (see response_init_ve).
         real(wp), allocatable, intent(inout) :: x(:,:,:)
         real(wp), allocatable :: y(:,:,:)
         integer :: k
         allocate(y(NLAM+NLAM_TOR, size(x,2), size(x,3)))
         !$omp parallel do default(shared) private(k) schedule(static)
         do k = 1, size(x,3)
            y(1:NLAM,:,k) = x(:,:,k)
            y(NLAM+1:,:,k) = 0.0_wp
         end do
         !$omp end parallel do
         call move_alloc(y, x)
      end subroutine widen
   end subroutine enable_toroidal

   subroutine ensure_commit_scratch_tor(self)
      !! The implicit commit's W endpoint drift and convergence buffer.
      type(response), intent(inout) :: self
      if (allocated(self%edWn_re)) return
      allocate(self%edWn_re(self%nr,self%nk), self%edWn_im(self%nr,self%nk))
      allocate(self%dWa_prev(self%nk))
      self%edWn_re = 0.0_wp;  self%edWn_im = 0.0_wp;  self%dWa_prev = (0.0_wp,0.0_wp)
   end subroutine ensure_commit_scratch_tor

   subroutine response_enable_lateral_visc_from_nodes(self, sht, visc_node, lid_depth, lid_log10max, &
                                                      log10_cap)
      !! Rung 6c — enable laterally-varying viscosity from a NODE-based ABSOLUTE
      !! log10(η) field on the Gauss grid, visc_node(nphi*nlat, nr) (as produced by
      !! vilma_read_visc_3d). Bridges node→element by the log10-mean of the two
      !! bracketing nodes (geometric mean of η), forms the per-element log10
      !! perturbation against the element's radial reference viscosity
      !! η_radial(e) = μ(e)/MkPerDt(e), and calls enable_lateral_visc. Elastic/
      !! fluid elements (MkPerDt=0) keep pert=0 — irrelevant, they stay memory-free.
      !!
      !! Lid rule (optional, all three arguments together): a Maxwell element lying
      !! wholly above lid_depth [m] below the surface whose log10 η exceeds
      !! lid_log10max is set to log10_cap. It acts per ELEMENT, after the node→element
      !! bridging: a node-level cap would change the node at the lid base, which the
      !! element below shares, and leave intermediate (slow) elements on both sides.
      type(response), intent(inout) :: self
      type(sht_grid),     intent(in)    :: sht
      real(wp),           intent(in)    :: visc_node(:,:)   !! (nphi*nlat, nr) log10(η)
      real(wp), optional, intent(in)    :: lid_depth, lid_log10max, log10_cap
      real(wp), allocatable :: pert(:,:,:)
      real(wp) :: elem_abs, logeta_ref, r_lid
      logical  :: in_lid
      integer  :: e, i, j, sp
      if (size(visc_node,1) /= sht%nphi*sht%nlat .or. size(visc_node,2) /= self%nr) &
         error stop 'enable_lateral_visc_from_nodes: visc_node must be (nphi*nlat, nr)'
      if (present(lid_depth) .neqv. (present(lid_log10max) .and. present(log10_cap))) &
         error stop 'enable_lateral_visc_from_nodes: lid_depth, lid_log10max, log10_cap go together'
      r_lid = huge(1.0_wp)                                   ! no element is above it: rule off
      if (present(lid_depth)) then
         if (lid_depth > 0.0_wp) r_lid = self%r(self%nr) - lid_depth
      end if
      allocate(pert(sht%nphi, sht%nlat, self%ne));  pert = 0.0_wp
      do e = 1, self%ne
         if (self%MkPerDt(e) == 0.0_wp) cycle             ! elastic/fluid: stay as-is
         logeta_ref = log10(self%mu(e)/self%MkPerDt(e))    ! log10 η_radial(e) (μ/MkPerDt)
         in_lid = self%r(e) >= r_lid - 1.0e-6_wp*self%r(self%nr)   ! element wholly above the lid base
         do j = 1, sht%nlat
            do i = 1, sht%nphi
               sp = i + (j-1)*sht%nphi
               elem_abs = 0.5_wp*(visc_node(sp,e) + visc_node(sp,e+1))   ! log10-mean of nodes
               if (in_lid) then
                  if (elem_abs > lid_log10max) elem_abs = log10_cap
               end if
               pert(i,j,e) = elem_abs - logeta_ref
            end do
         end do
      end do
      call response_enable_lateral_visc(self, sht, pert)
      deallocate(pert)
   end subroutine response_enable_lateral_visc_from_nodes

   function response_radial_rate(self) result(rate)
      !! The radial profile of the viscosity the response integrates, as the per-element
      !! Maxwell rate μ/η (ne). A 1-D element gives its scalar rate, which after a
      !! lateral-viscosity enable is already the lateral geometric mean; a genuinely
      !! 3-D element gives the lateral geometric mean of its grid rate (the unweighted
      !! grid mean of log10 η, as the 1-D collapse in enable_lateral_visc takes it).
      !! Elastic/fluid elements stay 0. This is what a radially symmetric consumer
      !! (the degree-2 rotation channels) should see.
      type(response), intent(in) :: self
      real(wp) :: rate(self%ne)
      integer  :: k, e
      rate = self%MkPerDt
      if (.not. self%lat_visc) return
      do k = 1, self%ne3d
         e = self%e3d(k)
         rate(e) = 10.0_wp**(sum(log10(self%MkPerDt3(:,:,e))) &
                             / real(size(self%MkPerDt3,1)*size(self%MkPerDt3,2), wp))
      end do
   end function response_radial_rate

   subroutine advance_memory_3d(self, sht, sigma_lm)
      !! Tensor-correct pseudo-spectral FE memory advance for laterally-varying
      !! viscosity (rung 6, general order). The Maxwell update τ⁺=(1−M)τ−2μM·ε is
      !! pointwise in PHYSICAL space, so per radial element and per radial shape-
      !! coefficient (A,B,C) the memory and strain TENSORS are reconstructed on the
      !! Gauss grid via their six dyadic components (vilma_tensor_sh; Martinec 2000
      !! B10/B11), advanced pointwise with the lateral field M(θ,φ), and projected
      !! back. With a uniform M the dyadic round trip is the identity ⇒ reduces to the
      !! 1-D advance.
      !!
      !! Parallel over elements: each element's dyadic transforms run on the calling
      !! thread's PRIVATE SHTns config (tsh%thread_cfg) — a single config is not safe
      !! for concurrent calls, but the per-thread pool (vilma_tensor_sh) makes the element
      !! loop embarrassingly parallel. Per-thread coeff/grid scratch is allocated
      !! inside the region; the memory writeback touches a distinct element per
      !! iteration, so there is no race.
      type(response), intent(inout) :: self
      type(sht_grid),     intent(in)    :: sht
      complex(wp),        intent(in)    :: sigma_lm(:)
      complex(wp), allocatable :: cma(:,:), cmb(:,:), cmc(:,:)   ! memory coeffs (nlam,nlm)
      complex(wp), allocatable :: cea(:,:), ceb(:,:), cec(:,:)   ! strain coeffs (nlam,nlm)
      complex(wp), allocatable :: cdel(:,:)                       ! analysed increment (nlam,nlm)
      real(wp),    allocatable :: dtau(:,:,:), deps(:,:,:)        ! (nphi,nlat,6)
      type(c_ptr) :: cfg
      integer  :: e, ei, k, lm, nc

      if (self%ne3d == 0) return        ! no genuinely-3-D element ⇒ all handled spectrally
      !$omp parallel default(shared) &
      !$omp   private(e, ei, k, lm, nc, cma, cmb, cmc, cea, ceb, cec, cdel, dtau, deps, cfg)
      ! Block width = the channels carried: TLAM_SPH, or TLAM once toroidal (the
      ! local channel orders of vilma_tensor_sh and vilma_viscoelastic agree, λ3,4 at 5,6).
      allocate(cma(self%nlam,sht%nlm), cmb(self%nlam,sht%nlm), cmc(self%nlam,sht%nlm))
      allocate(cea(self%nlam,sht%nlm), ceb(self%nlam,sht%nlm), cec(self%nlam,sht%nlm))
      allocate(cdel(self%nlam,sht%nlm))
      allocate(dtau(sht%nphi,sht%nlat,6), deps(sht%nphi,sht%nlat,6))
      cfg = tensor_sh_thread_cfg(self%tsh)                          ! this thread's private config
      !$omp do schedule(dynamic)
      do ei = 1, self%ne3d                                 ! only the genuinely-3-D elements
         e = self%e3d(ei)
         call gather_tensor_coeffs(self, sigma_lm, e, cma, cmb, cmc, cea, ceb, cec)
         call advance_shape_tensor(self, sht, e, cma, cea, dtau, deps, cdel, cfg)
         call advance_shape_tensor(self, sht, e, cmb, ceb, dtau, deps, cdel, cfg)
         call advance_shape_tensor(self, sht, e, cmc, cec, dtau, deps, cdel, cfg)
         do k = 1, self%nk                                 ! write updated memory back
            lm = self%k2lm(k)
            nc = size(cma,1)
            self%Are(1:nc,e,k) = real(cma(:,lm), wp);  self%Aim(1:nc,e,k) = aimag(cma(:,lm))
            self%Bre(1:nc,e,k) = real(cmb(:,lm), wp);  self%Bim(1:nc,e,k) = aimag(cmb(:,lm))
            self%Cre(1:nc,e,k) = real(cmc(:,lm), wp);  self%Cim(1:nc,e,k) = aimag(cmc(:,lm))
         end do
      end do
      !$omp end do
      deallocate(cma, cmb, cmc, cea, ceb, cec, cdel, dtau, deps)
      !$omp end parallel
   end subroutine advance_memory_3d

   subroutine gather_tensor_coeffs(self, sigma_lm, e, cma, cmb, cmc, cea, ceb, cec)
      !! Per element e, gather the memory shape-coeffs (Are/Aim …) and the current
      !! strain shape-coeffs (strain_coeffs of σ·xUn + drift, exactly as vilma_advance)
      !! into complex (nlam, nlm) blocks for the dyadic transform — all (l,m). With
      !! the toroidal channels carried, channels 5,6 take the λ3,4 memory and the
      !! strain of the nodal W drift (no load term: eq 84 does not force W).
      type(response), intent(in)  :: self
      complex(wp),        intent(in)  :: sigma_lm(:)
      integer,            intent(in)  :: e
      complex(wp),        intent(out) :: cma(:,:), cmb(:,:), cmc(:,:), cea(:,:), ceb(:,:), cec(:,:)
      real(wp) :: ar(NLAM), br(NLAM), cr(NLAM), ai(NLAM), bi(NLAM), ci(NLAM)
      real(wp) :: at(NLAM_TOR), bt(NLAM_TOR), ct(NLAM_TOR), ati(NLAM_TOR), bti(NLAM_TOR), cti(NLAM_TOR)
      real(wp) :: sre, sim, Ur, Ur1, Vr, Vr1, Ui, Ui1, Vi, Vi1
      integer  :: k, l, lm, lam
      cma = (0.0_wp,0.0_wp); cmb = (0.0_wp,0.0_wp); cmc = (0.0_wp,0.0_wp)
      cea = (0.0_wp,0.0_wp); ceb = (0.0_wp,0.0_wp); cec = (0.0_wp,0.0_wp)
      do k = 1, self%nk
         l  = self%kdeg(k);  lm = self%k2lm(k)
         sre = real(sigma_lm(lm), wp);  sim = aimag(sigma_lm(lm))
         Ur  = sre*self%xUn(e,  l) + self%dUn_re(e,  k);  Ui  = sim*self%xUn(e,  l) + self%dUn_im(e,  k)
         Ur1 = sre*self%xUn(e+1,l) + self%dUn_re(e+1,k);  Ui1 = sim*self%xUn(e+1,l) + self%dUn_im(e+1,k)
         Vr  = sre*self%xVn(e,  l) + self%dVn_re(e,  k);  Vi  = sim*self%xVn(e,  l) + self%dVn_im(e,  k)
         Vr1 = sre*self%xVn(e+1,l) + self%dVn_re(e+1,k);  Vi1 = sim*self%xVn(e+1,l) + self%dVn_im(e+1,k)
         call strain_coeffs(Ur, Ur1, Vr, Vr1, self%Jr(l), ar, br, cr)
         call strain_coeffs(Ui, Ui1, Vi, Vi1, self%Jr(l), ai, bi, ci)
         do lam = 1, NLAM
            cea(lam,lm) = cmplx(ar(lam), ai(lam), wp)
            ceb(lam,lm) = cmplx(br(lam), bi(lam), wp)
            cec(lam,lm) = cmplx(cr(lam), ci(lam), wp)
            cma(lam,lm) = cmplx(self%Are(lam,e,k), self%Aim(lam,e,k), wp)
            cmb(lam,lm) = cmplx(self%Bre(lam,e,k), self%Bim(lam,e,k), wp)
            cmc(lam,lm) = cmplx(self%Cre(lam,e,k), self%Cim(lam,e,k), wp)
         end do
         if (self%nlam == NLAM) cycle
         call strain_coeffs_tor(self%dWn_re(e,k), self%dWn_re(e+1,k), self%Jr(l), at, bt, ct)
         call strain_coeffs_tor(self%dWn_im(e,k), self%dWn_im(e+1,k), self%Jr(l), ati, bti, cti)
         do lam = 1, NLAM_TOR
            cea(NLAM+lam,lm) = cmplx(at(lam), ati(lam), wp)
            ceb(NLAM+lam,lm) = cmplx(bt(lam), bti(lam), wp)
            cec(NLAM+lam,lm) = cmplx(ct(lam), cti(lam), wp)
            cma(NLAM+lam,lm) = cmplx(self%Are(NLAM+lam,e,k), self%Aim(NLAM+lam,e,k), wp)
            cmb(NLAM+lam,lm) = cmplx(self%Bre(NLAM+lam,e,k), self%Bim(NLAM+lam,e,k), wp)
            cmc(NLAM+lam,lm) = cmplx(self%Cre(NLAM+lam,e,k), self%Cim(NLAM+lam,e,k), wp)
         end do
      end do
   end subroutine gather_tensor_coeffs

   subroutine advance_shape_tensor(self, sht, e, c, eps, dtau, deps, cdel, cfg)
      !! One radial shape-coefficient: reconstruct the memory τ and strain ε tensors
      !! on the grid (six dyadic components), and apply τ⁺=(1−M)τ−2μM·ε pointwise per
      !! component with the lateral field M(θ,φ)=Mk3(:,:,e). c is updated in place;
      !! dtau/deps/cdel are caller-provided scratch; cfg is the calling thread's
      !! private SHTns config (for the parallel element loop).
      !!
      !! THE UPDATE IS APPLIED AS AN INCREMENT, not as a replacement. Writing it as
      !!       τ⁺ = τ − M·(τ + 2μ·ε)
      !! and transforming only the increment Δ = −M·(τ + 2μ·ε) is algebraically the
      !! same, but numerically it is not, and the difference is not small.
      !!
      !! The dyadic analysis is exact to ~8e-14 RELATIVE TO THE FIELD IT ANALYSES. A
      !! step only changes τ by M·(τ + 2με), so analysing τ⁺ itself carries an error
      !! of 8e-14·|τ| into a quantity whose true size is M·|τ| — a relative error of
      !! 8e-14/M in the increment. Measured (tests/diag_tensor_grid.f90, LOG session
      !! 32k): 2e-12 at M = 4e-2, 7.9e-5 at M = 1e-9, 7.8e-1 at M = 1e-13. Block D's
      !! viscosity reaches 1e30 Pa s (cratonic lid and lithosphere), which with
      !! μ ≈ 5.7e10 and Δt ≈ 50 yr puts the stiffest Maxwell elements at M ≈ 1e-10:
      !! their memory advance was wrong by percent PER SUB-STEP, and the whole model
      !! drifted 0.50 m rms from the 1-D path on an identical laterally-uniform field
      !! — 57 % of the lateral-viscosity signal it exists to compute.
      !!
      !! Analysing Δ instead puts the same 8e-14 on Δ, so the relative error is 8e-14
      !! at every M. That matches the 1-D path (advance_memory), which does this same
      !! near-total cancellation in coefficient space per mode with no transform in
      !! between and is therefore exact. Cost is unchanged: two synths, one analysis.
      type(response), intent(in)    :: self
      type(sht_grid),     intent(in)    :: sht
      integer,            intent(in)    :: e
      complex(wp),        intent(inout) :: c(:,:)
      complex(wp),        intent(in)    :: eps(:,:)
      real(wp),           intent(inout) :: dtau(:,:,:), deps(:,:,:)
      complex(wp),        intent(inout) :: cdel(:,:)
      type(c_ptr),        intent(in)    :: cfg
      real(wp) :: twoMu
      integer  :: p
      twoMu = 2.0_wp*self%mu(e)
      call tensor_sh_synth(self%tsh, sht, c,   dtau, cfg)
      call tensor_sh_synth(self%tsh, sht, eps, deps, cfg)
      do p = 1, 6
         dtau(:,:,p) = -self%Mk3(:,:,e)*(dtau(:,:,p) + twoMu*deps(:,:,p))
      end do
      call tensor_sh_analysis(self%tsh, sht, dtau, cdel, cfg)
      c = c + cdel
   end subroutine advance_shape_tensor

   subroutine advance_memory_3d_trap(self, sht, sigma_lm)
      !! Trapezoidal (Crank–Nicolson) pseudo-spectral memory advance for laterally-
      !! varying viscosity — the 3D analogue of trapezoid_advance_all. One endpoint
      !! pass: reset to τ_n (the *0 snapshot) and advance per radial shape-coefficient
      !! with ε_n (σ_n·xUn + dUn) and ε_{n+1} (σ_{n+1}·xUn + edUn), applying the
      !! pointwise rule τ⁺ = [(1−M/2)τ_n − μM(ε_n+ε_{n+1})]/(1+M/2) on the Gauss grid
      !! with the lateral field M=Mk3(θ,φ). With uniform M the dyadic round trip is the
      !! identity ⇒ reduces to the 1-D trapezoidal advance per (l,m). Parallel over
      !! elements on per-thread configs, exactly like the FE advance_memory_3d.
      type(response), intent(inout) :: self
      type(sht_grid),     intent(in)    :: sht
      complex(wp),        intent(in)    :: sigma_lm(:)              ! σ_{n+1}
      complex(wp), allocatable :: cm0a(:,:), cm0b(:,:), cm0c(:,:)   ! τ_n coeffs (nlam,nlm)
      complex(wp), allocatable :: cna(:,:),  cnb(:,:),  cnc(:,:)    ! ε_n coeffs
      complex(wp), allocatable :: c1a(:,:),  c1b(:,:),  c1c(:,:)    ! ε_{n+1} coeffs
      complex(wp), allocatable :: cdel(:,:)                         ! analysed increment
      real(wp),    allocatable :: dt0(:,:,:), den(:,:,:), de1(:,:,:) ! (nphi,nlat,6)
      type(c_ptr) :: cfg
      integer  :: e, ei, k, lm, nc

      if (self%ne3d == 0) return        ! no genuinely-3-D element ⇒ all handled spectrally
      !$omp parallel default(shared) &
      !$omp   private(e, ei, k, lm, nc, cm0a, cm0b, cm0c, cna, cnb, cnc, c1a, c1b, c1c, cdel, dt0, den, de1, cfg)
      allocate(cm0a(self%nlam,sht%nlm), cm0b(self%nlam,sht%nlm), cm0c(self%nlam,sht%nlm))
      allocate(cna(self%nlam,sht%nlm),  cnb(self%nlam,sht%nlm),  cnc(self%nlam,sht%nlm))
      allocate(c1a(self%nlam,sht%nlm),  c1b(self%nlam,sht%nlm),  c1c(self%nlam,sht%nlm))
      allocate(cdel(self%nlam,sht%nlm))
      allocate(dt0(sht%nphi,sht%nlat,6), den(sht%nphi,sht%nlat,6), de1(sht%nphi,sht%nlat,6))
      cfg = tensor_sh_thread_cfg(self%tsh)
      !$omp do schedule(dynamic)
      do ei = 1, self%ne3d                                 ! only the genuinely-3-D elements
         e = self%e3d(ei)
         call gather_tensor_coeffs_trap(self, sigma_lm, e, cm0a, cm0b, cm0c, &
                                        cna, cnb, cnc, c1a, c1b, c1c)
         call advance_shape_tensor_trap(self, sht, e, cm0a, cna, c1a, dt0, den, de1, cdel, cfg)
         call advance_shape_tensor_trap(self, sht, e, cm0b, cnb, c1b, dt0, den, de1, cdel, cfg)
         call advance_shape_tensor_trap(self, sht, e, cm0c, cnc, c1c, dt0, den, de1, cdel, cfg)
         do k = 1, self%nk                                 ! write updated memory back
            lm = self%k2lm(k)
            nc = size(cm0a,1)
            self%Are(1:nc,e,k) = real(cm0a(:,lm), wp);  self%Aim(1:nc,e,k) = aimag(cm0a(:,lm))
            self%Bre(1:nc,e,k) = real(cm0b(:,lm), wp);  self%Bim(1:nc,e,k) = aimag(cm0b(:,lm))
            self%Cre(1:nc,e,k) = real(cm0c(:,lm), wp);  self%Cim(1:nc,e,k) = aimag(cm0c(:,lm))
         end do
      end do
      !$omp end do
      deallocate(cm0a, cm0b, cm0c, cna, cnb, cnc, c1a, c1b, c1c, cdel, dt0, den, de1)
      !$omp end parallel
   end subroutine advance_memory_3d_trap

   subroutine gather_tensor_coeffs_trap(self, sigma_lm, e, cm0a, cm0b, cm0c, &
                                        cna, cnb, cnc, c1a, c1b, c1c)
      !! Per element e, gather for the trapezoidal advance: τ_n (the *0 snapshot), the
      !! start strain ε_n (σ_n·xUn + dUn) and the endpoint strain ε_{n+1} (σ_{n+1}·xUn +
      !! edUn), each as complex (nlam,nlm) blocks. σ_n is sigma_n when primed, else the
      !! first-step fallback σ_{n+1} (matching trapezoid_advance_all). The toroidal
      !! channels, when carried, take W from dWn (ε_n) and edWn (ε_{n+1}).
      type(response), intent(in)  :: self
      complex(wp),        intent(in)  :: sigma_lm(:)              ! σ_{n+1}
      integer,            intent(in)  :: e
      complex(wp),        intent(out) :: cm0a(:,:), cm0b(:,:), cm0c(:,:)
      complex(wp),        intent(out) :: cna(:,:), cnb(:,:), cnc(:,:)
      complex(wp),        intent(out) :: c1a(:,:), c1b(:,:), c1c(:,:)
      real(wp) :: arn(NLAM), brn(NLAM), crn(NLAM), ain(NLAM), bin(NLAM), cin(NLAM)
      real(wp) :: ar1(NLAM), br1(NLAM), cr1(NLAM), ai1(NLAM), bi1(NLAM), ci1(NLAM)
      real(wp) :: atn(NLAM_TOR), btn(NLAM_TOR), ctn(NLAM_TOR), atin(NLAM_TOR), btin(NLAM_TOR), ctin(NLAM_TOR)
      real(wp) :: at1(NLAM_TOR), bt1(NLAM_TOR), ct1(NLAM_TOR), ati1(NLAM_TOR), bti1(NLAM_TOR), cti1(NLAM_TOR)
      real(wp) :: sre, sim, srn, sin
      real(wp) :: Urn, Urn1, Vrn, Vrn1, Uin, Uin1, Vin, Vin1
      real(wp) :: Ur1, Ur11, Vr1, Vr11, Ui1, Ui11, Vi1, Vi11
      integer  :: k, l, lm, lam
      cm0a = (0.0_wp,0.0_wp); cm0b = (0.0_wp,0.0_wp); cm0c = (0.0_wp,0.0_wp)
      cna  = (0.0_wp,0.0_wp); cnb  = (0.0_wp,0.0_wp); cnc  = (0.0_wp,0.0_wp)
      c1a  = (0.0_wp,0.0_wp); c1b  = (0.0_wp,0.0_wp); c1c  = (0.0_wp,0.0_wp)
      do k = 1, self%nk
         l  = self%kdeg(k);  lm = self%k2lm(k)
         sre = real(sigma_lm(lm), wp);  sim = aimag(sigma_lm(lm))          ! σ_{n+1}
         if (self%sigma_primed) then
            srn = real(self%sigma_n(lm), wp);  sin = aimag(self%sigma_n(lm))   ! σ_n
         else
            srn = sre;  sin = sim     ! first step: ε_n uses σ_{n+1} (load present at t=0)
         end if
         ! ε_n nodal (σ_n·xUn + begin_step drift dUn)
         Urn  = srn*self%xUn(e,  l) + self%dUn_re(e,  k);  Uin  = sin*self%xUn(e,  l) + self%dUn_im(e,  k)
         Urn1 = srn*self%xUn(e+1,l) + self%dUn_re(e+1,k);  Uin1 = sin*self%xUn(e+1,l) + self%dUn_im(e+1,k)
         Vrn  = srn*self%xVn(e,  l) + self%dVn_re(e,  k);  Vin  = sin*self%xVn(e,  l) + self%dVn_im(e,  k)
         Vrn1 = srn*self%xVn(e+1,l) + self%dVn_re(e+1,k);  Vin1 = sin*self%xVn(e+1,l) + self%dVn_im(e+1,k)
         ! ε_{n+1} nodal (σ_{n+1}·xUn + endpoint drift edUn)
         Ur1  = sre*self%xUn(e,  l) + self%edUn_re(e,  k);  Ui1  = sim*self%xUn(e,  l) + self%edUn_im(e,  k)
         Ur11 = sre*self%xUn(e+1,l) + self%edUn_re(e+1,k);  Ui11 = sim*self%xUn(e+1,l) + self%edUn_im(e+1,k)
         Vr1  = sre*self%xVn(e,  l) + self%edVn_re(e,  k);  Vi1  = sim*self%xVn(e,  l) + self%edVn_im(e,  k)
         Vr11 = sre*self%xVn(e+1,l) + self%edVn_re(e+1,k);  Vi11 = sim*self%xVn(e+1,l) + self%edVn_im(e+1,k)
         call strain_coeffs(Urn, Urn1, Vrn, Vrn1, self%Jr(l), arn, brn, crn)
         call strain_coeffs(Uin, Uin1, Vin, Vin1, self%Jr(l), ain, bin, cin)
         call strain_coeffs(Ur1, Ur11, Vr1, Vr11, self%Jr(l), ar1, br1, cr1)
         call strain_coeffs(Ui1, Ui11, Vi1, Vi11, self%Jr(l), ai1, bi1, ci1)
         do lam = 1, NLAM
            cna(lam,lm) = cmplx(arn(lam), ain(lam), wp)
            cnb(lam,lm) = cmplx(brn(lam), bin(lam), wp)
            cnc(lam,lm) = cmplx(crn(lam), cin(lam), wp)
            c1a(lam,lm) = cmplx(ar1(lam), ai1(lam), wp)
            c1b(lam,lm) = cmplx(br1(lam), bi1(lam), wp)
            c1c(lam,lm) = cmplx(cr1(lam), ci1(lam), wp)
            cm0a(lam,lm) = cmplx(self%Are0(lam,e,k), self%Aim0(lam,e,k), wp)
            cm0b(lam,lm) = cmplx(self%Bre0(lam,e,k), self%Bim0(lam,e,k), wp)
            cm0c(lam,lm) = cmplx(self%Cre0(lam,e,k), self%Cim0(lam,e,k), wp)
         end do
         if (self%nlam == NLAM) cycle
         call strain_coeffs_tor(self%dWn_re(e,k),  self%dWn_re(e+1,k),  self%Jr(l), atn, btn, ctn)
         call strain_coeffs_tor(self%dWn_im(e,k),  self%dWn_im(e+1,k),  self%Jr(l), atin, btin, ctin)
         call strain_coeffs_tor(self%edWn_re(e,k), self%edWn_re(e+1,k), self%Jr(l), at1, bt1, ct1)
         call strain_coeffs_tor(self%edWn_im(e,k), self%edWn_im(e+1,k), self%Jr(l), ati1, bti1, cti1)
         do lam = 1, NLAM_TOR
            cna(NLAM+lam,lm) = cmplx(atn(lam), atin(lam), wp)
            cnb(NLAM+lam,lm) = cmplx(btn(lam), btin(lam), wp)
            cnc(NLAM+lam,lm) = cmplx(ctn(lam), ctin(lam), wp)
            c1a(NLAM+lam,lm) = cmplx(at1(lam), ati1(lam), wp)
            c1b(NLAM+lam,lm) = cmplx(bt1(lam), bti1(lam), wp)
            c1c(NLAM+lam,lm) = cmplx(ct1(lam), cti1(lam), wp)
            cm0a(NLAM+lam,lm) = cmplx(self%Are0(NLAM+lam,e,k), self%Aim0(NLAM+lam,e,k), wp)
            cm0b(NLAM+lam,lm) = cmplx(self%Bre0(NLAM+lam,e,k), self%Bim0(NLAM+lam,e,k), wp)
            cm0c(NLAM+lam,lm) = cmplx(self%Cre0(NLAM+lam,e,k), self%Cim0(NLAM+lam,e,k), wp)
         end do
      end do
   end subroutine gather_tensor_coeffs_trap

   subroutine advance_shape_tensor_trap(self, sht, e, c0, eps_n, eps_1, dt0, den, de1, cdel, cfg)
      !! One radial shape-coefficient, trapezoidal: reconstruct τ_n, ε_n, ε_{n+1} on the
      !! grid (six dyadic components) and apply the pointwise Crank–Nicolson update
      !! τ⁺ = [(1−M/2)τ_n − μM(ε_n+ε_{n+1})]/(1+M/2) per component with M=Mk3(:,:,e).
      !! c0 holds τ_n on entry, τ_{n+1} on exit; dt0/den/de1/cdel are per-thread scratch.
      !!
      !! AS AN INCREMENT, for the reason spelled out in advance_shape_tensor:
      !!       τ⁺ − τ_n = −M·[τ_n + μ(ε_n + ε_{n+1})] / (1 + M/2)
      !! which is the same expression rearranged, and keeps the transform error
      !! proportional to what the step changes rather than to τ itself. Without it the
      !! stiff elements (M ~ 1e-10 for the 1e30 Pa s lid) advance on round-off.
      type(response), intent(in)    :: self
      type(sht_grid),     intent(in)    :: sht
      integer,            intent(in)    :: e
      complex(wp),        intent(inout) :: c0(:,:)
      complex(wp),        intent(in)    :: eps_n(:,:), eps_1(:,:)
      real(wp),           intent(inout) :: dt0(:,:,:), den(:,:,:), de1(:,:,:)
      complex(wp),        intent(inout) :: cdel(:,:)
      type(c_ptr),        intent(in)    :: cfg
      real(wp), dimension(size(dt0,1),size(dt0,2)) :: wtau, weps
      integer  :: p
      call tensor_sh_synth(self%tsh, sht, c0,    dt0, cfg)
      call tensor_sh_synth(self%tsh, sht, eps_n, den, cfg)
      call tensor_sh_synth(self%tsh, sht, eps_1, de1, cfg)
      wtau = self%Mk3(:,:,e)            / (1.0_wp + 0.5_wp*self%Mk3(:,:,e))
      weps = self%mu(e)*self%Mk3(:,:,e) / (1.0_wp + 0.5_wp*self%Mk3(:,:,e))
      do p = 1, 6
         dt0(:,:,p) = -wtau*dt0(:,:,p) - weps*(den(:,:,p) + de1(:,:,p))
      end do
      call tensor_sh_analysis(self%tsh, sht, dt0, cdel, cfg)
      c0 = c0 + cdel
   end subroutine advance_shape_tensor_trap

   subroutine solve_endpoint_drift(self, sht)
      !! solve_drift into the endpoint (ε_{n+1}) targets, W included when carried.
      type(response), intent(inout) :: self
      type(sht_grid),     intent(in)    :: sht
      if (self%nlam > NLAM) then
         call solve_drift(self, sht, self%edUn_re, self%edUn_im, self%edVn_re, self%edVn_im, &
                          self%edWn_re, self%edWn_im)
      else
         call solve_drift(self, sht, self%edUn_re, self%edUn_im, self%edVn_re, self%edVn_im)
      end if
   end subroutine solve_endpoint_drift

   subroutine keep_drift(self)
      !! Remember the current surface drift for the next convergence test.
      type(response), intent(inout) :: self
      self%dUa_prev = self%dUa
      if (self%nlam > NLAM) self%dWa_prev = self%dWa
   end subroutine keep_drift

   subroutine drift_change(self, cnorm, snorm)
      !! Change in the surface drift since keep_drift, and its size: the coupling
      !! iteration's convergence signal. W(a) joins U(a) once it is carried — the
      !! toroidal field is part of the fixed point, and U(a) alone could settle
      !! while W is still moving.
      type(response), intent(in)  :: self
      real(wp),       intent(out) :: cnorm, snorm
      cnorm = maxval(abs(self%dUa - self%dUa_prev))
      snorm = maxval(abs(self%dUa))
      if (self%nlam > NLAM) then
         cnorm = max(cnorm, maxval(abs(self%dWa - self%dWa_prev)))
         snorm = max(snorm, maxval(abs(self%dWa)))
      end if
   end subroutine drift_change

   subroutine snapshot_taun(self)
      !! Snapshot the entering memory τ_n into the *0 arrays so every trapezoid pass
      !! advances from τ_n (not compounding). ε_n nodal drift is in self%dUn_* (set by
      !! begin_step) and stays fixed; ε_{n+1} drift is re-solved into self%edUn_*.
      type(response), intent(inout) :: self
      self%Are0 = self%Are;  self%Aim0 = self%Aim
      self%Bre0 = self%Bre;  self%Bim0 = self%Bim
      self%Cre0 = self%Cre;  self%Cim0 = self%Cim
   end subroutine snapshot_taun

   subroutine trapezoid_advance_all(self, sht, sigma_lm)
      !! One trapezoid endpoint advance for all (l,m): reset memory to τ_n (the *0
      !! snapshot), then advance with the entering strain ε_n (σ·xUn + dUn) and the
      !! endpoint strain ε_{n+1} (σ·xUn + edUn). Reads self%edUn_*/edVn_* (the current
      !! τ_{n+1} drift estimate); writes self%Are…. Does not touch self%dUa. With
      !! laterally-varying viscosity the trapezoid factor M is a field, so the advance
      !! goes pseudo-spectral (advance_memory_3d_trap), exactly as FE uses advance_memory_3d.
      type(response), intent(inout) :: self
      type(sht_grid),     intent(in)    :: sht
      complex(wp),        intent(in)    :: sigma_lm(:)
      real(wp), allocatable :: Ure(:), Uim(:), Vre(:), Vim(:)
      real(wp), allocatable :: Ure_n(:), Uim_n(:), Vre_n(:), Vim_n(:)
      real(wp) :: sre, sim, srn, sin
      integer  :: k, l, lm, node
      integer(kind=8) :: pc0, pc1, prate           ! PROFILE: memory-advance wall-clock
      call system_clock(pc0, prate)

      ! Spectral trapezoid over the 1-D-effective elements FIRST (masked by active1d).
      ! The per-k reset writes τ_n to every element, including the 3-D ones, so the
      ! genuinely-3-D path MUST run afterwards (below) to land their τ_{n+1} last.
      !$omp parallel default(shared) &
      !$omp   private(k, l, lm, node, sre, sim, srn, sin, Ure, Uim, Vre, Vim, Ure_n, Uim_n, Vre_n, Vim_n)
      allocate(Ure(self%nr), Uim(self%nr), Vre(self%nr), Vim(self%nr), &
               Ure_n(self%nr), Uim_n(self%nr), Vre_n(self%nr), Vim_n(self%nr))
      !$omp do schedule(static)
      do k = 1, self%nk
         l  = self%kdeg(k)
         lm = self%k2lm(k)
         sre = real(sigma_lm(lm), wp);  sim = aimag(sigma_lm(lm))          ! σ_{n+1}
         if (self%sigma_primed) then
            srn = real(self%sigma_n(lm), wp);  sin = aimag(self%sigma_n(lm))   ! σ_n
         else
            srn = sre;  sin = sim     ! first step: ε_n uses σ_{n+1} (load present at t=0)
         end if
         do node = 1, self%nr
            Ure_n(node) = srn*self%xUn(node,l) + self%dUn_re(node,k)   ! ε_n (σ_n)
            Vre_n(node) = srn*self%xVn(node,l) + self%dVn_re(node,k)
            Uim_n(node) = sin*self%xUn(node,l) + self%dUn_im(node,k)
            Vim_n(node) = sin*self%xVn(node,l) + self%dVn_im(node,k)
            Ure(node)   = sre*self%xUn(node,l) + self%edUn_re(node,k)  ! ε_{n+1}
            Vre(node)   = sre*self%xVn(node,l) + self%edVn_re(node,k)
            Uim(node)   = sim*self%xUn(node,l) + self%edUn_im(node,k)
            Vim(node)   = sim*self%xVn(node,l) + self%edVn_im(node,k)
         end do
         ! reset to τ_n, then trapezoid-advance the 1-D-effective elements (ε_n, ε_{n+1})
         self%Are(:,:,k) = self%Are0(:,:,k);  self%Bre(:,:,k) = self%Bre0(:,:,k)
         self%Cre(:,:,k) = self%Cre0(:,:,k)
         call advance_memory(self%ne, self%mu, self%Mk, Ure, Vre, self%Jr(l), &
              self%Are(:,:,k), self%Bre(:,:,k), self%Cre(:,:,k), &
              scheme=SCHEME_TRAP, Un_prev=Ure_n, Vn_prev=Vre_n, active=self%active1d)
         self%Aim(:,:,k) = self%Aim0(:,:,k);  self%Bim(:,:,k) = self%Bim0(:,:,k)
         self%Cim(:,:,k) = self%Cim0(:,:,k)
         call advance_memory(self%ne, self%mu, self%Mk, Uim, Vim, self%Jr(l), &
              self%Aim(:,:,k), self%Bim(:,:,k), self%Cim(:,:,k), &
              scheme=SCHEME_TRAP, Un_prev=Uim_n, Vn_prev=Vim_n, active=self%active1d)
         if (self%nlam > NLAM) then     ! W: ε_n from τ_n's drift, ε_{n+1} from the endpoint's
            call advance_memory_tor(self%ne, self%mu, self%Mk, self%edWn_re(:,k), self%Jr(l), &
                 self%Are(:,:,k), self%Bre(:,:,k), self%Cre(:,:,k), &
                 scheme=SCHEME_TRAP, Wn_prev=self%dWn_re(:,k), active=self%active1d)
            call advance_memory_tor(self%ne, self%mu, self%Mk, self%edWn_im(:,k), self%Jr(l), &
                 self%Aim(:,:,k), self%Bim(:,:,k), self%Cim(:,:,k), &
                 scheme=SCHEME_TRAP, Wn_prev=self%dWn_im(:,k), active=self%active1d)
         end if
      end do
      !$omp end do
      deallocate(Ure, Uim, Vre, Vim, Ure_n, Uim_n, Vre_n, Vim_n)
      !$omp end parallel

      ! genuinely-3-D elements last, reading the intact τ_n snapshot (Are0 …)
      if (self%lat_visc) call advance_memory_3d_trap(self, sht, sigma_lm)
      call system_clock(pc1)
      self%t_mem = self%t_mem + real(pc1-pc0,wp)/prate;  self%n_mem = self%n_mem + 1
   end subroutine trapezoid_advance_all

   subroutine ve_response_prepare_endpoint(self, sht)
      !! Open a co-converging step (§3c 3b): snapshot τ_n and seed the endpoint drift
      !! ε_{n+1} with the entering τ_n drift (begin_step's dUn). The SLE driver then
      !! converges σ against the current report drift (dUa, = drift(τ_n) on entry) and
      !! calls advance_endpoint, which refreshes the report drift to τ_{n+1}.
      type(response), intent(inout) :: self
      type(sht_grid),     intent(in)    :: sht
      self%couple_pass = 0
      self%couple_done = .false.
      if (.not. scheme_is_implicit(self%scheme)) return    ! FE: nothing to snapshot
      call ensure_commit_scratch(self)
      call snapshot_taun(self)
      self%edUn_re = self%dUn_re;  self%edUn_im = self%dUn_im
      self%edVn_re = self%dVn_re;  self%edVn_im = self%dVn_im
      if (self%nlam > NLAM) then
         self%edWn_re = self%dWn_re;  self%edWn_im = self%dWn_im
      end if
      call keep_drift(self)
   end subroutine ve_response_prepare_endpoint

   subroutine ve_response_advance_endpoint(self, sht, sigma_lm)
      !! One co-convergence pass with the SLE-converged load σ (§3c 3b). FE: a single
      !! Maxwell update (1st-order; no co-iteration). TRAP: trapezoid-advance τ_n→τ_{n+1}
      !! using the current endpoint-drift estimate, THEN refresh the report drift dUa
      !! (and ε_{n+1}=edUn) from the new τ_{n+1} — so the driver's next σ-convergence,
      !! and the next pass's endpoint strain, see the advanced memory. Sets couple_done
      !! when the surface drift settles to couple_tol. Does NOT advance time.
      type(response), intent(inout) :: self
      type(sht_grid),     intent(in)    :: sht
      complex(wp),        intent(in)    :: sigma_lm(:)
      real(wp) :: cnorm, snorm

      if (.not. scheme_is_implicit(self%scheme)) then
         call vilma_advance(self, sht, sigma_lm)
         self%couple_pass = 1
         self%couple_done = .true.
         return
      end if

      call trapezoid_advance_all(self, sht, sigma_lm)   ! ε_{n+1} from the previous estimate
      self%sigma_next = sigma_lm                   ! stage σ_{n+1}; finalize commits it to σ_n
      ! refresh dUa + ε_{n+1} drift from the new τ_{n+1} (Are…), ready for the next pass
      call solve_endpoint_drift(self, sht)
      self%couple_pass = self%couple_pass + 1
      call drift_change(self, cnorm, snorm)
      call keep_drift(self)
      self%couple_done = (self%couple_pass >= 2 .and. &
                          cnorm <= self%couple_tol*max(snorm, tiny(1.0_wp)))
   end subroutine ve_response_advance_endpoint

   logical function ve_response_endpoint_converged(self) result(done)
      type(response), intent(in) :: self
      done = self%couple_done
   end function ve_response_endpoint_converged

   subroutine ve_response_finalize_step(self, sht)
      !! Close a co-converging step: the memory already holds the converged τ_{n+1}
      !! (advance_endpoint left it there), so only advance time.
      type(response), intent(inout) :: self
      type(sht_grid),     intent(in)    :: sht
      self%couple_iters_last = self%couple_pass
      if (allocated(self%sigma_next)) then          ! TRAP: commit σ_{n+1} as next σ_n
         self%sigma_n = self%sigma_next;  self%sigma_primed = .true.
      end if
      self%time = self%time + self%dt
   end subroutine ve_response_finalize_step

   subroutine response_set_dt(self, dt)
      !! Change the step size. Δt enters only through Mk = (μ/η)·Δt, so rescale Mk from
      !! the Δt-invariant rate MkPerDt — exact (no drift from repeated halving/restoring)
      !! and no operator re-factor (the band LU is Δt-independent). The adaptive
      !! controller uses this to try a step, halve for the fine sub-steps, and restore.
      type(response), intent(inout) :: self
      real(wp),           intent(in)    :: dt
      self%dt = dt
      ! Only RESP_VE carries a Maxwell factor. RESP_ELASTIC / RESP_NULL never
      ! allocate Mk/MkPerDt at all — so returning early for everything but VE is what keeps
      ! the elastic path off an unallocated array (vilma_timestep calls this
      ! unconditionally on the explicit path, which elastic reaches).
      if (self%kind /= RESP_VE) return
      self%Mk = self%MkPerDt * dt
      if (self%lat_visc) self%Mk3 = self%MkPerDt3 * dt
   end subroutine response_set_dt

   subroutine ensure_state_scratch(self)
      !! Lazily allocate the controller's state buffers A (τ_n) and B (τ_coarse).
      type(response), intent(inout) :: self
      if (allocated(self%Are_s)) return
      allocate(self%Are_s(self%nlam,self%ne,self%nk), self%Aim_s(self%nlam,self%ne,self%nk))
      allocate(self%Bre_s(self%nlam,self%ne,self%nk), self%Bim_s(self%nlam,self%ne,self%nk))
      allocate(self%Cre_s(self%nlam,self%ne,self%nk), self%Cim_s(self%nlam,self%ne,self%nk))
      allocate(self%Are_c(self%nlam,self%ne,self%nk), self%Aim_c(self%nlam,self%ne,self%nk))
      allocate(self%Bre_c(self%nlam,self%ne,self%nk), self%Bim_c(self%nlam,self%ne,self%nk))
      allocate(self%Cre_c(self%nlam,self%ne,self%nk), self%Cim_c(self%nlam,self%ne,self%nk))
      allocate(self%sigma_n_s(self%nlm))
   end subroutine ensure_state_scratch

   subroutine response_save_state(self)
      !! Snapshot the entering prognostic state (memory τ_n + time + σ_n) into buffer A.
      !! A rejected step or the fine sub-step path restores to this with restore_state.
      type(response), intent(inout) :: self
      integer :: k
      call ensure_state_scratch(self)
      ! Threaded on the same schedule(static) partition over k as every other loop
      ! that touches these arrays, so a thread copies the slices it already owns
      ! (and, on the first call, first-touches the destination pages onto its own
      ! NUMA domain). This snapshot is taken on EVERY sub-step and was the single
      ! largest phase of the driver's residual cost bucket -- ~61 % of it -- purely
      ! because it ran serially. Slices of the last dimension are contiguous, so
      ! these are plain memcpys with no temporaries.
      !$omp parallel do default(shared) private(k) schedule(static)
      do k = 1, self%nk
         self%Are_s(:,:,k) = self%Are(:,:,k)
         self%Aim_s(:,:,k) = self%Aim(:,:,k)
         self%Bre_s(:,:,k) = self%Bre(:,:,k)
         self%Bim_s(:,:,k) = self%Bim(:,:,k)
         self%Cre_s(:,:,k) = self%Cre(:,:,k)
         self%Cim_s(:,:,k) = self%Cim(:,:,k)
      end do
      self%time_s = self%time
      if (allocated(self%sigma_n)) self%sigma_n_s = self%sigma_n
      self%sigma_primed_s = self%sigma_primed
   end subroutine response_save_state

   subroutine response_restore_state(self)
      !! Restore the prognostic state saved by save_state (buffer A).
      type(response), intent(inout) :: self
      integer :: k
      !$omp parallel do default(shared) private(k) schedule(static)
      do k = 1, self%nk
         self%Are(:,:,k) = self%Are_s(:,:,k)
         self%Aim(:,:,k) = self%Aim_s(:,:,k)
         self%Bre(:,:,k) = self%Bre_s(:,:,k)
         self%Bim(:,:,k) = self%Bim_s(:,:,k)
         self%Cre(:,:,k) = self%Cre_s(:,:,k)
         self%Cim(:,:,k) = self%Cim_s(:,:,k)
      end do
      self%time = self%time_s
      if (allocated(self%sigma_n)) self%sigma_n = self%sigma_n_s
      self%sigma_primed = self%sigma_primed_s
   end subroutine response_restore_state

   subroutine response_stash_coarse(self)
      !! Snapshot the current memory (the coarse one-Δt τ_{n+1}) into buffer B for the
      !! step-doubling error estimate, to be compared against the fine result.
      type(response), intent(inout) :: self
      integer :: k
      call ensure_state_scratch(self)
      !$omp parallel do default(shared) private(k) schedule(static)
      do k = 1, self%nk
         self%Are_c(:,:,k) = self%Are(:,:,k)
         self%Aim_c(:,:,k) = self%Aim(:,:,k)
         self%Bre_c(:,:,k) = self%Bre(:,:,k)
         self%Bim_c(:,:,k) = self%Bim(:,:,k)
         self%Cre_c(:,:,k) = self%Cre(:,:,k)
         self%Cim_c(:,:,k) = self%Cim(:,:,k)
      end do
   end subroutine response_stash_coarse

   subroutine response_coarse_fine_error(self, err_inf, tau_inf)
      !! After the fine path, return the coarse↔fine memory difference ‖τ_fine−τ_coarse‖∞
      !! (buffer B is τ_coarse, self%Are… is τ_fine) and the memory magnitude ‖τ_fine‖∞,
      !! for the controller's scaled local-error estimate.
      type(response), intent(in)  :: self
      real(wp),           intent(out) :: err_inf, tau_inf
      integer  :: k, nb
      ! Explicit threaded loop for the same reason response_memory_norm uses one:
      ! the maxval(abs(A - B)) form this replaced built TWELVE full (NLAM,ne,nk)
      ! heap temporaries per call, one per term, and ran serially. Each slot block
      ! is scanned flat (see solve_drift).
      nb = self%nlam*self%ne
      err_inf = 0.0_wp;  tau_inf = 0.0_wp
      !$omp parallel do default(shared) private(k) schedule(static) &
      !$omp   reduction(max:err_inf,tau_inf)
      do k = 1, self%nk
         err_inf = max(err_inf, absdiffmax6(nb, &
                       self%Are(1,1,k), self%Aim(1,1,k), self%Bre(1,1,k), &
                       self%Bim(1,1,k), self%Cre(1,1,k), self%Cim(1,1,k), &
                       self%Are_c(1,1,k), self%Aim_c(1,1,k), self%Bre_c(1,1,k), &
                       self%Bim_c(1,1,k), self%Cre_c(1,1,k), self%Cim_c(1,1,k)))
         tau_inf = max(tau_inf, absmax6(nb, self%Are(1,1,k), self%Aim(1,1,k), &
                       self%Bre(1,1,k), self%Bim(1,1,k), self%Cre(1,1,k), self%Cim(1,1,k)))
      end do
   end subroutine response_coarse_fine_error

   real(wp) function response_max_rate(self) result(rate)
      !! Largest Maxwell rate μ/η over all memory-carrying elements (and, with
      !! lateral viscosity, over the 3D grid). Sets the explicit forward-Euler
      !! stability ceiling Δt ≤ cfl/rate (M = μΔt/η ≤ cfl). rate = 0 ⇒ no Maxwell
      !! memory (purely elastic) ⇒ the caller may take a single step.
      type(response), intent(in) :: self
      rate = 0.0_wp
      if (allocated(self%MkPerDt)) rate = maxval(self%MkPerDt)
      if (self%lat_visc .and. allocated(self%MkPerDt3)) &
         rate = max(rate, maxval(self%MkPerDt3))
   end function response_max_rate

   real(wp) function response_memory_norm(self) result(nrm)
      !! ∞-norm of the viscoelastic memory stress over all coefficients/elements —
      !! the reactive-guard observable for the explicit stepper (a non-finite or
      !! runaway value flags an unstable sub-step). Explicit loop to avoid the
      !! abs(slice) heap temporaries maxval would create on these (NLAM,ne,nk) arrays.
      type(response), intent(in) :: self
      integer  :: k, nb
      nrm = 0.0_wp
      if (.not. allocated(self%Are)) return
      nb = self%nlam*self%ne                   ! slot blocks scanned flat (see solve_drift)
      !$omp parallel do default(shared) private(k) reduction(max:nrm) schedule(static)
      do k = 1, self%nk
         nrm = max(nrm, absmax6(nb, self%Are(1,1,k), self%Aim(1,1,k), self%Bre(1,1,k), &
                                    self%Bim(1,1,k), self%Cre(1,1,k), self%Cim(1,1,k)))
      end do
      !$omp end parallel do
   end function response_memory_norm

   subroutine ensure_commit_scratch(self)
      !! Lazily allocate the implicit-commit scratch (the τ_n memory snapshot doubles
      !! the memory footprint, so it is only paid when a TRAP commit is first used).
      type(response), intent(inout) :: self
      if (allocated(self%Are0)) return
      allocate(self%Are0(self%nlam,self%ne,self%nk), self%Aim0(self%nlam,self%ne,self%nk))
      allocate(self%Bre0(self%nlam,self%ne,self%nk), self%Bim0(self%nlam,self%ne,self%nk))
      allocate(self%Cre0(self%nlam,self%ne,self%nk), self%Cim0(self%nlam,self%ne,self%nk))
      allocate(self%edUn_re(self%nr,self%nk), self%edUn_im(self%nr,self%nk))
      allocate(self%edVn_re(self%nr,self%nk), self%edVn_im(self%nr,self%nk))
      allocate(self%dUa_prev(self%nk))
      if (self%nlam > NLAM) call ensure_commit_scratch_tor(self)
      call ensure_sigma(self)
   end subroutine ensure_commit_scratch

   subroutine ensure_sigma(self)
      !! Lazily allocate the start-of-step load buffers (σ_n / σ_next). Separate from
      !! ensure_commit_scratch so prime_sigma can seed σ_0 before the first commit.
      type(response), intent(inout) :: self
      if (allocated(self%sigma_n)) return
      allocate(self%sigma_n(self%nlm), self%sigma_next(self%nlm))
      self%sigma_n = (0.0_wp, 0.0_wp)     ! σ at t=0; primed by prime_sigma or step 1
      self%sigma_primed = .false.
   end subroutine ensure_sigma

   subroutine response_prime_sigma(self, sigma_lm)
      !! Seed the start-of-step load σ_n with a known load (the elastic-consistent SLE
      !! load at t=0) and mark it tracked, so the trapezoidal ε_n on the FIRST step uses
      !! the true σ_0 rather than the σ_{n+1} proxy — making that step 2nd→3rd order.
      type(response), intent(inout) :: self
      complex(wp),        intent(in)    :: sigma_lm(:)
      call ensure_sigma(self)
      self%sigma_n = sigma_lm
      self%sigma_primed = .true.
   end subroutine response_prime_sigma

   subroutine ve_response_destroy(self)
      type(response), intent(inout) :: self
      integer :: l
      if (allocated(self%ops)) then
         do l = 1, size(self%ops);  call radial_operator_destroy(self%ops(l));  end do
         deallocate(self%ops)
      end if
      if (allocated(self%r))      deallocate(self%r)
      if (allocated(self%mu))     deallocate(self%mu)
      if (allocated(self%Mk3))      deallocate(self%Mk3)
      if (allocated(self%MkPerDt3)) deallocate(self%MkPerDt3)
      if (allocated(self%e3d))      deallocate(self%e3d)
      if (allocated(self%active1d)) deallocate(self%active1d)
      self%ne3d = 0
      if (self%lat_visc) call tensor_sh_destroy(self%tsh)
      self%lat_visc = .false.
      if (allocated(self%Mk))     deallocate(self%Mk)
      if (allocated(self%Jr))     deallocate(self%Jr)
      if (allocated(self%nrmc))   deallocate(self%nrmc)
      if (allocated(self%sa))     deallocate(self%sa)
      if (allocated(self%sb))     deallocate(self%sb)
      if (allocated(self%sc))     deallocate(self%sc)
      if (allocated(self%gu))     deallocate(self%gu)
      if (allocated(self%gn))     deallocate(self%gn)
      if (allocated(self%gv))     deallocate(self%gv)
      if (allocated(self%xUn))    deallocate(self%xUn)
      if (allocated(self%xVn))    deallocate(self%xVn)
      if (allocated(self%Are))    deallocate(self%Are)
      if (allocated(self%Aim))    deallocate(self%Aim)
      if (allocated(self%Bre))    deallocate(self%Bre)
      if (allocated(self%Bim))    deallocate(self%Bim)
      if (allocated(self%Cre))    deallocate(self%Cre)
      if (allocated(self%Cim))    deallocate(self%Cim)
      if (allocated(self%dUa))    deallocate(self%dUa)
      if (allocated(self%dFa))    deallocate(self%dFa)
      if (allocated(self%dVa))    deallocate(self%dVa)
      if (allocated(self%dUn_re)) deallocate(self%dUn_re)
      if (allocated(self%dUn_im)) deallocate(self%dUn_im)
      if (allocated(self%dVn_re)) deallocate(self%dVn_re)
      if (allocated(self%dVn_im)) deallocate(self%dVn_im)
      if (allocated(self%k2lm))   deallocate(self%k2lm)
      if (allocated(self%kdeg))   deallocate(self%kdeg)
      if (allocated(self%kbeg))   deallocate(self%kbeg)
      if (allocated(self%mnorm))  deallocate(self%mnorm)
      if (allocated(self%Are0))   deallocate(self%Are0)
      if (allocated(self%Aim0))   deallocate(self%Aim0)
      if (allocated(self%Bre0))   deallocate(self%Bre0)
      if (allocated(self%Bim0))   deallocate(self%Bim0)
      if (allocated(self%Cre0))   deallocate(self%Cre0)
      if (allocated(self%Cim0))   deallocate(self%Cim0)
      if (allocated(self%edUn_re)) deallocate(self%edUn_re)
      if (allocated(self%edUn_im)) deallocate(self%edUn_im)
      if (allocated(self%edVn_re)) deallocate(self%edVn_re)
      if (allocated(self%edVn_im)) deallocate(self%edVn_im)
      if (allocated(self%dUa_prev)) deallocate(self%dUa_prev)
      if (allocated(self%sigma_n))    deallocate(self%sigma_n)
      if (allocated(self%sigma_next)) deallocate(self%sigma_next)
      if (allocated(self%MkPerDt))    deallocate(self%MkPerDt)
      if (allocated(self%Are_s))      deallocate(self%Are_s, self%Aim_s, self%Bre_s, &
                                                 self%Bim_s, self%Cre_s, self%Cim_s)
      if (allocated(self%Are_c))      deallocate(self%Are_c, self%Aim_c, self%Bre_c, &
                                                 self%Bim_c, self%Cre_c, self%Cim_c)
      if (allocated(self%sigma_n_s))  deallocate(self%sigma_n_s)
      if (allocated(self%tops)) then
         do l = 1, size(self%tops);  call toroidal_operator_destroy(self%tops(l));  end do
         deallocate(self%tops)
      end if
      if (allocated(self%nrmt))     deallocate(self%nrmt, self%sat, self%sbt, self%sct)
      if (allocated(self%dWa))      deallocate(self%dWa)
      if (allocated(self%dWn_re))   deallocate(self%dWn_re, self%dWn_im)
      if (allocated(self%edWn_re))  deallocate(self%edWn_re, self%edWn_im)
      if (allocated(self%dWa_prev)) deallocate(self%dWa_prev)
      self%nlam = NLAM
      self%lmax = 0;  self%nlm = 0;  self%nk = 0
      self%sigma_primed = .false.;  self%sigma_primed_s = .false.
   end subroutine ve_response_destroy

end module vilma_response
