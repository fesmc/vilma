module vilma_rotation
   !! Rotational feedback / true polar wander (Spada et al. 2011 §2.1.1; the
   !! time-domain rotational theory of Martinec & Hagedoorn 2014, as in VILMA-v1).
   !!
   !! A surface load and the deformation it drives perturb the off-diagonal
   !! inertia (I₁₃, I₂₃); the linearized Liouville equation maps that to equatorial
   !! polar motion m = m₁ + i m₂; the shifted pole perturbs the centrifugal
   !! potential — a degree-2, order-1 potential — which deforms the Earth back
   !! (tidal response) and feeds the geoid / sea-level equation. We integrate the
   !! GIA (quasi-static) Liouville equation with the Chandler wobble neglected
   !! (Spada eq. 7), the regime the GIA timescale lives in:
   !!
   !!     [1 − k^T(t)/k_s] ∗ m(t) = Ψ_L(t),                              (Spada 7)
   !!     Ψ_L(t) = I(t)/(C−A),  I(t) = [δ(t)+k^L(t)] ∗ I_rigid(t),       (Spada 20)
   !!
   !! with k^T, k^L the degree-2 tidal / loading Love numbers and k_s ≡ k^T_f the
   !! secular (fluid) tidal Love number (Spada eq. 11). NO explicit Ω appears —
   !! the centrifugal scaling is absorbed into m (= ω/Ω) and k_s by construction.
   !!
   !! Method (composes with the time-domain Maxwell machinery, no normal modes,
   !! no convolution quadrature). Two compact degree-2 viscoelastic channels carry
   !! the convolutions as Maxwell memory (reusing the per-element kernel of
   !! vilma_viscoelastic):
   !!   - a LOADING channel: forced by the rigid inertia I_rigid, returns
   !!     I(t) = [1+k^L]∗I_rigid  ⇒  Ψ_L = I/(C−A);
   !!   - a TIDAL channel: forced by the centrifugal potential ∝ m, returns the
   !!     induced potential k^T∗m  ⇒  the rotational feedback.
   !! Each step the feedback makes the Liouville equation ALGEBRAIC in m (the affine
   !! begin/apply/commit structure of the field driver):
   !!
   !!     m_n = [ Ψ_L,n − dF_tidal/k_s ] / [ 1 − k^T_e/k_s ],
   !!
   !! where k^T_e is the elastic tidal Love number and dF_tidal is the tidal
   !! channel's frozen memory (the m-history). The memory is then advanced with the
   !! converged m_n. Rotation is purely degree 2, so this is a single complex
   !! coefficient per channel — self-contained, decoupled from the SLE field driver
   !! (the centrifugal potential is fed back into the SLE in a later step).
   !!
   !! 3-D ready: I_rigid is a direct Gauss-grid quadrature of the actual load
   !! (any field), so no axisymmetric assumption enters; only the (1+k^L)/k^T
   !! channels use the 1-D radial relaxation (laterally-varying η is rung 6).
   use vilma_precision,       only: wp
   use vilma_constants,       only: pi, grav_G
   use vilma_earth_structure, only: earth_n_layers, earth_gravity_at, earth_model, RHEOL_FLUID
   use vilma_radial_fe,       only: radial_operator_load_rhs, radial_operator_tidal_rhs, radial_operator_destroy, radial_operator_solve_vec, radial_operator_assemble, radial_mesh_build, radial_mesh, radial_operator, tidal_love, &
                                 idx_u, idx_v, idx_f, ndof_of
   use vilma_viscoelastic,    only: NLAM, ve_strain_constants, dissipative_rhs, &
                                 advance_memory, SCHEME_FE
   use vilma_sht,             only: sht_grid, sht_grid_surface_integral
   implicit none
   private

   public :: rotation_state
   public :: channel_init, channel_set_dt, channel_begin, channel_commit, channel_destroy, rotation_init, rotation_begin_step, rotation_solve_m, rotation_s_rot, rotation_commit, rotation_update, rotation_destroy
   public :: rotation_ne, rotation_get_memory, rotation_set_memory, ROT_NCOMP

   integer, parameter :: JROT = 2          !! rotation is purely degree 2
   integer, parameter :: ROT_NCOMP = 6     !! packed memory components per channel:
                                           !! [Are, Aim, Bre, Bim, Cre, Cim]

   type :: deg2_channel
      !! A single degree-2 viscoelastic response to a COMPLEX forcing coefficient,
      !! carrying Maxwell memory. `tidal` selects the forcing kind (external
      !! potential via tidal_rhs, vs. surface load via load_rhs); the surface
      !! perturbed-potential coefficient F(a) is the readout. Affine in the current
      !! forcing: F = Fe·coeff + dF, with dF frozen from the entering memory each
      !! step (begin_step) and the memory advanced with the converged coeff (commit).
      logical  :: tidal = .false.
      integer  :: nr = 0, ne = 0, ndof = 0
      real(wp) :: Jr = 0.0_wp, dt = 0.0_wp
      type(radial_operator) :: op                 !! degree-2 operator (assembled once)
      real(wp), allocatable :: r(:), mu(:), Mk(:), MkPerDt(:)
      real(wp) :: norm(NLAM), sa(4,NLAM), sb(4,NLAM), sc(4,NLAM)
      real(wp) :: Fe = 0.0_wp                      !! elastic surface-F per unit forcing
      real(wp) :: Ue = 0.0_wp                      !! elastic surface-U per unit forcing (→ h^T)
      real(wp), allocatable :: xUn(:), xVn(:)      !! unit-forcing nodal U,V
      ! per-element memory stress, real/imag (NLAM, ne)
      real(wp), allocatable :: Are(:,:), Aim(:,:), Bre(:,:), Bim(:,:), Cre(:,:), Cim(:,:)
      ! frozen drift from the entering memory τ_n (begin_step)
      complex(wp) :: dF = (0.0_wp, 0.0_wp)         !! surface F drift
      complex(wp) :: dU = (0.0_wp, 0.0_wp)         !! surface U drift (uplift)
      real(wp), allocatable :: dUn_re(:), dUn_im(:), dVn_re(:), dVn_im(:)  !! nodal ε_n drift
   end type deg2_channel

   type :: rotation_state
      logical     :: enabled = .false.       !! set from p%rotation by the coupling init
      complex(wp) :: m = (0.0_wp, 0.0_wp)     !! polar motion m₁ + i m₂ [rad]
      real(wp)    :: time = 0.0_wp            !! model time [s]
      ! physics constants (defaults: Spada 2011 Table 2; overridable for deep time)
      real(wp)    :: a       = 6.371e6_wp     !! Earth radius [m]
      real(wp)    :: g       = 9.81_wp        !! surface gravity [m s⁻²]
      real(wp)    :: CminusA = 2.63e35_wp     !! C − A [kg m²]
      real(wp)    :: Omega   = 7.292115e-5_wp !! mean rotation rate Ω [s⁻¹]
      real(wp)    :: k_s     = 0.0_wp         !! secular tidal Love number used (k_s)
      real(wp)    :: k_s_fluid = 0.0_wp       !! model relaxed limit k^T_f (Spada eq. 11 benchmark value)
      real(wp)    :: k_s_flat  = 0.0_wp       !! observed-flattening k_s = 3G(C−A)/(a⁵Ω²) (Adhikari/Mitrovica)
      real(wp)    :: kTe     = 0.0_wp         !! elastic tidal Love number k^T_e (degree 2)
      real(wp)    :: hTe     = 0.0_wp         !! elastic tidal Love number h^T_e (degree 2)
      real(wp)    :: dt_fe_max = huge(1.0_wp) !! forward-Euler stability ceiling, Δt <= min(η/μ):
                                              !! half the 2·min(η/μ) limit, i.e. the 0.5
                                              !! safety factor is already folded in
                                              !! (see rotation_init)
                                              !! (the channels use explicit FE; the driver
                                              !! sub-steps a coupling interval to respect this)
      complex(wp) :: cload   = (0.0_wp,0.0_wp)!! load-channel operator coefficient (set by solve_m, used by commit)
      type(deg2_channel) :: load_ch          !! (1+k^L)∗ channel
      type(deg2_channel) :: tidal_ch         !! k^T∗ channel
   end type rotation_state

contains

   ! === rotation_state ========================================================

   subroutine rotation_init(self, earth, sht, dt, k_s, CminusA)
      !! Build the two degree-2 channels, the elastic tidal Love number k^T_e, and
      !! the secular k_s = k^T_f (the relaxed tidal limit = elastic tidal solve of
      !! the model with every Maxwell layer fluidized). Pass k_s to override it with
      !! an observed-flattening value (Mitrovica et al. 2005) for deep-time runs,
      !! and CminusA to override C − A (which normalises Ψ_L and sets k_s_flat).
      type(rotation_state), intent(inout) :: self
      type(earth_model),     intent(in)    :: earth
      type(sht_grid),        intent(in)    :: sht
      real(wp),              intent(in)    :: dt
      real(wp), optional,    intent(in)    :: k_s
      real(wp), optional,    intent(in)    :: CminusA
      type(radial_mesh) :: mesh

      call rotation_destroy(self)
      call radial_mesh_build(mesh, earth)
      self%a = earth%r_earth
      self%g = earth_gravity_at(earth, earth%r_earth)
      call channel_init(self%load_ch, earth, mesh, dt, tidal=.false.)
      call channel_init(self%tidal_ch, earth, mesh, dt, tidal=.true.)
      ! elastic tidal Love numbers from the tidal channel's unit response (φ_t = 1):
      ! k^T = −F(a)/φ_t − 1, h^T = g U(a)/φ_t (tidal_love convention).
      self%kTe = -self%tidal_ch%Fe - 1.0_wp
      self%hTe =  earth_gravity_at(earth, earth%r_earth) * self%tidal_ch%Ue
      ! two secular Love numbers (Spada eq. 11 vs Adhikari/Mitrovica): the model
      ! relaxed limit k^T_f reproduces the Spada Test 3/2 benchmark; the observed-
      ! flattening closed form k_s = 3G(C−A)/(a⁵Ω²) avoids the lithosphere-thickness
      ! paradox and is the recommended deep-time value.
      if (present(CminusA)) self%CminusA = CminusA
      self%k_s_fluid = fluid_tidal_k(earth, mesh)
      self%k_s_flat  = 3.0_wp*grav_G*self%CminusA/(self%a**5*self%Omega**2)
      ! Forward-Euler stability ceiling: Mk = (μ/η)Δt < 2 ⇒ Δt < 2/max(μ/η), with a
      ! 0.5 safety factor. The driver sub-steps any coupling interval larger than this.
      if (maxval(self%load_ch%MkPerDt) > 0.0_wp) &
         self%dt_fe_max = 1.0_wp/maxval(self%load_ch%MkPerDt)
      if (present(k_s)) then
         self%k_s = k_s                       ! explicit override (e.g. observed flattening)
      else
         self%k_s = self%k_s_fluid            ! default: model fluid limit (benchmark)
      end if
      self%m = (0.0_wp, 0.0_wp);  self%time = 0.0_wp
   end subroutine rotation_init

   subroutine rotation_begin_step(self, sht, dt)
      !! Open a timestep: set Δt and freeze both channels' relaxation drift from the
      !! entering memory τ_n. The polar motion (solve_m) and the rotational SLE field
      !! (s_rot) are then AFFINE in the current load / m, so the caller may iterate the
      !! rotation ↔ SLE fixed point without advancing memory; commit closes the step.
      type(rotation_state), intent(inout) :: self
      type(sht_grid),        intent(in)    :: sht
      real(wp),              intent(in)    :: dt
      if (.not. self%enabled) return
      if (dt /= self%load_ch%dt) then
         call channel_set_dt(self%load_ch, dt);  call channel_set_dt(self%tidal_ch, dt)
      end if
      call channel_begin(self%load_ch)
      call channel_begin(self%tidal_ch)
   end subroutine rotation_begin_step

   subroutine rotation_solve_m(self, sht, load)
      !! Solve the algebraic (Chandler-neglected) Liouville equation for the polar
      !! motion under the surface mass load `load` [kg m⁻²], using the drift frozen by
      !! begin_step (pure — no memory advance, safe inside the fixed point). Sets
      !! self%m and self%cload (the load-channel coefficient commit will advance with).
      type(rotation_state), intent(inout) :: self
      type(sht_grid),        intent(in)    :: sht
      real(wp),              intent(in)    :: load(:,:)
      complex(wp) :: Irig, Itot, psiL
      real(wp)    :: scl
      if (.not. self%enabled) return
      Irig = inertia21(sht, load, self%a)
      ! LOADING: feed σ whose own degree-2 potential equals I_rigid (φ^L = 4πGaσ/(2j+1)),
      ! so −F = [1+k^L]∗I_rigid = I(t); Ψ_L = I/(C−A).
      scl        = real(2*JROT+1, wp)/(4.0_wp*pi*grav_G*self%a)
      self%cload = Irig*scl
      Itot       = -(self%load_ch%Fe*self%cload + self%load_ch%dF)
      psiL       = Itot/self%CminusA
      ! Liouville: m = Ψ_L + (1/k_s)(k^T_e m − dF_tidal) ⇒ solve for m.
      self%m = (psiL - self%tidal_ch%dF/self%k_s)/(1.0_wp - self%kTe/self%k_s)
   end subroutine rotation_solve_m

   subroutine rotation_s_rot(self, sht, srot)
      !! Build the rotational-feedback contribution to relative sea level on the Gauss
      !! grid, s_rot = N_rot − u_rot, from the current self%m (call after solve_m). The
      !! centrifugal potential Λ = −Ω²a² sinθcosθ (m₁cosφ + m₂sinφ) is a degree-2 order-1
      !! field; the sea surface and solid respond with the tidal Love numbers (Adhikari
      !! et al. 2016, eq. 8): N_rot = (1+k^T)Λ/g, u_rot = h^T Λ/g. The sign follows from
      !! Λ = ½Ω²a² sin²θ' about the displaced pole ẑ + m₁x̂ + m₂ŷ, where cosθ' = cosθ +
      !! sinθ(m₁cosφ + m₂sinφ): a point the pole moves toward comes closer to the axis,
      !! so sea level falls there. The VE (1+k^T),h^T are the tidal channel's affine
      !! response to m, which is linear, so the sign is applied to the assembled field:
      !! total potential coeff = m + P_ind with P_ind = k^T_e m − dF_tidal, uplift coeff
      !! C_u = U_e m + dU_tidal (so g·C_u = h^T∗m).
      type(rotation_state), intent(inout) :: self
      type(sht_grid),        intent(in)    :: sht
      real(wp),              intent(out)   :: srot(:,:)
      complex(wp) :: cN, cU
      real(wp)    :: kN, ku, gam, cphi, sphi
      integer     :: il, ip
      if (.not. self%enabled) then
         srot = 0.0_wp;  return
      end if
      cN = self%m + (self%kTe*self%m - self%tidal_ch%dF)     ! (1+k^T)∗m total potential coeff
      cU = self%tidal_ch%Ue*self%m + self%tidal_ch%dU        ! uplift coeff (g·cU = h^T∗m)
      ! N_rot = −(Ω²a²/g)·γ·[Re(cN)cosφ+Im(cN)sinφ]; u_rot = −Ω²a²·γ·[Re(cU)cosφ+Im(cU)sinφ]
      kN = self%Omega**2 * self%a**2 / self%g
      ku = self%Omega**2 * self%a**2
      do il = 1, sht%nlat
         gam = sin(sht%colat(il))*cos(sht%colat(il))
         do ip = 1, sht%nphi
            cphi = cos(sht%lon(ip));  sphi = sin(sht%lon(ip))
            srot(ip,il) = -gam*( kN*(real(cN,wp)*cphi + aimag(cN)*sphi) &
                               -  ku*(real(cU,wp)*cphi + aimag(cU)*sphi) )
         end do
      end do
   end subroutine rotation_s_rot

   subroutine rotation_commit(self, sht)
      !! Close the step: advance both channels' Maxwell memory with the converged
      !! state (loading with self%cload, tidal with self%m) and advance time.
      type(rotation_state), intent(inout) :: self
      type(sht_grid),        intent(in)    :: sht
      if (.not. self%enabled) return
      call channel_commit(self%load_ch, self%cload)
      call channel_commit(self%tidal_ch, self%m)
      self%time = self%time + self%load_ch%dt
   end subroutine rotation_commit

   subroutine rotation_update(self, sht, load, dt)
      !! Standalone (no SLE feedback) one-step advance of the polar motion under the
      !! surface mass load `load` [kg m⁻²]: begin_step + solve_m + commit. Reports m at
      !! the entry time (first call ⇒ elastic m₀), then advances both channels' memory.
      !! The SLE-coupled driver instead calls begin_step / solve_m / s_rot / commit so
      !! it can iterate the rotation ↔ sea-level fixed point before committing.
      type(rotation_state), intent(inout) :: self
      type(sht_grid),        intent(in)    :: sht
      real(wp),              intent(in)    :: load(:,:)
      real(wp),              intent(in)    :: dt
      if (.not. self%enabled) return
      call rotation_begin_step(self, sht, dt)
      call rotation_solve_m(self, sht, load)
      call rotation_commit(self, sht)
   end subroutine rotation_update

   subroutine rotation_destroy(self)
      type(rotation_state), intent(inout) :: self
      call channel_destroy(self%load_ch)
      call channel_destroy(self%tidal_ch)
      self%m = (0.0_wp, 0.0_wp);  self%time = 0.0_wp
      self%k_s = 0.0_wp;  self%kTe = 0.0_wp;  self%hTe = 0.0_wp
   end subroutine rotation_destroy

   ! === restart serialization =================================================
   ! The prognostic state of the rotation solver is the polar motion m plus both
   ! channels' per-element Maxwell memory stress; everything else (the operators,
   ! elastic Love numbers, secular constants) is rebuilt deterministically by
   ! rotation_init. The drift fields are intra-step (frozen by begin_step from the
   ! memory each step), so they are not persisted. These accessors keep the channel
   ! internals private to this module — vilma_io moves only opaque packed arrays.

   pure integer function rotation_ne(self) result(ne)
      !! Maxwell elements per degree-2 channel (both channels share the mesh); 0 if
      !! the channels were never allocated (rotation not initialised).
      type(rotation_state), intent(in) :: self
      ne = 0
      if (allocated(self%load_ch%Are)) ne = size(self%load_ch%Are, 2)
   end function rotation_ne

   subroutine rotation_get_memory(self, m, load_mem, tidal_mem)
      !! Serialize the prognostic state: the polar motion m and both channels'
      !! memory stress, packed along the last axis (ROT_NCOMP = 6) as
      !! [Are, Aim, Bre, Bim, Cre, Cim]. load_mem/tidal_mem are (NLAM, rotation_ne, 6).
      type(rotation_state), intent(in)  :: self
      complex(wp),          intent(out) :: m
      real(wp),             intent(out) :: load_mem(:,:,:), tidal_mem(:,:,:)
      m = self%m
      call pack_channel(self%load_ch,  load_mem)
      call pack_channel(self%tidal_ch, tidal_mem)
   end subroutine rotation_get_memory

   subroutine rotation_set_memory(self, m, load_mem, tidal_mem)
      !! Inverse of rotation_get_memory: restore m and both channels' memory stress.
      type(rotation_state), intent(inout) :: self
      complex(wp),          intent(in)    :: m
      real(wp),             intent(in)    :: load_mem(:,:,:), tidal_mem(:,:,:)
      self%m = m
      call unpack_channel(self%load_ch,  load_mem)
      call unpack_channel(self%tidal_ch, tidal_mem)
   end subroutine rotation_set_memory

   subroutine pack_channel(ch, mem)
      type(deg2_channel), intent(in)  :: ch
      real(wp),           intent(out) :: mem(:,:,:)   !! (NLAM, ne, ROT_NCOMP)
      mem(:,:,1) = ch%Are;  mem(:,:,2) = ch%Aim
      mem(:,:,3) = ch%Bre;  mem(:,:,4) = ch%Bim
      mem(:,:,5) = ch%Cre;  mem(:,:,6) = ch%Cim
   end subroutine pack_channel

   subroutine unpack_channel(ch, mem)
      type(deg2_channel), intent(inout) :: ch
      real(wp),           intent(in)    :: mem(:,:,:)
      ch%Are = mem(:,:,1);  ch%Aim = mem(:,:,2)
      ch%Bre = mem(:,:,3);  ch%Bim = mem(:,:,4)
      ch%Cre = mem(:,:,5);  ch%Cim = mem(:,:,6)
   end subroutine unpack_channel

   ! === degree-2 inertia from the load (3-D-ready grid quadrature) =============

   complex(wp) function inertia21(sht, load, a) result(I21)
      !! Off-diagonal inertia perturbation I₁₃ + i I₂₃ of a surface mass load on the
      !! Gauss grid: I₁₃ = −a⁴∫σ sinθcosθ cosφ dΩ, I₂₃ = −a⁴∫σ sinθcosθ sinφ dΩ
      !! (= −a⁴∫σ sinθcosθ e^{iφ} dΩ, the degree-2 order-1 mass moment). Direct
      !! quadrature — no spherical-harmonic normalization enters, so it is exact for
      !! any (3-D) load field.
      type(sht_grid), intent(in) :: sht
      real(wp),       intent(in) :: load(:,:)   !! (nphi, nlat) [kg m⁻²]
      real(wp),       intent(in) :: a
      real(wp), allocatable :: w13(:,:), w23(:,:)
      real(wp) :: st, ct, sc2
      integer  :: il, ip
      allocate(w13(sht%nphi, sht%nlat), w23(sht%nphi, sht%nlat))
      do il = 1, sht%nlat
         st  = sin(sht%colat(il));  ct = cos(sht%colat(il))
         sc2 = st*ct                                   ! sinθ cosθ
         do ip = 1, sht%nphi
            w13(ip,il) = load(ip,il)*sc2*cos(sht%lon(ip))
            w23(ip,il) = load(ip,il)*sc2*sin(sht%lon(ip))
         end do
      end do
      I21 = cmplx(-a**4*sht_grid_surface_integral(sht, w13), &
                  -a**4*sht_grid_surface_integral(sht, w23), wp)
   end function inertia21

   ! === secular (fluid) tidal Love number =====================================

   real(wp) function fluid_tidal_k(earth, mesh) result(k_s)
      !! k_s = k^T_f: the relaxed (t→∞) degree-2 tidal Love number = the ELASTIC
      !! tidal solve of the model with every MAXWELL (viscous) layer fluidized (μ=0).
      !! The elastic lithosphere (RHEOL_ELASTIC, η→∞) is kept elastic and the inviscid
      !! core (RHEOL_FLUID) is unchanged — same construction as the loading fluid limit
      !! in test_benchmark_love. This is Spada eq. 11's secular Love number, and the
      !! rotational secular slope is pathologically sensitive to it (the lithosphere-
      !! thickness paradox, Mitrovica et al. 2005): fluidizing the lithosphere by
      !! mistake inflates k^T_f and badly under-drives the late-time polar motion.
      use vilma_earth_structure, only: RHEOL_MAXWELL
      type(earth_model), intent(in) :: earth
      type(radial_mesh), intent(in) :: mesh
      type(earth_model)     :: ef
      type(radial_operator) :: op
      real(wp), allocatable :: x(:)
      real(wp) :: h, l, ua, va, fa
      integer  :: lay
      ef = earth
      do lay = 1, earth_n_layers(ef)
         if (ef%layers(lay)%rheology == RHEOL_MAXWELL) then
            ef%layers(lay)%mu = 0.0_wp;  ef%layers(lay)%rheology = RHEOL_FLUID
         end if
      end do
      call radial_operator_assemble(op, ef, mesh, JROT)
      allocate(x(op%ndof))
      call radial_operator_solve_vec(op, radial_operator_tidal_rhs(op, 1.0_wp), x)
      ua = x(idx_u(mesh%nr));  va = x(idx_v(mesh%nr));  fa = x(idx_f(mesh%nr))
      call tidal_love(ef, JROT, 1.0_wp, ua, va, fa, h, l, k_s)
      call radial_operator_destroy(op)
   end function fluid_tidal_k

   ! === deg2_channel ==========================================================

   subroutine channel_init(self, earth, mesh, dt, tidal)
      !! Assemble the degree-2 operator, the unit-forcing response (Fe + nodal U,V),
      !! the per-element Maxwell factors, and zero the memory. `tidal` picks the
      !! forcing kind (tidal_rhs vs load_rhs).
      type(deg2_channel), intent(inout) :: self
      type(earth_model),   intent(in)    :: earth
      type(radial_mesh),   intent(in)    :: mesh
      real(wp),            intent(in)    :: dt
      logical,             intent(in)    :: tidal
      real(wp), allocatable :: x(:)
      real(wp) :: eta_e
      integer  :: e, lay, node

      call channel_destroy(self)
      self%tidal = tidal
      self%nr = mesh%nr;  self%ne = mesh%ne;  self%ndof = ndof_of(mesh%nr)
      self%Jr = real(JROT, wp)*real(JROT+1, wp);  self%dt = dt

      allocate(self%r(self%nr));  self%r = mesh%r
      allocate(self%mu(self%ne), self%Mk(self%ne), self%MkPerDt(self%ne))
      do e = 1, self%ne
         lay = mesh%elem_layer(e)
         self%mu(e) = earth%layers(lay)%mu
         eta_e      = earth%layers(lay)%eta
         if (eta_e > 0.0_wp) then
            self%MkPerDt(e) = self%mu(e)/eta_e
         else
            self%MkPerDt(e) = 0.0_wp
         end if
      end do
      self%Mk = self%MkPerDt*dt

      call ve_strain_constants(self%Jr, self%norm, self%sa, self%sb, self%sc)
      call radial_operator_assemble(self%op, earth, mesh, JROT)

      ! unit-forcing response: Fe (surface F) + nodal U,V for the memory forcing
      allocate(x(self%ndof), self%xUn(self%nr), self%xVn(self%nr))
      if (tidal) then
         call radial_operator_solve_vec(self%op, radial_operator_tidal_rhs(self%op, 1.0_wp), x)
      else
         call radial_operator_solve_vec(self%op, radial_operator_load_rhs(self%op, 1.0_wp), x)
      end if
      self%Fe = x(idx_f(self%nr))
      self%Ue = x(idx_u(self%nr))
      do node = 1, self%nr
         self%xUn(node) = x(idx_u(node));  self%xVn(node) = x(idx_v(node))
      end do

      allocate(self%Are(NLAM,self%ne), self%Aim(NLAM,self%ne))
      allocate(self%Bre(NLAM,self%ne), self%Bim(NLAM,self%ne))
      allocate(self%Cre(NLAM,self%ne), self%Cim(NLAM,self%ne))
      self%Are = 0.0_wp; self%Aim = 0.0_wp; self%Bre = 0.0_wp
      self%Bim = 0.0_wp; self%Cre = 0.0_wp; self%Cim = 0.0_wp
      allocate(self%dUn_re(self%nr), self%dUn_im(self%nr), &
               self%dVn_re(self%nr), self%dVn_im(self%nr))
      self%dUn_re = 0.0_wp; self%dUn_im = 0.0_wp
      self%dVn_re = 0.0_wp; self%dVn_im = 0.0_wp
      self%dF = (0.0_wp, 0.0_wp)
   end subroutine channel_init

   subroutine channel_set_dt(self, dt)
      !! Rescale the Maxwell factor for a new Δt (Mk = (μ/η)·Δt); no re-factor.
      type(deg2_channel), intent(inout) :: self
      real(wp),            intent(in)    :: dt
      self%dt = dt;  self%Mk = self%MkPerDt*dt
   end subroutine channel_set_dt

   subroutine channel_begin(self)
      !! Freeze the drift from the entering memory τ_n: solve the load-free memory
      !! forcing −∫τ^V:δε for the real and imaginary parts, storing the surface F
      !! drift (self%dF) and the nodal strain drift (ε_n term for the commit).
      type(deg2_channel), intent(inout) :: self
      real(wp), allocatable :: fre(:), fim(:), xre(:), xim(:)
      integer :: node
      allocate(fre(self%ndof), fim(self%ndof), xre(self%ndof), xim(self%ndof))
      fre = 0.0_wp;  fim = 0.0_wp
      call dissipative_rhs(self%ne, self%r, self%sa, self%sb, self%sc, self%norm, &
                           self%Are, self%Bre, self%Cre, fre)
      call dissipative_rhs(self%ne, self%r, self%sa, self%sb, self%sc, self%norm, &
                           self%Aim, self%Bim, self%Cim, fim)
      call radial_operator_solve_vec(self%op, fre, xre)
      call radial_operator_solve_vec(self%op, fim, xim)
      self%dF = cmplx(xre(idx_f(self%nr)), xim(idx_f(self%nr)), wp)
      self%dU = cmplx(xre(idx_u(self%nr)), xim(idx_u(self%nr)), wp)
      do node = 1, self%nr
         self%dUn_re(node) = xre(idx_u(node));  self%dUn_im(node) = xim(idx_u(node))
         self%dVn_re(node) = xre(idx_v(node));  self%dVn_im(node) = xim(idx_v(node))
      end do
   end subroutine channel_begin

   subroutine channel_commit(self, coeff)
      !! Advance the Maxwell memory (forward-Euler) with the converged forcing
      !! coefficient: total nodal strain = coeff·(unit response) + drift(τ_n).
      type(deg2_channel), intent(inout) :: self
      complex(wp),         intent(in)    :: coeff
      real(wp), allocatable :: Ure(:), Uim(:), Vre(:), Vim(:)
      real(wp) :: cr, ci
      integer  :: node
      allocate(Ure(self%nr), Uim(self%nr), Vre(self%nr), Vim(self%nr))
      cr = real(coeff, wp);  ci = aimag(coeff)
      do node = 1, self%nr
         Ure(node) = cr*self%xUn(node) + self%dUn_re(node)
         Uim(node) = ci*self%xUn(node) + self%dUn_im(node)
         Vre(node) = cr*self%xVn(node) + self%dVn_re(node)
         Vim(node) = ci*self%xVn(node) + self%dVn_im(node)
      end do
      call advance_memory(self%ne, self%mu, self%Mk, Ure, Vre, self%Jr, &
                          self%Are, self%Bre, self%Cre)
      call advance_memory(self%ne, self%mu, self%Mk, Uim, Vim, self%Jr, &
                          self%Aim, self%Bim, self%Cim)
   end subroutine channel_commit

   subroutine channel_destroy(self)
      type(deg2_channel), intent(inout) :: self
      call radial_operator_destroy(self%op)
      if (allocated(self%r))      deallocate(self%r)
      if (allocated(self%mu))     deallocate(self%mu)
      if (allocated(self%Mk))     deallocate(self%Mk)
      if (allocated(self%MkPerDt))deallocate(self%MkPerDt)
      if (allocated(self%xUn))    deallocate(self%xUn)
      if (allocated(self%xVn))    deallocate(self%xVn)
      if (allocated(self%Are))    deallocate(self%Are, self%Aim, self%Bre, &
                                             self%Bim, self%Cre, self%Cim)
      if (allocated(self%dUn_re)) deallocate(self%dUn_re, self%dUn_im, &
                                             self%dVn_re, self%dVn_im)
      self%Fe = 0.0_wp;  self%dF = (0.0_wp, 0.0_wp)
      self%nr = 0;  self%ne = 0;  self%ndof = 0
   end subroutine channel_destroy

end module vilma_rotation
