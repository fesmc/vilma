module vilma_sht
   !! Thin Fortran wrapper around the SHTns spherical-harmonic transform library.
   !!
   !! Horizontal transforms are the performance kernel of the spectral-finite-
   !! element solver (Martinec 2000): every time step maps fields between the
   !! Gauss-Legendre spatial grid (where lateral viscosity and the sea-level
   !! equation live) and the spectral coefficients (where the radial solves are
   !! per-degree and decoupled). This module isolates the SHTns C API behind a
   !! small derived type so the rest of VILMA2 never touches iso_c_binding.
   !!
   !! Convention: fully-normalized real spherical harmonics, no Condon-Shortley
   !! phase (SHT_ORTHONORMAL + SHT_NO_CS_PHASE), on a Gauss grid with a
   !! phi-contiguous spatial layout. Spectral arrays hold m >= 0 only (the field
   !! is real). See doc/design.md for the rationale and the degree-1 / frame
   !! caveats that the load Love numbers will later depend on.
   use, intrinsic :: iso_c_binding
   use vilma_precision, only: wp
   implicit none
   private

   ! SHTns Fortran 2003 interface (installed by `make install` into the SHTns
   ! prefix include dir; found via -I$(SHTNSROOT)/include).
   include 'shtns.f03'

   real(wp), parameter :: pi = acos(-1.0_wp)

   public :: sht_grid

   type :: sht_grid
      !! One configured SHTns transform (a spectral truncation + matching grid).
      type(c_ptr) :: cfg = c_null_ptr   !! opaque SHTns configuration handle
      integer :: lmax = 0               !! max spherical-harmonic degree
      integer :: mmax = 0               !! max order (in units of mres)
      integer :: mres = 1              !! order step
      integer :: nlat = 0              !! # latitudes (Gauss nodes)
      integer :: nphi = 0              !! # longitudes
      integer :: nlm  = 0              !! # spectral coefficients (m >= 0)
      real(c_double) :: eps = 1.0e-10_c_double  !! polar-optimisation threshold (for clone_cfg)
      ! Physical grid geometry, cached at init. Spatial fields are (nphi, nlat).
      real(wp), allocatable :: colat(:)    !! colatitude of each latitude row [rad] (nlat)
      real(wp), allocatable :: lon(:)      !! longitude of each column [rad] (nphi)
      real(wp), allocatable :: gauss_w(:)  !! Gauss quadrature weights [-] (nlat), sum = 2
   end type sht_grid

   public :: sht_free_cfg   !! release a config from clone_cfg
   public :: sht_grid_init, sht_grid_destroy, sht_grid_synthesis, sht_grid_analysis, sht_grid_sph_synthesis, sht_grid_sph_analysis, sht_grid_tor_synthesis, sht_grid_sphtor_synthesis, sht_grid_sphtor_analysis, sht_grid_eval_point, sht_grid_eval_point_horiz, sht_grid_lmidx, sht_grid_surface_integral, sht_grid_clone_cfg

contains

   subroutine sht_grid_init(self, lmax, nlat, nphi, mmax, mres, eps_polar)
      !! Create the SHTns config and its Gauss grid.
      !!
      !! Grid defaults to the smallest that resolves degree `lmax` exactly with
      !! Gauss quadrature (nlat >= lmax+1, nphi >= 2*mmax+1). For quadratic
      !! products (the sea-level equation's ocean-function multiply) the caller
      !! should pass a de-aliased grid: nlat ~ 3*lmax/2, nphi ~ 3*lmax.
      type(sht_grid), intent(inout) :: self
      integer,  intent(in)           :: lmax
      integer,  intent(in), optional :: nlat, nphi, mmax, mres
      real(wp), intent(in), optional :: eps_polar

      type(shtns_info), pointer :: info
      real(c_double) :: eps
      integer :: mmax_, mres_, nlat_, nphi_

      ! SHTns aborts the process ("nlat is too small!") whenever nlat < 4*VSIZE2,
      ! where VSIZE2 is the width of the SIMD vector its *build* uses -- a
      ! property of the machine, not of the transform. It is 1 on a scalar build
      ! and up to 8 on AVX-512 (shtns_simd.h), so 32 clears the ceiling for every
      ! build of SHTns 3.7.5. Below that the on-the-fly kernels would index past
      ! the end of the hemisphere arrays, which is why SHTns refuses rather than
      ! degrading.
      !
      ! The floor is invisible at production resolutions -- nlat = lmax+2 already
      ! exceeds it for lmax >= 30 -- and costs only a slightly finer grid for the
      ! small-lmax unit tests. Gauss quadrature stays exact on the finer grid, so
      ! padding changes no result; it only avoids an abort that otherwise depends
      ! on which machine the code was compiled for.
      integer, parameter :: NLAT_MIN = 32

      mres_ = 1;          if (present(mres)) mres_ = mres
      mmax_ = lmax/mres_; if (present(mmax)) mmax_ = mmax
      eps   = 1.0e-10_c_double
      if (present(eps_polar)) eps = real(eps_polar, c_double)

      nlat_ = lmax + 2;       if (present(nlat)) nlat_ = nlat
      nphi_ = 2*(mmax_ + 1);  if (present(nphi)) nphi_ = nphi
      nlat_ = max(nlat_, NLAT_MIN)    ! see NLAT_MIN above; applies to an explicit
                                      ! nlat too, since SHTns would otherwise abort

      self%lmax = lmax
      self%mmax = mmax_
      self%mres = mres_
      self%eps  = eps

      self%cfg = create_cfg(lmax, mmax_, mres_, nlat_, nphi_, eps)

      ! Read back the realized grid sizes from the SHTns info struct.
      call c_f_pointer(self%cfg, info)
      self%nlat = info%nlat
      self%nphi = info%nphi
      self%nlm  = info%nlm

      call cache_geometry(self, info)
   end subroutine sht_grid_init

   subroutine cache_geometry(self, info)
      !! Cache the physical grid coordinates and quadrature weights from SHTns.
      type(sht_grid),  intent(inout) :: self
      type(shtns_info), intent(in)    :: info
      real(c_double), pointer     :: cos_theta(:)
      real(c_double), allocatable :: wts_half(:)
      integer :: j, nh

      ! Colatitudes from the cos(theta) array SHTns exposes on its struct.
      call c_f_pointer(info%ct, cos_theta, [self%nlat])
      self%colat = acos(real(cos_theta, wp))

      ! Uniform longitudes.
      allocate(self%lon(self%nphi))
      do j = 1, self%nphi
         self%lon(j) = 2.0_wp*pi*real(j-1, wp)/real(self%nphi, wp)
      end do

      ! Gauss weights: SHTns returns one hemisphere (nlat_2); mirror to full grid.
      nh = info%nlat_2
      allocate(wts_half(nh))
      call shtns_gauss_wts(self%cfg, wts_half)
      allocate(self%gauss_w(self%nlat))
      self%gauss_w(1:nh) = real(wts_half, wp)
      self%gauss_w(self%nlat:self%nlat-nh+1:-1) = real(wts_half, wp)
      deallocate(wts_half)
   end subroutine cache_geometry

   function create_cfg(lmax, mmax, mres, nlat, nphi, eps) result(cfg)
      !! Create one SHTns config + Gauss grid with this project's fixed conventions
      !! (orthonormal, no Condon-Shortley phase; Gauss grid, phi-contiguous). Used
      !! both for the primary config and for the per-thread clones (clone_cfg).
      integer,        intent(in) :: lmax, mmax, mres, nlat, nphi
      real(c_double), intent(in) :: eps
      type(c_ptr) :: cfg
      integer :: norm, layout
      norm   = SHT_ORTHONORMAL + SHT_NO_CS_PHASE
      ! SHT_QUICK_INIT (not SHT_GAUSS): skip SHTns's per-config algorithm auto-tuning,
      ! which costs ~17 s/config — prohibitive for the per-thread tensor-SH pool (8
      ! identical configs ⇒ ~136 s of init). Quick-init builds the same Gauss grid with
      ! a default (untuned) algorithm in ~0 s; transforms run ~20% slower than the tuned
      ! optimum. For long production runs, add SHT_LOAD_SAVE_CFG to SHT_GAUSS to recover
      ! tuned transforms with cached (one-time) init.
      layout = SHT_QUICK_INIT + SHT_PHI_CONTIGUOUS
      cfg = shtns_create(lmax, mmax, mres, norm)
      call shtns_set_grid(cfg, layout, eps, nlat, nphi)
   end function create_cfg

   function sht_grid_clone_cfg(self) result(cfg)
      !! Return a fresh, independent SHTns config identical to self%cfg. Distinct
      !! configs can run transforms concurrently (one config is NOT safe for
      !! concurrent calls), so a thread pool of these enables OpenMP over the
      !! element loop in the tensor-SH memory advance. Creation is NOT thread-safe
      !! (FFTW planning) — build the pool serially before any parallel region.
      type(sht_grid), intent(in) :: self
      type(c_ptr) :: cfg
      cfg = create_cfg(self%lmax, self%mmax, self%mres, self%nlat, self%nphi, self%eps)
   end function sht_grid_clone_cfg

   subroutine sht_free_cfg(cfg)
      !! Release a config obtained from clone_cfg.
      type(c_ptr), intent(inout) :: cfg
      if (c_associated(cfg)) call shtns_destroy(cfg)
      cfg = c_null_ptr
   end subroutine sht_free_cfg

   subroutine sht_grid_destroy(self)
      !! Release the SHTns configuration and cached geometry.
      type(sht_grid), intent(inout) :: self
      if (c_associated(self%cfg)) call shtns_destroy(self%cfg)
      self%cfg = c_null_ptr
      self%lmax = 0; self%mmax = 0; self%nlat = 0; self%nphi = 0; self%nlm = 0
      if (allocated(self%colat))   deallocate(self%colat)
      if (allocated(self%lon))     deallocate(self%lon)
      if (allocated(self%gauss_w)) deallocate(self%gauss_w)
   end subroutine sht_grid_destroy

   real(wp) function sht_grid_surface_integral(self, f) result(total)
      !! Surface integral ∫ f dΩ over the unit sphere by Gauss-Legendre
      !! quadrature in latitude and the uniform rule in longitude. Used for the
      !! mass/area integrals the sea-level equation depends on; for f≡1 it
      !! returns 4π.
      type(sht_grid), intent(in) :: self
      real(wp),        intent(in) :: f(:,:)   !! (nphi, nlat)
      real(wp) :: dlon, row
      integer  :: i, j
      dlon  = 2.0_wp*pi/real(self%nphi, wp)
      total = 0.0_wp
      do j = 1, self%nlat
         row = 0.0_wp
         do i = 1, self%nphi
            row = row + f(i, j)
         end do
         total = total + self%gauss_w(j)*row*dlon
      end do
   end function sht_grid_surface_integral

   subroutine sht_grid_synthesis(self, slm, sh, cfg)
      !! Spectral coefficients -> spatial field, shape (nphi, nlat). Pass `cfg` (a
      !! clone_cfg handle) to run on a thread-local config instead of self%cfg.
      type(sht_grid), intent(in)  :: self
      complex(wp),     intent(in)  :: slm(:)   !! length nlm
      real(wp),        intent(out) :: sh(:,:)  !! (nphi, nlat)
      type(c_ptr), intent(in), optional :: cfg
      type(c_ptr) :: c
      c = self%cfg;  if (present(cfg)) c = cfg
      call SH_to_spat(c, slm, sh)
   end subroutine sht_grid_synthesis

   subroutine sht_grid_analysis(self, sh, slm, cfg)
      !! Spatial field -> spectral coefficients.
      !! NB: SHTns overwrites the input spatial array, hence intent(inout).
      type(sht_grid), intent(in)    :: self
      real(wp),        intent(inout) :: sh(:,:)  !! (nphi, nlat)
      complex(wp),     intent(out)   :: slm(:)   !! length nlm
      type(c_ptr), intent(in), optional :: cfg
      type(c_ptr) :: c
      c = self%cfg;  if (present(cfg)) c = cfg
      call spat_to_SH(c, sh, slm)
   end subroutine sht_grid_analysis

   subroutine sht_grid_sph_synthesis(self, slm, vth, vph, cfg)
      !! Spheroidal (gradient) vector synthesis: from the spectral coefficients slm
      !! of a scalar potential, return the surface-gradient field on the grid,
      !!   vth = ∂_θ(Σ slm Y_lm),   vph = (1/sinθ) ∂_φ(Σ slm Y_lm),
      !! i.e. the E_lm and F_lm tensor-harmonic building blocks (Martinec 2000 B11).
      !! Used to bootstrap the tensor-SH dyadic basis (vilma_tensor_sh).
      type(sht_grid), intent(in)  :: self
      complex(wp),     intent(in)  :: slm(:)      !! length nlm
      real(wp),        intent(out) :: vth(:,:), vph(:,:)  !! (nphi, nlat)
      type(c_ptr), intent(in), optional :: cfg
      type(c_ptr) :: c
      c = self%cfg;  if (present(cfg)) c = cfg
      call SHsph_to_spat(c, slm, vth, vph)
   end subroutine sht_grid_sph_synthesis

   subroutine sht_grid_sph_analysis(self, vth, vph, slm, cfg)
      !! Spheroidal vector analysis — the inverse/adjoint of sph_synthesis: from a
      !! horizontal field (vth,vph) recover the spheroidal potential coefficients slm
      !! (the toroidal part is discarded; sphtor_analysis keeps it). NB: SHTns
      !! overwrites the inputs.
      type(sht_grid), intent(in)    :: self
      real(wp),        intent(inout) :: vth(:,:), vph(:,:)  !! (nphi, nlat)
      complex(wp),     intent(out)   :: slm(:)              !! length nlm
      type(c_ptr), intent(in), optional :: cfg
      complex(wp) :: tlm(self%nlm)
      type(c_ptr) :: c
      c = self%cfg;  if (present(cfg)) c = cfg
      call spat_to_SHsphtor(c, vth, vph, slm, tlm)
   end subroutine sht_grid_sph_analysis

   subroutine sht_grid_tor_synthesis(self, tlm, vth, vph, cfg)
      !! Toroidal vector synthesis: from the coefficients tlm of a scalar, return
      !! the rotated surface gradient on the grid,
      !!   (vth, vph) = e_r × ∇₁(Σ tlm Y_lm) = (−F, E),
      !! with E = ∂_θ, F = (1/sinθ)∂_φ as in sph_synthesis. This is Martinec's
      !! S⁽⁰⁾ direction, the one the toroidal displacement W multiplies, and the
      !! building block of the toroidal tensor harmonics Z³, Z⁴ (vilma_tensor_sh).
      !!
      !! SHTns' own toroidal field is ∇×(T e_r) = −e_r × ∇T = (F, −E), the
      !! opposite orientation, so its output is negated here. test_sht pins the
      !! sign against (−F, E) formed from sph_synthesis of the same coefficients.
      type(sht_grid), intent(in)  :: self
      complex(wp),     intent(in)  :: tlm(:)             !! length nlm
      real(wp),        intent(out) :: vth(:,:), vph(:,:) !! (nphi, nlat)
      type(c_ptr), intent(in), optional :: cfg
      type(c_ptr) :: c
      c = self%cfg;  if (present(cfg)) c = cfg
      call SHtor_to_spat(c, tlm, vth, vph)
      vth = -vth;  vph = -vph
   end subroutine sht_grid_tor_synthesis

   subroutine sht_grid_sphtor_synthesis(self, slm, tlm, vth, vph, cfg)
      !! Both vector syntheses in one transform: (vth, vph) = ∇₁S + e_r × ∇₁T, in
      !! this module's toroidal orientation (see tor_synthesis). Costs about one
      !! sph_synthesis, against two for the separate calls.
      type(sht_grid), intent(in)  :: self
      complex(wp),     intent(in)  :: slm(:), tlm(:)     !! length nlm
      real(wp),        intent(out) :: vth(:,:), vph(:,:) !! (nphi, nlat)
      type(c_ptr), intent(in), optional :: cfg
      complex(wp) :: tneg(self%nlm)
      type(c_ptr) :: c
      c = self%cfg;  if (present(cfg)) c = cfg
      tneg = -tlm                                  ! SHTns orientation is −e_r×∇
      call SHsphtor_to_spat(c, slm, tneg, vth, vph)
   end subroutine sht_grid_sphtor_synthesis

   subroutine sht_grid_sphtor_analysis(self, vth, vph, slm, tlm, cfg)
      !! Full horizontal vector analysis, the inverse of the pair above: a field
      !! (vth,vph) = ∇₁S + e_r × ∇₁T is split into its spheroidal slm and toroidal
      !! tlm coefficients, both in this module's orientation (tlm is the negative
      !! of SHTns' own toroidal output; see tor_synthesis). NB: SHTns overwrites
      !! the inputs.
      type(sht_grid), intent(in)    :: self
      real(wp),        intent(inout) :: vth(:,:), vph(:,:)  !! (nphi, nlat)
      complex(wp),     intent(out)   :: slm(:), tlm(:)      !! length nlm
      type(c_ptr), intent(in), optional :: cfg
      type(c_ptr) :: c
      c = self%cfg;  if (present(cfg)) c = cfg
      call spat_to_SHsphtor(c, vth, vph, slm, tlm)
      tlm = -tlm
   end subroutine sht_grid_sphtor_analysis

   subroutine sht_grid_eval_point(self, f_lm, colat, lon, val)
      !! Evaluate a scalar field (spectral coefficients f_lm) at an ARBITRARY point
      !! (colat, lon) [rad] -- not restricted to a grid node. Uses SHTns'
      !! SHqst_to_point with zero spheroidal/toroidal parts, so the normalization
      !! matches analysis/synthesis exactly. Used to sample model fields along the
      !! Martinec-2018 benchmark profiles (circles of constant lon or lat).
      type(sht_grid), intent(in) :: self
      complex(wp),     intent(in) :: f_lm(:)    !! length nlm
      real(wp),        intent(in) :: colat, lon !! [rad]
      real(wp),        intent(out):: val
      complex(wp) :: q(self%nlm), s(self%nlm), t(self%nlm)
      real(wp)    :: vr(1), vt(1), vp(1)
      q = f_lm;  s = (0.0_wp, 0.0_wp);  t = (0.0_wp, 0.0_wp)
      call SHqst_to_point(self%cfg, q, s, t, cos(colat), lon, vr, vt, vp)
      val = vr(1)
   end subroutine sht_grid_eval_point

   subroutine sht_grid_eval_point_horiz(self, s_lm, colat, lon, vth, vph, t_lm)
      !! Evaluate the HORIZONTAL field of a spheroidal scalar (coefficients s_lm)
      !! at an arbitrary point (colat, lon) [rad]: ∇₁(Σ s_lm Y_lm) =
      !! (∂_θ, (1/sinθ)∂_φ)(Σ s_lm Y_lm), via SHqst_to_point with q = 0 and
      !! spheroidal = s_lm. Returns vth = θ-component, vph = φ-component. Used for
      !! the horizontal-displacement (u_θ, u_φ) columns of the Martinec benchmark;
      !! s_lm = response%horizontal output (the per-degree V(a) coefficients).
      !! With t_lm, the toroidal part e_r × ∇₁(Σ t_lm Y_lm) is added, in this
      !! module's orientation (see tor_synthesis; SHTns' own is the negative).
      type(sht_grid), intent(in) :: self
      complex(wp),     intent(in) :: s_lm(:)    !! length nlm (spheroidal V(a))
      real(wp),        intent(in) :: colat, lon !! [rad]
      real(wp),        intent(out):: vth, vph   !! θ-, φ-components
      complex(wp),     intent(in), optional :: t_lm(:)  !! length nlm (toroidal W(a))
      complex(wp) :: q(self%nlm), s(self%nlm), t(self%nlm)
      real(wp)    :: vr(1), vt(1), vp(1)
      q = (0.0_wp, 0.0_wp);  s = s_lm;  t = (0.0_wp, 0.0_wp)
      if (present(t_lm)) t = -t_lm                 ! SHTns orientation is −e_r×∇
      call SHqst_to_point(self%cfg, q, s, t, cos(colat), lon, vr, vt, vp)
      vth = vt(1);  vph = vp(1)
   end subroutine sht_grid_eval_point_horiz

   integer function sht_grid_lmidx(self, l, m) result(lm)
      !! 1-based index into the spectral array for harmonic (l, m).
      type(sht_grid), intent(in) :: self
      integer,         intent(in) :: l, m
      lm = shtns_lmidx(self%cfg, l, m)
   end function sht_grid_lmidx

end module vilma_sht
