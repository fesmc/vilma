program dump_reference
   !! Reference dumper for the Julia port (vilma.jl). Runs the validated
   !! Fortran model through the same call paths as the unit tests and writes
   !! machine-readable NetCDF reference data (ncio), one file per topic.
   !!
   !! Usage:  dump_reference.x [outdir] [item ...]
   !!   outdir  default ./reference (run from the VILMA root)
   !!   items   any of: sht radial ve_degree response sle disc rotation coupling
   !!           martinec visc3d martinec_sle   (default: all)
   !!
   !! All values SI (m, s, kg m^-2, Pa) except explicit time axes in years. The
   !! variable documentation lives in <outdir>/README.md.
   use vilma_precision,       only: wp
   use vilma_constants,       only: pi, grav_G, sec_per_year, kyr, rho_ice, rho_water, rad2deg
   use vilma_earth_structure, only: earth_model, earth_layer, build_M3L70V01, earth_gravity_at, &
                                 earth_n_layers, RHEOL_ELASTIC, RHEOL_MAXWELL, RHEOL_FLUID
   use vilma_radial_fe,       only: radial_mesh, radial_mesh_build, radial_operator, &
                                 radial_operator_assemble, radial_operator_solve_vec, &
                                 radial_operator_load_rhs, radial_operator_tidal_rhs, &
                                 radial_operator_destroy, build_dense_operator, uniq_weight, &
                                 shell_Rk, loading_love, tidal_love, idx_u, idx_v, idx_f, ndof_of, &
                                 radial_fe_finalize
   use vilma_viscoelastic,    only: ve_degree, ve_init, ve_step, ve_destroy, NLAM, SCHEME_FE, SCHEME_TRAP
   use vilma_sht,             only: sht_grid, sht_grid_init, sht_grid_destroy, sht_grid_lmidx, &
                                 sht_grid_synthesis, sht_grid_analysis, sht_grid_sph_synthesis, &
                                 sht_grid_surface_integral, sht_grid_eval_point, &
                                 sht_grid_eval_point_horiz
   use vilma_field,           only: spherical_cap, exp_basin
   use vilma_tensor_sh,       only: tensor_sh, TLAM_SPH, tensor_sh_init, tensor_sh_synth, tensor_sh_analysis, &
                                 tensor_sh_destroy
   use vilma_response,        only: response, response_init_elastic, response_init_ve, &
                                 response_begin_step, response_apply, response_horizontal, &
                                 response_commit_step, response_set_dt, response_destroy, &
                                 response_enable_lateral_visc
   use vilma_sle,             only: sle_solver, sle_result, sle_solve, sle_result_grids
   use vilma_rotation,        only: rotation_state, rotation_init, rotation_update, rotation_destroy
   use vilma_params,          only: vilma_param_class
   use vilma_coupling,        only: solid_earth, solid_earth_init, solid_earth_update, solid_earth_finalize
   use vilma_io,              only: vilma_restart_write, vilma_io_set_table
   use ncio
   implicit none

   character(len=*), parameter :: VILMA_ROOT = "./"   ! run from the VILMA root
   character(len=*), parameter :: EARTH_NAME = "M3-L70-V01"
   integer,  parameter :: LMAX32 = 32
   real(wp), parameter :: DT100 = 100.0_wp*sec_per_year     ! 100 yr in s

   ! --- Martinec-2018 SLE benchmark (item `martinec_sle`, mirrors tests/test_benchmark_sle.f90)
   real(wp), parameter :: MDEG = pi/180.0_wp
   real(wp), parameter :: MS_ICE_COLAT = 25.0_wp*MDEG, MS_ICE_LON = 75.0_wp*MDEG
   real(wp), parameter :: MS_ALPHA = 10.0_wp*MDEG, MS_SIGB = 26.0_wp*MDEG
   real(wp), parameter :: MS_GGIAPY = 9.815_wp              ! giapy reference gravity of col5
   real(wp), parameter :: MS_DT = 0.02_wp*kyr               ! 20 yr
   integer,  parameter :: MS_NGROW = 500, MS_NHOLD = 250
   integer,  parameter :: MS_MAXSPIN = 12
   real(wp), parameter :: MS_SPINTOL = 1.0_wp               ! F1: mean|present - B2| [m]

   character(len=512) :: outdir, arg
   logical :: want(11), all_items
   integer :: nargs, i

   outdir = "./reference"
   nargs = command_argument_count()
   if (nargs >= 1) call get_command_argument(1, outdir)
   want = .false.;  all_items = (nargs < 2)
   do i = 2, nargs
      call get_command_argument(i, arg)
      select case (trim(arg))
      case ("sht");       want(1)  = .true.
      case ("radial");    want(2)  = .true.
      case ("ve_degree"); want(3)  = .true.
      case ("response");  want(4)  = .true.
      case ("sle");       want(5)  = .true.
      case ("disc");      want(6)  = .true.
      case ("rotation");  want(7)  = .true.
      case ("coupling");  want(8)  = .true.
      case ("martinec");  want(9)  = .true.
      case ("visc3d");    want(10) = .true.
      case ("martinec_sle"); want(11) = .true.
      case default
         write(*,'(3a)') ' unknown item "', trim(arg), '"';  error stop 1
      end select
   end do
   if (all_items) want = .true.
   call execute_command_line("mkdir -p '"//trim(outdir)//"'")
   call vilma_io_set_table(VILMA_ROOT//"input/vilma-variables.md")

   if (want(1))  call dump_sht(trim(outdir)//"/sht.nc")
   if (want(2))  call dump_radial(trim(outdir)//"/radial.nc")
   if (want(3))  call dump_ve_degree(trim(outdir)//"/ve_degree.nc")
   if (want(4))  call dump_response(trim(outdir)//"/response.nc")
   if (want(5))  call dump_sle(trim(outdir)//"/sle.nc")
   if (want(6))  call dump_disc(trim(outdir)//"/disc.nc")
   if (want(7))  call dump_rotation(trim(outdir)//"/rotation.nc")
   if (want(8))  call dump_coupling(trim(outdir)//"/coupling.nc", trim(outdir))
   if (want(9))  call dump_martinec(trim(outdir)//"/martinec_A.nc")
   if (want(10)) call dump_visc3d(trim(outdir)//"/visc3d.nc")
   if (want(11)) call dump_martinec_sle(trim(outdir)//"/martinec_sle.nc")

   call radial_fe_finalize()
   write(*,'(a)') ' dump_reference: done'

contains

   ! =====================================================================
   ! helpers
   ! =====================================================================

   subroutine stamp(f, item)
      !! Global provenance attributes common to every file.
      character(len=*), intent(in) :: f, item
      call nc_write_attr(f, "item", item)
      call nc_write_attr(f, "source", "VILMA Fortran tests/dump_reference.f90")
      call nc_write_attr(f, "earth_model", EARTH_NAME)
      call nc_write_attr(f, "grav_G", grav_G)
      call nc_write_attr(f, "sec_per_year", sec_per_year)
      call nc_write_attr(f, "rho_ice", rho_ice)
      call nc_write_attr(f, "rho_water", rho_water)
      call nc_write_attr(f, "units_convention", "SI (m, s, kg m-2, Pa); time axes in years")
   end subroutine stamp

   subroutine w_one(f)
      character(len=*), intent(in) :: f
      call nc_write_dim(f, "one", x=1, dx=1, nx=1, units="1")
   end subroutine w_one

   subroutine ws(f, name, val, units, ln)
      !! scalar real on the "one" dimension
      character(len=*), intent(in) :: f, name, units, ln
      real(wp),         intent(in) :: val
      call nc_write(f, name, val, dim1="one", units=units, long_name=ln)
   end subroutine ws

   subroutine wsi(f, name, val, ln)
      character(len=*), intent(in) :: f, name, ln
      integer,          intent(in) :: val
      call nc_write(f, name, val, dim1="one", units="1", long_name=ln)
   end subroutine wsi

   subroutine wc1(f, name, z, d1, units, ln)
      !! complex 1-D array -> name_re / name_im
      character(len=*), intent(in) :: f, name, d1, units, ln
      complex(wp),      intent(in) :: z(:)
      real(wp), allocatable :: t(:)
      allocate(t(size(z)))
      t = real(z, wp);  call nc_write(f, trim(name)//"_re", t, dim1=d1, units=units, long_name=trim(ln)//" (real part)")
      t = aimag(z);     call nc_write(f, trim(name)//"_im", t, dim1=d1, units=units, long_name=trim(ln)//" (imag part)")
   end subroutine wc1

   subroutine wc2(f, name, z, d1, d2, units, ln)
      character(len=*), intent(in) :: f, name, d1, d2, units, ln
      complex(wp),      intent(in) :: z(:,:)
      real(wp), allocatable :: t(:,:)
      allocate(t(size(z,1), size(z,2)))
      t = real(z, wp);  call nc_write(f, trim(name)//"_re", t, dim1=d1, dim2=d2, units=units, long_name=trim(ln)//" (real part)")
      t = aimag(z);     call nc_write(f, trim(name)//"_im", t, dim1=d1, dim2=d2, units=units, long_name=trim(ln)//" (imag part)")
   end subroutine wc2

   subroutine w_lm_index(f, g)
      !! SHTns coefficient ordering: degree l(lm) and order m(lm) (m-major).
      character(len=*), intent(in) :: f
      type(sht_grid),   intent(in) :: g
      integer, allocatable :: ld(:), md(:)
      integer :: l, m, lm
      allocate(ld(g%nlm), md(g%nlm));  ld = -1;  md = -1
      do m = 0, g%mmax*g%mres, g%mres
         do l = m, g%lmax
            lm = sht_grid_lmidx(g, l, m);  ld(lm) = l;  md(lm) = m
         end do
      end do
      call nc_write(f, "l", ld, dim1="lm", units="1", long_name="degree l of coefficient lm (SHTns order)")
      call nc_write(f, "m", md, dim1="lm", units="1", long_name="order m of coefficient lm (SHTns order)")
   end subroutine w_lm_index

   subroutine w_grid(f, g, lonname, colatname)
      !! Gauss-grid coordinate dimensions in radians (+ weights).
      character(len=*), intent(in) :: f, lonname, colatname
      type(sht_grid),   intent(in) :: g
      call nc_write_dim(f, lonname,   x=g%lon,   units="radians", long_name="longitude of column i")
      call nc_write_dim(f, colatname, x=g%colat, units="radians", long_name="colatitude of Gauss row j")
      call nc_write(f, "gauss_w_"//colatname, g%gauss_w, dim1=colatname, units="1", &
           long_name="Gauss-Legendre weight per latitude row (sum = 2)")
   end subroutine w_grid

   subroutine grid_dims_std(f, g)
      !! standard lmax=32 field grid: lon(nphi), colat(nlat), lm(nlm)
      character(len=*), intent(in) :: f
      type(sht_grid),   intent(in) :: g
      call w_grid(f, g, "lon", "colat")
      call nc_write_dim(f, "lm", x=1, dx=1, nx=g%nlm, units="1")
      call w_lm_index(f, g)
      call nc_write_attr(f, "lmax", g%lmax);  call nc_write_attr(f, "mmax", g%mmax)
      call nc_write_attr(f, "nlat", g%nlat);  call nc_write_attr(f, "nphi", g%nphi)
      call nc_write_attr(f, "nlm",  g%nlm)
      call nc_write_attr(f, "sht_norm", "SHT_ORTHONORMAL + SHT_NO_CS_PHASE, real field, m>=0 stored")
   end subroutine grid_dims_std

   real(wp) function frand(s) result(r)
      !! LCG from test_tensor_sh (deterministic pseudo-random in [-1,1])
      integer, intent(inout) :: s
      s = mod(1103515245*s + 12345, 2147483647)
      r = 2.0_wp*real(s,wp)/2147483647.0_wp - 1.0_wp
   end function frand

   subroutine fluidise(e)
      !! Relax every Maxwell layer to an inviscid fluid (mu=0): the t->inf limit.
      type(earth_model), intent(inout) :: e
      integer :: k
      do k = 1, earth_n_layers(e)
         if (e%layers(k)%rheology == RHEOL_MAXWELL) then
            e%layers(k)%mu = 0.0_wp;  e%layers(k)%rheology = RHEOL_FLUID
         end if
      end do
   end subroutine fluidise

   ! =====================================================================
   ! 1. sht.nc
   ! =====================================================================

   subroutine dump_sht(f)
      character(len=*), intent(in) :: f
      type(sht_grid)  :: g
      type(tensor_sh) :: tsh
      real(wp),    allocatable :: fld(:,:), tmp(:,:), vth(:,:), vph(:,:), dyad(:,:,:)
      complex(wp), allocatable :: slm0(:), flm(:), rt(:), c(:,:), c2(:,:)
      integer  :: l, m, lm, i, j, lam, seed
      real(wp) :: th, ph, integ

      write(*,'(a)') ' [sht] lmax=32 Gauss grid 64x128, SHT conventions'
      call sht_grid_init(g, LMAX32, nlat=2*LMAX32, nphi=4*LMAX32)
      call nc_create(f, overwrite=.true.)
      call stamp(f, "sht")
      call w_one(f)
      call grid_dims_std(f, g)
      call nc_write_dim(f, "lam",    x=1, dx=1, nx=TLAM_SPH, units="1")
      call nc_write_dim(f, "comp",   x=1, dx=1, nx=6,    units="1")
      call nc_write_dim(f, "degree", x=0, dx=1, nx=LMAX32+1, units="1")
      call nc_write_attr(f, "tensor_lam_map", "lam 1..4 = Martinec lambda 1,2,5,6")
      call nc_write_attr(f, "dyad_comp_map", "comp 1..6 = rr, r-theta, r-phi, theta-theta, theta-phi, phi-phi")

      allocate(fld(g%nphi,g%nlat), tmp(g%nphi,g%nlat), vth(g%nphi,g%nlat), vph(g%nphi,g%nlat))
      allocate(slm0(g%nlm), flm(g%nlm), rt(g%nlm))

      ! analytic real test field on the grid
      do j = 1, g%nlat
         th = g%colat(j)
         do i = 1, g%nphi
            ph = g%lon(i)
            fld(i,j) = 3.0_wp + 2.0_wp*cos(th) + sin(th)*cos(ph) + sin(th)**2*cos(2.0_wp*ph) &
                     + 0.5_wp*sin(th)*cos(th)*sin(ph) + 0.25_wp*cos(3.0_wp*th)
         end do
      end do
      call nc_write(f, "f_test", fld, dim1="lon", dim2="colat", units="1", &
           long_name="analytic test field 3+2cos(th)+sin(th)cos(ph)+sin^2(th)cos(2ph)+0.5sin(th)cos(th)sin(ph)+0.25cos(3th)")
      tmp = fld
      call sht_grid_analysis(g, tmp, flm)
      call wc1(f, "f_test_lm", flm, "lm", "1", "spat_to_SH coefficients of f_test")
      integ = sht_grid_surface_integral(g, fld)
      call ws(f, "f_test_integral", integ, "1", "sht_grid_surface_integral(f_test) = int f dOmega")

      ! deterministic coefficient vector (test_sht), synthesis, round trip
      slm0 = (0.0_wp, 0.0_wp)
      do l = 0, LMAX32
         do m = 0, l
            lm = sht_grid_lmidx(g, l, m)
            if (m == 0) then
               slm0(lm) = cmplx(1.0_wp/real(l+1, wp), 0.0_wp, wp)
            else
               slm0(lm) = cmplx(1.0_wp/real(l+1, wp), 0.1_wp*real(m, wp)/real(l+1, wp), wp)
            end if
         end do
      end do
      call wc1(f, "c_test", slm0, "lm", "1", "test coefficient vector: 1/(l+1) + i 0.1 m/(l+1)")
      call sht_grid_synthesis(g, slm0, fld)
      call nc_write(f, "c_test_grid", fld, dim1="lon", dim2="colat", units="1", &
           long_name="SH_to_spat synthesis of c_test")
      tmp = fld
      call sht_grid_analysis(g, tmp, rt)
      call wc1(f, "c_test_rt", rt, "lm", "1", "spat_to_SH(SH_to_spat(c_test)) round trip")
      call sht_grid_sph_synthesis(g, slm0, vth, vph)
      call nc_write(f, "c_test_vth", vth, dim1="lon", dim2="colat", units="1", &
           long_name="SHsph_to_spat theta component of c_test: d/dtheta sum c Y")
      call nc_write(f, "c_test_vph", vph, dim1="lon", dim2="colat", units="1", &
           long_name="SHsph_to_spat phi component of c_test: (1/sin theta) d/dphi sum c Y")

      ! tensor-SH dyadic synthesis of a deterministic coefficient set (test_tensor_sh)
      call tensor_sh_init(tsh, g)
      allocate(c(TLAM_SPH,g%nlm), c2(TLAM_SPH,g%nlm), dyad(g%nphi,g%nlat,6))
      seed = 1;  c = (0.0_wp, 0.0_wp)
      do m = 0, LMAX32
         do l = m, LMAX32
            lm = sht_grid_lmidx(g, l, m)
            do lam = 1, TLAM_SPH
               if (lam == 2 .and. l < 1) cycle
               if (lam == 3 .and. l < 1) cycle
               if (lam == 4 .and. l < 2) cycle
               if (m == 0) then
                  c(lam,lm) = cmplx(frand(seed), 0.0_wp, wp)
               else
                  c(lam,lm) = cmplx(frand(seed), frand(seed), wp)
               end if
            end do
         end do
      end do
      call wc2(f, "tensor_c", c, "lam", "lm", "1", "tensor-SH coefficients (lam, lm), LCG seed=1 as test_tensor_sh")
      call tensor_sh_synth(tsh, g, c, dyad)
      call nc_write(f, "tensor_dyad", dyad, dim1="lon", dim2="colat", dim3="comp", units="1", &
           long_name="tensor_sh_synth of tensor_c: six dyadic grid components")
      call tensor_sh_analysis(tsh, g, dyad, c2)
      call wc2(f, "tensor_c_rt", c2, "lam", "lm", "1", "tensor_sh_analysis(tensor_sh_synth(tensor_c))")
      call nc_write(f, "n6", tsh%n6, dim1="degree", units="1", &
           long_name="calibrated spin-2 channel norm n6(l) (adjoint S6* S6 diagonal)")
      call tensor_sh_destroy(tsh)
      call sht_grid_destroy(g)
   end subroutine dump_sht

   ! =====================================================================
   ! 2. radial.nc
   ! =====================================================================

   subroutine dump_radial(f)
      character(len=*), intent(in) :: f
      integer, parameter :: NJ = 8, JMAX = 64
      integer, parameter :: jl(NJ) = [1, 2, 3, 5, 10, 20, 32, 64]
      type(earth_model)     :: e, ef
      type(radial_mesh)     :: mesh
      type(radial_operator) :: op
      real(wp), allocatable :: A(:,:), w(:), x(:), b(:), Rk(:)
      real(wp), allocatable :: rho_e(:), mu_e(:), eta_e(:), lrhs(:,:), lsol(:,:), trhs(:,:), tsol(:,:)
      real(wp), allocatable :: he(:), le(:), ke(:), hf(:), lf(:), kf(:)
      real(wp), allocatable :: rhe(:), rle(:), rke(:), rhf(:), rlf(:), rkf(:)
      integer,  allocatable :: rheol_e(:), rows(:), cols(:)
      real(wp), allocatable :: vals(:)
      integer  :: nr, ne, nd, ij, j, i, k, nnz, p, lay
      real(wp) :: ua, va, fa, h, l, kk, g
      character(len=8) :: js
      logical :: okread

      write(*,'(a)') ' [radial] mesh, dense operators, RHS/solutions, Love numbers'
      e = build_M3L70V01();  call radial_mesh_build(mesh, e)
      nr = mesh%nr;  ne = mesh%ne;  nd = ndof_of(nr)
      g = earth_gravity_at(e, e%r_earth)

      call nc_create(f, overwrite=.true.)
      call stamp(f, "radial")
      call w_one(f)
      call nc_write_dim(f, "node",   x=1, dx=1, nx=nr, units="1")
      call nc_write_dim(f, "elem",   x=1, dx=1, nx=ne, units="1")
      call nc_write_dim(f, "row",    x=1, dx=1, nx=nd, units="1")
      call nc_write_dim(f, "col",    x=1, dx=1, nx=nd, units="1")
      call nc_write_dim(f, "jlist",  x=jl, units="1", long_name="degrees with operator dumps")
      call nc_write_dim(f, "degree", x=1, dx=1, nx=JMAX, units="1")
      call nc_write_dim(f, "layer",  x=1, dx=1, nx=earth_n_layers(e), units="1")
      call nc_write_attr(f, "nr", nr);  call nc_write_attr(f, "ne", ne);  call nc_write_attr(f, "ndof", nd)
      call nc_write_attr(f, "dof_layout", "node-interleaved: U_k=4(k-1)+1, V_k=4(k-1)+2, F_k=4(k-1)+3, Pi_e=4e; ndof=4nr-1")
      call ws(f, "r_earth", e%r_earth, "m", "surface radius a")
      call ws(f, "r_core",  e%r_core,  "m", "core-mantle boundary radius")
      call ws(f, "g_surf",  g, "m s-2", "earth_gravity_at(a) = G M / a^2")

      ! layers (surface-first) + mesh
      block
         real(wp), allocatable :: lrb(:), lrt(:), lrho(:), lmu(:), leta(:)
         integer,  allocatable :: lrh(:)
         allocate(lrb(size(e%layers)), lrt(size(e%layers)), lrho(size(e%layers)), &
                  lmu(size(e%layers)), leta(size(e%layers)), lrh(size(e%layers)))
         do k = 1, earth_n_layers(e)
            lrb(k) = e%layers(k)%r_bot;  lrt(k) = e%layers(k)%r_top;  lrho(k) = e%layers(k)%rho
            lmu(k) = e%layers(k)%mu;     leta(k) = e%layers(k)%eta;   lrh(k)  = e%layers(k)%rheology
         end do
         call nc_write(f, "layer_r_bot", lrb, dim1="layer", units="m", long_name="layer inner radius (surface-first)")
         call nc_write(f, "layer_r_top", lrt, dim1="layer", units="m", long_name="layer outer radius")
         call nc_write(f, "layer_rho",   lrho, dim1="layer", units="kg m-3", long_name="layer density")
         call nc_write(f, "layer_mu",    lmu, dim1="layer", units="Pa", long_name="layer shear modulus")
         call nc_write(f, "layer_eta",   leta, dim1="layer", units="Pa s", long_name="layer viscosity (huge(1d0) = elastic, 0 = fluid)")
         call nc_write(f, "layer_rheology", lrh, dim1="layer", units="1", long_name="0=elastic 1=Maxwell 2=fluid")
      end block
      call nc_write(f, "r", mesh%r, dim1="node", units="m", long_name="radial mesh node radii (ascending, r(1)=0, r(nr)=a)")
      call nc_write(f, "elem_layer", mesh%elem_layer, dim1="elem", units="1", long_name="earth layer index of each element")
      allocate(rho_e(ne), mu_e(ne), eta_e(ne), rheol_e(ne))
      do i = 1, ne
         lay = mesh%elem_layer(i)
         rho_e(i) = e%layers(lay)%rho;  mu_e(i) = e%layers(lay)%mu
         eta_e(i) = e%layers(lay)%eta;  rheol_e(i) = e%layers(lay)%rheology
      end do
      call nc_write(f, "elem_rho", rho_e, dim1="elem", units="kg m-3", long_name="element density")
      call nc_write(f, "elem_mu",  mu_e,  dim1="elem", units="Pa",     long_name="element shear modulus")
      call nc_write(f, "elem_eta", eta_e, dim1="elem", units="Pa s",   long_name="element viscosity (huge(1d0)=elastic, 0=fluid)")
      call nc_write(f, "elem_rheology", rheol_e, dim1="elem", units="1", long_name="0=elastic 1=Maxwell 2=fluid")
      Rk = shell_Rk(e, mesh)
      call nc_write(f, "Rk", Rk, dim1="elem", units="kg m-3 m3", long_name="shell_Rk: accumulated density-jump moment (Martinec eq 77)")
      w = uniq_weight(mesh)
      call nc_write(f, "w_uniq", w, dim1="row", units="m3", long_name="degree-1 KKT border / E_uniq weight vector uniq_weight(mesh) (ndof)")

      ! per-degree operators, RHS, solutions
      allocate(lrhs(nd,NJ), lsol(nd,NJ), trhs(nd,NJ), tsol(nd,NJ), x(nd))
      do ij = 1, NJ
         j = jl(ij);  write(js,'(i0)') j
         A = build_dense_operator(e, mesh, j, with_uniq=.false.)
         if (j <= 2) then
            call nc_write(f, "A_dense_j"//trim(js), A, dim1="row", dim2="col", units="mixed SI", &
                 long_name="build_dense_operator(j="//trim(js)//", with_uniq=.false.), before equilibration")
         end if
         nnz = count(A /= 0.0_wp)
         allocate(rows(nnz), cols(nnz), vals(nnz));  p = 0
         do k = 1, nd
            do i = 1, nd
               if (A(i,k) == 0.0_wp) cycle
               p = p + 1;  rows(p) = i;  cols(p) = k;  vals(p) = A(i,k)
            end do
         end do
         call nc_write_dim(f, "nnz_j"//trim(js), x=1, dx=1, nx=nnz, units="1")
         call nc_write(f, "A_row_j"//trim(js), rows, dim1="nnz_j"//trim(js), units="1", long_name="COO row index (1-based) of A(j="//trim(js)//")")
         call nc_write(f, "A_col_j"//trim(js), cols, dim1="nnz_j"//trim(js), units="1", long_name="COO col index (1-based) of A(j="//trim(js)//")")
         call nc_write(f, "A_val_j"//trim(js), vals, dim1="nnz_j"//trim(js), units="mixed SI", long_name="COO value of A(j="//trim(js)//") (column-major order)")
         deallocate(rows, cols, vals)

         call radial_operator_assemble(op, e, mesh, j)
         b = radial_operator_load_rhs(op, 1.0_wp);   lrhs(:,ij) = b
         call radial_operator_solve_vec(op, b, x);   lsol(:,ij) = x
         b = radial_operator_tidal_rhs(op, 1.0_wp);  trhs(:,ij) = b
         call radial_operator_solve_vec(op, b, x);   tsol(:,ij) = x
         call radial_operator_destroy(op)
      end do
      call nc_write(f, "load_rhs", lrhs, dim1="row", dim2="jlist", units="mixed SI", &
           long_name="radial_operator_load_rhs(sigma=1 kg m-2) per degree in jlist")
      call nc_write(f, "load_sol", lsol, dim1="row", dim2="jlist", units="mixed SI (U,V,F nodal; Pi elements)", &
           long_name="solution x of A x = load_rhs (radial_operator_solve_vec; j=1 via KKT border)")
      call nc_write(f, "tidal_rhs", trhs, dim1="row", dim2="jlist", units="mixed SI", &
           long_name="radial_operator_tidal_rhs(phi_t=1 m2 s-2) per degree in jlist")
      call nc_write(f, "tidal_sol", tsol, dim1="row", dim2="jlist", units="mixed SI", &
           long_name="solution x of A x = tidal_rhs")

      ! Love numbers elastic + fluid, j=1..64
      allocate(he(JMAX), le(JMAX), ke(JMAX), hf(JMAX), lf(JMAX), kf(JMAX))
      ef = e;  call fluidise(ef)
      do j = 1, JMAX
         call radial_operator_assemble(op, e, mesh, j)
         call radial_operator_solve_vec(op, radial_operator_load_rhs(op, 1.0_wp), x)
         ua = x(idx_u(nr));  va = x(idx_v(nr));  fa = x(idx_f(nr))
         call loading_love(e, j, 1.0_wp, ua, va, fa, h, l, kk)
         he(j) = h;  le(j) = l;  ke(j) = kk
         call radial_operator_destroy(op)
         call radial_operator_assemble(op, ef, mesh, j)
         call radial_operator_solve_vec(op, radial_operator_load_rhs(op, 1.0_wp), x)
         ua = x(idx_u(nr));  va = x(idx_v(nr));  fa = x(idx_f(nr))
         call loading_love(ef, j, 1.0_wp, ua, va, fa, h, l, kk)
         hf(j) = h;  lf(j) = l;  kf(j) = kk
         call radial_operator_destroy(op)
      end do
      call nc_write(f, "h_el", he, dim1="degree", units="1", long_name="elastic loading Love number h (loading_love)")
      call nc_write(f, "l_el", le, dim1="degree", units="1", long_name="elastic loading Love number l = g V(a)/phi_L")
      call nc_write(f, "k_el", ke, dim1="degree", units="1", long_name="elastic loading Love number k = -F(a)/phi_L - 1")
      call nc_write(f, "h_fl", hf, dim1="degree", units="1", long_name="fluid-limit loading h (Maxwell layers fluidised, mu=0)")
      call nc_write(f, "l_fl", lf, dim1="degree", units="1", long_name="fluid-limit loading l")
      call nc_write(f, "k_fl", kf, dim1="degree", units="1", long_name="fluid-limit loading k")

      ! tidal degree-2 Love numbers elastic + fluid
      call radial_operator_assemble(op, e, mesh, 2)
      call radial_operator_solve_vec(op, radial_operator_tidal_rhs(op, 1.0_wp), x)
      call tidal_love(e, 2, 1.0_wp, x(idx_u(nr)), x(idx_v(nr)), x(idx_f(nr)), h, l, kk)
      call ws(f, "hT2_el", h,  "1", "elastic degree-2 tidal Love number h^T (tidal_love)")
      call ws(f, "lT2_el", l,  "1", "elastic degree-2 tidal Love number l^T")
      call ws(f, "kT2_el", kk, "1", "elastic degree-2 tidal Love number k^T")
      call radial_operator_destroy(op)
      call radial_operator_assemble(op, ef, mesh, 2)
      call radial_operator_solve_vec(op, radial_operator_tidal_rhs(op, 1.0_wp), x)
      call tidal_love(ef, 2, 1.0_wp, x(idx_u(nr)), x(idx_v(nr)), x(idx_f(nr)), h, l, kk)
      call ws(f, "hT2_fl", h,  "1", "fluid-limit degree-2 tidal h^T")
      call ws(f, "lT2_fl", l,  "1", "fluid-limit degree-2 tidal l^T")
      call ws(f, "kT2_fl", kk, "1", "fluid-limit degree-2 tidal k^T = k_s (secular)")
      call radial_operator_destroy(op)

      ! benchmark normal-mode Love table (giapy mod_M3-L70-V01), degrees 1..64
      allocate(rhe(256), rle(256), rke(256), rhf(256), rlf(256), rkf(256))
      rhe = 0;  rle = 0;  rke = 0;  rhf = 0;  rlf = 0;  rkf = 0
      call read_love_ref(VILMA_ROOT//"data/benchmarks/love_M3-L70-V01/mod_M3-L70-V01", &
                         rhe, rle, rke, rhf, rlf, rkf, okread)
      if (okread) then
         call nc_write(f, "h_el_ref", rhe(1:JMAX), dim1="degree", units="1", long_name="benchmark table elastic h (TABOO/ALMA normal modes)")
         call nc_write(f, "l_el_ref", rle(1:JMAX), dim1="degree", units="1", long_name="benchmark table elastic l")
         call nc_write(f, "k_el_ref", rke(1:JMAX), dim1="degree", units="1", long_name="benchmark table elastic k")
         call nc_write(f, "h_fl_ref", rhf(1:JMAX), dim1="degree", units="1", long_name="benchmark table fluid h")
         call nc_write(f, "l_fl_ref", rlf(1:JMAX), dim1="degree", units="1", long_name="benchmark table fluid l")
         call nc_write(f, "k_fl_ref", rkf(1:JMAX), dim1="degree", units="1", long_name="benchmark table fluid k")
      else
         write(*,'(a)') '   WARNING: could not read the benchmark Love table; *_ref omitted'
      end if
   end subroutine dump_radial

   subroutine read_love_ref(fname, he, le, ke, hf, lf, kf, okread)
      !! Parser copied from test_benchmark_love.
      character(len=*), intent(in)  :: fname
      real(wp),         intent(out) :: he(:), le(:), ke(:), hf(:), lf(:), kf(:)
      logical,          intent(out) :: okread
      integer  :: u, i, ni, nm, mm, ios
      real(wp) :: kv, hv, lv, dum
      okread = .false.
      he = 0; le = 0; ke = 0; hf = 0; lf = 0; kf = 0
      open(newunit=u, file=fname, status='old', action='read', iostat=ios)
      if (ios /= 0) return
      do i = 1, 5
         read(u,*,iostat=ios)
         if (ios /= 0) then;  close(u);  return;  end if
      end do
      do
         read(u,*,iostat=ios) ni, nm, kv, hv, lv
         if (ios /= 0) exit
         if (ni < 1 .or. ni > size(he)) exit
         ke(ni) = kv;  he(ni) = hv;  le(ni) = lv
         do mm = 1, nm
            read(u,*,iostat=ios) dum
            if (ios /= 0) then;  close(u);  return;  end if
         end do
         read(u,*,iostat=ios) ni, nm, kv, hv, lv
         if (ios /= 0) then;  close(u);  return;  end if
         kf(ni) = kv;  hf(ni) = hv;  lf(ni) = lv
      end do
      close(u)
      okread = .true.
   end subroutine read_love_ref

   ! =====================================================================
   ! 3. ve_degree.nc
   ! =====================================================================

   subroutine dump_ve_degree(f)
      character(len=*), intent(in) :: f
      integer, parameter :: NSTEP = 200, NMEM = 10
      integer, parameter :: jl(2) = [2, 10]
      type(earth_model) :: e
      type(radial_mesh) :: mesh
      type(ve_degree)   :: vd
      real(wp), allocatable :: tyr(:), U(:), V(:), Fa(:)
      real(wp) :: t
      integer  :: ij, j, is, isch
      character(len=8)  :: js
      character(len=8)  :: sn(2) = [character(len=8) :: "fe", "trap"]

      write(*,'(a)') ' [ve_degree] 1-D Maxwell stepper j=2,10; fe + trap; dt=100 yr, 200 steps'
      e = build_M3L70V01();  call radial_mesh_build(mesh, e)
      call nc_create(f, overwrite=.true.)
      call stamp(f, "ve_degree")
      call w_one(f)
      allocate(tyr(NSTEP), U(NSTEP), V(NSTEP), Fa(NSTEP))
      do is = 1, NSTEP;  tyr(is) = real(is-1, wp)*100.0_wp;  end do
      call nc_write_dim(f, "time", x=tyr, units="years", long_name="report time of each ve_step (state before advancing)")
      call nc_write_dim(f, "lam",  x=1, dx=1, nx=NLAM,    units="1")
      call nc_write_dim(f, "elem", x=1, dx=1, nx=mesh%ne, units="1")
      call nc_write_dim(f, "node", x=1, dx=1, nx=mesh%nr, units="1")
      call ws(f, "dt", DT100, "s", "time step (100 yr)")
      call ws(f, "sigma", 1.0_wp, "kg m-2", "held unit load coefficient")
      call nc_write_attr(f, "trap_settings", "scheme=SCHEME_TRAP, max_couple_iter=50, couple_tol=1e-10")
      call nc_write_attr(f, "fe_settings", "scheme=SCHEME_FE, max_couple_iter=1")

      do ij = 1, 2
         j = jl(ij);  write(js,'(i0)') j
         do isch = 1, 2
            call ve_init(vd, e, mesh, j, DT100)
            ! ve_init does not reset the scheme knobs: set them explicitly for BOTH runs
            if (isch == 2) then
               vd%scheme = SCHEME_TRAP;  vd%max_couple_iter = 50;  vd%couple_tol = 1.0e-10_wp
            else
               vd%scheme = SCHEME_FE;    vd%max_couple_iter = 1
            end if
            do is = 1, NSTEP
               call ve_step(vd, 1.0_wp, t, U(is), V(is), Fa(is))
               if (abs(t - tyr(is)*sec_per_year) > 1.0e-6_wp*sec_per_year) &
                  error stop 'dump_ve_degree: time axis mismatch'
            end do
            call nc_write(f, "U_j"//trim(js)//"_"//trim(sn(isch)), U,  dim1="time", units="m", &
                 long_name="surface radial displacement U(a,t), degree "//trim(js)//", scheme "//trim(sn(isch)))
            call nc_write(f, "V_j"//trim(js)//"_"//trim(sn(isch)), V,  dim1="time", units="m", &
                 long_name="surface spheroidal V(a,t)")
            call nc_write(f, "F_j"//trim(js)//"_"//trim(sn(isch)), Fa, dim1="time", units="m2 s-2", &
                 long_name="surface potential coefficient F(a,t) (Martinec phi_1; N=-F/g)")
            call ve_destroy(vd)
         end do
         ! memory arrays after NMEM fe steps
         call ve_init(vd, e, mesh, j, DT100)
         vd%scheme = SCHEME_FE;  vd%max_couple_iter = 1
         do is = 1, NMEM
            call ve_step(vd, 1.0_wp, t, U(1), V(1), Fa(1))
         end do
         call nc_write(f, "Am_j"//trim(js), vd%Am, dim1="lam", dim2="elem", units="Pa m", &
              long_name="memory coefficient A (eq 109) after 10 fe steps, degree "//trim(js))
         call nc_write(f, "Bm_j"//trim(js), vd%Bm, dim1="lam", dim2="elem", units="Pa m", &
              long_name="memory coefficient B after 10 fe steps")
         call nc_write(f, "Cm_j"//trim(js), vd%Cm, dim1="lam", dim2="elem", units="Pa m", &
              long_name="memory coefficient C after 10 fe steps")
         call nc_write(f, "Un_j"//trim(js), vd%Un, dim1="node", units="m", &
              long_name="nodal U reported at step 10 (strain used for the 10th memory advance)")
         call nc_write(f, "Vn_j"//trim(js), vd%Vn, dim1="node", units="m", &
              long_name="nodal V reported at step 10")
         call nc_write(f, "Mk_j"//trim(js), vd%Mk, dim1="elem", units="1", long_name="element Maxwell factor M = mu dt / eta")
         call nc_write(f, "mu_j"//trim(js), vd%mu, dim1="elem", units="Pa", long_name="element shear modulus")
         call ve_destroy(vd)
      end do
   end subroutine dump_ve_degree

   ! =====================================================================
   ! shared held-load spectrum for items 4, 10
   ! =====================================================================

   subroutine held_load(g, slm)
      type(sht_grid), intent(in)  :: g
      complex(wp),    intent(out) :: slm(:)
      slm = (0.0_wp, 0.0_wp)
      slm(sht_grid_lmidx(g, 2, 1)) = cmplx(1.0_wp,  0.5_wp, wp)
      slm(sht_grid_lmidx(g, 5, 3)) = cmplx(0.7_wp, -0.3_wp, wp)
   end subroutine held_load

   ! =====================================================================
   ! 4. response.nc
   ! =====================================================================

   subroutine dump_response(f)
      character(len=*), intent(in) :: f
      integer, parameter :: NSTEP = 20
      type(sht_grid)    :: g
      type(earth_model) :: e
      type(response)    :: ve
      complex(wp), allocatable :: slm(:), ulm(:), nlm(:), vlm(:), us(:,:), ns(:,:), vs(:,:)
      real(wp), allocatable :: tyr(:)
      integer :: is

      write(*,'(a)') ' [response] ve_response field driver lmax=32, held (2,1)+(5,3), 20 x 100 yr'
      call sht_grid_init(g, LMAX32, nlat=2*LMAX32, nphi=4*LMAX32)
      e = build_M3L70V01()
      call response_init_ve(ve, e, g, DT100)
      allocate(slm(g%nlm), ulm(g%nlm), nlm(g%nlm), vlm(g%nlm))
      allocate(us(g%nlm,NSTEP), ns(g%nlm,NSTEP), vs(g%nlm,NSTEP), tyr(NSTEP))
      call held_load(g, slm)
      do is = 1, NSTEP
         tyr(is) = ve%time/sec_per_year
         call response_begin_step(ve, g)
         call response_apply(ve, g, slm, ulm, nlm)
         call response_horizontal(ve, g, slm, vlm)
         us(:,is) = ulm;  ns(:,is) = nlm;  vs(:,is) = vlm
         call response_commit_step(ve, g, slm)
      end do

      call nc_create(f, overwrite=.true.)
      call stamp(f, "response")
      call w_one(f)
      call grid_dims_std(f, g)
      call nc_write_dim(f, "time", x=tyr, units="years", long_name="report time (state before commit)")
      call nc_write_dim(f, "degree", x=0, dx=1, nx=LMAX32+1, units="1")
      call nc_write_dim(f, "degree1", x=1, dx=1, nx=LMAX32, units="1")
      call nc_write_dim(f, "node", x=1, dx=1, nx=ve%nr, units="1")
      call ws(f, "dt", DT100, "s", "time step")
      call ws(f, "g_surf", ve%g, "m s-2", "surface gravity used for N = -F/g")
      call nc_write_attr(f, "scheme", "fe (SCHEME_FE), skip_tol=1e-4 (default)")
      call wc1(f, "sigma", slm, "lm", "kg m-2", "held load: (2,1)=1+0.5i, (5,3)=0.7-0.3i")
      call wc2(f, "u_lm", us, "lm", "time", "m", "radial displacement coefficients (response_apply)")
      call wc2(f, "N_lm", ns, "lm", "time", "m", "geoid coefficients N = -F/g (N_1 = 0, CM frame)")
      call wc2(f, "V_lm", vs, "lm", "time", "m", "spheroidal horizontal V(a) coefficients (response_horizontal)")
      call nc_write(f, "gu", ve%gu, dim1="degree", units="m / (kg m-2)", long_name="elastic gain U(a) per unit sigma_l")
      call nc_write(f, "gn", ve%gn, dim1="degree", units="m / (kg m-2)", long_name="elastic gain N(a) per unit sigma_l (gn(1)=0)")
      call nc_write(f, "gv", ve%gv, dim1="degree", units="m / (kg m-2)", long_name="elastic gain V(a) per unit sigma_l")
      call nc_write_dim(f, "elem", x=1, dx=1, nx=ve%ne, units="1")
      call nc_write(f, "Mk", ve%Mk, dim1="elem", units="1", long_name="element Maxwell factor M = mu dt/eta")
      call nc_write(f, "xUn", ve%xUn, dim1="node", dim2="degree1", units="m / (kg m-2)", long_name="unit-load nodal U per degree l=1..lmax")
      call nc_write(f, "xVn", ve%xVn, dim1="node", dim2="degree1", units="m / (kg m-2)", long_name="unit-load nodal V per degree l=1..lmax")
      call response_destroy(ve);  call sht_grid_destroy(g)
   end subroutine dump_response

   ! =====================================================================
   ! 5. sle.nc
   ! =====================================================================

   subroutine dump_sle(f)
      character(len=*), intent(in) :: f
      integer, parameter :: NSTEP = 10
      type(sht_grid)    :: g
      type(earth_model) :: e
      type(response)    :: el, ve
      type(sle_solver)  :: sle
      type(sle_result)  :: res
      real(wp), allocatable :: topo0(:,:), d_ice(:,:), ice(:,:), S(:,:), C(:,:), ug(:,:), Ng(:,:)
      real(wp), allocatable :: rsl(:,:,:), Cs(:,:,:), us(:,:,:), Ns(:,:,:)
      real(wp), allocatable :: esl(:), mres(:), resid(:), tyr(:)
      integer,  allocatable :: ninner(:), nouter(:)
      integer  :: i, j, is
      real(wp) :: th

      write(*,'(a)') ' [sle] sea-level equation lmax=32, elastic + ve (10 x 100 yr)'
      call sht_grid_init(g, LMAX32, nlat=2*LMAX32, nphi=4*LMAX32)
      e = build_M3L70V01()
      allocate(topo0(g%nphi,g%nlat), d_ice(g%nphi,g%nlat), ice(g%nphi,g%nlat), S(g%nphi,g%nlat), C(g%nphi,g%nlat))
      do j = 1, g%nlat
         th = g%colat(j)
         do i = 1, g%nphi
            topo0(i,j) = merge(500.0_wp, -4000.0_wp, th < 60.0_wp*pi/180.0_wp)
            d_ice(i,j) = merge(2000.0_wp, 0.0_wp,    th < 40.0_wp*pi/180.0_wp)
         end do
      end do
      ice = d_ice

      call nc_create(f, overwrite=.true.)
      call stamp(f, "sle")
      call w_one(f)
      call grid_dims_std(f, g)
      allocate(tyr(NSTEP))
      do is = 1, NSTEP;  tyr(is) = real(is-1, wp)*100.0_wp;  end do
      call nc_write_dim(f, "time", x=tyr, units="years", long_name="report time of each ve SLE step")
      call nc_write_attr(f, "sle_settings", "n_outer=3 n_inner=20 tol=1e-7 subgrid=.true. fixed_ocean=.false. warm_start=.false. (sle_solver defaults)")
      call nc_write_attr(f, "fields", "topo0: +500 m for colat<60 deg else -4000 m; d_ice = ice = 2000 m for colat<40 deg (test_sle_ve)")
      call nc_write(f, "topo0", topo0, dim1="lon", dim2="colat", units="m", long_name="reference topography (solid surface vs sea surface)")
      call nc_write(f, "d_ice", d_ice, dim1="lon", dim2="colat", units="m", long_name="grounded ice thickness change (load); ice = d_ice")

      ! elastic (stateless: one solve)
      call response_init_elastic(el, e, lmax=LMAX32)
      call sle_solve(sle, g, el, d_ice, ice, topo0, S, C, res)
      call nc_write(f, "rsl_el", S, dim1="lon", dim2="colat", units="m", long_name="elastic SLE: relative sea level change rsl = N - u + dphi (full field)")
      call nc_write(f, "C_el",   C, dim1="lon", dim2="colat", units="1", long_name="elastic SLE: converged ocean function")
      allocate(ug(g%nphi,g%nlat), Ng(g%nphi,g%nlat))
      call sle_result_grids(res, g, u=ug, N=Ng)
      call nc_write(f, "u_el",   ug, dim1="lon", dim2="colat", units="m", long_name="elastic SLE: converged solid uplift")
      call nc_write(f, "N_el",   Ng, dim1="lon", dim2="colat", units="m", long_name="elastic SLE: converged geoid rise")
      call ws(f, "esl_el", res%esl, "m", "elastic SLE: eustatic offset dphi (uniform sea-surface shift)")
      call ws(f, "mass_resid_el", res%mass_resid, "1", "elastic SLE: relative ocean-mass residual")
      call ws(f, "resid_el", res%resid, "m", "elastic SLE: last inner max|dS|")
      call ws(f, "ocean_frac_el", res%ocean_frac, "1", "elastic SLE: int C dOmega / 4pi")
      call wsi(f, "n_inner_el", res%n_inner_last, "elastic SLE: inner iterations in last outer pass")
      call wsi(f, "n_outer_el", res%n_outer_done, "elastic SLE: outer (coastline) passes")
      call response_destroy(el)

      ! viscoelastic, 10 steps
      call response_init_ve(ve, e, g, DT100)
      allocate(rsl(g%nphi,g%nlat,NSTEP), Cs(g%nphi,g%nlat,NSTEP), us(g%nphi,g%nlat,NSTEP), Ns(g%nphi,g%nlat,NSTEP))
      allocate(esl(NSTEP), mres(NSTEP), resid(NSTEP), ninner(NSTEP), nouter(NSTEP))
      do is = 1, NSTEP
         if (abs(ve%time - tyr(is)*sec_per_year) > 1.0e-6_wp*sec_per_year) error stop 'dump_sle: time mismatch'
         call sle_solve(sle, g, ve, d_ice, ice, topo0, S, C, res)
         call sle_result_grids(res, g, u=ug, N=Ng)
         rsl(:,:,is) = S;  Cs(:,:,is) = C;  us(:,:,is) = ug;  Ns(:,:,is) = Ng
         esl(is) = res%esl;  mres(is) = res%mass_resid;  resid(is) = res%resid
         ninner(is) = res%n_inner_last;  nouter(is) = res%n_outer_done
      end do
      call nc_write(f, "rsl_ve", rsl, dim1="lon", dim2="colat", dim3="time", units="m", long_name="ve SLE: rsl = N - u + dphi per step")
      call nc_write(f, "C_ve",   Cs,  dim1="lon", dim2="colat", dim3="time", units="1", long_name="ve SLE: converged ocean function per step")
      call nc_write(f, "u_ve",   us,  dim1="lon", dim2="colat", dim3="time", units="m", long_name="ve SLE: converged solid uplift per step")
      call nc_write(f, "N_ve",   Ns,  dim1="lon", dim2="colat", dim3="time", units="m", long_name="ve SLE: converged geoid rise per step")
      call nc_write(f, "esl_ve", esl, dim1="time", units="m", long_name="ve SLE: eustatic offset dphi per step")
      call nc_write(f, "mass_resid_ve", mres, dim1="time", units="1", long_name="ve SLE: relative ocean-mass residual per step")
      call nc_write(f, "resid_ve", resid, dim1="time", units="m", long_name="ve SLE: last inner max|dS| per step")
      call nc_write(f, "n_inner_ve", ninner, dim1="time", units="1", long_name="ve SLE: inner iterations (last outer pass)")
      call nc_write(f, "n_outer_ve", nouter, dim1="time", units="1", long_name="ve SLE: outer passes")
      call ws(f, "dt", DT100, "s", "ve time step")
      call nc_write_attr(f, "ve_scheme", "fe (SCHEME_FE); each sle_solve = begin_step + fixed point + advance_endpoint + finalize (time += dt)")
      call response_destroy(ve);  call sht_grid_destroy(g)
   end subroutine dump_sle

   ! =====================================================================
   ! Legendre helpers (disc / martinec)
   ! =====================================================================

   subroutine disc_coeffs(ca, sig0, sig)
      !! sig_n = sig0/2 [P_{n-1}(ca) - P_{n+1}(ca)], n = 1..size(sig)
      real(wp), intent(in)  :: ca, sig0
      real(wp), intent(out) :: sig(:)
      real(wp) :: Pnm1, Pn, Pnp1
      integer  :: n
      Pnm1 = 1.0_wp;  Pn = ca
      do n = 1, size(sig)
         Pnp1   = (real(2*n+1,wp)*ca*Pn - real(n,wp)*Pnm1)/real(n+1,wp)
         sig(n) = 0.5_wp*sig0*(Pnm1 - Pnp1)
         Pnm1 = Pn;  Pn = Pnp1
      end do
   end subroutine disc_coeffs

   real(wp) function legsum(coef, x) result(s)
      !! Sum_{n=1}^{N} coef(n) P_n(x)
      real(wp), intent(in) :: coef(:), x
      real(wp) :: pm2, pm1, pc
      integer  :: n
      pm2 = 1.0_wp;  pm1 = x;  s = 0.0_wp
      do n = 1, size(coef)
         if (n == 1) then
            pc = x
         else
            pc = (real(2*n-1,wp)*x*pm1 - real(n-1,wp)*pm2)/real(n,wp)
            pm2 = pm1;  pm1 = pc
         end if
         s = s + coef(n)*pc
      end do
   end function legsum

   subroutine read_mat(fname, x, okr)
      character(*), intent(in)  :: fname
      real(wp),     intent(out) :: x(:,:)
      logical,      intent(out) :: okr
      integer :: u, i, ios
      okr = .false.;  x = 0.0_wp
      open(newunit=u, file=fname, status='old', action='read', iostat=ios)
      if (ios /= 0) return
      do i = 1, size(x,1)
         read(u,*,iostat=ios) x(i,:)
         if (ios /= 0) then;  close(u);  return;  end if
      end do
      close(u);  okr = .true.
   end subroutine read_mat

   ! =====================================================================
   ! 6. disc.nc
   ! =====================================================================

   subroutine dump_disc(f)
      character(len=*), intent(in) :: f
      integer,  parameter :: NTH = 201, NT = 6, NMAX_EL = 256, NMAX_VE = 128, NTV = 5
      integer,  parameter :: tsteps(NTV) = [0, 50, 100, 250, 500]     ! @ dt = 20 yr
      real(wp), parameter :: ALPHA_DEG = 10.0_wp, H_ICE = 1000.0_wp
      real(wp), parameter :: tref_kyr(NT) = [0.0_wp, 1.0_wp, 2.0_wp, 5.0_wp, 10.0_wp, 100.0_wp]
      type(earth_model) :: em
      type(radial_mesh) :: mm
      type(ve_degree)   :: ve
      type(response)    :: el
      real(wp) :: uref(NTH,NT), nref(NTH,NT), sig(NMAX_EL), theta(NTH)
      real(wp) :: cu(NMAX_EL), cn(NMAX_EL), uel(NTH), nel(NTH)
      real(wp) :: Un(NMAX_VE,NTV), Nn(NMAX_VE,NTV), Vn(NMAX_VE,NTV), uve(NTH,NTV), nve(NTH,NTV), tve(NTV)
      real(wp) :: g, sig0, ca, dt, t1, ua, va, fa, x
      integer  :: n, i, it, istep
      logical  :: okr

      write(*,'(a)') ' [disc] Spada 2011 disc: elastic profile (NMAX=256) + VE profiles (NMAX=128, dt=20 yr)'
      call read_mat(VILMA_ROOT//'data/benchmarks/disc_spada2011/u_disc.txt', uref, okr)
      if (.not. okr) error stop 'dump_disc: cannot read u_disc.txt'
      call read_mat(VILMA_ROOT//'data/benchmarks/disc_spada2011/n_disc.txt', nref, okr)
      if (.not. okr) error stop 'dump_disc: cannot read n_disc.txt'

      em = build_M3L70V01();  g = earth_gravity_at(em, em%r_earth)
      sig0 = rho_ice*H_ICE;  ca = cos(ALPHA_DEG*pi/180.0_wp)
      call disc_coeffs(ca, sig0, sig)
      do i = 1, NTH;  theta(i) = real(i-1,wp)*0.1_wp;  end do

      ! elastic profile
      call response_init_elastic(el, em, lmax=NMAX_EL)
      do n = 1, NMAX_EL
         cu(n) = sig(n)*el%ugain(n);  cn(n) = sig(n)*el%ngain(n)
      end do
      do i = 1, NTH
         x = cos(theta(i)*pi/180.0_wp)
         uel(i) = legsum(cu, x);  nel(i) = legsum(cn, x)
      end do

      ! viscoelastic per-degree transient
      dt = 0.02_wp*kyr
      call radial_mesh_build(mm, em)
      Un = 0.0_wp;  Nn = 0.0_wp;  Vn = 0.0_wp
      do n = 2, NMAX_VE
         call ve_init(ve, em, mm, n, dt)
         it = 1
         do istep = 0, tsteps(NTV)
            call ve_step(ve, 1.0_wp, t1, ua, va, fa)
            if (it <= NTV) then
               if (istep == tsteps(it)) then
                  Un(n,it) = ua;  Vn(n,it) = va;  Nn(n,it) = -fa/g;  tve(it) = t1/sec_per_year;  it = it + 1
               end if
            end if
         end do
         call ve_destroy(ve)
      end do
      do it = 1, NTV
         do i = 1, NTH
            x = cos(theta(i)*pi/180.0_wp)
            uve(i,it) = legsum(sig(1:NMAX_VE)*Un(:,it), x)
            nve(i,it) = legsum(sig(1:NMAX_VE)*Nn(:,it), x)
         end do
      end do

      call nc_create(f, overwrite=.true.)
      call stamp(f, "disc")
      call w_one(f)
      call nc_write_dim(f, "theta", x=theta, units="degrees", long_name="colatitude from disc centre")
      call nc_write_dim(f, "tref",  x=tref_kyr, units="kyr", long_name="benchmark times")
      call nc_write_dim(f, "tve",   x=tve, units="years", long_name="VE report times (0,1,2,5,10 kyr)")
      call nc_write_dim(f, "degree_el", x=1, dx=1, nx=NMAX_EL, units="1")
      call nc_write_dim(f, "degree_ve", x=1, dx=1, nx=NMAX_VE, units="1")
      call nc_write_dim(f, "degree0", x=0, dx=1, nx=NMAX_EL+1, units="1")
      call nc_write_attr(f, "load", "uniform disc radius 10 deg, 1000 m ice (rho_ice=931), sigma_n = sigma0/2 [P_{n-1}-P_{n+1}](cos alpha)")
      call nc_write_attr(f, "reference", "data/benchmarks/disc_spada2011 (u_disc.txt, n_disc.txt; ascending theta)")
      call ws(f, "sigma0", sig0, "kg m-2", "disc surface mass density rho_ice*H")
      call ws(f, "g_surf", g, "m s-2", "surface gravity")
      call ws(f, "dt_ve", dt, "s", "VE time step (20 yr)")
      call nc_write(f, "u_ref", uref, dim1="theta", dim2="tref", units="m", long_name="Spada 2011 reference vertical displacement")
      call nc_write(f, "N_ref", nref, dim1="theta", dim2="tref", units="m", long_name="Spada 2011 reference geoid")
      call nc_write(f, "sig_n", sig, dim1="degree_el", units="kg m-2", long_name="disc load Legendre coefficients sigma_n")
      call nc_write(f, "ugain", el%ugain, dim1="degree0", units="m / (kg m-2)", long_name="elastic gain U(a)/sigma_l (response_init_elastic lmax=256)")
      call nc_write(f, "ngain", el%ngain, dim1="degree0", units="m / (kg m-2)", long_name="elastic gain N(a)/sigma_l (ngain(1)=0)")
      call nc_write(f, "u_el", uel, dim1="theta", units="m", long_name="model elastic u(theta) = sum sig_n ugain_n P_n, n=1..256")
      call nc_write(f, "N_el", nel, dim1="theta", units="m", long_name="model elastic N(theta) = sum sig_n ngain_n P_n")
      call nc_write(f, "U_n_ve", Un, dim1="degree_ve", dim2="tve", units="m", long_name="per-degree U(a,t) per unit sigma (ve_degree, fe, dt=20 yr; n=1 skipped)")
      call nc_write(f, "V_n_ve", Vn, dim1="degree_ve", dim2="tve", units="m", long_name="per-degree V(a,t) per unit sigma")
      call nc_write(f, "N_n_ve", Nn, dim1="degree_ve", dim2="tve", units="m", long_name="per-degree N(a,t) = -F/g per unit sigma")
      call nc_write(f, "u_ve", uve, dim1="theta", dim2="tve", units="m", long_name="model VE u(theta,t) = sum_{n=2}^{128} sig_n U_n(t) P_n")
      call nc_write(f, "N_ve", nve, dim1="theta", dim2="tve", units="m", long_name="model VE N(theta,t) = sum_{n=2}^{128} sig_n N_n(t) P_n")
      call response_destroy(el)
   end subroutine dump_disc

   ! =====================================================================
   ! 7. rotation.nc
   ! =====================================================================

   subroutine dump_rotation(f)
      character(len=*), intent(in) :: f
      real(wp), parameter :: deg = acos(-1.0_wp)/180.0_wp
      real(wp), parameter :: yr_sp = 3.15576e7_wp                ! Spada Table 2 sec/yr
      real(wp), parameter :: alpha = 10.0_wp*deg, thetac = 25.0_wp*deg, lambdac = 75.0_wp*deg
      real(wp), parameter :: t_ref(6) = [0.0_wp,1.0_wp,2.0_wp,5.0_wp,10.0_wp,20.0_wp]
      real(wp), parameter :: m_cap(6) = [0.0132_wp,0.0161_wp,0.0180_wp,0.0211_wp,0.0240_wp,0.0292_wp]
      real(wp), parameter :: m_dsc(6) = [0.0131_wp,0.0160_wp,0.0179_wp,0.0210_wp,0.0239_wp,0.0290_wp]
      real(wp), parameter :: G_cap(2) = [-0.541e-4_wp, -0.202e-3_wp]
      real(wp), parameter :: G_dsc(2) = [-0.539e-4_wp, -0.201e-3_wp]
      integer,  parameter :: NSTEP = 801
      type(earth_model)    :: earth
      type(sht_grid)       :: g
      type(rotation_state) :: rot
      real(wp), allocatable :: load(:,:), tyr(:), m1(:), m2(:), mabs(:)
      real(wp) :: dt
      integer  :: ic, is
      character(len=4) :: nm(2) = [character(len=4) :: "cap", "disc"]
      real(wp) :: hh(2) = [1.5e3_wp, 1.0e3_wp]

      write(*,'(a)') ' [rotation] Spada Test 3/2 polar motion, lmax=128 grid, dt=25 yr, 20 kyr'
      earth = build_M3L70V01()
      call sht_grid_init(g, 128, nlat=256, nphi=512)
      allocate(load(g%nphi, g%nlat), tyr(NSTEP), m1(NSTEP), m2(NSTEP), mabs(NSTEP))
      dt = 25.0_wp*yr_sp

      call nc_create(f, overwrite=.true.)
      call stamp(f, "rotation")
      call w_one(f)
      call w_grid(f, g, "lon", "colat")
      call nc_write_dim(f, "tref", x=t_ref, units="kyr", long_name="Table 14 times")
      call nc_write_dim(f, "reim", x=1, dx=1, nx=2, units="1")
      call nc_write_attr(f, "grid", "lmax=128 nlat=256 nphi=512 (test_rotation)")
      call nc_write_attr(f, "sec_per_year_used", yr_sp)
      call nc_write_attr(f, "loads", "cap: sigma=rho_i h sqrt((cos g - cos a)/(1-cos a)), h=1500 m; disc: sigma=rho_i h, h=1000 m; alpha=10 deg, centroid colat 25 deg lon 75 deg; rho_i=931")
      call ws(f, "dt", dt, "s", "rotation channel step (25 yr Spada)")
      call nc_write(f, "m_cap_ref", m_cap, dim1="tref", units="degrees", long_name="Table 14 |m| cap, Cw=0 (Gs)")
      call nc_write(f, "m_disc_ref", m_dsc, dim1="tref", units="degrees", long_name="Table 14 |m| disc, Cw=0 (Gs)")
      call nc_write(f, "G_cap_ref", G_cap, dim1="reim", units="1", long_name="published geometrical factor G (cap), re/im")
      call nc_write(f, "G_disc_ref", G_dsc, dim1="reim", units="1", long_name="published geometrical factor G (disc), re/im")

      do ic = 1, 2
         call build_rot_load(g, nm(ic), hh(ic), load)
         call nc_write(f, "load_"//trim(nm(ic)), load, dim1="lon", dim2="colat", units="kg m-2", &
              long_name="surface mass load, "//trim(nm(ic)))
         call rotation_init(rot, earth, g, dt)
         rot%enabled = .true.
         if (ic == 1) then
            call ws(f, "kTe", rot%kTe, "1", "elastic degree-2 tidal Love number k^T_e (rotation_init)")
            call ws(f, "hTe", rot%hTe, "1", "elastic degree-2 tidal Love number h^T_e")
            call ws(f, "k_s", rot%k_s, "1", "secular Love number used (= k_s_fluid)")
            call ws(f, "k_s_fluid", rot%k_s_fluid, "1", "model relaxed tidal k^T_f")
            call ws(f, "k_s_flat", rot%k_s_flat, "1", "observed-flattening k_s = 3G(C-A)/(a^5 Omega^2)")
            call ws(f, "CminusA", rot%CminusA, "kg m2", "C - A")
            call ws(f, "Omega", rot%Omega, "rad s-1", "rotation rate")
            call ws(f, "dt_fe_max", rot%dt_fe_max, "s", "forward-Euler stability ceiling of the channels")
         end if
         call nc_write(f, "Irig_"//trim(nm(ic)), [real(inertia21(g, load, earth%r_earth), wp), &
              aimag(inertia21(g, load, earth%r_earth))], dim1="reim", units="kg m2", &
              long_name="rigid (2,1) inertia I13 + i I23 = -a^4 int sigma sin cos e^{i phi} dOmega")
         do is = 1, NSTEP
            tyr(is) = rot%time/yr_sp
            call rotation_update(rot, g, load, dt)
            m1(is) = real(rot%m, wp);  m2(is) = aimag(rot%m);  mabs(is) = abs(rot%m)
         end do
         if (ic == 1) call nc_write_dim(f, "time", x=tyr, units="years", long_name="report time (m at entry time, before advancing)")
         call nc_write(f, "m1_"//trim(nm(ic)), m1, dim1="time", units="rad", long_name="polar motion m1(t), "//trim(nm(ic)))
         call nc_write(f, "m2_"//trim(nm(ic)), m2, dim1="time", units="rad", long_name="polar motion m2(t)")
         call nc_write(f, "mabs_"//trim(nm(ic)), mabs, dim1="time", units="rad", long_name="|m|(t)")
         call rotation_destroy(rot)
      end do
      call sht_grid_destroy(g)
   end subroutine dump_rotation

   subroutine build_rot_load(g, name, h, load)
      type(sht_grid),   intent(in)  :: g
      character(len=*), intent(in)  :: name
      real(wp),         intent(in)  :: h
      real(wp),         intent(out) :: load(:,:)
      real(wp), parameter :: deg = acos(-1.0_wp)/180.0_wp
      real(wp), parameter :: alpha = 10.0_wp*deg, thetac = 25.0_wp*deg, lambdac = 75.0_wp*deg
      real(wp), parameter :: rho_i = 931.0_wp
      real(wp) :: ca, cg, frac
      integer  :: il, ip
      ca = cos(alpha);  load = 0.0_wp
      do il = 1, g%nlat
         do ip = 1, g%nphi
            cg = cos(thetac)*cos(g%colat(il)) + sin(thetac)*sin(g%colat(il))*cos(g%lon(ip) - lambdac)
            if (cg >= ca) then
               if (name == 'cap') then
                  frac = (cg - ca)/(1.0_wp - ca)
                  load(ip,il) = rho_i*h*sqrt(max(frac, 0.0_wp))
               else
                  load(ip,il) = rho_i*h
               end if
            end if
         end do
      end do
   end subroutine build_rot_load

   complex(wp) function inertia21(g, load, a) result(I21)
      type(sht_grid), intent(in) :: g
      real(wp),       intent(in) :: load(:,:), a
      real(wp), allocatable :: w13(:,:), w23(:,:)
      integer  :: il, ip
      allocate(w13(g%nphi,g%nlat), w23(g%nphi,g%nlat))
      do il = 1, g%nlat
         do ip = 1, g%nphi
            w13(ip,il) = load(ip,il)*sin(g%colat(il))*cos(g%colat(il))*cos(g%lon(ip))
            w23(ip,il) = load(ip,il)*sin(g%colat(il))*cos(g%colat(il))*sin(g%lon(ip))
         end do
      end do
      I21 = cmplx(-a**4*sht_grid_surface_integral(g, w13), -a**4*sht_grid_surface_integral(g, w23), wp)
   end function inertia21

   ! =====================================================================
   ! 8. coupling.nc (+ coupling_restart.nc)
   ! =====================================================================

   subroutine dump_coupling(f, outdir)
      character(len=*), intent(in) :: f, outdir
      integer,  parameter :: NINT = 5
      real(wp), parameter :: DT_YR = 100.0_wp
      type(sht_grid), target :: g
      type(vilma_param_class)   :: p
      type(solid_earth)      :: se
      real(wp), allocatable  :: z_bed_eq(:,:), h_ice_eq(:,:), h_ice(:,:)
      real(wp), allocatable  :: rsl(:,:,:), zb(:,:,:), Cs(:,:,:), bsl(:), tyr(:), wm(:)
      integer  :: i, j, is, isch
      real(wp) :: thd
      character(len=8) :: sn(2) = [character(len=8) :: "fe", "trap"]

      write(*,'(a)') ' [coupling] solid_earth init/update lmax=32, 5 x 100 yr, fe + trap (+ restart)'
      call sht_grid_init(g, LMAX32, nlat=2*LMAX32, nphi=4*LMAX32)
      allocate(z_bed_eq(g%nphi,g%nlat), h_ice_eq(g%nphi,g%nlat), h_ice(g%nphi,g%nlat))
      do j = 1, g%nlat
         thd = g%colat(j)*rad2deg
         do i = 1, g%nphi
            z_bed_eq(i,j) = merge(500.0_wp, -4000.0_wp, thd < 50.0_wp)
            h_ice(i,j)    = merge(2000.0_wp, 0.0_wp,     thd < 30.0_wp)
         end do
      end do
      h_ice_eq = 0.0_wp

      call nc_create(f, overwrite=.true.)
      call stamp(f, "coupling")
      call w_one(f)
      call w_grid(f, g, "lon", "colat")
      allocate(tyr(NINT), bsl(NINT), wm(NINT))
      do is = 1, NINT;  tyr(is) = real(is, wp)*DT_YR;  end do
      call nc_write_dim(f, "time", x=tyr, units="years", long_name="model time after each solid_earth_update interval")
      call nc_write_attr(f, "params", "vilma_param_class defaults except lmax=32 nlat=64 nphi=128 rotation=.false. scheme=fe|trap; earth_response=ve; sle defaults; cfl=1; rtol=1e-4")
      call nc_write_attr(f, "fields", "z_bed_eq: +500 m colat<50 deg else -4000 m; h_ice_eq=0; h_ice: 2000 m colat<30 deg held (test_coupling)")
      call nc_write(f, "z_bed_eq", z_bed_eq, dim1="lon", dim2="colat", units="m", long_name="reference (relaxed) bedrock")
      call nc_write(f, "h_ice_eq", h_ice_eq, dim1="lon", dim2="colat", units="m", long_name="reference grounded ice")
      call nc_write(f, "h_ice",    h_ice,    dim1="lon", dim2="colat", units="m", long_name="forcing grounded ice (held over all intervals)")
      allocate(rsl(g%nphi,g%nlat,NINT), zb(g%nphi,g%nlat,NINT), Cs(g%nphi,g%nlat,NINT))

      do isch = 1, 2
         p%lmax = LMAX32;  p%nlat = 2*LMAX32;  p%nphi = 4*LMAX32
         p%rotation = .false.
         p%scheme = sn(isch)
         se%par = p
         call solid_earth_init(se, z_bed_eq, h_ice_eq)
         do is = 1, NINT
            call solid_earth_update(se, h_ice, DT_YR)
            rsl(:,:,is) = se%rsl;  zb(:,:,is) = se%z_bed;  Cs(:,:,is) = se%gg%C
            bsl(is) = se%bsl;  wm(is) = se%worst_mass_resid
            if (abs(se%time - tyr(is)) > 1.0e-9_wp) error stop 'dump_coupling: time mismatch'
         end do
         call nc_write(f, "rsl_"//trim(sn(isch)), rsl, dim1="lon", dim2="colat", dim3="time", units="m", &
              long_name="relative sea level change after each interval, scheme "//trim(sn(isch)))
         call nc_write(f, "z_bed_"//trim(sn(isch)), zb, dim1="lon", dim2="colat", dim3="time", units="m", &
              long_name="bedrock z_bed = z_bed_eq - rsl")
         call nc_write(f, "C_"//trim(sn(isch)), Cs, dim1="lon", dim2="colat", dim3="time", units="1", &
              long_name="ocean function")
         call nc_write(f, "bsl_"//trim(sn(isch)), bsl, dim1="time", units="m", long_name="barystatic sea level vs reference")
         call nc_write(f, "worst_mass_resid_"//trim(sn(isch)), wm, dim1="time", units="1", &
              long_name="worst SLE mass residual over the interval")
         if (isch == 1) then
            call vilma_restart_write(se, se%time, filename="coupling_restart.nc", folder=outdir, init=.true.)
         end if
         call solid_earth_finalize(se)
      end do
      call sht_grid_destroy(g)
   end subroutine dump_coupling

   ! =====================================================================
   ! 9. martinec_A.nc
   ! =====================================================================

   subroutine dump_martinec(f)
      character(len=*), intent(in) :: f
      integer,  parameter :: NROW = 721, NMAX = 128, NQ = 8000, NSTEP = 500
      real(wp), parameter :: ALPHA_DEG = 10.0_wp, H0 = 1500.0_wp
      type(earth_model) :: em
      type(radial_mesh) :: mm
      type(ve_degree)   :: ve
      real(wp) :: Un(NMAX), Vn(NMAX), Nn(NMAX), sig(NMAX)
      real(wp) :: colat(NROW), uref(NROW), href(NROW), gref(NROW), us(NROW), vs(NROW), ns(NROW)
      real(wp) :: alpha, ca, g, dt, t1, ua, va, fa
      real(wp) :: th, x, st, pn, pnm1, pnp1, dpn, d, dd, w, s, cd
      integer  :: n, istep, i, iq, u, ios

      write(*,'(a)') ' [martinec] Martinec 2018 case A: cap load, 10 kyr (NMAX=128, dt=20 yr)'
      open(newunit=u, file=VILMA_ROOT//'data/benchmarks/sle_martinec2018/A_fig10_SBK.dat', status='old', action='read', iostat=ios)
      if (ios /= 0) error stop 'dump_martinec: cannot read A_fig10_SBK.dat'
      read(u,*)
      do i = 1, NROW
         read(u,*,iostat=ios) colat(i), uref(i), href(i), gref(i)
         if (ios /= 0) error stop 'dump_martinec: short read'
      end do
      close(u)

      em = build_M3L70V01();  call radial_mesh_build(mm, em)
      g  = earth_gravity_at(em, em%r_earth)
      alpha = ALPHA_DEG*pi/180.0_wp;  ca = cos(alpha)
      ! cap coefficients (midpoint quadrature as in the test)
      sig = 0.0_wp;  dd = alpha/real(NQ,wp)
      do iq = 1, NQ
         d  = (real(iq,wp)-0.5_wp)*dd;  cd = cos(d)
         s  = rho_ice*H0*sqrt(max((cd-ca)/(1.0_wp-ca), 0.0_wp))
         w  = s*sin(d)*dd
         pnm1 = 1.0_wp;  pn = cd
         do n = 1, NMAX
            sig(n) = sig(n) + w*pn
            pnp1 = (real(2*n+1,wp)*cd*pn - real(n,wp)*pnm1)/real(n+1,wp)
            pnm1 = pn;  pn = pnp1
         end do
      end do
      do n = 1, NMAX;  sig(n) = sig(n)*real(2*n+1,wp)*0.5_wp;  end do

      dt = 0.02_wp*kyr
      Un = 0.0_wp;  Vn = 0.0_wp;  Nn = 0.0_wp
      do n = 2, NMAX
         call ve_init(ve, em, mm, n, dt)
         do istep = 1, NSTEP
            call ve_step(ve, 1.0_wp, t1, ua, va, fa)
         end do
         Un(n) = ua;  Vn(n) = va;  Nn(n) = -fa/g
         call ve_destroy(ve)
      end do
      do i = 1, NROW
         th = colat(i)*pi/180.0_wp;  x = cos(th);  st = sin(th)
         pnm1 = 1.0_wp;  pn = x
         us(i) = 0.0_wp;  vs(i) = 0.0_wp;  ns(i) = 0.0_wp
         do n = 1, NMAX
            if (st > 1.0e-12_wp) then;  dpn = real(n,wp)*(x*pn - pnm1)/st;  else;  dpn = 0.0_wp;  end if
            if (n >= 2) then
               us(i) = us(i) + sig(n)*Un(n)*pn
               vs(i) = vs(i) + sig(n)*Vn(n)*dpn
               ns(i) = ns(i) + sig(n)*Nn(n)*pn
            end if
            pnp1 = (real(2*n+1,wp)*x*pn - real(n,wp)*pnm1)/real(n+1,wp)
            pnm1 = pn;  pn = pnp1
         end do
      end do

      call nc_create(f, overwrite=.true.)
      call stamp(f, "martinec_A")
      call w_one(f)
      call nc_write_dim(f, "colat", x=colat, units="degrees", long_name="reference colatitude (180 -> 0, 0.25 deg)")
      call nc_write_dim(f, "degree", x=1, dx=1, nx=NMAX, units="1")
      call nc_write_attr(f, "load", "spherical cap alpha=10 deg, h0=1500 m, rho=931, sigma=rho h0 sqrt((cos d - cos a)/(1-cos a)); Heaviside held 10 kyr; degrees 0,1 dropped")
      call nc_write_attr(f, "reference", "data/benchmarks/sle_martinec2018/A_fig10_SBK.dat (giapy j256)")
      call ws(f, "dt", dt, "s", "VE time step (20 yr)");  call wsi(f, "nstep", NSTEP, "steps (t = 10 kyr after last advance)")
      call ws(f, "g_surf", g, "m s-2", "surface gravity")
      call nc_write(f, "sig_n", sig, dim1="degree", units="kg m-2", long_name="cap load Legendre coefficients")
      call nc_write(f, "U_n", Un, dim1="degree", units="m", long_name="per-degree U(a) per unit sigma reported at step 500 (t=9.98 kyr entry; see README)")
      call nc_write(f, "V_n", Vn, dim1="degree", units="m", long_name="per-degree V(a) per unit sigma")
      call nc_write(f, "N_n", Nn, dim1="degree", units="m", long_name="per-degree N = -F/g per unit sigma")
      call nc_write(f, "u_model", us, dim1="colat", units="m", long_name="model uplift sum sig_n U_n P_n")
      call nc_write(f, "h_model", vs, dim1="colat", units="m", long_name="model horizontal (theta) sum sig_n V_n dP_n/dtheta")
      call nc_write(f, "N_model", ns, dim1="colat", units="m", long_name="model geoid sum sig_n N_n P_n")
      call nc_write(f, "u_ref", uref, dim1="colat", units="m", long_name="reference uplift")
      call nc_write(f, "h_ref", href, dim1="colat", units="m", long_name="reference horizontal displacement")
      call nc_write(f, "N_ref", gref, dim1="colat", units="m", long_name="reference geoid")
   end subroutine dump_martinec

   ! =====================================================================
   ! 10. visc3d.nc
   ! =====================================================================

   subroutine dump_visc3d(f)
      character(len=*), intent(in) :: f
      integer, parameter :: NSTEP = 10
      integer(8) :: tic3, toc3, trate3
      type(sht_grid)    :: g
      type(earth_model) :: e
      type(radial_mesh) :: mesh
      type(response)    :: ve3d, veu, ve1d
      real(wp),    allocatable :: pert(:,:,:), tyr(:), pert3(:,:,:), mk3(:,:,:)
      complex(wp), allocatable :: slm(:), ulm(:), nlm(:)
      complex(wp), allocatable :: u3(:,:), n3(:,:), uu(:,:), nu(:,:), u1(:,:), n1(:,:)
      integer :: i, j, ie, is, k

      write(*,'(a)') ' [visc3d] lateral viscosity lmax=32 (96x96, mmax=32), held (2,1), 10 x 100 yr'
      call sht_grid_init(g, LMAX32, nlat=3*LMAX32, nphi=3*LMAX32, mmax=LMAX32)
      e = build_M3L70V01();  call radial_mesh_build(mesh, e)
      allocate(slm(g%nlm), ulm(g%nlm), nlm(g%nlm), tyr(NSTEP))
      allocate(u3(g%nlm,NSTEP), n3(g%nlm,NSTEP), uu(g%nlm,NSTEP), nu(g%nlm,NSTEP), u1(g%nlm,NSTEP), n1(g%nlm,NSTEP))
      slm = (0.0_wp, 0.0_wp)
      slm(sht_grid_lmidx(g, 2, 1)) = cmplx(1000.0_wp, 500.0_wp, wp)

      ! (A) analytic lateral anomaly on the upper-mantle layer (elem_layer == 2)
      call response_init_ve(ve3d, e, g, DT100)
      allocate(pert(g%nphi, g%nlat, ve3d%ne));  pert = 0.0_wp
      do ie = 1, ve3d%ne
         if (mesh%elem_layer(ie) /= 2) cycle
         do j = 1, g%nlat
            do i = 1, g%nphi
               pert(i,j,ie) = 0.5_wp*sin(g%colat(j))*cos(g%lon(i))
            end do
         end do
      end do
      call response_enable_lateral_visc(ve3d, g, pert)
      call system_clock(tic3, trate3)
      call drive3d(ve3d, g, slm, NSTEP, u3, n3, tyr)
      call system_clock(toc3)
      write(*,'(a,f9.4,a)') '  visc3d timing run A 3-D 35 elements: ', real(toc3-tic3,wp)/real(trate3,wp)/real(NSTEP,wp)*1.0e3_wp, ' ms/step'
      ! (B) uniform (zero) perturbation forced through the 3-D kernel
      call response_init_ve(veu, e, g, DT100)
      veu%visc3d_tol = -1.0_wp
      pert = 0.0_wp
      call response_enable_lateral_visc(veu, g, pert)
      call system_clock(tic3, trate3)
      call drive3d(veu, g, slm, NSTEP, uu, nu, tyr)
      call system_clock(toc3)
      write(*,'(a,f9.4,a)') '  visc3d timing run B all Maxwell elements 3-D: ', real(toc3-tic3,wp)/real(trate3,wp)/real(NSTEP,wp)*1.0e3_wp, ' ms/step'
      ! (C) plain 1-D
      call response_init_ve(ve1d, e, g, DT100)
      call system_clock(tic3, trate3)
      call drive3d(ve1d, g, slm, NSTEP, u1, n1, tyr)
      call system_clock(toc3)
      write(*,'(a,f9.4,a)') '  visc3d timing run C 1-D: ', real(toc3-tic3,wp)/real(trate3,wp)/real(NSTEP,wp)*1.0e3_wp, ' ms/step'

      call nc_create(f, overwrite=.true.)
      call stamp(f, "visc3d")
      call w_one(f)
      call grid_dims_std(f, g)
      call nc_write_dim(f, "time", x=tyr, units="years", long_name="report time (state before commit)")
      call nc_write_dim(f, "elem", x=1, dx=1, nx=ve3d%ne, units="1")
      call nc_write_dim(f, "e3d",  x=1, dx=1, nx=ve3d%ne3d, units="1")
      call nc_write_dim(f, "node", x=1, dx=1, nx=ve3d%nr, units="1")
      call nc_write_attr(f, "pert_formula", "log10(eta) perturbation = 0.5*sin(colat)*cos(lon) on elements with elem_layer==2 (upper mantle, 5951-6301 km), 0 elsewhere; eta_eff = eta*10^pert")
      call nc_write_attr(f, "run_A", "ve3d: enable_lateral_visc(pert), visc3d_tol=1e-3 (default), scheme fe")
      call nc_write_attr(f, "run_B", "veu: pert=0 with visc3d_tol=-1 (every Maxwell element forced through the pseudo-spectral tensor advance)")
      call nc_write_attr(f, "run_C", "ve1d: plain 1-D ve_response")
      call ws(f, "dt", DT100, "s", "time step")
      call nc_write(f, "elem_layer", mesh%elem_layer, dim1="elem", units="1", long_name="earth layer per element")
      call nc_write(f, "r_node", ve3d%r, dim1="node", units="m", long_name="node radii")
      call nc_write(f, "e3d", ve3d%e3d, dim1="e3d", units="1", long_name="element indices flagged genuinely 3-D (run A)")
      call nc_write(f, "MkPerDt_1d", ve1d%MkPerDt, dim1="elem", units="s-1", long_name="1-D rate mu/eta per element")
      call nc_write(f, "MkPerDt_A", ve3d%MkPerDt, dim1="elem", units="s-1", long_name="run A scalar rate per element (lateral mean for 1-D-effective elements)")
      allocate(pert3(g%nphi,g%nlat,ve3d%ne3d), mk3(g%nphi,g%nlat,ve3d%ne3d))
      do k = 1, ve3d%ne3d
         ie = ve3d%e3d(k)
         mk3(:,:,k) = ve3d%Mk3(:,:,ie)
         do j = 1, g%nlat
            do i = 1, g%nphi
               pert3(i,j,k) = 0.5_wp*sin(g%colat(j))*cos(g%lon(i))
            end do
         end do
      end do
      call nc_write(f, "pert_3d", pert3, dim1="lon", dim2="colat", dim3="e3d", units="dex", long_name="log10 eta perturbation on the 3-D elements (run A)")
      call nc_write(f, "Mk3", mk3, dim1="lon", dim2="colat", dim3="e3d", units="1", long_name="run A lateral Maxwell factor M = mu dt / eta_eff on the 3-D elements")
      call wc1(f, "sigma", slm, "lm", "kg m-2", "held load (2,1) = 1000 + 500 i")
      call wc2(f, "u_lm_3d", u3, "lm", "time", "m", "run A u coefficients")
      call wc2(f, "N_lm_3d", n3, "lm", "time", "m", "run A N coefficients")
      call wc2(f, "u_lm_uniform3d", uu, "lm", "time", "m", "run B u coefficients (uniform field through 3-D kernel)")
      call wc2(f, "N_lm_uniform3d", nu, "lm", "time", "m", "run B N coefficients")
      call wc2(f, "u_lm_1d", u1, "lm", "time", "m", "run C u coefficients (1-D)")
      call wc2(f, "N_lm_1d", n1, "lm", "time", "m", "run C N coefficients")
      call response_destroy(ve3d);  call response_destroy(veu);  call response_destroy(ve1d)
      call sht_grid_destroy(g)

   end subroutine dump_visc3d

   subroutine drive3d(ve, g, slm, nstep, us, ns, tyr)
      !! begin/apply/commit a held load nstep times, recording u, N per step.
      type(response), intent(inout) :: ve
      type(sht_grid), intent(in)    :: g
      complex(wp),    intent(in)    :: slm(:)
      integer,        intent(in)    :: nstep
      complex(wp),    intent(out)   :: us(:,:), ns(:,:)
      real(wp),       intent(out)   :: tyr(:)
      complex(wp), allocatable :: ulm(:), nlm(:)
      integer :: is
      allocate(ulm(g%nlm), nlm(g%nlm))
      do is = 1, nstep
         tyr(is) = ve%time/sec_per_year
         call response_begin_step(ve, g)
         call response_apply(ve, g, slm, ulm, nlm)
         us(:,is) = ulm;  ns(:,is) = nlm
         call response_commit_step(ve, g, slm)
      end do
   end subroutine drive3d


   ! =====================================================================
   ! 11. martinec_sle.nc -- Martinec 2018 SLE cases C2/D3/E2/F1
   !     (mirrors tests/test_benchmark_sle.f90 exactly, at lmax 64 and 128)
   ! =====================================================================

   subroutine msle_history(sht, em, resp, sle, h0, ngrow, nsteps, heav, topo0, rsl, C, ice_now, res)
      !! One time history on the given reference topography `topo0`, from a relaxed, ice-free
      !! state (the response is destroyed and re-initialised, zeroing the Maxwell memory).
      !! T0 (heav): full ice every step. T1: linear growth over `ngrow` steps, then held.
      type(sht_grid),    intent(in)    :: sht
      type(earth_model), intent(in)    :: em
      type(response),    intent(inout) :: resp
      type(sle_solver),  intent(inout) :: sle
      real(wp),          intent(in)    :: h0
      integer,           intent(in)    :: ngrow, nsteps
      logical,           intent(in)    :: heav
      real(wp),          intent(in)    :: topo0(:,:)
      real(wp),          intent(inout) :: rsl(:,:), C(:,:), ice_now(:,:)
      type(sle_result),  intent(inout) :: res
      integer  :: istep
      real(wp) :: frac, alpha, h
      call response_destroy(resp)
      call response_init_ve(resp, em, sht, MS_DT)
      rsl = 0.0_wp
      do istep = 1, nsteps
         if (heav) then
            frac = 1.0_wp
         else
            frac = min(1.0_wp, real(istep,wp)/real(ngrow,wp))
         end if
         alpha = frac*MS_ALPHA;  h = frac*h0
         call spherical_cap(sht, MS_ICE_COLAT, MS_ICE_LON, alpha, h, ice_now)
         call sle_solve(sle, sht, resp, ice_now, ice_now, topo0, rsl, C, res)
      end do
   end subroutine msle_history


   subroutine msle_run(sht, em, icase, h0, bcol, blon, bmx, bz, ngrow, nsteps, heav, fixo, spn, &
                       topo0, rsl, C, u_lm, N_lm, rsl_lm, v_lm, esl, ocfrac, mresid, nspin)
      !! Full case: build the exponential basin, (F1) iterate the paleotopography, run the
      !! history and leave the converged spectral fields of the final step.
      type(sht_grid),    intent(in)  :: sht
      type(earth_model), intent(in)  :: em
      integer,           intent(in)  :: icase, ngrow, nsteps
      real(wp),          intent(in)  :: h0, bcol, blon, bmx, bz
      logical,           intent(in)  :: heav, fixo, spn
      real(wp),          intent(out) :: topo0(:,:), rsl(:,:), C(:,:)
      complex(wp),       intent(out) :: u_lm(:), N_lm(:), rsl_lm(:), v_lm(:)
      real(wp),          intent(out) :: esl, ocfrac, mresid
      integer,           intent(out) :: nspin
      type(response)   :: resp
      type(sle_solver) :: sle
      type(sle_result) :: res
      real(wp), allocatable :: basin(:,:), ice_now(:,:), tmp(:,:), dtopo(:,:), C0(:,:)
      complex(wp), allocatable :: load_lm(:)
      real(wp) :: spinerr
      integer  :: spin

      allocate(basin(sht%nphi,sht%nlat), ice_now(sht%nphi,sht%nlat), tmp(sht%nphi,sht%nlat), &
               dtopo(sht%nphi,sht%nlat), C0(sht%nphi,sht%nlat), load_lm(sht%nlm))
      call response_init_ve(resp, em, sht, MS_DT)
      sle%n_inner     = 20                 ! n_outer left at the default 3
      sle%fixed_ocean = fixo
      sle%subgrid     = .true.             ! the library default (giapy is a subgrid coast)
      sle%warm_start  = .true.             ! production path; rsl is zeroed per history

      call exp_basin(sht, bcol, blon, bmx, bz, MS_SIGB, basin)
      topo0 = basin

      nspin = 1
      if (spn) then
         do spin = 1, MS_MAXSPIN
            call msle_history(sht, em, resp, sle, h0, ngrow, nsteps, heav, topo0, rsl, C, ice_now, res)
            dtopo   = (topo0 - rsl) - basin
            spinerr = sht_grid_surface_integral(sht, abs(dtopo))/(16.0_wp*atan(1.0_wp))
            nspin   = spin
            write(*,'(a,i0,a,i2,a,f10.4,a)') '   [', icase, '] spinup ', spin, &
                 ': mean|present - B2| = ', spinerr, ' m'
            if (spinerr < MS_SPINTOL) exit
            topo0 = topo0 - 0.5_wp*dtopo
         end do
      else
         call msle_history(sht, em, resp, sle, h0, ngrow, nsteps, heav, topo0, rsl, C, ice_now, res)
      end if

      esl    = res%esl
      ocfrac = res%ocean_frac
      mresid = res%mass_resid
      ! converged surface load INCLUDING the subgrid sloping-coast correction (the load the
      ! solver actually applied): w_corr = -rho_water*topo0*(C - C0), with the reference ocean
      ! function C0 of the ice-free reference state (d_ice = ice, so ice_ref = 0).
      where (topo0 < 0.0_wp)
         C0 = 1.0_wp
      elsewhere
         C0 = 0.0_wp
      end where
      tmp = rho_ice*ice_now*(1.0_wp - C) + rho_water*(C*rsl) - rho_water*topo0*(C - C0)
      call sht_grid_analysis(sht, tmp, load_lm)
      call response_horizontal(resp, sht, load_lm, v_lm)
      u_lm = res%u_lm;  N_lm = res%N_lm
      tmp = rsl;    call sht_grid_analysis(sht, tmp, rsl_lm)
      call response_destroy(resp)
      deallocate(basin, ice_now, tmp, dtopo, C0, load_lm)
   end subroutine msle_run


   subroutine msle_read_fig(cname, fname, nrow, c1, ref)
      !! Read a 7-column *_SBK.dat profile: col1; col2 u; col3 v_th; col4 v_ph; col5 -N g;
      !! col6 sea surface; col7 SLE. Stored as ref(:,1:6) = u, v_th, v_ph, N, ss, sle.
      character(*), intent(in)  :: cname, fname
      integer,      intent(in)  :: nrow
      real(wp),     intent(out) :: c1(:), ref(:,:)
      character(256) :: line
      real(wp) :: v(7)
      integer  :: u, ios, n
      open(newunit=u, file=VILMA_ROOT//'data/benchmarks/sle_martinec2018/'//cname//'_'//fname// &
           '_SBK.dat', status='old', action='read', iostat=ios)
      if (ios /= 0) error stop 'dump_martinec_sle: cannot read SBK profile'
      n = 0
      do
         read(u,'(a)',iostat=ios) line
         if (ios /= 0) exit
         if (len_trim(line) == 0) cycle
         if (line(1:1) == '#') cycle
         read(line,*,iostat=ios) v
         if (ios /= 0) cycle
         n = n + 1
         if (n > nrow) exit
         c1(n)    = v(1)
         ref(n,1) = v(2);  ref(n,2) = v(3);  ref(n,3) = v(4)
         ref(n,4) = -v(5)/MS_GGIAPY
         ref(n,5) = v(6);  ref(n,6) = v(7)
      end do
      close(u)
      if (n /= nrow) error stop 'dump_martinec_sle: unexpected SBK row count'
   end subroutine msle_read_fig

   subroutine msle_sample(sht, u_lm, N_lm, rsl_lm, v_lm, esl, ptyp, pf, nrow, c1, mdl)
      !! Evaluate the converged model curves at the profile points. ptyp 'Z': col1 = colatitude
      !! at the fixed longitude pf; 'Y': col1 = 180 + lon, sampled at lon = col1, fixed colat pf.
      type(sht_grid), intent(in)  :: sht
      complex(wp),    intent(in)  :: u_lm(:), N_lm(:), rsl_lm(:), v_lm(:)
      real(wp),       intent(in)  :: esl, pf, c1(:)
      character(1),   intent(in)  :: ptyp
      integer,        intent(in)  :: nrow
      real(wp),       intent(out) :: mdl(:,:)
      real(wp) :: colat, lon, um, nm, slem, vthm, vphm
      integer  :: i
      do i = 1, nrow
         if (ptyp == 'Z') then
            colat = c1(i)*MDEG;  lon = pf
         else
            colat = pf;          lon = c1(i)*MDEG
         end if
         call sht_grid_eval_point(sht, u_lm,   colat, lon, um)
         call sht_grid_eval_point(sht, N_lm,   colat, lon, nm)
         call sht_grid_eval_point(sht, rsl_lm, colat, lon, slem)
         call sht_grid_eval_point_horiz(sht, v_lm, colat, lon, vthm, vphm)
         mdl(i,1) = um;    mdl(i,2) = vthm;  mdl(i,3) = vphm
         mdl(i,4) = nm;    mdl(i,5) = nm + esl;  mdl(i,6) = slem
      end do
   end subroutine msle_sample

   subroutine msle_errors(nrow, mdl, ref, err)
      !! Peak-normalised max errors, the horizontal pair normalised by the peak vector magnitude.
      integer,  intent(in)  :: nrow
      real(wp), intent(in)  :: mdl(:,:), ref(:,:)
      real(wp), intent(out) :: err(:)
      real(wp) :: pk(6)
      integer  :: q
      do q = 1, 6
         pk(q) = maxval(abs(ref(1:nrow,q)))
      end do
      pk(2) = maxval(sqrt(ref(1:nrow,2)**2 + ref(1:nrow,3)**2));  pk(3) = pk(2)
      do q = 1, 6
         err(q) = maxval(abs(mdl(1:nrow,q) - ref(1:nrow,q)))/pk(q)
      end do
   end subroutine msle_errors


   subroutine dump_martinec_sle(f)
      character(len=*), intent(in) :: f
      integer,  parameter :: NC = 4, NFG = 4, NRZ = 361, NRY = 720, NQ = 6, NRES = 2
      integer,  parameter :: LMS(NRES) = [64, 128]
      character(2), parameter :: CN(NC)  = ['C2', 'D3', 'E2', 'F1']
      character(5), parameter :: FN(NFG) = ['fig10', 'fig11', 'fig12', 'fig13']
      character(1), parameter :: PT(NFG) = ['Z', 'Y', 'Z', 'Y']
      character(3), parameter :: QN(NQ)  = ['u  ', 'vth', 'vph', 'N  ', 'ss ', 'sle']
      character(3), parameter :: QU(NQ)  = ['m  ', 'm  ', 'm  ', 'm  ', 'm  ', 'm  ']
      ! Martinec 2018 Table 4 (SBK file letter = paper letter + 1), as tests/test_benchmark_sle.f90
      real(wp), parameter :: H0C(NC) = [1500.0_wp, 500.0_wp, 500.0_wp, 500.0_wp]
      real(wp), parameter :: BCL(NC) = [100.0_wp, 35.0_wp, 35.0_wp, 35.0_wp]*MDEG
      real(wp), parameter :: BLN(NC) = [320.0_wp, 25.0_wp, 25.0_wp, 25.0_wp]*MDEG
      real(wp), parameter :: BMX(NC) = [760.0_wp, 3800.0_wp, 3800.0_wp, 3800.0_wp]
      real(wp), parameter :: BZ(NC)  = [1200.0_wp, 6000.0_wp, 6000.0_wp, 6000.0_wp]
      logical,  parameter :: HEAV(NC) = [.true., .false., .false., .false.]
      logical,  parameter :: FIXO(NC) = [.true., .true., .false., .false.]
      logical,  parameter :: SPN(NC)  = [.false., .false., .false., .true.]
      integer,  parameter :: NST(NC) = [MS_NGROW, MS_NGROW + MS_NHOLD, MS_NGROW + MS_NHOLD, &
                                        MS_NGROW + MS_NHOLD]
      real(wp), parameter :: PFX(NFG,NC) = reshape( &
         [75.0_wp,  25.0_wp, -40.0_wp, 100.0_wp, &
          75.0_wp,  25.0_wp,  25.0_wp,  35.0_wp, &
          75.0_wp,  25.0_wp,  25.0_wp,  35.0_wp, &
          75.0_wp,  25.0_wp,  25.0_wp,  35.0_wp], [NFG, NC])*MDEG

      real(wp), allocatable :: c1Z(:,:,:), refZ(:,:,:,:), mdlZ(:,:,:,:,:)
      real(wp), allocatable :: c1Y(:,:,:), refY(:,:,:,:), mdlY(:,:,:,:,:)
      real(wp), allocatable :: rtZ(:,:), rtY(:,:), err(:,:,:,:)
      real(wp), allocatable :: topo0(:,:), rsl(:,:), C(:,:)
      complex(wp), allocatable :: u_lm(:,:), N_lm(:,:), rsl_lm(:,:), v_lm(:,:)
      real(wp), allocatable :: esl(:), ocfrac(:), mresid(:)
      integer,  allocatable :: nspin(:)
      type(sht_grid)    :: sht
      type(earth_model) :: em
      character(8) :: sfx
      integer :: ic, k, ir, q, iz, iy
      integer :: tic, toc, trate
      integer :: lmv(NRES), nlatv(NRES), nphiv(NRES), nlmv(NRES)
      real(wp), allocatable :: wall(:,:)

      write(*,'(a)') ' [martinec_sle] Martinec 2018 SLE cases C2/D3/E2/F1 at lmax 64 and 128'
      allocate(c1Z(NRZ,2,NC), refZ(NRZ,2,NC,NQ), mdlZ(NRZ,2,NC,NQ,NRES))
      allocate(c1Y(NRY,2,NC), refY(NRY,2,NC,NQ), mdlY(NRY,2,NC,NQ,NRES))
      allocate(rtZ(NRZ,NQ), rtY(NRY,NQ), err(NQ,NFG,NC,NRES))
      allocate(esl(NC), ocfrac(NC), mresid(NC), nspin(NC), wall(NC,NRES))

      do ic = 1, NC
         do k = 1, NFG
            if (PT(k) == 'Z') then
               iz = (k + 1)/2                                   ! fig10 -> 1, fig12 -> 2
               call msle_read_fig(CN(ic), FN(k), NRZ, c1Z(:,iz,ic), rtZ)
               do q = 1, NQ;  refZ(:,iz,ic,q) = rtZ(:,q);  end do
            else
               iy = k/2                                          ! fig11 -> 1, fig13 -> 2
               call msle_read_fig(CN(ic), FN(k), NRY, c1Y(:,iy,ic), rtY)
               do q = 1, NQ;  refY(:,iy,ic,q) = rtY(:,q);  end do
            end if
         end do
      end do

      call nc_create(f, overwrite=.true.)
      call stamp(f, "martinec_sle")
      call w_one(f)
      call nc_write_dim(f, "case", x=1, dx=1, nx=NC, units="1", long_name="C2, D3, E2, F1")
      call nc_write_dim(f, "fig",  x=1, dx=1, nx=NFG, units="1", long_name="fig10..fig13")
      call nc_write_dim(f, "figZ", x=1, dx=1, nx=2, units="1", long_name="fig10, fig12 (col1 = colatitude)")
      call nc_write_dim(f, "figY", x=1, dx=1, nx=2, units="1", long_name="fig11, fig13 (col1 = 180 + longitude)")
      call nc_write_dim(f, "rowZ", x=1, dx=1, nx=NRZ, units="1")
      call nc_write_dim(f, "rowY", x=1, dx=1, nx=NRY, units="1")
      call nc_write_dim(f, "qty",  x=1, dx=1, nx=NQ, units="1", long_name="u, vth, vph, N, ss, sle")
      call nc_write_attr(f, "cases", "C2 D3 E2 F1")
      call nc_write_attr(f, "figs", "fig10 fig11 fig12 fig13")
      call nc_write_attr(f, "qtys", "u vth vph N ss sle")
      call nc_write_attr(f, "profile_type", "Z Y Z Y (Z: col1 = colatitude at fixed lon; Y: col1 = 180+lon, sampled at lon = col1, fixed colat)")
      call nc_write_attr(f, "reference", "data/benchmarks/sle_martinec2018/<case>_<fig>_SBK.dat (Martinec et al. 2018, giapy)")
      call nc_write_attr(f, "settings", "sle_solver: n_outer=3, n_inner=20, tol=1e-7, subgrid=.true., warm_start=.true.; response_init_ve (scheme fe), M3-L70-V01; rsl zeroed at the start of each history")
      call nc_write_attr(f, "tolerances", "TOL_U=0.30 TOL_N=0.12 TOL_SS=0.08 TOL_SLE=0.12 TOL_H=0.35 (peak-normalised, tests/test_benchmark_sle.f90)")
      call ws(f, "dt", MS_DT, "s", "time step (20 yr)")
      call wsi(f, "ngrow", MS_NGROW, "T1 ramp steps (10 kyr)")
      call wsi(f, "nhold", MS_NHOLD, "T1 hold steps (5 kyr)")
      call ws(f, "ice_colat", MS_ICE_COLAT, "radians", "ice-cap centre colatitude")
      call ws(f, "ice_lon", MS_ICE_LON, "radians", "ice-cap centre longitude")
      call ws(f, "alpha_max", MS_ALPHA, "radians", "ice-cap angular radius at full load")
      call ws(f, "sigma_basin", MS_SIGB, "radians", "exponential-basin decay rate")
      call ws(f, "g_giapy", MS_GGIAPY, "m s-2", "gravity used to convert SBK col5 to the geoid N = -col5/g")
      call nc_write(f, "h0", H0C, dim1="case", units="m", long_name="ice-cap height at full load")
      call nc_write(f, "basin_colat", BCL, dim1="case", units="radians", long_name="basin centre colatitude")
      call nc_write(f, "basin_lon", BLN, dim1="case", units="radians", long_name="basin centre longitude")
      call nc_write(f, "basin_bmax", BMX, dim1="case", units="m", long_name="basin bmax")
      call nc_write(f, "basin_b0", BZ, dim1="case", units="m", long_name="basin b0")
      call nc_write(f, "nsteps", NST, dim1="case", units="1", long_name="steps in the time history")
      call nc_write(f, "heaviside", merge(1, 0, HEAV), dim1="case", units="1", long_name="T0 Heaviside load history (else T1 ramp then hold)")
      call nc_write(f, "fixed_ocean", merge(1, 0, FIXO), dim1="case", units="1", long_name="SLE1 fixed ocean geometry (else SLE2 migrating)")
      call nc_write(f, "spinup", merge(1, 0, SPN), dim1="case", units="1", long_name="F1 paleotopography fixed-point spin-up")
      call nc_write(f, "pfix", PFX, dim1="fig", dim2="case", units="radians", long_name="fixed longitude (Z profiles) or colatitude (Y profiles)")
      call nc_write(f, "c1_Z", c1Z, dim1="rowZ", dim2="figZ", dim3="case", units="degrees", long_name="SBK col1: colatitude")
      call nc_write(f, "c1_Y", c1Y, dim1="rowY", dim2="figY", dim3="case", units="degrees", long_name="SBK col1: 180 + longitude")
      do q = 1, NQ
         call nc_write(f, trim(QN(q))//"_Z_ref", refZ(:,:,:,q), dim1="rowZ", dim2="figZ", dim3="case", &
              units=trim(QU(q)), long_name="SBK reference "//trim(QN(q))//" on the Z profiles")
         call nc_write(f, trim(QN(q))//"_Y_ref", refY(:,:,:,q), dim1="rowY", dim2="figY", dim3="case", &
              units=trim(QU(q)), long_name="SBK reference "//trim(QN(q))//" on the Y profiles")
      end do

      ! --- run the four cases at each resolution and sample the profiles -------------------
      call system_clock(count_rate=trate)
      do ir = 1, NRES
         call sht_grid_init(sht, LMS(ir), nlat=2*LMS(ir), nphi=4*LMS(ir))
         em = build_M3L70V01()
         lmv(ir) = LMS(ir);  nlatv(ir) = sht%nlat;  nphiv(ir) = sht%nphi;  nlmv(ir) = sht%nlm
         write(*,'(a,i0,a,i0,a,i0,a,i0)') '   lmax=', LMS(ir), '  grid ', sht%nphi, ' x ', &
              sht%nlat, ', nlm=', sht%nlm
         allocate(topo0(sht%nphi,sht%nlat), rsl(sht%nphi,sht%nlat), C(sht%nphi,sht%nlat))
         allocate(u_lm(sht%nlm,NC), N_lm(sht%nlm,NC), rsl_lm(sht%nlm,NC), v_lm(sht%nlm,NC))
         do ic = 1, NC
            call system_clock(tic)
            call msle_run(sht, em, ic, H0C(ic), BCL(ic), BLN(ic), BMX(ic), BZ(ic), MS_NGROW, &
                          NST(ic), HEAV(ic), FIXO(ic), SPN(ic), topo0, rsl, C, &
                          u_lm(:,ic), N_lm(:,ic), rsl_lm(:,ic), v_lm(:,ic), &
                          esl(ic), ocfrac(ic), mresid(ic), nspin(ic))
            call system_clock(toc)
            wall(ic,ir) = real(toc - tic, wp)/real(trate, wp)
            write(*,'(3a,i0,a,i0,a,i2,a,f9.4,a,f7.4,a,es9.2,a,f8.1,a)') '   ', CN(ic), &
                 ' lmax=', LMS(ir), ', ', NST(ic), ' steps, nspin=', nspin(ic), ':  esl=', esl(ic), &
                 '  ocean_frac=', ocfrac(ic), '  mass_resid=', mresid(ic), '  [', wall(ic,ir), ' s]'
            do k = 1, NFG
               if (PT(k) == 'Z') then
                  iz = (k + 1)/2
                  call msle_sample(sht, u_lm(:,ic), N_lm(:,ic), rsl_lm(:,ic), v_lm(:,ic), esl(ic), &
                                   PT(k), PFX(k,ic), NRZ, c1Z(:,iz,ic), rtZ)
                  do q = 1, NQ;  mdlZ(:,iz,ic,q,ir) = rtZ(:,q);  end do
                  call msle_errors(NRZ, rtZ, refZ(:,iz,ic,:), err(:,k,ic,ir))
               else
                  iy = k/2
                  call msle_sample(sht, u_lm(:,ic), N_lm(:,ic), rsl_lm(:,ic), v_lm(:,ic), esl(ic), &
                                   PT(k), PFX(k,ic), NRY, c1Y(:,iy,ic), rtY)
                  do q = 1, NQ;  mdlY(:,iy,ic,q,ir) = rtY(:,q);  end do
                  call msle_errors(NRY, rtY, refY(:,iy,ic,:), err(:,k,ic,ir))
               end if
               write(*,'(6x,3a,6f8.2)') CN(ic), ' ', FN(k), 100*err(:,k,ic,ir)
            end do
         end do
         write(sfx,'(a,i0)') "_lm", LMS(ir)
         do q = 1, NQ
            call nc_write(f, trim(QN(q))//"_Z"//trim(sfx), mdlZ(:,:,:,q,ir), dim1="rowZ", &
                 dim2="figZ", dim3="case", units=trim(QU(q)), &
                 long_name="Fortran model "//trim(QN(q))//" on the Z profiles")
            call nc_write(f, trim(QN(q))//"_Y"//trim(sfx), mdlY(:,:,:,q,ir), dim1="rowY", &
                 dim2="figY", dim3="case", units=trim(QU(q)), &
                 long_name="Fortran model "//trim(QN(q))//" on the Y profiles")
         end do
         call nc_write(f, "err"//trim(sfx), err(:,:,:,ir), dim1="qty", dim2="fig", dim3="case", &
              units="1", long_name="peak-normalised max error of the Fortran model vs SBK")
         call nc_write(f, "esl"//trim(sfx), esl, dim1="case", units="m", &
              long_name="eustatic offset dphi of the final step")
         call nc_write(f, "ocean_frac"//trim(sfx), ocfrac, dim1="case", units="1", &
              long_name="migrated ocean fraction int C dOmega / 4pi")
         call nc_write(f, "mass_resid"//trim(sfx), mresid, dim1="case", units="1", &
              long_name="relative ocean-mass residual of the final step")
         call nc_write(f, "nspin"//trim(sfx), nspin, dim1="case", units="1", &
              long_name="paleotopography spin-up passes (1 without spin-up)")
         call nc_write(f, "wall"//trim(sfx), wall(:,ir), dim1="case", units="s", &
              long_name="Fortran wall time of the case")
         deallocate(topo0, rsl, C, u_lm, N_lm, rsl_lm, v_lm)
         call sht_grid_destroy(sht)
      end do
      call nc_write_dim(f, "res", x=1, dx=1, nx=NRES, units="1", long_name="lmax 64, 128")
      call nc_write(f, "lmax", lmv, dim1="res", units="1", long_name="truncation degree")
      call nc_write(f, "nlat", nlatv, dim1="res", units="1", long_name="Gauss latitudes")
      call nc_write(f, "nphi", nphiv, dim1="res", units="1", long_name="longitudes")
      call nc_write(f, "nlm", nlmv, dim1="res", units="1", long_name="spectral coefficients")

      deallocate(wall,c1Z, refZ, mdlZ, c1Y, refY, mdlY, rtZ, rtY, err, esl, ocfrac, mresid, nspin)
   end subroutine dump_martinec_sle

end program dump_reference
