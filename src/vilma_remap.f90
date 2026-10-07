module vilma_remap
   !! Bidirectional conservative/bilinear remapping between a host lon-lat grid and
   !! the model's SHTns Gauss-Legendre grid, built on the fesm-utils `coords` library
   !! (in-package SCRIP-style weights, gen="coords", no CDO).
   !!
   !! Two legs, mirroring CLIMBER-X's VILMA1 coupling:
   !!   host lon-lat --(conservative)--> Gauss   [remap_to_gauss]   mass-bearing fields (ice)
   !!   Gauss --(bilinear)--> host lon-lat        [remap_to_ll]      smooth fields (rsl)
   !!
   !! Build the map pair once with remap_init, then apply per field / time slice. The
   !! coupling layer (vilma_coupling) holds one of these when the host grid differs from
   !! the model Gauss grid, and drives h_ice in / rsl out through it; the standalone
   !! driver (vilma_drive) uses it to remap lon-lat forcing onto the Gauss grid.
   !!
   !! Caching: map_init writes each weight set to a SCRIP(-superset) NetCDF cache under
   !! `fldr` (default "maps") and reloads it on the next run, keyed by the grid names —
   !! which we tag with the grid dimensions so a resolution change invalidates the
   !! cache. The conservative weight build (polygon clipping) is the costly step; the
   !! cache lets a restart skip it. (Same-dimension/different-axes host grids would
   !! collide on the name — clear `fldr` if the host grid changes at fixed resolution.)
   !!
   !! Mass note: `coords` derives target cell boundaries from axis midpoints, so the
   !! Gauss-cell areas it uses are not exactly the SHTns Gauss quadrature weights. The
   !! conservative map is conservative w.r.t. its own areas, but when SHTns re-integrates
   !! the remapped field its quadrature total differs by O(1/nlat^2). For a mass-bearing
   !! field (ice) pass conserve_mass=.true.: a single global factor rescales the result
   !! so its SHTns surface integral equals the source area-integral exactly (the whole
   !! sphere is 4*pi sr, used to convert the conserved coords-area total to steradians
   !! without needing the planet radius). Geometry fields (bed) are remapped as-is.
   use vilma_precision, only: wp
   use vilma_sht,       only: sht_grid, sht_grid_surface_integral
   use coords,       only: grid_class, grid_init, map_class, map_init, map_field
   implicit none
   private

   public :: remap_ll_gauss
   public :: remap_init, remap_to_gauss, remap_to_ll

   real(wp), parameter :: RAD2DEG = 57.295779513082323_wp
   real(wp), parameter :: FOURPI  = 12.566370614359172_wp
   integer,  parameter :: BIL_NEIGHBORS = 8   !! neighbour pool for the bilinear leg

   type :: remap_ll_gauss
      !! A precomputed map pair between a host lon-lat grid and the model's SHTns
      !! Gauss grid. Build once with remap_init; remap_to_gauss / remap_to_ll apply it.
      type(grid_class) :: ll                  !! host lon-lat grid
      type(grid_class) :: gauss               !! Gauss target (matches sht), lat ascending
      type(map_class)  :: to_gauss            !! ll -> Gauss, conservative
      type(map_class)  :: to_ll               !! Gauss -> ll, bilinear
      integer :: nlon = 0, nlat_ll = 0        !! host dimensions
      integer :: nphi = 0, nlat = 0           !! Gauss (SHTns) dimensions
   end type remap_ll_gauss

contains

   subroutine remap_init(self, sht, lon_src, lat_src, fldr)
      !! Build (or load from cache) the conservative (ll->Gauss) and bilinear
      !! (Gauss->ll) maps for source cell-centre axes lon_src(:), lat_src(:) [degrees,
      !! ascending] against the SHTns Gauss grid. The Gauss grid is built with ASCENDING
      !! latitude (south first), which coords expects; the apply routines flip into/out
      !! of the SHTns (north-first) row order. fldr is the weight-cache directory.
      type(remap_ll_gauss), intent(out) :: self
      type(sht_grid),       intent(in)  :: sht
      real(wp),             intent(in)  :: lon_src(:)   !! source longitudes [deg]
      real(wp),             intent(in)  :: lat_src(:)   !! source latitudes  [deg]
      character(len=*), optional, intent(in) :: fldr    !! weight-cache dir (default "maps")
      real(wp), allocatable :: lon_t(:), lat_t(:)
      character(len=:), allocatable :: cache
      character(len=64) :: ll_name, gauss_name
      integer :: j

      self%nlon = size(lon_src);  self%nlat_ll = size(lat_src)
      self%nphi = sht%nphi;       self%nlat    = sht%nlat
      cache = "maps";  if (present(fldr)) cache = trim(fldr)
      call execute_command_line("mkdir -p '"//cache//"'")

      allocate(lon_t(self%nphi), lat_t(self%nlat))
      do j = 1, self%nphi
         lon_t(j) = sht%lon(j)*RAD2DEG                  ! SHTns longitudes [0,360)
      end do
      do j = 1, self%nlat
         ! SHTns colat is ascending (north -> south), so lat descends; reverse it so
         ! the target axis is ascending (south -> north), as coords expects.
         lat_t(j) = 90.0_wp - sht%colat(self%nlat - j + 1)*RAD2DEG
      end do

      ! tag grid names with their dimensions so the weight cache is resolution-specific
      write(ll_name,    '(a,i0,a,i0)') "ll_",    self%nlon, "x", self%nlat_ll
      write(gauss_name, '(a,i0,a,i0)') "gauss_", self%nphi, "x", self%nlat
      call grid_init(self%ll,    name=trim(ll_name),    mtype="latlon", units="degrees", &
                     lon180=.true.,  x=lon_src, y=lat_src)
      call grid_init(self%gauss, name=trim(gauss_name), mtype="latlon", units="degrees", &
                     lon180=.false., x=lon_t,   y=lat_t)

      ! in-package weights (gen="coords"), auto-cached under fldr and reloaded if present
      call map_init(self%to_gauss, self%ll, self%gauss, method="con", gen="coords", fldr=cache)
      call map_init(self%to_ll, self%gauss, self%ll, method="bilinear", &
                    max_neighbors=BIL_NEIGHBORS, gen="coords", fldr=cache)
   end subroutine remap_init

   subroutine remap_to_gauss(self, sht, f_ll, f_gauss, conserve_mass)
      !! Conservatively remap f_ll(nlon, nlat_ll) onto f_gauss(nphi, nlat) in the SHTns
      !! spatial layout. conserve_mass (default .false.): rescale f_gauss so its SHTns
      !! surface integral equals the source area-integral (use for ice thickness; leave
      !! off for geometry such as bed).
      type(remap_ll_gauss), intent(in)  :: self
      type(sht_grid),       intent(in)  :: sht
      real(wp),             intent(in)  :: f_ll(:,:)    !! (nlon, nlat_ll)
      real(wp),             intent(out) :: f_gauss(:,:) !! (nphi, nlat)
      logical, optional,    intent(in)  :: conserve_mass
      real(wp), allocatable :: vt(:,:)
      logical,  allocatable :: m2(:,:)
      real(wp) :: src_int_sr, gauss_int
      integer  :: j
      logical  :: do_mass

      if (size(f_ll,1) /= self%nlon .or. size(f_ll,2) /= self%nlat_ll) &
         error stop 'vilma_remap: source field shape /= map source grid'
      if (size(f_gauss,1) /= self%nphi .or. size(f_gauss,2) /= self%nlat) &
         error stop 'vilma_remap: target field shape /= SHTns grid'

      allocate(vt(self%nphi, self%nlat), m2(self%nphi, self%nlat))
      call map_field(self%to_gauss, "f", f_ll, vt, stat="mean", mask2=m2)
      where (.not. m2) vt = 0.0_wp                       ! uncovered cells (global src: none)

      do j = 1, self%nlat                                 ! ascending-lat -> SHTns north-first
         f_gauss(:, j) = vt(:, self%nlat - j + 1)
      end do

      do_mass = .false.;  if (present(conserve_mass)) do_mass = conserve_mass
      if (do_mass) then
         ! conserved source total as a steradian integral: the map preserves
         ! sum(f_ll*src_area); convert to sr via 4*pi/total_area (full-sphere grid).
         src_int_sr = sum(f_ll * self%ll%area) / sum(self%ll%area) * FOURPI
         gauss_int  = sht_grid_surface_integral(sht, f_gauss)
         if (abs(gauss_int) > tiny(1.0_wp)) f_gauss = f_gauss * (src_int_sr/gauss_int)
      end if
   end subroutine remap_to_gauss

   subroutine remap_to_ll(self, f_gauss, f_ll)
      !! Bilinearly remap f_gauss(nphi, nlat) in the SHTns spatial layout onto the host
      !! lon-lat grid f_ll(nlon, nlat_ll). Used for smooth fields (relative sea level)
      !! returned to the host; the host reconstructs z_bed on its own high-resolution
      !! bed as z_bed_eq - f_ll, so no mass rescale is applied here.
      type(remap_ll_gauss), intent(in)  :: self
      real(wp),             intent(in)  :: f_gauss(:,:) !! (nphi, nlat), SHTns north-first
      real(wp),             intent(out) :: f_ll(:,:)    !! (nlon, nlat_ll)
      real(wp), allocatable :: g_asc(:,:)
      logical,  allocatable :: m2(:,:)
      integer :: j

      if (size(f_gauss,1) /= self%nphi .or. size(f_gauss,2) /= self%nlat) &
         error stop 'vilma_remap: source field shape /= SHTns grid'
      if (size(f_ll,1) /= self%nlon .or. size(f_ll,2) /= self%nlat_ll) &
         error stop 'vilma_remap: target field shape /= map target grid'

      allocate(g_asc(self%nphi, self%nlat), m2(self%nlon, self%nlat_ll))
      do j = 1, self%nlat                                 ! SHTns north-first -> ascending-lat
         g_asc(:, j) = f_gauss(:, self%nlat - j + 1)
      end do

      call map_field(self%to_ll, "f", g_asc, f_ll, method="bilinear", mask2=m2)
      where (.not. m2) f_ll = 0.0_wp                      ! uncovered host cells (full-sphere Gauss src: none)
   end subroutine remap_to_ll

end module vilma_remap
