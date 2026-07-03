# astrovodka.sources.sph_source — SPHSource

Abstracts an HI-emitting object from SPH particle arrays. Supports legacy distance mode and new cosmological redshift mode.

## SPHSource

### `__init__(self, *, distance=None, redshift=None, cosmology=None, vpeculiar=0 km/s, rotation=None, L_coords=None, ra=0 deg, dec=0 deg, h=0.7, T_g=None, mHI_g=None, xyz_g=None, vxyz_g=None, hsm_g=None, coordinate_axis=None, coordinate_frame=ICRS())`

| Name | Type | Description |
|------|------|-------------|
| `distance` | `Quantity[Mpc] \| None` | Legacy mode: physical distance (use with `h`) |
| `redshift` | `float \| None` | Cosmological redshift mode (use with `cosmology`) |
| `cosmology` | `Cosmology \| None` | Astropy cosmology (default: Planck15) |
| `vpeculiar` | `Quantity[km/s]` | Peculiar velocity (default: 0) |
| `rotation` | `Rotation \| None` | scipy Rotation |
| `L_coords` | `L_coords \| None` | Angular momentum alignment |
| `ra` / `dec` | `Quantity[deg]` | Source centroid (default: 0) |
| `h` | `float` | Hubble constant (legacy mode, default: 0.7) |
| `T_g` | `Quantity[K] \| None` | Particle temperatures |
| `mHI_g` | `Quantity[M_sun]` | HI mass per particle |
| `xyz_g` | `Quantity[kpc]` | Positions. Shape (3,N) or (N,3). x = LOS |
| `vxyz_g` | `Quantity[km/s]` | Velocities. Must match xyz_g |
| `hsm_g` | `Quantity[kpc] \| None` | SPH smoothing lengths (**FWHM**) |
| `coordinate_axis` | `int \| None` | Auto-detected |
| `coordinate_frame` | `BaseRADecFrame` | Frame (default: ICRS) |

> **Modes:**
> - **Redshift** (new): set `redshift` + optional `cosmology`. Computes D<sub>L</sub>, D<sub>A</sub>, vsys=H(z)×D<sub>L</sub>/(1+z).
> - **Distance** (legacy): set `distance` + `h`. Computes vsys=h×100×D + vpeculiar.
> - Cannot set both (raises ValueError).
> - Smoothing lengths MUST be FWHM (not kernel radius).

---

## Key Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `vsys` | `Quantity[km/s]` | Systemic velocity. Use as `DataCube.spectral_centre` |
| `distance` | `Quantity[Mpc]` | Luminosity distance D<sub>L</sub> |
| `angular_diameter_distance` | `Quantity[Mpc]` | D<sub>A</sub> |
| `vhubble` | `Quantity[km/s]` | Systemic Hubble velocity |
| `npart` | `int` | Particles after pruning |
| `input_mass` | `Quantity[M_sun]` | Total HI mass before pruning |

## Methods

### `rotate(self, rotation=None, *, L_coords=None)`

Apply rotation. Returns `Rotation`.

### `translate(self, translation_vector: Quantity[kpc])`

Translate positions.

### `boost(self, boost_vector: Quantity[km/s])`

Add velocity offset.

### `preview(self, max_points=5000, fig=1, lim=None, vlim=None, point_scaling='auto', title='', save=None)`

3-panel figure: sky + 2 PV diagrams. Returns `Figure`.

### `save_current_rotation(self, fname: str)`

Save 3×3 rotation matrix to text.

### `save_current_affine_transformations(self, fname: str)`

Save (2,4,4) affine transforms to `.npy`.

### `load_affine_transformations(self, fname: str)`

Load and apply affine transforms from `.npy`.
