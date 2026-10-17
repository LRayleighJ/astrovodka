# astrovodka.L_coords — Angular Momentum Coordinates

A `NamedTuple` specifying source orientation based on angular momentum.

MARTINI 2.1.16 and later ship the same type as `martini.L_coords`; on 2.1.9–2.1.15 astrovodka defines it locally.

## `L_coords(incl=0 deg, pa=270 deg, az_rot=0 deg)`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `incl` | `Quantity[deg]` | 0° | Inclination (0°=face-on, 90°=edge-on) |
| `pa` | `Quantity[deg]` | 270° | Position angle on sky (E of N, receding side) |
| `az_rot` | `Quantity[deg]` | 0° | Azimuthal rotation around L axis |

**Algorithm:** The angular momentum vector of the central 1/3 HI particles is aligned to x-axis (face-on). Then: `az_rot` around L → `incl` tilt (y-axis) → `pa` sky rotation (x-axis). All rotations right-handed and extrinsic.

**Usage:** `source.rotate(L_coords=L_coords(incl=60*U.deg, pa=270*U.deg))`
