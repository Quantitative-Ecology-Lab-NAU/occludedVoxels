# occludedVoxels

Voxel-based occlusion and visibility analysis for terrestrial LiDAR (TLS) point clouds. This repo classifies voxels inside a quadrat (or other AOI) as **populated**, **empty-visible**, or **occluded** using explicit ray tracing from a scanner origin.

## What this does

Given a voxelized point cloud and a scanner location, `occludedVoxels`:

1. Builds a set of **populated (occupied) voxels** from LiDAR returns.
2. Enumerates the full **analysis voxel grid** within the AOI.
3. Casts a ray from the scanner to each voxel center.
4. Determines whether each voxel is:
   - **Populated (1)**: contains ≥1 LiDAR return
   - **Empty visible (0)**: contains no returns and is reachable by an unobstructed ray
   - **Occluded (2)**: contains no returns but the ray intersects populated voxels before reaching it
5. Summarizes counts (and optional volumes) of populated, empty-visible, and occluded voxels.

This partitions the AOI volume into observed structure, confirmed empty space, and unobservable (masked) space.
<img width="1498" height="673" alt="image" src="https://github.com/user-attachments/assets/c158c0a4-a11f-451d-9933-d1ef8db29859" />
## How voxel states are determined

Let:
- `P` = set of populated voxels (from LiDAR returns)
- `S` = scanner origin
- `v` = target voxel center
- `R(S→v)` = ray from scanner to target

Classification:

- **Populated (1)** if `v ∈ P`
- **Empty visible (0)** if `v ∉ P` and the ray reaches `v` without intersecting any populated voxel:
  - `R(S→v) ∩ P = ∅`
- **Occluded (2)** if `v ∉ P` and the ray intersects populated voxels before reaching `v`:
  - `R(S→v) ∩ P ≠ ∅`

All populated voxels can act as occluders, including vegetation between the scanner and the AOI.
