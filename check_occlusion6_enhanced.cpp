#include <Rcpp.h>
#include <omp.h>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <atomic>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// Type alias for voxel keys
using VoxelKey = std::tuple<int, int, int>;

// Spatial hash function for fast voxel lookup
struct VoxelHash {
  std::size_t operator()(const VoxelKey& k) const {
    auto h1 = std::hash<int>{}(std::get<0>(k));
    auto h2 = std::hash<int>{}(std::get<1>(k));
    auto h3 = std::hash<int>{}(std::get<2>(k));
    return h1 ^ (h2 << 1) ^ (h3 << 2);
  }
};

// Convert continuous coordinates to voxel grid indices
inline VoxelKey get_voxel_key(double x, double y, double z, double voxel_size) {
  return std::make_tuple(
    static_cast<int>(floor(x / voxel_size)),
    static_cast<int>(floor(y / voxel_size)),
    static_cast<int>(floor(z / voxel_size))
  );
}

// ============================================================
// Flat 3D boolean grid for O(1) obstacle lookups.
// Replaces unordered_set for ~5-10x faster ray tracing.
// Memory: NX * NY * NZ bytes (e.g., 500x500x150 = 37 MB).
// ============================================================
struct VoxelGrid3D {
  std::vector<uint8_t> data;
  int ox, oy, oz;   // origin offsets (minimum voxel index)
  int nx, ny, nz;   // grid dimensions
  
  void init(int min_x, int min_y, int min_z,
            int max_x, int max_y, int max_z) {
    ox = min_x; oy = min_y; oz = min_z;
    nx = max_x - min_x + 1;
    ny = max_y - min_y + 1;
    nz = max_z - min_z + 1;
    data.assign((size_t)nx * ny * nz, 0);
  }
  
  inline void set(int ix, int iy, int iz) {
    int lx = ix - ox, ly = iy - oy, lz = iz - oz;
    if (lx >= 0 && lx < nx && ly >= 0 && ly < ny && lz >= 0 && lz < nz)
      data[(size_t)lx * ny * nz + (size_t)ly * nz + lz] = 1;
  }
  
  inline bool get(int ix, int iy, int iz) const {
    int lx = ix - ox, ly = iy - oy, lz = iz - oz;
    if (lx < 0 || lx >= nx || ly < 0 || ly >= ny || lz < 0 || lz >= nz)
      return false;
    return data[(size_t)lx * ny * nz + (size_t)ly * nz + lz] != 0;
  }
};

// ============================================================
// DDA (Amanatides-Woo) 3D ray traversal using flat grid.
// FAT-RAY variant: at each step, checks the traversed voxel
// PLUS its two face-adjacent neighbors perpendicular to the
// axis just stepped.  This closes the 1 cm inter-voxel gaps
// that allow rays to thread between adjacent obstacle voxels
// at oblique angles, while adding only 4 extra array lookups
// per step (still O(1) each on the flat grid).
//
// "last_axis" tracks which axis was just stepped so we know
// which plane to check perpendicular neighbors in.
// ============================================================
inline bool is_ray_occluded_dda_fat(
    double start_x, double start_y, double start_z,
    double end_x, double end_y, double end_z,
    double voxel_size,
    const VoxelGrid3D& grid) {
  
  int ix = static_cast<int>(floor(start_x / voxel_size));
  int iy = static_cast<int>(floor(start_y / voxel_size));
  int iz = static_cast<int>(floor(start_z / voxel_size));
  
  int ix_end = static_cast<int>(floor(end_x / voxel_size));
  int iy_end = static_cast<int>(floor(end_y / voxel_size));
  int iz_end = static_cast<int>(floor(end_z / voxel_size));
  
  if (ix == ix_end && iy == iy_end && iz == iz_end) return false;
  
  double dx = end_x - start_x;
  double dy = end_y - start_y;
  double dz = end_z - start_z;
  
  int step_x = (dx >= 0) ? 1 : -1;
  int step_y = (dy >= 0) ? 1 : -1;
  int step_z = (dz >= 0) ? 1 : -1;
  
  double tDeltaX, tDeltaY, tDeltaZ;
  double tMaxX, tMaxY, tMaxZ;
  
  if (dx != 0.0) {
    tDeltaX = fabs(voxel_size / dx);
    tMaxX = (dx >= 0) ? ((ix + 1) * voxel_size - start_x) / dx
                       : (ix * voxel_size - start_x) / dx;
  } else { tDeltaX = 1e30; tMaxX = 1e30; }
  
  if (dy != 0.0) {
    tDeltaY = fabs(voxel_size / dy);
    tMaxY = (dy >= 0) ? ((iy + 1) * voxel_size - start_y) / dy
                       : (iy * voxel_size - start_y) / dy;
  } else { tDeltaY = 1e30; tMaxY = 1e30; }
  
  if (dz != 0.0) {
    tDeltaZ = fabs(voxel_size / dz);
    tMaxZ = (dz >= 0) ? ((iz + 1) * voxel_size - start_z) / dz
                       : (iz * voxel_size - start_z) / dz;
  } else { tDeltaZ = 1e30; tMaxZ = 1e30; }
  
  int max_steps = abs(ix_end - ix) + abs(iy_end - iy) + abs(iz_end - iz) + 2;
  int last_axis = -1; // 0=x, 1=y, 2=z
  
  for (int step = 0; step < max_steps; ++step) {
    if (tMaxX < tMaxY) {
      if (tMaxX < tMaxZ) { ix += step_x; tMaxX += tDeltaX; last_axis = 0; }
      else               { iz += step_z; tMaxZ += tDeltaZ; last_axis = 2; }
    } else {
      if (tMaxY < tMaxZ) { iy += step_y; tMaxY += tDeltaY; last_axis = 1; }
      else               { iz += step_z; tMaxZ += tDeltaZ; last_axis = 2; }
    }
    
    if (ix == ix_end && iy == iy_end && iz == iz_end) break;
    
    // Check the exact voxel
    if (grid.get(ix, iy, iz)) return true;
    
    // Check 4 face-adjacent neighbors perpendicular to the step axis.
    // This closes gaps at oblique angles without changing the obstacle set.
    if (last_axis == 0) {       // stepped along X → check ±Y and ±Z
      if (grid.get(ix, iy+1, iz) || grid.get(ix, iy-1, iz) ||
          grid.get(ix, iy, iz+1) || grid.get(ix, iy, iz-1)) return true;
    } else if (last_axis == 1) { // stepped along Y → check ±X and ±Z
      if (grid.get(ix+1, iy, iz) || grid.get(ix-1, iy, iz) ||
          grid.get(ix, iy, iz+1) || grid.get(ix, iy, iz-1)) return true;
    } else {                    // stepped along Z → check ±X and ±Y
      if (grid.get(ix+1, iy, iz) || grid.get(ix-1, iy, iz) ||
          grid.get(ix, iy+1, iz) || grid.get(ix, iy-1, iz)) return true;
    }
  }
  
  return false;
}

// ============================================================
// DDA using hash set (kept for the legacy check_occlusion fn)
// ============================================================
inline bool is_ray_occluded_dda(
    double start_x, double start_y, double start_z,
    double end_x, double end_y, double end_z,
    double voxel_size,
    const std::unordered_set<VoxelKey, VoxelHash>& obstacles) {
  
  int ix = static_cast<int>(floor(start_x / voxel_size));
  int iy = static_cast<int>(floor(start_y / voxel_size));
  int iz = static_cast<int>(floor(start_z / voxel_size));
  
  int ix_end = static_cast<int>(floor(end_x / voxel_size));
  int iy_end = static_cast<int>(floor(end_y / voxel_size));
  int iz_end = static_cast<int>(floor(end_z / voxel_size));
  
  if (ix == ix_end && iy == iy_end && iz == iz_end) return false;
  
  double dx = end_x - start_x;
  double dy = end_y - start_y;
  double dz = end_z - start_z;
  
  int step_x = (dx >= 0) ? 1 : -1;
  int step_y = (dy >= 0) ? 1 : -1;
  int step_z = (dz >= 0) ? 1 : -1;
  
  double tDeltaX, tDeltaY, tDeltaZ;
  double tMaxX, tMaxY, tMaxZ;
  
  if (dx != 0.0) {
    tDeltaX = fabs(voxel_size / dx);
    tMaxX = (dx >= 0) ? ((ix + 1) * voxel_size - start_x) / dx
                       : (ix * voxel_size - start_x) / dx;
  } else { tDeltaX = 1e30; tMaxX = 1e30; }
  
  if (dy != 0.0) {
    tDeltaY = fabs(voxel_size / dy);
    tMaxY = (dy >= 0) ? ((iy + 1) * voxel_size - start_y) / dy
                       : (iy * voxel_size - start_y) / dy;
  } else { tDeltaY = 1e30; tMaxY = 1e30; }
  
  if (dz != 0.0) {
    tDeltaZ = fabs(voxel_size / dz);
    tMaxZ = (dz >= 0) ? ((iz + 1) * voxel_size - start_z) / dz
                       : (iz * voxel_size - start_z) / dz;
  } else { tDeltaZ = 1e30; tMaxZ = 1e30; }
  
  int max_steps = abs(ix_end - ix) + abs(iy_end - iy) + abs(iz_end - iz) + 2;
  
  for (int step = 0; step < max_steps; ++step) {
    if (tMaxX < tMaxY) {
      if (tMaxX < tMaxZ) { ix += step_x; tMaxX += tDeltaX; }
      else               { iz += step_z; tMaxZ += tDeltaZ; }
    } else {
      if (tMaxY < tMaxZ) { iy += step_y; tMaxY += tDeltaY; }
      else               { iz += step_z; tMaxZ += tDeltaZ; }
    }
    
    if (ix == ix_end && iy == iy_end && iz == iz_end) break;
    if (obstacles.count(std::make_tuple(ix, iy, iz)) > 0) return true;
  }
  
  return false;
}

// ============================================================
// Original fixed-step ray marching (kept for coarse viz grid)
// ============================================================
inline bool is_ray_occluded(
    double start_x, double start_y, double start_z,
    double end_x, double end_y, double end_z,
    double voxel_size,
    const std::unordered_map<VoxelKey, int, VoxelHash>& voxel_map) {
  
  double dx = end_x - start_x;
  double dy = end_y - start_y;
  double dz = end_z - start_z;
  double distance = sqrt(dx*dx + dy*dy + dz*dz);
  
  if (distance < voxel_size * 0.1) return false;
  
  dx /= distance;
  dy /= distance;
  dz /= distance;
  
  int num_steps = static_cast<int>(distance / (voxel_size * 0.5)) + 1;
  double step_size = distance / num_steps;
  
  for (int step = 1; step < num_steps; ++step) {
    double t = step * step_size;
    double check_x = start_x + t * dx;
    double check_y = start_y + t * dy;
    double check_z = start_z + t * dz;
    
    auto key = get_voxel_key(check_x, check_y, check_z, voxel_size);
    auto it = voxel_map.find(key);
    if (it != voxel_map.end() && it->second > 0) return true;
  }
  
  return false;
}

// ============================================================
// Original occlusion check (used for coarse visualization grid)
// 
// NOTE: All Rcpp/R objects are copied into plain C++ containers
// BEFORE the OpenMP parallel region, and results are copied back
// AFTER. This avoids R API calls (and R's stack-limit check)
// inside OpenMP threads, which would intermittently crash R with
// "C stack usage ... is too close to the limit" because OpenMP
// thread stacks are at different addresses than R's main stack.
// ============================================================
// [[Rcpp::export]]
IntegerVector check_occlusion(NumericMatrix voxels, NumericVector scanner_location, 
                               double voxel_size, int ncores) {
  int n = voxels.nrow();
  
  // --- Copy Rcpp data into plain C++ containers (R-API-free) ---
  std::vector<double> vox_x(n), vox_y(n), vox_z(n), vox_v(n);
  for (int i = 0; i < n; ++i) {
    vox_x[i] = voxels(i, 0);
    vox_y[i] = voxels(i, 1);
    vox_z[i] = voxels(i, 2);
    vox_v[i] = voxels(i, 3);
  }
  double scan_x = scanner_location[0];
  double scan_y = scanner_location[1];
  double scan_z = scanner_location[2];
  std::vector<int> occ(n, 0);  // plain C++ result vector
  
  // Build voxel map (single-threaded, uses plain C++ data)
  std::unordered_map<VoxelKey, int, VoxelHash> voxel_map;
  voxel_map.reserve(n);
  for (int i = 0; i < n; ++i) {
    auto key = get_voxel_key(vox_x[i], vox_y[i], vox_z[i], voxel_size);
    voxel_map[key] = static_cast<int>(vox_v[i]);
  }
  
  omp_set_num_threads(ncores > 1 ? ncores : 1);
  
  int progress_interval = n / 20;
  if (progress_interval < 10000) progress_interval = 10000;
  std::atomic<int> completed(0);
  
  Rcout << "  [Viz grid] Processing " << n << " voxels with " << ncores << " cores...\n";
  Rcout << "  Progress: 0%";
  Rcout.flush();
  
  // --- OpenMP region: only plain C++ data is touched ---
#pragma omp parallel for schedule(dynamic, 1024) if(ncores > 1)
  for (int i = 0; i < n; ++i) {
    if (vox_v[i] > 0) {
      occ[i] = 1;
    } else {
      bool ray_blocked = is_ray_occluded(
        scan_x, scan_y, scan_z,
        vox_x[i], vox_y[i], vox_z[i],
        voxel_size, voxel_map);
      occ[i] = ray_blocked ? 2 : 0;
    }
    
    int current = ++completed;
    if (current % progress_interval == 0) {
#pragma omp critical
      {
        // printf is thread-safe; Rcout is NOT safe inside OpenMP
        Rprintf("\r  Progress: %d%%", (current * 100) / n);
      }
    }
  }
  
  Rcout << "\r  Progress: 100% - Complete!\n";
  Rcout.flush();
  
  // --- Copy results back to Rcpp IntegerVector ---
  IntegerVector occluded(n);
  for (int i = 0; i < n; ++i) occluded[i] = occ[i];
  return occluded;
}

// ============================================================
// ACCURATE occlusion check for single-return TLS.
//
// Uses backward DDA (scanner → target voxel) with a fat-ray
// neighbor check to close 1cm inter-voxel gaps.
//
// At each DDA step the traversed voxel + its 4 face-adjacent
// neighbors perpendicular to the step axis are checked.
// This prevents oblique rays from threading through diagonal
// gaps between adjacent obstacle voxels, ensuring contiguous
// shadow regions behind vegetation surfaces extend properly.
//
// Obstacle lookup uses a flat 3D boolean array (O(1) per check,
// cache-friendly) rather than an unordered_set.
// ============================================================
// [[Rcpp::export]]
IntegerVector check_occlusion_accurate(
    NumericMatrix analysis_grid,    // N x 4: X, Y, Z, V1(0/1)
    NumericMatrix obstacles,        // M x 3: X, Y, Z  (all populated voxels)
    NumericVector scanner_location,
    double voxel_size,
    int ncores) {
  
  int n = analysis_grid.nrow();
  int m = obstacles.nrow();
  
  // --- Copy Rcpp data into plain C++ containers ---
  std::vector<double> grid_x(n), grid_y(n), grid_z(n), grid_v(n);
  for (int i = 0; i < n; ++i) {
    grid_x[i] = analysis_grid(i, 0);
    grid_y[i] = analysis_grid(i, 1);
    grid_z[i] = analysis_grid(i, 2);
    grid_v[i] = analysis_grid(i, 3);
  }
  double scan_x = scanner_location[0];
  double scan_y = scanner_location[1];
  double scan_z = scanner_location[2];
  std::vector<int> res(n, 0);
  
  // --- Compute bounding box over all coordinates ---
  int min_ix = static_cast<int>(floor(scan_x / voxel_size));
  int max_ix = min_ix;
  int min_iy = static_cast<int>(floor(scan_y / voxel_size));
  int max_iy = min_iy;
  int min_iz = static_cast<int>(floor(scan_z / voxel_size));
  int max_iz = min_iz;
  
  for (int i = 0; i < m; ++i) {
    int ix = static_cast<int>(floor(obstacles(i, 0) / voxel_size));
    int iy = static_cast<int>(floor(obstacles(i, 1) / voxel_size));
    int iz = static_cast<int>(floor(obstacles(i, 2) / voxel_size));
    if (ix < min_ix) min_ix = ix; if (ix > max_ix) max_ix = ix;
    if (iy < min_iy) min_iy = iy; if (iy > max_iy) max_iy = iy;
    if (iz < min_iz) min_iz = iz; if (iz > max_iz) max_iz = iz;
  }
  for (int i = 0; i < n; ++i) {
    int ix = static_cast<int>(floor(grid_x[i] / voxel_size));
    int iy = static_cast<int>(floor(grid_y[i] / voxel_size));
    int iz = static_cast<int>(floor(grid_z[i] / voxel_size));
    if (ix < min_ix) min_ix = ix; if (ix > max_ix) max_ix = ix;
    if (iy < min_iy) min_iy = iy; if (iy > max_iy) max_iy = iy;
    if (iz < min_iz) min_iz = iz; if (iz > max_iz) max_iz = iz;
  }
  
  // Extra margin for fat-ray neighbor checks (2 voxels)
  min_ix -= 2; min_iy -= 2; min_iz -= 2;
  max_ix += 2; max_iy += 2; max_iz += 2;
  
  // --- Build flat obstacle grid ---
  VoxelGrid3D obstacle_grid;
  obstacle_grid.init(min_ix, min_iy, min_iz, max_ix, max_iy, max_iz);
  for (int i = 0; i < m; ++i) {
    int ix = static_cast<int>(floor(obstacles(i, 0) / voxel_size));
    int iy = static_cast<int>(floor(obstacles(i, 1) / voxel_size));
    int iz = static_cast<int>(floor(obstacles(i, 2) / voxel_size));
    obstacle_grid.set(ix, iy, iz);
  }
  
  size_t grid_bytes = (size_t)obstacle_grid.nx * obstacle_grid.ny * obstacle_grid.nz;
  Rcout << "       Obstacle grid: " << obstacle_grid.nx << " x "
        << obstacle_grid.ny << " x " << obstacle_grid.nz
        << " (" << (grid_bytes / (1024*1024)) << " MB)\n";
  Rcout.flush();
  
  omp_set_num_threads(ncores > 1 ? ncores : 1);
  
  // --- Backward fat-ray DDA for each analysis voxel ---
#pragma omp parallel for schedule(dynamic, 1024) if(ncores > 1)
  for (int i = 0; i < n; ++i) {
    if (grid_v[i] > 0) {
      res[i] = 1;  // Populated: scanner got a return here
    } else {
      bool blocked = is_ray_occluded_dda_fat(
        scan_x, scan_y, scan_z,
        grid_x[i], grid_y[i], grid_z[i],
        voxel_size, obstacle_grid);
      res[i] = blocked ? 2 : 0;
    }
  }
  
  IntegerVector result(n);
  for (int i = 0; i < n; ++i) result[i] = res[i];
  return result;
}