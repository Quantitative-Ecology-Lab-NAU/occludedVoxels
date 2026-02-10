#include <Rcpp.h>
#include <omp.h>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <atomic>
#include <cmath>
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
// DDA (Amanatides-Woo) 3D ray traversal
// Visits EVERY voxel the ray passes through - no gaps possible
// More accurate and faster than fixed-step marching
// ============================================================
inline bool is_ray_occluded_dda(
    double start_x, double start_y, double start_z,
    double end_x, double end_y, double end_z,
    double voxel_size,
    const std::unordered_set<VoxelKey, VoxelHash>& obstacles) {
  
  // Current voxel indices (scanner position)
  int ix = static_cast<int>(floor(start_x / voxel_size));
  int iy = static_cast<int>(floor(start_y / voxel_size));
  int iz = static_cast<int>(floor(start_z / voxel_size));
  
  // Target voxel indices
  int ix_end = static_cast<int>(floor(end_x / voxel_size));
  int iy_end = static_cast<int>(floor(end_y / voxel_size));
  int iz_end = static_cast<int>(floor(end_z / voxel_size));
  
  // If start and end are in the same voxel, nothing can block the ray
  if (ix == ix_end && iy == iy_end && iz == iz_end) return false;
  
  // Ray direction
  double dx = end_x - start_x;
  double dy = end_y - start_y;
  double dz = end_z - start_z;
  
  // Step direction (+1 or -1 along each axis)
  int step_x = (dx >= 0) ? 1 : -1;
  int step_y = (dy >= 0) ? 1 : -1;
  int step_z = (dz >= 0) ? 1 : -1;
  
  // Parameter t increment for one full voxel step in each axis
  double tDeltaX, tDeltaY, tDeltaZ;
  double tMaxX, tMaxY, tMaxZ;
  
  if (dx != 0.0) {
    tDeltaX = fabs(voxel_size / dx);
    tMaxX = (dx >= 0) ? ((ix + 1) * voxel_size - start_x) / dx
                       : (ix * voxel_size - start_x) / dx;
  } else {
    tDeltaX = 1e30;
    tMaxX = 1e30;
  }
  
  if (dy != 0.0) {
    tDeltaY = fabs(voxel_size / dy);
    tMaxY = (dy >= 0) ? ((iy + 1) * voxel_size - start_y) / dy
                       : (iy * voxel_size - start_y) / dy;
  } else {
    tDeltaY = 1e30;
    tMaxY = 1e30;
  }
  
  if (dz != 0.0) {
    tDeltaZ = fabs(voxel_size / dz);
    tMaxZ = (dz >= 0) ? ((iz + 1) * voxel_size - start_z) / dz
                       : (iz * voxel_size - start_z) / dz;
  } else {
    tDeltaZ = 1e30;
    tMaxZ = 1e30;
  }
  
  // Maximum traversal steps (Manhattan distance + safety margin)
  int max_steps = abs(ix_end - ix) + abs(iy_end - iy) + abs(iz_end - iz) + 2;
  
  for (int step = 0; step < max_steps; ++step) {
    // Advance to next voxel boundary (Amanatides-Woo algorithm)
    if (tMaxX < tMaxY) {
      if (tMaxX < tMaxZ) {
        ix += step_x;
        tMaxX += tDeltaX;
      } else {
        iz += step_z;
        tMaxZ += tDeltaZ;
      }
    } else {
      if (tMaxY < tMaxZ) {
        iy += step_y;
        tMaxY += tDeltaY;
      } else {
        iz += step_z;
        tMaxZ += tDeltaZ;
      }
    }
    
    // Stop when we reach the target voxel (don't check target itself)
    if (ix == ix_end && iy == iy_end && iz == iz_end) break;
    
    // Check if this voxel contains an obstacle (populated vegetation)
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
// ACCURATE occlusion check with separate obstacle map
// Uses Amanatides-Woo DDA ray tracing through ALL populated
// 1cm voxels for exact voxel-level classification.
// 
// Key differences from check_occlusion():
//   1. Obstacles are separate from the analysis grid - rays are
//      blocked by ALL populated voxels (including transect vegetation
//      between scanner and quadrat), not just those in the AOI
//   2. Uses true DDA traversal that visits every voxel the ray
//      passes through, with zero gaps
//   3. Analysis at 1cm = no scaling/multiplication needed
// ============================================================
// [[Rcpp::export]]
IntegerVector check_occlusion_accurate(
    NumericMatrix analysis_grid,    // N x 4: X, Y, Z, V1(0/1) - voxels to classify
    NumericMatrix obstacles,        // M x 3: X, Y, Z - ALL populated voxels that block rays
    NumericVector scanner_location, 
    double voxel_size,              // voxel size (same for analysis and obstacles)
    int ncores) {
  
  int n = analysis_grid.nrow();
  int m = obstacles.nrow();
  
  // --- Copy Rcpp data into plain C++ containers (R-API-free) ---
  // This prevents R's stack-limit check from firing inside OpenMP threads
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
  std::vector<int> res(n, 0);  // plain C++ result vector
  
  // Build obstacle hash set from ALL populated voxels (entire clip area)
  // This includes vegetation between scanner and quadrat that blocks rays
  std::unordered_set<VoxelKey, VoxelHash> obstacle_set;
  obstacle_set.reserve(m);
  
  for (int i = 0; i < m; ++i) {
    auto key = get_voxel_key(obstacles(i, 0), obstacles(i, 1), obstacles(i, 2), voxel_size);
    obstacle_set.insert(key);
  }
  
  omp_set_num_threads(ncores > 1 ? ncores : 1);
  
  // --- OpenMP region: only plain C++ data is touched ---
  // NO R API calls (Rprintf, Rcout, etc.) inside the parallel region
  // to avoid "C stack usage too close to the limit" on OpenMP thread stacks
#pragma omp parallel for schedule(dynamic, 1024) if(ncores > 1)
  for (int i = 0; i < n; ++i) {
    if (grid_v[i] > 0) {
      res[i] = 1;  // Populated (has actual point cloud data)
    } else {
      // DDA ray trace from scanner through ALL vegetation obstacles
      bool blocked = is_ray_occluded_dda(
        scan_x, scan_y, scan_z,
        grid_x[i], grid_y[i], grid_z[i],
        voxel_size, obstacle_set);
      res[i] = blocked ? 2 : 0;  // 2 = occluded, 0 = empty visible
    }
  }
  
  // --- Copy results back to Rcpp IntegerVector ---
  IntegerVector result(n);
  for (int i = 0; i < n; ++i) result[i] = res[i];
  return result;
}