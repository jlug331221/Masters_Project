#include "../Header_Files/Headers.h"
#include "../Header_Files/Defs.h"
#include "../Header_Files/kdtree.h"
#include "../Header_Files/helping_procedures.h"

int kdtree(int dim, int ndata, double *data, int k,
           int *cluster_size, int *cluster_start, double **cluster_bdry,
           double **cluster_centroid, int *cluster_assign)
{
  int treeDepth = 0;
  int currDepth = 0, dLeft = 0, dRight = ndata, kLeft = 0, kRight = k;

  treeDepth = (int) floor(log2(k));

  kdtree_helper(dim, ndata, data, treeDepth, currDepth, dLeft, dRight,
                kLeft, kRight, cluster_size, cluster_start, cluster_bdry,
                cluster_centroid, cluster_assign);

  return 0;
}

void kdtree_helper(int dim, int ndata, double *data, int treeDepth,
                   int currDepth, int dLeft, int dRight, int kLeft,
                   int kRight, int *cluster_size, int *cluster_start,
                   double **cluster_bdry, double **cluster_centroid,
                   int *cluster_assign)
{
  int midCluster, midData;

  if(currDepth != treeDepth) {
    bipartition(dim, dLeft, dRight, data, &cluster_size[kLeft], &cluster_start[kLeft],
                &cluster_bdry[kLeft], &cluster_centroid[kLeft],
                cluster_assign);
  }

  if(currDepth >= treeDepth) {
    // Assign clusters
    int i;
    for(i = dLeft; i < dRight; i++) {
      cluster_assign[i] = kLeft;
    }

    return;
  }

  midCluster = (kLeft + kRight) / 2;
  midData = cluster_start[kLeft] + cluster_size[kLeft]; //midData = cluster_start[kLeft];

  // Recursive call on left
  kdtree_helper(dim, ndata, data, treeDepth, currDepth + 1, dLeft, midData,
                kLeft, midCluster, cluster_size, cluster_start, cluster_bdry,
                cluster_centroid, cluster_assign);

  // Recursive call on right
  kdtree_helper(dim, ndata, data, treeDepth, currDepth + 1, midData, dRight,
                midCluster, kRight, cluster_size, cluster_start, cluster_bdry,
                cluster_centroid, cluster_assign);
}

int bipartition(int dim, int io, int im, double *data,
                int cluster_size[2], int cluster_start[2],
                double *cluster_bdry[2], double *cluster_centroid[2],
                int *cluster_assign)
{
  int i, j, partitionDim = -1;
  double mean = 0.0, maxVar = 0.0, currVar = 0.0, partitionMean = 0.0;

  // Find means of each dimension and determine the dimension with max variance
  for(i = 0; i < dim; i++) {
    mean = calc_mean(i, dim, io * dim, im * dim, data);

    currVar = calc_variance(i, dim, io * dim, im * dim, data, mean);

    if(currVar > maxVar) {
      maxVar = currVar;
      partitionDim = i;
      partitionMean = mean;
    }
  }

  // Arrange data around partition mean (**Hoare Partition**)
  // Both 'i' and 'j' are points in the data set and NOT indices.
  i = io; // start 'i' at beginning point
  j = im - 1; // start j at end - 1 because 'im' is non-inclusive
  int leftPartitionIndex = i * dim + partitionDim,
      rightPartitionIndex = im * dim - dim + partitionDim;
  while(1) {
    while(data[leftPartitionIndex] < partitionMean && i < im) {
      i++;
      leftPartitionIndex += dim;
    }

    while(data[rightPartitionIndex] > partitionMean && j > io) {
      j--;
      rightPartitionIndex -= dim;
    }

    if(i >= j) { break; }

    kdtree_swap_points(dim, data, i, j);
    //kdtree_swap_labels(labels, i, j);
  }

  // Set cluster start
  cluster_start[0] = io;
  cluster_start[1] = j+1;

  // Set cluster size
  cluster_size[0] = cluster_start[1] - cluster_start[0];
  cluster_size[1] = im - cluster_start[1];

  int bdryMinMaxIndex = 0;
  // Find boundaries of each cluster
  for(i = 0; i < dim; i++) {
    cluster_bdry[0][bdryMinMaxIndex] = find_boundary_min(i, dim, data,
                                                         cluster_start[0] * dim,
                                                         (cluster_start[0] + cluster_size[0]) * dim);
    cluster_bdry[0][bdryMinMaxIndex + 1] = find_boundary_max(i, dim, data,
                                                             cluster_start[0] * dim,
                                                             (cluster_start[0] + cluster_size[0]) * dim);
    cluster_bdry[1][bdryMinMaxIndex] = find_boundary_min(i, dim, data,
                                                         cluster_start[1] * dim,
                                                         (cluster_start[1] + cluster_size[1]) * dim);
    cluster_bdry[1][bdryMinMaxIndex + 1] = find_boundary_max(i, dim, data,
                                                             cluster_start[1]*dim,
                                                             (cluster_start[1] + cluster_size[1]) * dim);

    bdryMinMaxIndex += 2;
  }

  return 0;
}

void kdtree_search_clusters_for_approx_neighbors(int dim, int test_size, int k,
                                                 double *train_feature_data, double *test_feature_data,
                                                 int *train_non_feature_data, int *test_non_feature_data,
                                                 int *cluster_size, int *cluster_start, double **cluster_bdry)
{
  int i, j, a, b, c = 0, dim_min_max_index, closest_cluster = -1;
  double distance, min_cluster_distance = (double) INT_MAX, closest_neighbor_dist = 0.0,
         total_closest_neighbor_dist = 0.0, *query = malloc(dim * sizeof(double)),
         pts_searched_in_closest_cluster = 0.0, total_pts_searched = 0.0;

  double *cluster_distances = malloc(k * sizeof(double));
  for(i = 0; i < k; i++) {
    cluster_distances[i] = -1.0;
  }

  for(a = 0; a < test_size; a++) {
    for(b = a * dim; b < a * dim + dim; b++) {
      query[c] = test_feature_data[b];
      c++;
    }
    c = 0;

    // Begin searching for the nearest cluster (cluster with shortest distance to query)
    for(i = 0; i < k; i++) {
      distance = 0.0;
      dim_min_max_index = 0;
      for(j = 0; j < dim; j++) {
        if(query[j] < cluster_bdry[i][dim_min_max_index]) {
          distance += (query[j] - cluster_bdry[i][dim_min_max_index]) *
                      (query[j] - cluster_bdry[i][dim_min_max_index]);
        }
        if(query[j] > cluster_bdry[i][dim_min_max_index + 1]) {
          distance += (query[j] - cluster_bdry[i][dim_min_max_index + 1]) *
                      (query[j] - cluster_bdry[i][dim_min_max_index + 1]);
        }
        dim_min_max_index += 2;
      }
      cluster_distances[i] = sqrt(distance);

      if(cluster_distances[i] < min_cluster_distance) {
        min_cluster_distance = cluster_distances[i]; closest_cluster = i;
      }
    }

    // Calculate distances for each point in closest_cluster
    closest_neighbor_dist = kdtree_search_points_in_cluster(dim, query, train_feature_data, closest_cluster,
                                                                cluster_start, cluster_size,
                                                                &pts_searched_in_closest_cluster);

    total_closest_neighbor_dist += closest_neighbor_dist;
    total_pts_searched += pts_searched_in_closest_cluster;

    pts_searched_in_closest_cluster = 0.0; closest_neighbor_dist = 0.0;
  }

  free(cluster_distances); free(query);

  print_search_results(test_size, total_closest_neighbor_dist, total_pts_searched);
}

double kdtree_search_points_in_cluster(int dim, double *query, double *train_feature_data, int closest_cluster,
                                       int *cluster_start, int *cluster_size,
                                       double *pts_searched_in_closest_cluster)
{
  int i, j, dataPointIndex, cluster_end = cluster_start[closest_cluster] + cluster_size[closest_cluster];
  double distance, min_point_distance = (double) INT_MAX;

  for(i = cluster_start[closest_cluster]; i < cluster_end; i++) {
    distance = 0.0;
    for(j = 0; j < dim; j++) {
      dataPointIndex = i * dim + j;
      if(query[j] < train_feature_data[dataPointIndex]) {
        distance += (query[j] - train_feature_data[dataPointIndex]) *
                    (query[j] - train_feature_data[dataPointIndex]);
      }

      if(query[j] > train_feature_data[dataPointIndex]) {
        distance += (query[j] - train_feature_data[dataPointIndex]) *
                    (query[j] - train_feature_data[dataPointIndex]);
      }
    }

    distance = sqrt(distance);
    if(distance < min_point_distance) {
      min_point_distance = distance;
    }

    *pts_searched_in_closest_cluster = *pts_searched_in_closest_cluster + 1;
  }

  return min_point_distance;
}

void kdtree_swap_points(int dim, double *data, int point1, int point2)
{
  int i;
  double tmp;
  int point1Index = point1 * dim, point2Index = point2 * dim;
  for(i = point1Index; i < point1Index + dim; i++) {
    tmp = data[point2Index];
    data[point2Index] = data[i];
    data[i] = tmp;

    point2Index++;
  }
}

void kdtree_swap_labels(int *labels, int point1, int point2)
{
  int tmp;
  tmp = labels[point2];
  labels[point2] = labels[point1];
  labels[point1] = tmp;
}

double find_boundary_min(int currDim, int totalDim, double *data, int start, int end)
{
  int i;
  double min = (double) INT_MAX;
  for(i = start+currDim; i < end; i+=totalDim) {
    if(data[i] < min) {
      min = data[i];
    }
  }

  return min;
}

double find_boundary_max(int currDim, int totalDim, double *data, int start, int end)
{
  int i;
  double max = (double) INT_MIN;
  for(i = start+currDim; i < end; i+=totalDim) {
    if(data[i] > max) {
      max = data[i];
    }
  }

  return max;
}

double calc_mean(int currDim, int totalDim, int start, int end, double *data)
{
  int i, ndata = 0;
  double sum = 0.0;
  for(i = currDim+start; i < end; i+=totalDim) {
    sum += data[i];
    ndata++;
  }

  return sum / ndata;
}

double calc_variance(int currDim, int totalDim, int start, int end, double *data,
                     double dimMean)
{
  int i = 0, ndata = 0;
  double sumOfSquares = 0.0;

  // First, store sum of squares for each data item for variance calculation
  for(i = currDim+start; i < end; i+=totalDim) {
    sumOfSquares += data[i] * data[i];
    ndata++;
  }

  // Then calculate and return the variance
  return (sumOfSquares / ndata) - (dimMean * dimMean);
}
