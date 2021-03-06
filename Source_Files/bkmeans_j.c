#include "../Header_Files/Headers.h"
#include "../Header_Files/Defs.h"
#include "../Header_Files/bkmeans_j.h"
#include "../Header_Files/helping_procedures.h"

int bisecting_kmeans(int dim, int ndata, double *data, int k,
                     int *cluster_size, int *cluster_start,
                     double *cluster_radius, double **cluster_centroid,
                     int *cluster_assign)
{
  int i, num_iterations = 0, curr_cluster_count = 1,
      cluster_x = 0, cluster_y = 1, cluster_with_min_SSE = -1;

  double *cluster_sse = malloc(k * sizeof(double));
  // Initialize cluster_sse
  for(i = 0; i < k; i++) {
    cluster_sse[i] = 0.0;
  }

  // Assign all points to cluster_x; initially, cluster_x == 0
  for(i = 0; i < ndata; i++) {
    cluster_assign[i] = cluster_x;
  }

  // Bisect cluster_x -> cluster with max SSE
  // Initially, there is only one cluster, so cluster_x is still cluster with max SSE
  cluster_size[cluster_x] = ndata;

  while(curr_cluster_count != k) {
    set_new_cluster_centroids(dim, ndata, data, cluster_x, cluster_y,
                              cluster_centroid, cluster_assign, cluster_size);

    num_iterations += two_kmeans(dim, ndata, data, cluster_x, cluster_y,
                                 cluster_size, cluster_centroid, cluster_assign);

    cluster_with_min_SSE = get_max_SSE(dim, ndata, data, k, cluster_centroid,
                                       cluster_sse, cluster_assign);

    cluster_x = cluster_with_min_SSE;
    cluster_y = cluster_y + 1;

    curr_cluster_count++;
  }

  num_iterations += kmeans_bkm(dim, ndata, data, k,
                               cluster_size, cluster_start, cluster_radius,
                               cluster_centroid, cluster_assign);

  return num_iterations;
}

int two_kmeans(int dim, int ndata, double *data, int cluster_x, int cluster_y,
               int *cluster_size, double **cluster_centroid, int *cluster_assign)
{
  int i, num_iterations = 0, thresh_hold = 2;
  bool notFinishedClustering = true;

  int *prev_cluster_assign = malloc(ndata * sizeof(int));
  for(i = 0; i < ndata; i++) {
    prev_cluster_assign[i] = -1;
  }

  while(notFinishedClustering && num_iterations < thresh_hold) {
    num_iterations++;

    // Assign data points to cluster_x or cluster_y
    for(i = 0; i < ndata; i++) {
      // Only check points that are in cluster_x -> cluster with max SSE
      if(cluster_assign[i] == cluster_x || cluster_assign[i] == cluster_y) {
        two_assign_pt_to_cluster(dim, cluster_x, cluster_y, i, data, cluster_centroid,
                                 cluster_assign, prev_cluster_assign);
      }
    }

    // Reset cluster sizes because of new data point cluster assignments
    cluster_size[cluster_x] = 0; cluster_size[cluster_y] = 0;

    // Set cluster sizes for cluster_x and cluster_y
    for(i = 0; i < ndata; i++) {
      if(cluster_assign[i] == cluster_x) {
        cluster_size[cluster_x]++;
      }

      if(cluster_assign[i] == cluster_y) {
        cluster_size[cluster_y]++;
      }
    }

    if(! cluster_assignments_changed_bkm(ndata, cluster_assign, prev_cluster_assign)) {
      notFinishedClustering = false;
    }
    else {
      if(num_iterations < thresh_hold) {
        set_prev_cluster_assignments_bkm(ndata, cluster_assign, prev_cluster_assign);
      }

      two_update_cluster_centroids(dim, cluster_x, cluster_y, ndata, data, cluster_centroid,
                                   cluster_size, cluster_assign);
    }
  }

  return num_iterations;
}

int kmeans_bkm(int dim, int ndata, double *data, int k,
               int *cluster_size, int *cluster_start,
               double *cluster_radius, double **cluster_centroid,
               int *cluster_assign)
{
  int i, j, num_iterations = 0, thresh_hold = 1, *prev_cluster_assign = NULL;

  if(thresh_hold > 1) {
    prev_cluster_assign = malloc(ndata * sizeof(double));
    for (i = 0; i < ndata; i++) {
      prev_cluster_assign[i] = -1;
    }
  }

  while (num_iterations < thresh_hold) {
    num_iterations++;

    // Assign data point to clusters and update centroids
    for (i = 0; i < ndata; i++) {
      assign_pt_to_cluster_bkm(dim, k, i, data, cluster_centroid, cluster_assign);
    }

    // Reset cluster sizes because of new data point cluster assignments
    for (i = 0; i < k; i++) {
      cluster_size[i] = 0;
    }

    // Set cluster sizes
    for (i = 0; i < ndata; i++) {
      for (j = 0; j < k; j++) {
        if (cluster_assign[i] == j) {
          cluster_size[j]++;
        }
      }
    }

    if(thresh_hold > 1) {
      set_prev_cluster_assignments_bkm(ndata, cluster_assign, prev_cluster_assign);
    }

    update_cluster_centroids_bkm(dim, k, ndata, data, cluster_centroid,
                                 cluster_size, cluster_assign);
  }

  printf("\nSorting the data...\n");
  quick_sort_data_bkm(dim, 0, ndata, data, cluster_assign);

  set_cluster_start_bkm(k, cluster_size, cluster_start);

  for (i = 0; i < k; i++) {
    set_cluster_radius_bkm(dim, i, cluster_size, cluster_start, data,
                           cluster_centroid, cluster_radius);
  }

  return num_iterations;
}

void set_new_cluster_centroids(int dim, int ndata, double *data, int cluster_x, int cluster_y,
                               double **cluster_centroid, int *cluster_assign, int *cluster_size)
{
  int i, j, centroid_pt_x = -1, centroid_pt_y = -1;
  double distance, maxDistance = (double) INT_MIN;

  int *cluster_x_pts = malloc(cluster_size[cluster_x] * sizeof(int));

  j = 0;
  for(i = 0; i < ndata; i++) {
    if(cluster_assign[i] == cluster_x) {
      cluster_x_pts[j] = i;
      j++;
    }
  }

  // centroid_pt_x is randomly chosen from the data set of cluster_x
  centroid_pt_x = cluster_x_pts[cluster_size[cluster_x] + (rand() / (RAND_MAX / (0 - cluster_size[cluster_x])))];

  for(j = 0; j < dim; j++) {
    cluster_centroid[cluster_x][j] = data[centroid_pt_x * dim + j];
  }

  // Second centroid is the furthest away point from centroid_pt_x in cluster_x
  for(i = 0; i < ndata; i++) {
    if(cluster_assign[i] == cluster_x) {
      distance = 0.0;
      for(j = 0; j < dim; j++) {
        distance += (cluster_centroid[cluster_x][j] - data[i * dim + j]) *
                    (cluster_centroid[cluster_x][j] - data[i * dim + j]);
      }

      distance = sqrt(distance);
      if(distance > maxDistance) {
        maxDistance = distance;
        centroid_pt_y = i;
      }
    }
  }

  for(j = 0; j < dim; j++) {
    cluster_centroid[cluster_y][j] = data[centroid_pt_y * dim + j];
  }

  free(cluster_x_pts);
}

int get_max_SSE(int dim, int ndata, double *data, int k, double **cluster_centroid,
                double *cluster_sse, int *cluster_assign)
{
  int i, cluster_with_max_SSE = -1;
  double max_SSE = (double) INT_MIN;

  // Calculate SSE for each cluster
  for(i = 0; i < k; i++) {
    cluster_sse[i] = calc_SSE(dim, ndata, data, i, cluster_centroid, cluster_assign);

    if(cluster_sse[i] > max_SSE) {
      max_SSE = cluster_sse[i];
      cluster_with_max_SSE = i;
    }
  }

  return cluster_with_max_SSE;
}


double calc_SSE(int dim, int ndata, double *data, int curr_cluster,
                double **cluster_centroid, int *cluster_assign)
{
  int i, j;
  double sum_squares_error = 0.0;

  for(i = 0; i < ndata; i++) {
    if(cluster_assign[i] == curr_cluster) {
      for(j = 0; j < dim; j++) {
        sum_squares_error += (data[i * dim + j] - cluster_centroid[curr_cluster][j]) *
                             (data[i * dim + j] - cluster_centroid[curr_cluster][j]);
      }
    }
  }

  return sum_squares_error;
}

void assign_pt_to_cluster_bkm(int dim, int totalClusters, int data_pt, double *data,
                              double **cluster_centroid, int *cluster_assign)
{
  int i, j, closestCluster = -1;
  double distance, minDistToCluster = (double) INT_MAX;
  for (i = 0; i < totalClusters; i++) {
    distance = 0.0;
    for (j = 0; j < dim; j++) {
      distance += (cluster_centroid[i][j] - data[data_pt * dim + j]) *
                  (cluster_centroid[i][j] - data[data_pt * dim + j]);
    }

    distance = sqrt(distance);
    if (distance < minDistToCluster) {
      minDistToCluster = distance;
      closestCluster = i;
    }
  }

  cluster_assign[data_pt] = closestCluster;
}

void two_assign_pt_to_cluster(int dim, int cluster_x, int cluster_y, int data_pt, double *data,
                             double **cluster_centroid, int *cluster_assign,
                             int *prev_cluster_assign)
{
  int j;
  double distance_to_cluster_x = 0.0, distance_to_cluster_y = 0.0;

  if(cluster_assign[data_pt] != prev_cluster_assign[data_pt]) {
    for(j = 0; j < dim; j++) {
      distance_to_cluster_x += (cluster_centroid[cluster_x][j] - data[data_pt * dim + j]) *
                               (cluster_centroid[cluster_x][j] - data[data_pt * dim + j]);
    }
    distance_to_cluster_x = sqrt(distance_to_cluster_x);

    for(j = 0; j < dim; j++) {
      distance_to_cluster_y += (cluster_centroid[cluster_y][j] - data[data_pt * dim + j]) *
                               (cluster_centroid[cluster_y][j] - data[data_pt * dim + j]);
    }
    distance_to_cluster_y = sqrt(distance_to_cluster_y);

    if(distance_to_cluster_x < distance_to_cluster_y) {
      cluster_assign[data_pt] = cluster_x;
    }
    else {
      cluster_assign[data_pt] = cluster_y;
    }
  }
}

void two_update_cluster_centroids(int dim, int cluster_x, int cluster_y, int ndata, double *data,
                                  double **cluster_centroid, int *cluster_size,
                                  int *cluster_assign)
{
  int j;

  // Only calculate when cluster_x size is > 0
  if(cluster_size[cluster_x] > 0) {
    for(j = 0; j < dim; j++) {
      cluster_centroid[cluster_x][j] = calc_centroid_bkm(dim, j, cluster_x, ndata, data,
                                                         cluster_size, cluster_assign);
    }
  }

  // Only calculate when cluster_y size is > 0
  if(cluster_size[cluster_y] > 0) {
    for(j = 0; j < dim; j++) {
      cluster_centroid[cluster_y][j] = calc_centroid_bkm(dim, j, cluster_y, ndata, data,
                                                         cluster_size, cluster_assign);
    }
  }
}

void update_cluster_centroids_bkm(int dim, int totalClusters, int ndata, double *data,
                                  double **cluster_centroid, int *cluster_size,
                                  int *cluster_assign)
{
  int i, j;
  for (i = 0; i < totalClusters; i++) {
    if (cluster_size[i] > 0) { // Only calculate when cluster size is > 0
      for (j = 0; j < dim; j++) {
        cluster_centroid[i][j] = calc_centroid_bkm(dim, j, i, ndata, data,
                                                  cluster_size, cluster_assign);
      }
    }
  }
}

double calc_centroid_bkm(int dim, int currDim, int currCluster, int ndata,
                         double *data, int *cluster_size, int *cluster_assign)
{
  int i;
  double sum = 0.0;

  for(i = 0; i < ndata; i++) {
    if(cluster_assign[i] == currCluster) {
      sum += data[i * dim + currDim];
    }
  }

  return sum / cluster_size[currCluster];
}

bool cluster_assignments_changed_bkm(int ndata, int *cluster_assign, int *prev_cluster_assign)
{
  int i;

  for(i = 0; i < ndata; i++) {
    if(cluster_assign[i] != prev_cluster_assign[i]) {
      return true;
    }
  }

  return false;
}

void set_prev_cluster_assignments_bkm(int ndata, int *cluster_assign, int *prev_cluster_assign)
{
  int i;

  for(i = 0; i < ndata; i++) {
    prev_cluster_assign[i] = cluster_assign[i];
  }
}

void quick_sort_data_bkm(int dim, int lo, int hi, double *data, int *cluster_assign)
{
  if(lo < hi) {
    int pivot = partition_bkm(dim, lo, hi, data, cluster_assign);
    quick_sort_data_bkm(dim, lo, pivot, data, cluster_assign);
    quick_sort_data_bkm(dim, pivot+1, hi, data, cluster_assign);
  }
}

int partition_bkm(int dim, int lo, int hi, double *data, int *cluster_assign)
{
  int p = cluster_assign[lo];
  int i = lo - 1;
  int j = hi + 1;

  while(true) {
    do {
      i++;
    } while(cluster_assign[i] < p && i < hi);

    do {
      j--;
    } while(cluster_assign[j] > p && j > lo);

    if(i >= j) {
      return j;
    }

    swap_cluster_assign_bkm(cluster_assign, i, j);
    swap_points_bkm(dim, i, j, data);

    //swap_labels_bkm(labels, i, j);
  }
}

void swap_cluster_assign_bkm(int *cluster_assign, int i, int j)
{
  int tmp;
  tmp = cluster_assign[j];
  cluster_assign[j] = cluster_assign[i];
  cluster_assign[i] = tmp;
}

void swap_points_bkm(int dim, int point1, int point2, double *data)
{
  int i, point1Index = point1 * dim, point2Index = point2 * dim;
  double temp;

  for(i = point1Index; i < point1Index + dim; i++) {
    temp = data[point2Index];
    data[point2Index] = data[i];
    data[i] = temp;

    point2Index++;
  }
}

void swap_labels_bkm(int *labels, int label1, int label2)
{
  int tmp;
  tmp = labels[label2];
  labels[label2] = labels[label1];
  labels[label1] = tmp;
}

void set_cluster_start_bkm(int totalClusters, int *cluster_size, int *cluster_start)
{
  int i, nextStartingPoint = cluster_size[0];

  cluster_start[0] = 0;
  for(i = 1; i < totalClusters; i++) {
    cluster_start[i] = nextStartingPoint;
    nextStartingPoint += cluster_size[i];
  }
}

void set_cluster_radius_bkm(int dim, int currCluster, int *cluster_size,
                           int *cluster_start, double *data,
                           double **cluster_centroid, double *cluster_radius)
{
  int i, j, end = cluster_size[currCluster] + cluster_start[currCluster];
  double distance, maxDistance = (double) INT_MIN;

  for(i = cluster_start[currCluster]; i < end; i++) {
    distance = 0.0;
    for(j = 0; j < dim; j++) {
      distance += (cluster_centroid[currCluster][j] - data[i * dim + j]) *
                  (cluster_centroid[currCluster][j] - data[i * dim + j]);
    }

    distance = sqrt(distance);
    if(distance > maxDistance) {
      maxDistance = distance;
    }
  }

  cluster_radius[currCluster] = maxDistance;
}

void bkmeans_search_clsuters_for_approx_neighbors(int dim, int test_size, int k,
                                                  double *train_feature_data, double *test_feature_data,
                                                  int *train_non_feature_data, int *test_non_feature_data,
                                                  int *cluster_size, int *cluster_start,
                                                  double *cluster_radius, double **cluster_centroid)
{
  int i, j, a, b, c = 0, closest_cluster = -1;
  double distance, min_cluster_distance = (double) INT_MAX, closest_neighbor_dist = 0.0,
         total_closest_neighbor_dist = 0.0, *query = malloc(dim * sizeof(double)),
         pts_searched_in_closest_cluster = 0.0, total_pts_searched = 0.0;;

  double *cluster_distances = malloc(k * sizeof(double));
  for(i = 0; i < k; i++) {
    cluster_distances[i] = 0.0;
  }

  for(a = 0; a < test_size; a++) {
    for(b = a * dim; b < a * dim + dim; b++) {
      query[c] = test_feature_data[b];
      c++;
    }
    c = 0;

    // Find closest cluster to the query point
    for(i = 0; i < k; i++) {
      distance = 0.0;
      for(j = 0; j < dim; j++) {
        distance += (query[j] - cluster_centroid[i][j]) *
                    (query[j] - cluster_centroid[i][j]);
      }
      cluster_distances[i] = sqrt(distance);

      if(cluster_distances[i] < min_cluster_distance) {
        min_cluster_distance = cluster_distances[i]; closest_cluster = i;
      }
    }

    // Search points in the closest cluster
    closest_neighbor_dist = search_points_in_cluster_bkm(dim, query, train_feature_data,
                                                             closest_cluster, cluster_start, cluster_size,
                                                             &pts_searched_in_closest_cluster);

    total_closest_neighbor_dist += closest_neighbor_dist;
    total_pts_searched += pts_searched_in_closest_cluster;

    pts_searched_in_closest_cluster = 0.0; closest_neighbor_dist = 0.0;
  }

  free(cluster_distances); free(query);

  print_search_results(test_size, total_closest_neighbor_dist, total_pts_searched);
}

double search_points_in_cluster_bkm(int dim, double *query, double *train_feature_data,
                                    int closest_cluster, int *cluster_start, int *cluster_size,
                                    double *pts_searched_in_closest_cluster)
{
  int i, j, end = cluster_size[closest_cluster] + cluster_start[closest_cluster];
  double distance, min_point_distance = (double) INT_MAX;

  for(i = cluster_start[closest_cluster]; i < end; i++) {
    distance = 0.0;
    for(j = 0; j < dim; j++) {
      distance += (query[j] - train_feature_data[i * dim + j]) * (query[j] - train_feature_data[i * dim + j]);
    }
    distance = sqrt(distance);

    if(distance < min_point_distance) {
      min_point_distance = distance;
    }

    *pts_searched_in_closest_cluster = *pts_searched_in_closest_cluster + 1;
  }

  return min_point_distance;
}