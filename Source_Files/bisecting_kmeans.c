#include "../Header_Files/Headers.h"
#include "../Header_Files/bisecting_kmeans.h"

int bisecting_kmeans(int dim, int ndata, double *data, int k,
                     int *cluster_size, int *cluster_start,
                     double *cluster_radius, double **cluster_centroid,
                     int *cluster_assign)
{
  int i, num_iterations = 0, curr_cluster_count = 1,
      cluster_x = 0, cluster_y = 1, cluster_with_max_SSE = -1;

  double *cluster_sse = malloc(k * sizeof(double));
  // Initialize cluster_sse
  for(i = 0; i < k; i++) {
    cluster_sse[i] = 0.0;
  }

  // Assign all points to cluster_x; initially, cluster_x == 0
  for(i = 0; i < ndata; i++) {
    cluster_assign[i] = cluster_x;
  }

  cluster_size[cluster_x] = ndata;

  while(curr_cluster_count != k) {
    // Bisect cluster_x -> cluster with max SSE
    // Initially, there is only one cluster, so cluster_x is still cluster with max SSE
    set_new_cluster_centroids(dim, ndata, data, cluster_x, cluster_y,
                              cluster_centroid, cluster_assign, cluster_size);

    num_iterations += two_kmeans(dim, ndata, data, cluster_x, cluster_y,
                                 cluster_size, cluster_start, cluster_radius,
                                 cluster_centroid, cluster_assign);

    cluster_with_max_SSE = get_max_SSE(dim, ndata, data, k, cluster_centroid,
                                       cluster_sse, cluster_assign);

    cluster_x = cluster_with_max_SSE;
    cluster_y = cluster_y + 1;

    curr_cluster_count++;
  }

  num_iterations += kmeans_bkm(dim, ndata, data, k, cluster_size,
                               cluster_start, cluster_radius,
                               cluster_centroid, cluster_assign);

  return num_iterations;
}

int two_kmeans(int dim, int ndata, double *data, int cluster_x, int cluster_y,
               int *cluster_size, int *cluster_start, double *cluster_radius,
               double **cluster_centroid, int *cluster_assign)
{
  int i, numIterations = 0;
  bool notFinishedClustering = true;

  int *prev_cluster_assign = malloc(ndata * sizeof(double));
  for(i = 0; i < ndata; i++) {
    prev_cluster_assign[i] = -1;
  }

  while(notFinishedClustering) {
    numIterations++;

    // Assign data points to cluster_x or cluster_y
    for(i = 0; i < ndata; i++) {
      // Only check points that are in cluster_x -> cluster with max SSE
      if(cluster_assign[i] == cluster_x || cluster_assign[i] == cluster_y) {
        two_assignPtToCluster(dim, cluster_x, cluster_y, i, data, cluster_centroid,
                              cluster_assign);
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

    if(! clusterAssignmentsChanged_bkm(ndata, cluster_assign, prev_cluster_assign)) {
      notFinishedClustering = false;
    }
    else {
      setPrevClusterAssignments_bkm(ndata, cluster_assign, prev_cluster_assign);

      two_updateClusterCentroids(dim, cluster_x, cluster_y, ndata, data, cluster_centroid,
                                 cluster_size, cluster_assign);
    }
  }

  return numIterations;
}

int kmeans_bkm(int dim, int ndata, double *data, int k,
               int *cluster_size, int *cluster_start,
               double *cluster_radius, double **cluster_centroid,
               int *cluster_assign)
{
  int i, j, numIterations = 0;

  int *prev_cluster_assign = malloc(ndata * sizeof(double));
  for (i = 0; i < ndata; i++) {
    prev_cluster_assign[i] = -1;
  }

  int thresh_hold = 2;
  while (numIterations < thresh_hold) {
    numIterations++;

    // Assign data point to clusters and update centroids
    for (i = 0; i < ndata; i++) {
      assignPtToCluster_bkm(dim, k, i, data, cluster_centroid, cluster_assign);
    }

    // Reset cluster sizes because of new data point cluster assignments
    for (i = 0; i < k; i++) {
      cluster_size[i] = 0;
    }

    // Assign cluster sizes
    for (i = 0; i < ndata; i++) {
      for (j = 0; j < k; j++) {
        if (cluster_assign[i] == j) {
          cluster_size[j]++;
        }
      }
    }

    setPrevClusterAssignments_bkm(ndata, cluster_assign, prev_cluster_assign);

    updateClusterCentroids_bkm(dim, k, ndata, data, cluster_centroid,
                               cluster_size, cluster_assign);
  }

  quick_sort_data_bkm(dim, 0, ndata, data, cluster_assign);

  setClusterStart_bkm(k, cluster_size, cluster_start);

  for (i = 0; i < k; i++) {
    setClusterRadius_bkm(dim, i, cluster_size, cluster_start, data,
                         cluster_centroid, cluster_radius);
  }

  return numIterations;
}

void set_new_cluster_centroids(int dim, int ndata, double *data, int cluster_x, int cluster_y,
                               double **cluster_centroid, int *cluster_assign, int *cluster_size)
{
  int i, j, centroid_pt_x = -1, centroid_pt_y = -1;
  double distance, maxDistance = (double) INT_MIN;

  int *cluster_x_pts = malloc(cluster_size[cluster_x] * sizeof(int));

  // centroid_pt_x is randomly chosen from the data set of cluster_x
  j = 0;
  for(i = 0; i < ndata; i++) {
    if(cluster_assign[i] == cluster_x) {
      cluster_x_pts[j] = i;
      j++;
    }
  }

//  srand(time(NULL));
  centroid_pt_x = cluster_x_pts[cluster_size[cluster_x] + (rand() / (RAND_MAX / (0 - cluster_size[cluster_x])))];

  for(j = 0; j < dim; j++) {
    cluster_centroid[cluster_x][j] = data[centroid_pt_x * dim + j];
  }

//  printf("Centroid point for cluster_x = %d : (", centroid_pt_x+1);
//  for(i = 0; i < dim; i++) {
//    if(i == dim - 1) { printf("%lf)\n", cluster_centroid[cluster_x][i]); }
//    else { printf("%lf, ", cluster_centroid[cluster_x][i]); }
//  }

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

//  printf("Centroid point for cluster_y = %d : (", centroid_pt_y+1);
//  for(i = 0; i < dim; i++) {
//    if(i == dim - 1) { printf("%lf)\n", cluster_centroid[cluster_y][i]); }
//    else { printf("%lf, ", cluster_centroid[cluster_y][i]); }
//  }

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

void assignPtToCluster_bkm(int dim, int totalClusters, int data_pt, double *data,
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

void two_assignPtToCluster(int dim, int cluster_x, int cluster_y, int data_pt, double *data,
                           double **cluster_centroid, int *cluster_assign)
{
  int j;
  double distance_to_cluster_x = 0.0, distance_to_cluster_y = 0.0;

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

void two_updateClusterCentroids(int dim, int cluster_x, int cluster_y, int ndata, double *data,
                                double **cluster_centroid, int *cluster_size,
                                int *cluster_assign)
{
  int j;

  // Only calculate when cluster_x size is > 0
  if(cluster_size[cluster_x] > 0) {
    for(j = 0; j < dim; j++) {
      cluster_centroid[cluster_x][j] = calcCentroid_bkm(dim, j, cluster_x, ndata, data,
                                                        cluster_size, cluster_assign);
    }
  }

  // Only calculate when cluster_y size is > 0
  if(cluster_size[cluster_y] > 0) {
    for(j = 0; j < dim; j++) {
      cluster_centroid[cluster_y][j] = calcCentroid_bkm(dim, j, cluster_y, ndata, data,
                                                        cluster_size, cluster_assign);
    }
  }
}

void updateClusterCentroids_bkm(int dim, int totalClusters, int ndata, double *data,
                                double **cluster_centroid, int *cluster_size,
                                int *cluster_assign)
{
  int i, j;
  for (i = 0; i < totalClusters; i++) {
    if (cluster_size[i] > 0) { // Only calculate when cluster size is > 0
      for (j = 0; j < dim; j++) {
        cluster_centroid[i][j] = calcCentroid_bkm(dim, j, i, ndata, data,
                                                  cluster_size, cluster_assign);
      }
    }
  }
}

double calcCentroid_bkm(int dim, int currDim, int currCluster, int ndata,
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

bool clusterAssignmentsChanged_bkm(int ndata, int *cluster_assign, int *prev_cluster_assign)
{
  int i;

  for(i = 0; i < ndata; i++) {
    if(cluster_assign[i] != prev_cluster_assign[i]) {
      return true;
    }
  }

  return false;
}

void setPrevClusterAssignments_bkm(int ndata, int *cluster_assign, int *prev_cluster_assign)
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

void setClusterStart_bkm(int totalClusters, int *cluster_size, int *cluster_start)
{
  int i, nextStartingPoint = cluster_size[0];

  cluster_start[0] = 0;
  for(i = 1; i < totalClusters; i++) {
    cluster_start[i] = nextStartingPoint;
    nextStartingPoint += cluster_size[i];
  }
}

void setClusterRadius_bkm(int dim, int currCluster, int *cluster_size,
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

void search_clusters_bkm(int dim, int ndata, double *data, int k, int *cluster_size,
                         int *cluster_start, double *cluster_radius,
                         double **cluster_centroid, double *query,
                         double *result_pt)
{
  int i, j, closestClusterIndex = -1, numPointsChecked = 0,
      clusterWithShorterDistanceIndex = -1;
  double distance, minClusterDistance = (double) INT_MAX,
      minCurrPointDist = (double) INT_MAX;
  double *current_pt = malloc(dim * sizeof(double));
  bool notFinishedSearching = false;

  double *cluster_distances = malloc(k * sizeof(double));
  // Initialize cluster distances to 0.0
  for(i = 0; i < k; i++) {
    cluster_distances[i] = 0.0;
  }

  // Find closest cluster to the query point
  for(i = 0; i < k; i++) {
    distance = 0.0;
    for(j = 0; j < dim; j++) {
      distance += (cluster_centroid[i][j] - query[j]) *
                  (cluster_centroid[i][j] - query[j]);
    }
    cluster_distances[i] = sqrt(distance);

    if(cluster_distances[i] < minClusterDistance) {
      minClusterDistance = cluster_distances[i];
      closestClusterIndex = i;
    }
  }

  // Search points in the closest cluster
  minCurrPointDist = searchPointsInCluster_bkm(dim, query, current_pt, data,
                                               closestClusterIndex, cluster_start,
                                               cluster_size, &numPointsChecked);

  for(i = 0; i < dim; i++) {
    result_pt[i] = current_pt[i];
  }

  // Set the cluster with current min distance to -1, so it is not checked
  // again.
  cluster_distances[closestClusterIndex] = -1.0;

  printf("\nFirst check of other clusters...\n");
  // Check other clusters
  clusterWithShorterDistanceIndex = checkOtherClusters_bkm(k, cluster_distances,
                                                           cluster_radius, closestClusterIndex,
                                                           minCurrPointDist);

  if(clusterWithShorterDistanceIndex != closestClusterIndex) {
    notFinishedSearching = true;
    cluster_distances[clusterWithShorterDistanceIndex] = -1.0;
  }
  else {
    printf("\n**Found closest point within the first checked cluster.**");
  }

  double minNextPointDistance = (double) INT_MAX;
  while(notFinishedSearching) {
    minNextPointDistance = searchPointsInCluster_bkm(dim, query, current_pt, data,
                                                     clusterWithShorterDistanceIndex, cluster_start,
                                                     cluster_size, &numPointsChecked);

    if(minNextPointDistance < minCurrPointDist) {
      minCurrPointDist = minNextPointDistance;
      for(i = 0; i < dim; i++) {
        result_pt[i] = current_pt[i];
      }
    }

    // Continue to check other clusters
    clusterWithShorterDistanceIndex = checkOtherClusters_bkm(k, cluster_distances,
                                                             cluster_radius, closestClusterIndex,
                                                             minCurrPointDist);

    if(clusterWithShorterDistanceIndex == closestClusterIndex) {
      notFinishedSearching = false;
    }
    else {
      cluster_distances[clusterWithShorterDistanceIndex] = -1.0;
      closestClusterIndex = clusterWithShorterDistanceIndex;
    }
  }

  printf("\n\nKMeans cluster search final result point = ( ");
  for(i = 0; i < dim; i++) {
    if(i == dim - 1) { printf("%lf ", result_pt[i]); }
    else { printf("%lf, ", result_pt[i]); }
  }
  printf("),\nwith distance = %lf\n", minCurrPointDist);
  printf("Number of points checked = %d\n", numPointsChecked);

  free(cluster_distances);
  free(current_pt);
}

double searchPointsInCluster_bkm(int dim, double *query, double *current_pt,
                                 double *data, int closestClusterIndex,
                                 int *cluster_start, int *cluster_size,
                                 int *numPointsChecked)
{
  int i, j, end = cluster_size[closestClusterIndex] + cluster_start[closestClusterIndex], closestPointIndex = -1;
  double distance, minPointDistance = (double) INT_MAX;

  for(i = cluster_start[closestClusterIndex]; i < end; i++) {
    distance = 0.0;
    for(j = 0; j < dim; j++) {
      distance += (query[j] - data[i * dim + j]) * (query[j] - data[i * dim + j]);
    }
    *numPointsChecked = *numPointsChecked + 1;
    distance = sqrt(distance);

    if(distance < minPointDistance) {
      minPointDistance = distance;
      closestPointIndex = i;
    }
  }

  j = 0;
  for(i = closestPointIndex * dim; i < closestPointIndex * dim + dim; i++) {
    current_pt[j] = data[i];
    j++;
  }

  return minPointDistance;
}

int checkOtherClusters_bkm(int totalClusters, double *cluster_distances,
                           double *cluster_radius, int closestClusterIndex,
                           double minCurrPointDist)
{
  int i;
  for(i = 0; i < totalClusters; i++) {
    if(cluster_distances[i] != -1.0) {
      if(minCurrPointDist > cluster_distances[i] - cluster_radius[i]) {
        return i;
      }
    }
  }

  return closestClusterIndex;
}

