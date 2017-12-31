#include "../Header_Files/Headers.h"
#include "../Header_Files/Defs.h"
#include "../Header_Files/kmeans.h"

int kmeans(int dim, int ndata, double *data, int *labels,
           int k, int *cluster_size, int *cluster_start,
           double *cluster_radius, double **cluster_centroid,
           int *cluster_assign, int thresh_hold)
{
  int i, j, numIterations = 0;
  bool notFinishedClustering = true;

  int *prev_cluster_assign = malloc(ndata * sizeof(double));
  for (i = 0; i < ndata; i++) {
    prev_cluster_assign[i] = -1;
  }

  setInitialClusterCentroids(dim, ndata, data, k, cluster_centroid);

  printf("\nDone setting initial centroids.\n");

  while (notFinishedClustering && numIterations < thresh_hold) {
    numIterations++;

    // Assign data point to clusters and update centroids
    for (i = 0; i < ndata; i++) {
      assignPtToCluster(dim, k, i, data, cluster_centroid, cluster_assign);
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

    if (! clusterAssignmentsChanged(ndata, cluster_assign, prev_cluster_assign)) {
      notFinishedClustering = false;
    } else {
      setPrevClusterAssignments(ndata, cluster_assign, prev_cluster_assign);

      updateClusterCentroids(dim, k, ndata, data, cluster_centroid,
                             cluster_size, cluster_assign);
    }
  }

  printf("\nSorting...\n\n");

  quick_sort_data(dim, 0, ndata, data, labels, cluster_assign);

  setClusterStart(k, cluster_size, cluster_start);

  for (i = 0; i < k; i++) {
    setClusterRadius(dim, i, cluster_size, cluster_start, data,
                     cluster_centroid, cluster_radius);
  }

  return numIterations;
}

void setInitialClusterCentroids(int dim, int ndata, double *data,
                                int totalClusters, double **cluster_centroid)
{
  int i, j;
  double distance, maxDistance = (double) INT_MIN;

  // First centroid is randomly chosen from the data set
  int firstCentroid_pt = rand() % ndata;
  for (i = 0; i < dim; i++) {
    cluster_centroid[0][i] = data[firstCentroid_pt * dim + i];
  }

  int secondCentroid_pt = -1;
  // Second centroid is the furthest away point from first centroid
  for (i = 0; i < ndata; i++) {
    distance = 0.0;
    for (j = 0; j < dim; j++) {
      distance += (cluster_centroid[0][j] - data[i * dim + j]) *
                  (cluster_centroid[0][j] - data[i * dim + j]);
    }

    distance = sqrt(distance);
    if (distance > maxDistance) {
      maxDistance = distance;
      secondCentroid_pt = i;
    }
  }

  for (i = 0; i < dim; i++) {
    cluster_centroid[1][i] = data[secondCentroid_pt * dim + i];
  }

  // Find the rest of the initial k-2 centroids; the first two are already
  // chosen.
  for (i = 2; i < totalClusters; i++) {
    findNextInitialCentroid(dim, i, ndata, data, cluster_centroid);
  }
}

void assignPtToCluster(int dim, int totalClusters, int data_pt, double *data,
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

void updateClusterCentroids(int dim, int totalClusters, int ndata, double *data,
                            double **cluster_centroid, int *cluster_size,
                            int *cluster_assign)
{
  int i, j;
  for (i = 0; i < totalClusters; i++) {
    if (cluster_size[i] > 0) { // Only calculate when cluster size is > 0
      for (j = 0; j < dim; j++) {
        cluster_centroid[i][j] = calcCentroid(dim, j, i, ndata, data,
                                              cluster_size, cluster_assign);
      }
    }
  }
}

double calcCentroid(int dim, int currDim, int currCluster, int ndata,
                    double *data, int *cluster_size, int *cluster_assign)
{
  int i;
  double sum = 0.0;
  for (i = 0; i < ndata; i++) {
    if (cluster_assign[i] == currCluster) {
      sum += data[i * dim + currDim];
    }
  }

  return sum / cluster_size[currCluster];
}

bool clusterAssignmentsChanged(int ndata, int *cluster_assign,
                               int *prev_cluster_assign)
{
  int i;
  for (i = 0; i < ndata; i++) {
    if (cluster_assign[i] != prev_cluster_assign[i]) {
      return true;
    }
  }

  return false;
}

void setPrevClusterAssignments(int ndata, int *cluster_assign,
                               int *prev_cluster_assign)
{
  int i;
  for (i = 0; i < ndata; i++) {
    prev_cluster_assign[i] = cluster_assign[i];
  }
}

void findNextInitialCentroid(int dim, int totalClusters, int ndata,
                             double *data, double **cluster_centroid)
{
  int i, j, currCluster = 0; // start at the first cluster
  double distance;
  double *minDistances = malloc(ndata * sizeof(double));
  for (i = 0; i < ndata; i++) {
    minDistances[i] = (double) INT_MAX;
  }

  // Calculate distance from each point to each cluster
  while (currCluster < totalClusters) {
    for (i = 0; i < ndata; i++) {
      distance = 0.0;
      for (j = 0; j < dim; j++) {
        distance += (cluster_centroid[currCluster][j] - data[i * dim + j]) *
                    (cluster_centroid[currCluster][j] - data[i * dim + j]);
      }

      distance = sqrt(distance);
      if (distance < minDistances[i]) {
        minDistances[i] = distance;
      }
    }

    currCluster++;
  }

  int nextClusterPt = findMaxMinimumPtDistance(ndata, minDistances);

  for (i = 0; i < dim; i++) {
    cluster_centroid[currCluster][i] = data[nextClusterPt * dim + i];
  }

  free(minDistances);
}

int findMaxMinimumPtDistance(int ndata, double *min_distances)
{
  int i, maxDistanceIndex = -1;
  double maxDistance = (double) INT_MIN;

  for (i = 0; i < ndata; i++) {
    if (min_distances[i] > maxDistance) {
      maxDistance = min_distances[i];
      maxDistanceIndex = i;
    }
  }

  return maxDistanceIndex;
}

void quick_sort_data(int dim, int lo, int hi, double *data, int *labels, int *cluster_assign)
{
  if(lo < hi) {
    int pivot = partition(dim, lo, hi, data, labels, cluster_assign);
    quick_sort_data(dim, lo, pivot, data, labels, cluster_assign);
    quick_sort_data(dim, pivot+1, hi, data, labels, cluster_assign);
  }
}

int partition(int dim, int lo, int hi, double *data, int *labels, int *cluster_assign)
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

    swap_cluster_assign(cluster_assign, i, j);
    kmeans_swap_points(dim, i, j, data);
    kmeans_swap_labels(labels, i, j);
  }
}

void swap_cluster_assign(int *cluster_assign, int i, int j)
{
  int tmp;
  tmp = cluster_assign[j];
  cluster_assign[j] = cluster_assign[i];
  cluster_assign[i] = tmp;
}

void kmeans_swap_points(int dim, int point1, int point2, double *data)
{

  int i;
  double temp;
  int point1Index = point1 * dim, point2Index = point2 * dim;
  for (i = point1Index; i < point1Index + dim; i++) {
    temp = data[point2Index];
    data[point2Index] = data[i];
    data[i] = temp;

    point2Index++;
  }
}

void kmeans_swap_labels(int *labels, int point1, int point2)
{
  int tmp;
  tmp = labels[point2];
  labels[point2] = labels[point1];
  labels[point1] = tmp;
}

void setClusterStart(int totalClusters, int *cluster_size, int *cluster_start)
{
  int i, nextStartingPoint = cluster_size[0];

  cluster_start[0] = 0;
  for (i = 1; i < totalClusters; i++) {
    cluster_start[i] = nextStartingPoint;
    nextStartingPoint += cluster_size[i];
  }
}

void setClusterRadius(int dim, int currCluster, int *cluster_size,
                      int *cluster_start, double *data,
                      double **cluster_centroid, double *cluster_radius)
{
  int i, j, end = cluster_size[currCluster] + cluster_start[currCluster];
  double distance, maxDistance = (double) INT_MIN;

  for (i = cluster_start[currCluster]; i < end; i++) {
    distance = 0.0;
    for (j = 0; j < dim; j++) {
      distance += (cluster_centroid[currCluster][j] - data[i * dim + j]) *
                  (cluster_centroid[currCluster][j] - data[i * dim + j]);
    }

    distance = sqrt(distance);
    if (distance > maxDistance) {
      maxDistance = distance;
    }
  }

  cluster_radius[currCluster] = maxDistance;
}

void search_clusters(int dim, int ndata, double *train_features, int *train_labels, int *test_labels,
                     int k, int query_index, int *cluster_size, int *cluster_start, double *cluster_radius,
                     double **cluster_centroid, double *query, int *correct_labeling_count)
{
  int i, j, closestClusterIndex = -1, closest_apprx_neighbor = -1, correct_label = 0;
  double distance, minClusterDistance = (double) INT_MAX;

  double *cluster_distances = malloc(k * sizeof(double));
  // Initialize cluster distances to 0.0
  for (i = 0; i < k; i++) {
    cluster_distances[i] = 0.0;
  }

  // Find closest cluster to the query point
  for (i = 0; i < k; i++) {
    distance = 0.0;
    for (j = 0; j < dim; j++) {
      distance += (cluster_centroid[i][j] - query[j]) *
                  (cluster_centroid[i][j] - query[j]);
    }
    cluster_distances[i] = sqrt(distance);

    if (cluster_distances[i] < minClusterDistance) {
      minClusterDistance = cluster_distances[i];
      closestClusterIndex = i;
    }
  }

  // Search points in the closest cluster
  closest_apprx_neighbor = kmeans_searchPointsInCluster(dim, query, train_features, closestClusterIndex,
                                                        cluster_start, cluster_size);

  if(train_labels[closest_apprx_neighbor] == test_labels[query_index]) {
    correct_label++; *correct_labeling_count += correct_label;
  }

  free(cluster_distances);
}

int kmeans_searchPointsInCluster(int dim, double *query, double *data, int closestClusterIndex,
                                 int *cluster_start, int *cluster_size)
{
  int i, j,
      end = cluster_size[closestClusterIndex] +
            cluster_start[closestClusterIndex],
      closest_apprx_point_index = -1;
  double distance, minPointDistance = (double) INT_MAX;

  for (i = cluster_start[closestClusterIndex]; i < end; i++) {
    distance = 0.0;
    for (j = 0; j < dim; j++) {
      distance += (query[j] - data[i * dim + j]) * (query[j] - data[i * dim + j]);
    }
    distance = sqrt(distance);

    if (distance < minPointDistance) {
      minPointDistance = distance;
      closest_apprx_point_index = i;
    }
  }

  return closest_apprx_point_index;
}
