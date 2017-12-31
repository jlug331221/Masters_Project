#include "../Header_Files/Headers.h"
#include "../Header_Files/Defs.h"
#include "../Header_Files/kdtree.h"

int kdtree(int dim, int ndata, double *data, int *labels, int k,
           int *cluster_size, int *cluster_start, double **cluster_bdry,
           double **cluster_centroid, int *cluster_assign)
{
  int treeDepth = 0;
  int currDepth = 0, dLeft = 0, dRight = ndata, kLeft = 0, kRight = k;

  treeDepth = (int) floor(log2(k));

  kdtreeHelper(dim, ndata, data, labels, treeDepth, currDepth, dLeft, dRight,
               kLeft, kRight, cluster_size, cluster_start, cluster_bdry,
               cluster_centroid, cluster_assign);

  return 0;
}

void kdtreeHelper(int dim, int ndata, double *data, int *labels, int treeDepth,
                  int currDepth, int dLeft, int dRight, int kLeft,
                  int kRight, int *cluster_size, int *cluster_start,
                  double **cluster_bdry, double **cluster_centroid,
                  int *cluster_assign)
{
  int midCluster, midData;

  if(currDepth != treeDepth) {
    bipartition(dim, dLeft, dRight, data, labels, &cluster_size[kLeft], &cluster_start[kLeft],
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
  kdtreeHelper(dim, ndata, data, labels, treeDepth, currDepth + 1, dLeft, midData,
               kLeft, midCluster, cluster_size, cluster_start, cluster_bdry,
               cluster_centroid, cluster_assign);

  // Recursive call on right
  kdtreeHelper(dim, ndata, data, labels, treeDepth, currDepth + 1, midData, dRight,
               midCluster, kRight, cluster_size, cluster_start, cluster_bdry,
               cluster_centroid, cluster_assign);
}

int bipartition(int dim, int io, int im, double *data, int *labels,
                int cluster_size[2], int cluster_start[2],
                double *cluster_bdry[2], double *cluster_centroid[2],
                int *cluster_assign)
{
  int i, j, partitionDim = -1;
  double mean = 0.0, maxVar = 0.0, currVar = 0.0, partitionMean = 0.0;

  // Find means of each dimension and determine the dimension with max variance
  for(i = 0; i < dim; i++) {
    mean = calcMean(i, dim, io * dim, im * dim, data);

    currVar = calcVariance(i, dim, io * dim, im * dim, data, mean);

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
    kdtree_swap_labels(labels, i, j);
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
    cluster_bdry[0][bdryMinMaxIndex] = findBoundaryMin(i, dim, data,
                                                       cluster_start[0]*dim, (cluster_start[0]+cluster_size[0])*dim);
    cluster_bdry[0][bdryMinMaxIndex+1] = findBoundaryMax(i, dim, data,
                                                         cluster_start[0]*dim, (cluster_start[0]+cluster_size[0])*dim);
    cluster_bdry[1][bdryMinMaxIndex] = findBoundaryMin(i, dim, data,
                                                       cluster_start[1]*dim, (cluster_start[1]+cluster_size[1])*dim);
    cluster_bdry[1][bdryMinMaxIndex+1] = findBoundaryMax(i, dim, data,
                                                         cluster_start[1]*dim, (cluster_start[1]+cluster_size[1])*dim);

    bdryMinMaxIndex += 2;
  }

  return 0;
}

void search_kdtree(int dim, int ndata, double *train_features, int *train_labels, int *test_labels,
                   int k, int query_index, int *cluster_size, int *cluster_start, double **cluster_bdry,
                   double *query, int *correct_labeling_count)
{
  int i, j, dimMinMaxIndex, closest_cluster = -1, correct_label = 0, apprx_closest_point = -1;
  double distance, min = (double) INT_MAX;

  double *cluster_distances = malloc(k * sizeof(double));
  // Initialize cluster distances to -1.0
  for(i = 0; i < k; i++) {
    cluster_distances[i] = -1.0;
  }

  // Begin searching for the cluster with the shortest distance
  for(i = 0; i < k; i++) {
    distance = 0.0;
    dimMinMaxIndex = 0;
    for(j = 0; j < dim; j++) {
      if(query[j] < cluster_bdry[i][dimMinMaxIndex]) {
        distance += (query[j] - cluster_bdry[i][dimMinMaxIndex]) *
                    (query[j] - cluster_bdry[i][dimMinMaxIndex]);
      }
      if(query[j] > cluster_bdry[i][dimMinMaxIndex+1]) {
        distance += (query[j] - cluster_bdry[i][dimMinMaxIndex+1]) *
                    (query[j] - cluster_bdry[i][dimMinMaxIndex+1]);
      }
      dimMinMaxIndex += 2;
    }
    cluster_distances[i] = sqrt(distance);

    if(cluster_distances[i] < min) {
      min = cluster_distances[i]; closest_cluster = i;
    }
  }

  // Calculate distances for each point in cluster with minimum distance
  apprx_closest_point = kdtree_searchPointsInCluster(dim, query, train_features, closest_cluster,
                                                     cluster_start, cluster_size);

  if(train_labels[apprx_closest_point] == test_labels[query_index]) {
    correct_label++; *correct_labeling_count += correct_label;
  }

  free(cluster_distances);
}

int kdtree_searchPointsInCluster(int dim, double *query, double *data, int closest_cluster,
                                 int *cluster_start, int *cluster_size)
{
  int i, j, dataPointIndex,
      apprx_closest_point = cluster_start[closest_cluster],
      cluster_end = cluster_start[closest_cluster] +
                    cluster_size[closest_cluster];
  double distance, minPointDistance = (double) INT_MAX;

  for(i = cluster_start[closest_cluster]; i < cluster_end; i++) {
    distance = 0.0;
    for(j = 0; j < dim; j++) {
      dataPointIndex = i * dim + j;
      if(query[j] < data[dataPointIndex]) {
        distance += (query[j] - data[dataPointIndex]) *
                    (query[j] - data[dataPointIndex]);
      }

      if(query[j] > data[dataPointIndex]) {
        distance += (query[j] - data[dataPointIndex]) *
                    (query[j] - data[dataPointIndex]);
      }
    }

    distance = sqrt(distance);
    if(distance < minPointDistance) {
      minPointDistance = distance;
      apprx_closest_point = i;
    }
  }

  return apprx_closest_point;
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

double findBoundaryMin(int currDim, int totalDim, double *data, int start, int end)
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

double findBoundaryMax(int currDim, int totalDim, double *data, int start, int end)
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

double calcMean(int currDim, int totalDim, int start, int end, double *data)
{
  int i, ndata = 0;
  double sum = 0.0;
  for(i = currDim+start; i < end; i+=totalDim) {
    sum += data[i];
    ndata++;
  }

  return sum / ndata;
}

double calcVariance(int currDim, int totalDim, int start, int end, double *data,
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
