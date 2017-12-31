#ifndef MNIST_KDTREE_KMEANS_LSH_KMEANS_H
#define MNIST_KDTREE_KMEANS_LSH_KMEANS_H

/**
 * Performs clustering of data using K means algorithm.
 */
int kmeans(int dim, int ndata, double *data, int *labels,
           int k, int *cluster_size, int *cluster_start,
           double *cluster_radius, double **cluster_centroid,
           int *cluster_assign, int thresh_hold);

/**
 * Set the initial cluster centroids. Initially, all centroids are points from
 * data.
 */
void setInitialClusterCentroids(int dim, int ndata, double *data,
                                int totalClusters, double **cluster_centroid);

/**
 * Assign data_pt to closest cluster centroid.
 */
void assignPtToCluster(int dim, int totalClusters, int data_pt, double *data,
                       double **cluster_centroid, int *cluster_assign);

/**
 * Updates cluster centroids by calculating the mean of each cluster.
 */
void updateClusterCentroids(int dim, int totalClusters, int ndata, double *data,
                            double **cluster_centroid, int *cluster_size,
                            int *cluster_assign);

/**
 * Returns the mean of currDim points in data.
 */
double calcCentroid(int dim, int currDim, int currCluster, int ndata,
                    double *data, int *cluster_size, int *cluster_assign);

/**
 * Returns true if the current cluster assignments are different from the previous cluster
 * assignments and false otherwise.
 */
bool clusterAssignmentsChanged(int ndata, int *cluster_assign,
                               int *prev_cluster_assign);

/**
 * Sets the previous cluster assignments with the current cluster assignments.
 */
void setPrevClusterAssignments(int ndata, int *cluster_assign,
                               int *prev_cluster_assign);

/**
 * Finds the next centroid point by choosing the data point in data with the maximum
 * minimum distance.
 */
void findNextInitialCentroid(int dim, int totalClusters, int ndata,
                             double *data, double **cluster_centroid);

/**
 * Returns the index of the maximum minimum distance in minDistances.
 */
int findMaxMinimumPtDistance(int ndata, double *min_distances);

/**
 * Perform quicksort on cluster_assign, swapping labels and data accordingly.
 */
void quick_sort_data(int dim, int lo, int hi, double *data, int *labels, int *cluster_assign);

/**
 * Hoare partition - used for quick_sort_data.
 */
int partition(int dim, int lo, int hi, double *data, int *labels, int *cluster_assign);

/**
 * Swap cluster[i] and cluster[j].
 */
void swap_cluster_assign(int *cluster_assign, int i, int j);

/**
 * Swap data[point1] and data[point2].
 */
void kmeans_swap_points(int dim, int point1, int point2, double *data);

/**
 * Swap labels[point1] and labels[point2].
 */
void kmeans_swap_labels(int *labels, int point1, int point2);

/**
 * Set the starting point for each cluster.
 */
void setClusterStart(int totalClusters, int *cluster_size, int *cluster_start);

/**
 * Sets the radius of currCluster.
 */
void setClusterRadius(int dim, int currCluster, int *cluster_size,
                      int *cluster_start, double *data,
                      double **cluster_centroid, double *cluster_radius);

/**
 * Perform a search for query in k clusters.
 */
void search_clusters(int dim, int ndata, double *train_features, int *train_labels, int *test_labels,
                     int k, int query_index, int *cluster_size, int *cluster_start, double *cluster_radius,
                     double **cluster_centroid, double *query, int *correct_labeling_count);

/**
 * Search data points within closestCluster. Return closest point index.
 */
int kmeans_searchPointsInCluster(int dim, double *query, double *data, int closestClusterIndex,
                                 int *cluster_start, int *cluster_size);

#endif //MNIST_KDTREE_KMEANS_LSH_KMEANS_H
