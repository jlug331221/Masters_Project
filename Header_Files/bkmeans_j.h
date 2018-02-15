#ifndef MASTERS_PROJECT_BISECTING_KMEANS_H
#define MASTERS_PROJECT_BISECTING_KMEANS_H

/**
 * Performs bisecting kmeans clustering on data.
 */
int bisecting_kmeans(int dim, int ndata, double *data, int k,
                     int *cluster_size, int *cluster_start,
                     double *cluster_radius, double **cluster_centroid,
                     int *cluster_assign);

/**
 * Performs clustering of data using K-means algorithm for cluster_x and cluster_y.
 */
int two_kmeans(int dim, int ndata, double *data, int cluster_x, int cluster_y,
               int *cluster_size, double **cluster_centroid, int *cluster_assign);

/**
 * Applies K-means clustering on data. This is the original K-means algorithm.
 */
int kmeans_bkm(int dim, int ndata, double *data, int k,
               int *cluster_size, int *cluster_start,
               double *cluster_radius, double **cluster_centroid,
               int *cluster_assign);

/**
 * Set the new cluster centroids. Points are taken from the data set of cluster_x.
 */
void set_new_cluster_centroids(int dim, int ndata, double *data, int cluster_x, int cluster_y,
                               double **cluster_centroid, int *cluster_assign, int *cluster_size);

/**
 * Returns the cluster with the max SSE. Stores subsequent SSE calculations for each cluster in
 * cluster_sse.
 */
int get_max_SSE(int dim, int ndata, double *data, int k, double **cluster_centroid,
                double *cluster_sse, int *cluster_assign);

/**
 * Returns the calculated sum of squares error for curr_cluster.
 */
double calc_SSE(int dim, int ndata, double *data, int curr_cluster,
                double **cluster_centroid, int *cluster_assign);

/**
 * Assign data_pt to the closest cluster_centroid.
 */
void assignPtToCluster_bkm(int dim, int totalClusters, int data_pt, double *data,
                           double **cluster_centroid, int *cluster_assign);

/**
 * Assign data_pt to closest cluster centroid (either cluster_x or cluster_y) if
 * cluster_assign[data_pt] != prev_cluster_assign[data_pt].
 */
void two_assignPtToCluster(int dim, int cluster_x, int cluster_y, int data_pt, double *data,
                           double **cluster_centroid, int *cluster_assign,
                           int *prev_cluster_assign);

/**
 * Updates cluster centroids by calculating the mean of cluster_x and cluster_y.
 */
void two_updateClusterCentroids(int dim, int cluster_x, int cluster_y, int ndata, double *data,
                                double **cluster_centroid, int *cluster_size,
                                int *cluster_assign);

/**
 * Update cluster centroids using the calculated mean.
 */
void updateClusterCentroids_bkm(int dim, int totalClusters, int ndata, double *data,
                                double **cluster_centroid, int *cluster_size,
                                int *cluster_assign);

/**
 * Returns the mean (centroid) of currDim points in data.
 */
double calcCentroid_bkm(int dim, int currDim, int currCluster, int ndata,
                        double *data, int *cluster_size, int *cluster_assign);

/**
 * Returns true if the current cluster assignments are different from the
 * previous cluster assignments and false otherwise.
 */
bool clusterAssignmentsChanged_bkm(int ndata, int *cluster_assign, int *prev_cluster_assign);

/**
 * Sets prev_cluster_assignments with assignments in cluster_assign.
 */
void setPrevClusterAssignments_bkm(int ndata, int *cluster_assign, int *prev_cluster_assign);

/**
 * Perform quicksort on cluster_assign, swapping data accordingly.
 */
void quick_sort_data_bkm(int dim, int lo, int hi, double *data, int *cluster_assign);

/**
 * Hoare partition - used for quick_sort_data_bkm.
 */
int partition_bkm(int dim, int lo, int hi, double *data, int *cluster_assign);

/**
 * Swap cluster[i] and cluster[j].
 */
void swap_cluster_assign_bkm(int *cluster_assign, int i, int j);

/**
 * Swap point1 with point2 in data.
 */
void swap_points_bkm(int dim, int point1, int point2, double *data);

/**
 * Swap label1 with label2 in labels.
 */
void swap_labels_bkm(int *labels, int label1, int label2);

/**
 * Set the starting point for each cluster.
 */
void setClusterStart_bkm(int totalClusters, int *cluster_size, int *cluster_start);

/**
 * Sets the radius of cluster_centroid[currCluster]. Radius is defined as the point
 * that is furthers from the cluster_centroid[currCluster].
 */
void setClusterRadius_bkm(int dim, int currCluster, int *cluster_size,
                          int *cluster_start, double *data,
                          double **cluster_centroid, double *cluster_radius);

/**
 * Perform a search for query in K clusters and return the number of points
 * searched.
 */
void search_clusters_bkm(int dim, int ndata, double *data, int *train_labels, int *test_labels,
                         int k, int query_index, int *cluster_size, int *cluster_start,
                         double *cluster_radius, double **cluster_centroid, double *query,
                         int *correct_labeling_count);

/**
 * Search data points within closestCluster. Returns closest neighbor point to query[dim].
 */
int searchPointsInCluster_bkm(int dim, double *query, double *data, int closestClusterIndex,
                              int *cluster_start, int *cluster_size);

/**
 * Check other cluster distances for a distance shorter than minCurrPointDist.
 * If there is such a cluster with a shorter distance, return the index of that
 * cluster. Otherwise, return closestClusterIndex.
 */
int checkOtherClusters_bkm(int totalClusters, double *cluster_distances,
                           double *cluster_radius, int closestClusterIndex,
                           double minCurrPointDist);

#endif //MASTERS_PROJECT_BISECTING_KMEANS_H
