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
void assign_pt_to_cluster_bkm(int dim, int totalClusters, int data_pt, double *data,
                              double **cluster_centroid, int *cluster_assign);

/**
 * Assign data_pt to closest cluster centroid (either cluster_x or cluster_y) if
 * cluster_assign[data_pt] != prev_cluster_assign[data_pt].
 */
void two_assign_pt_to_cluster(int dim, int cluster_x, int cluster_y, int data_pt, double *data,
                              double **cluster_centroid, int *cluster_assign,
                              int *prev_cluster_assign);

/**
 * Updates cluster centroids by calculating the mean of cluster_x and cluster_y.
 */
void two_update_cluster_centroids(int dim, int cluster_x, int cluster_y, int ndata, double *data,
                                  double **cluster_centroid, int *cluster_size,
                                  int *cluster_assign);

/**
 * Update cluster centroids using the calculated mean.
 */
void update_cluster_centroids_bkm(int dim, int totalClusters, int ndata, double *data,
                                  double **cluster_centroid, int *cluster_size,
                                  int *cluster_assign);

/**
 * Returns the mean (centroid) of currDim points in data.
 */
double calc_centroid_bkm(int dim, int currDim, int currCluster, int ndata,
                         double *data, int *cluster_size, int *cluster_assign);

/**
 * Returns true if the current cluster assignments are different from the
 * previous cluster assignments and false otherwise.
 */
bool cluster_assignments_changed_bkm(int ndata, int *cluster_assign, int *prev_cluster_assign);

/**
 * Sets prev_cluster_assignments with assignments in cluster_assign.
 */
void set_prev_cluster_assignments_bkm(int ndata, int *cluster_assign, int *prev_cluster_assign);

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
void set_cluster_start_bkm(int totalClusters, int *cluster_size, int *cluster_start);

/**
 * Sets the radius of cluster_centroid[currCluster]. Radius is defined as the point
 * that is furthers from the cluster_centroid[currCluster].
 */
void set_cluster_radius_bkm(int dim, int currCluster, int *cluster_size,
                           int *cluster_start, double *data,
                           double **cluster_centroid, double *cluster_radius);

/**
 * Searches k clusters for test_feature_data[test_size] query points.
 */
void bkmeans_search_clsuters_for_approx_neighbors(int dim, int test_size, int k,
                                                  double *train_feature_data, double *test_feature_data,
                                                  int *train_non_feature_data, int *test_non_feature_data,
                                                  int *cluster_size, int *cluster_start,
                                                  double *cluster_radius, double **cluster_centroid);

/**
 * Search data points within closest_cluster and return the closest neighbor point distance to query.
 */
double search_points_in_cluster_bkm(int dim, double *query, double *train_feature_data,
                                    int closest_cluster, int *cluster_start, int *cluster_size,
                                    double *pts_searched_in_closest_cluster);

#endif //MASTERS_PROJECT_BISECTING_KMEANS_H
