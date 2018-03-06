#ifndef KDTREE_H
#define KDTREE_H

/**
 * Builds a kd-tree of k clusters given the initial data set data.
 */
int kdtree(int dim, int ndata, double *data, int k,
           int *cluster_size, int *cluster_start, double **cluster_bdry,
           double **cluster_centroid, int *cluster_assign);

/**
 * Used for calling bipartition so as not to modify the given function
 * signature for kdtree.
 *
 * Partitions data of dim dimensions and assigns clusters.
 */
void kdtree_helper(int dim, int ndata, double *data, int treeDepth,
                   int currDepth, int dLeft, int dRight, int kLeft,
                   int kRight, int *cluster_size, int *cluster_start,
                   double **cluster_bdry, double **cluster_centroid,
                   int *cluster_assign);

/**
 * Partitions data using mean (centroid) of the highest dimension in data.
 */
int bipartition(int dim, int io, int im, double *data,
                int cluster_size[2], int cluster_start[2],
                double *cluster_bdry[2], double *cluster_centroid[2],
                int *cluster_assign);

/**
 * Searches k clusters for test_feature_data[test_size] query points.
 */
void kdtree_search_clusters_for_approx_neighbors(int dim, int test_size, int k,
                                                 double *train_feature_data, double *test_feature_data,
                                                 int *train_non_feature_data, int *test_non_feature_data,
                                                 int *cluster_size, int *cluster_start, double **cluster_bdry);

/**
 * Search closest_cluster for the nearest point (min distance) to query. Keeps track of the points searched
 * in closest_cluster (pts_searched_in_closest_cluster).
 *
 * Returns the min distance to the point closest to query.
 */
double kdtree_search_points_in_cluster(int dim, double *query, double *train_feature_data, int closest_cluster,
                                       int *cluster_start, int *cluster_size,
                                       double *pts_searched_in_closest_cluster);
/**
 * Swap data[point1] and data[point2] of dimensions dim.
 */
void kdtree_swap_points(int dim, double *data, int point1, int point2);

/**
 * Swap labels[point1] and labels[point2].
 */
void kdtree_swap_labels(int *labels, int point1, int point2);

/**
 * Returns min boundary in data from start to end.
 */
double find_boundary_min(int currDim, int totalDim, double *data, int start, int end);

/**
 * Returns max boundary in data from start to end.
 */
double find_boundary_max(int currDim, int totalDim, double *data, int start, int end);

/**
 * Returns mean (centroid) of dimension dim in data.
 */
double calc_mean(int currDim, int totalDim, int start, int end, double *data);

/**
 * Returns the variance of dimension currDim in data.
 */
double calc_variance(int currDim, int totalDim, int start, int end, double *data,
                     double dimMean);

#endif //KDTREE_H
