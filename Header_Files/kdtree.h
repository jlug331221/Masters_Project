#ifndef KDTREE_H
#define KDTREE_H

/**
 * Builds a kd-tree of k clusters given the initial data set data.
 */
int kdtree(int dim, int ndata, double *data, int *labels, int k,
           int *cluster_size, int *cluster_start, double **cluster_bdry,
           double **cluster_centroid, int *cluster_assign);

/**
 * Used for calling bipartition so as not to modify the given function
 * signature for kdtree.
 *
 * Partitions data of dim dimensions and assigns clusters.
 */
void kdtreeHelper(int dim, int ndata, double *data, int *labels, int treeDepth,
                  int currDepth, int dLeft, int dRight, int kLeft,
                  int kRight, int *cluster_size, int *cluster_start,
                  double **cluster_bdry, double **cluster_centroid,
                  int *cluster_assign);

/**
 * Partitions data using mean (centroid) of the highest dimension in data.
 */
int bipartition(int dim, int io, int im, double *data, int *labels,
                int cluster_size[2], int cluster_start[2],
                double *cluster_bdry[2], double *cluster_centroid[2],
                int *cluster_assign);

/**
 * Searches the kdtree for point query.
 */
void search_kdtree(int dim, int ndata, double *train_features, int *train_labels, int *test_labels,
                   int k, int query_index, int *cluster_size, int *cluster_start, double **cluster_bdry,
                   double *query, int *correct_labeling_count);

/**
 * Returns the approximate closest point to query in closest_cluster.
 */
int kdtree_searchPointsInCluster(int dim, double *query, double *data, int closest_cluster,
                                 int *cluster_start, int *cluster_size);
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
double findBoundaryMin(int currDim, int totalDim, double *data, int start, int end);

/**
 * Returns max boundary in data from start to end.
 */
double findBoundaryMax(int currDim, int totalDim, double *data, int start, int end);

/**
 * Returns mean (centroid) of dimension dim in data.
 */
double calcMean(int currDim, int totalDim, int start, int end, double *data);

/**
 * Returns the variance of dimension currDim in data.
 */
double calcVariance(int currDim, int totalDim, int start, int end, double *data,
                    double dimMean);

#endif //KDTREE_H
