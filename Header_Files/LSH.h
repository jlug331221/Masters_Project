#ifndef MNIST_KDTREE_KMEANS_LSH_LSH_H
#define MNIST_KDTREE_KMEANS_LSH_LSH_H

/**
 * Structure used for the data points in a cluster.
 */
typedef struct data_pt {
    int data_pt;
    struct data_pt *next;
} data_pt;

/**
 * Structure of a cluster. Each cluster has a hash and associated data points.
 */
typedef struct cluster {
    int *cluster_hash;
    data_pt *data_pts;
    struct cluster *next;
} cluster;

/**
 * Hash ndata points of data and form/return clusters accordingly.
 */
cluster* LSH(int dim, int ndata, double *data,
             int m, double **r, double *b, double w,
             int *num_clusters, int **H);

/**
 * Hash data_pt of dimensions dim using Gaussian distribution r. Store the result in H[data_pt][i]
 * where 0 <= i < m.
 */
void hash_pt(int dim, double *data_pt_vector, int data_pt, int m, double **r, double *b, double w, int **H);

/**
 * Return the hash of q_pt_vector using Gaussian distribution r.
 */
int* hash_q_pt(int dim, double *q_pt_vector, int m, double **r, double *b, double w);

/**
 * Return the dot product of vector_a and vector_b.
 */
double dot_product(int dim, const double *vector_a, const double *vector_b);

/**
 * Form and return the clusters from hash table H. Get the cluster count and store in num_clusters.
 */
cluster* form_clusters(int ndata, int m, int **H, int *num_clusters);

/**
 * Search hash table H a determine if H[d_pt] matches a cluster hash in clusters. If there is a match, add d_pt to
 * the cluster. If there is no match of hashes in clusters, create a new cluster and add d_pt to the cluster.
 */
cluster* compare_against_cluster_hash(int d_pt, int m, int **H, cluster *clusters);

/**
 * Search clusters for a matching hash of q_pt_hash. If there is a match, calculate the distance of the closest
 * neighbor and determine if the train label equals the test label. Increment correct_labeling_count if the
 * labels match.
 */
void search_clusters_for_apprx_neighbors(int dim, double *train_features, int *train_labels,
                                         int *test_labels, double *q_pt, int *q_pt_hash, int query_index,
                                         cluster *clusters, int m, int *correct_labeling_count);

/**
 * Calculate and return the distance from q_pt to neighbor_data_pt.
 */
double calc_dist_to_neighbor(int dim, double *data, double *q_pt, int neighbor_data_pt);

/**
 * Return random Gaussian distribution value. Taken from a method described by Abramowitz and Stegun.
 */
double gauss_rand();

/**
 * Print the cluster information of clusters. Used for debugging purposes only.
 */
void print_clusters_info(int m, cluster *clusters);

/**
 * Store the cluster count of clusters in num_clusters.
 */
void get_cluster_count(cluster *clusters, int *num_clusters);

/**
 * Read the binary data set. In this case, reading the features and labels.
 */
void read_binary_dataset(char *path, int size, int *labels, double *features);

/**
 * Normalize ndata*feature_dimensions values in data.
 */
void normalize_data(double *data, int feature_dimensions, int ndata);

#endif //MNIST_KDTREE_KMEANS_LSH_LSH_H
