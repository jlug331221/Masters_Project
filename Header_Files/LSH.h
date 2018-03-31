#ifndef LSH_H
#define LSH_H

#include "../Header_Files/LSH_cluster_ADT.h"

/**
 * Hash ndata points of data and form clusters accordingly.
 */
Tree LSH(int dim, int ndata, const double *data,
             int m, double **r, double *b, double w);

/**
 * Hash data_pt_vector of dim dimensions using Gaussian distribution r. Store the result in pt_hash[i],
 * where 0 <= i < m.
 */
void hash_pt(int dim, double *data_pt_vector, int m, double **r, const double *b, double w, int *pt_hash);

/**
 * Add pt to clusters using pt_hash. Create new cluster when pt_hash does not match
 * any of the hashes in clusters. Otherwise, add data_pt to matching cluster hash.
 */
Tree add_pt_to_cluster(Tree clusters, int pt, int *pt_hash, int m);

/**
 * Return the dot product of vector_a and vector_b.
 */
double dot_product(int dim, const double *vector_a, const double *vector_b);

/**
 * Return the hash of q_pt_vector using Gaussian distribution r.
 */
int* hash_q_pt(int dim, double *q_pt_vector, int m, double **r, const double *b, double w);

/**
 * Search clusters using test_feature_data. Calculate the distance of the closest neighbor and keep
 * an average of all distances to neighbors.
 */
void LSH_search_clusters_for_approx_neighbors(Tree clusters, int dim, int test_size,
                                              int m, double w, double *b, double **r,
                                              double *train_feature_data, double *test_feature_data,
                                              int *train_non_feature_data, int *test_non_feature_data);

/**
 * Calculate and return the distance from q_pt to neighbor_data_pt.
 */
double calc_dist_to_neighbor(int dim, double *data, double *q_pt, int neighbor_data_pt);

/**
 * Return random normal (Gaussian) distribution value taken from the Box-Muller transformation. This
 * generator is pseudo-random.
 */
double gauss_rand();

#endif //LSH_H
