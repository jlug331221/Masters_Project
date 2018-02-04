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
 * Search clusters for a matching hash of q_pt_hash. If there is a match, calculate the distance of the closest
 * neighbor and determine if the training label equals the test label.
 *
 * Increment correct_labeling_count if the labels match and return the number of points searched.
 */
int search_clusters_for_apprx_neighbors(int dim, double *train_features, int *train_labels,
                                         int *test_labels, double *q_pt, int *q_pt_hash,
                                         int test_query_index, Tree clusters, int m,
                                         int *correct_labeling_count);

/**
 * Calculate and return the distance from q_pt to neighbor_data_pt.
 */
double calc_dist_to_neighbor(int dim, double *data, double *q_pt, int neighbor_data_pt);

/**
 * Return random Gaussian distribution value. Taken from a method described by Abramowitz and Stegun.
 */
double gauss_rand();

#endif //LSH_H
