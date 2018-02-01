#ifndef LSH_H
#define LSH_H

#include "../Header_Files/LSH_cluster_ADT.h"

/**
 * Hash ndata points of data and form clusters accordingly.
 */
cluster* LSH(int dim, int ndata, double *data,
             int m, double **r, double *b, double w,
             int *num_clusters);

/**
 * Hash data_pt_vector of dim dimensions using Gaussian distribution r. Store the result in pt_hash[i],
 * where 0 <= i < m.
 */
void hash_pt(int dim, double *data_pt_vector, int m, double **r, double *b, double w, int *pt_hash);

/**
 * Return the hash of q_pt_vector using Gaussian distribution r.
 */
int* hash_q_pt(int dim, double *q_pt_vector, int m, double **r, double *b, double w);

/**
 * Return the dot product of vector_a and vector_b.
 */
double dot_product(int dim, const double *vector_a, const double *vector_b);

/**
 * Add pt to clusters using pt_hash. Create new cluster and increment num_clusters when pt_hash does not match
 * any of the hashes in clusters. Otherwise, add data_pt to matching cluster hash.
 */
cluster* add_pt_to_cluster(cluster* clusters, int pt, int *pt_hash, int m, int *num_clusters);

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

#endif //LSH_H
