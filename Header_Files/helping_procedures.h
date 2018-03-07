#ifndef HELPING_PROCEDURES_H
#define HELPING_PROCEDURES_H

#include "LSH_cluster_ADT.h"

void print_execution_error_message();

void print_execution_time(clock_t begin, clock_t end, char message[50]);

/**
 * Fetch the training and testing data sets.
 */
void fetch_datasets(int data_set, double **train_feature_data, double **test_feature_data,
                    int **train_non_feature_data, int **test_non_feature_data);

/**
 * Perform clustering_algorithm on data_set and normalize the data set if normalize_data == 'y'.
 *
 * clustering_algorithm:
 * 1 = LSH
 * 2 = Kdtree
 * 3 = BKmeans
 *
 * Then, perform search queries against the clustered data if desired by the user.
 */
void cluster_and_search(int clustering_algorithm, int data_set, char normalize_data,
                        double *train_feature_data, double *test_feature_data,
                        int *train_non_feature_data, int *test_non_feature_data);

Tree execute_LSH(int data_set, int dim, int train_size, double *train_feature_data);

void execute_kdtree(int dim, int k, int train_size, double *train_feature_data);

void execute_bkmeans_j(int data_set, int dim, int k,
                       int train_size, double *train_feature_data,
                       int test_size, double *test_data,
                       int *train_non_feature_data, int *test_non_feature_data);

void execute_kdtree_median(int *train_labels, double *train_features,
                           int *test_labels, double *test_features);

void execute_bkmeans_z(int *train_labels, double *train_features,
                       int *test_labels, double *test_features);

/**
 * Read the MNIST binary data set.
 */
void read_MNIST_binary_dataset(char *file_path, int size,
                               int *non_feature_data, double *feature_data);

/**
 * Read the BIO protein homology binary data set.
 */
void read_BIO_binary_dataset(char *file_path, int size,
                             int *non_feature_data, double *feature_data);

/**
 * Read the HIGGS binary data set.
 */
void read_HIGGS_binary_dataset(char *file_path, int size,
                               int *non_feature_data, double *feature_data);

/**
 * Normalize (ndata * feature_dimensions) of the feature values in data.
 */
void normalize_data_values(double *data, int feature_dimensions, int ndata,
                           int feature_min_value, int feature_max_value);

/**
 * Ask the user if they would like to perform search queries against the clustered data.
 *
 * If so, perform searches against data_set using test_feature_data.
 */
void perform_search_queries(int clustering_algorithm, int data_set,
                            double *train_feature_data, double *test_feature_data,
                            int *train_non_feature_data, int *test_non_feature_data);

/**
 * Perform searches against clustered data.
 */
void search_clusters(int clustering_algorithm, int data_set,
                     double *train_feature_data, double *test_feature_data,
                     int *train_non_feature_data, int *test_non_feature_data);

/**
 * Print the query search results.
 */
void print_search_results(int test_size, double total_closest_neighbor_distance, double total_pts_searched);

/**
 * Generate a random (double) number between M and N inclusive.
 *
 * M <= rand <= N
 */
double randMToN(double M, double N);

/**
 * Write the results of cluster assign for each point in data[dim * ndata] to disk.
 * This is used exclusively for debugging purposes.
 */
void writeResults(int dim, int ndata, double *data, int *cluster_assign);

#endif //HELPING_PROCEDURES_H
