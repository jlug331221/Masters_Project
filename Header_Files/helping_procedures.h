#ifndef HELPING_PROCEDURES_H
#define HELPING_PROCEDURES_H

void execute_LSH(int data_set, int dim,
                 int train_size, double *train_data,
                 int test_size, double *test_data,
                 int *train_non_feature_data, int *test_non_feature_data);

void execute_kdtree(int data_set, int dim, int k,
                    int train_size, double *train_data,
                    int test_size, double *test_data,
                    int *train_non_feature_data, int *test_non_feature_data);

void execute_bkmeans_j(int data_set, int dim, int k,
                       int train_size, double *train_data,
                       int test_size, double *test_data,
                       int *train_non_feature_data, int *test_non_feature_data);

void execute_kdtree_median(int *train_labels, double *train_features, int *test_labels, double *test_features);

void execute_bkmeans_z(int *train_labels, double *train_features, int *test_labels, double *test_features);

/**
 * Read the MNIST binary data set.
 */
void read_MNIST_binary_dataset(char *file_path, int size, int *non_feature_data, double *feature_data);

/**
 * Read the BIO protein homology binary data set.
 */
void read_BIO_binary_dataset(char *file_path, int size, int *non_feature_data, double *feature_data);

/**
 * Read the HIGGS binary data set.
 */
void read_HIGGS_binary_dataset(char *file_path, int size, int *non_feature_data, double *feature_data);

/**
 * Normalize (ndata * feature_dimensions) of MNIST data values.
 */
void normalize_MNIST_data(double *data, int feature_dimensions, int ndata);

/**
 * Normalize (ndata * feature_dimensions) of BIO data values.
 */
void normalize_BIO_data(double *data, int feature_dimensions, int ndata);

/**
 * Normalize (ndata * feature_dimensions) of HIGGS data values.
 */
void normalize_HIGGS_data(double *data, int feature_dimensions, int ndata);

/**
 * Ask the user if they would like to perform search queries against the clustered data.
 *
 * If so, perform searches against data_set using test_feature_data.
 */
void perform_queries(int clustering_algorithm, int data_set,
                     double *train_feature_data, double *test_feature_data,
                     int *train_non_feature_data, int *test_non_feature_data);

/**
 * Perform searches against clustered data.
 */
void search_clusters(int clustering_algorithm, int data_set,
                     double *train_feature_data, double *test_feature_data,
                     int *train_non_feature_data, int *test_non_feature_data);

#endif //HELPING_PROCEDURES_H
