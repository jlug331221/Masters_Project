#ifndef HELPING_PROCEDURES_H
#define HELPING_PROCEDURES_H

void execute_LSH(int dim, int train_size, double *train_data, int data_set);

void execute_kdtree(int *train_labels, double *train_features, int *test_labels, double *test_features);

void execute_bkmeans_j(int *train_labels, double * train_features, int *test_labels, double *test_features);

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

#endif //HELPING_PROCEDURES_H
