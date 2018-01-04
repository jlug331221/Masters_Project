#ifndef HELPING_PROCEDURES_H
#define HELPING_PROCEDURES_H

void execute_kdtree(int *train_labels, double *train_features, int *test_labels, double *test_features);

void execute_kmeans(int *train_labels, double *train_features, int *test_labels, double *test_features);

void execute_LSH(int *train_labels, double *train_features, int *test_labels, double *test_features);

/**
 * Read the binary data set. In this case, reading the features and labels.
 */
void read_binary_dataset(char *path, int size, int *labels, double *features);

/**
 * Normalize ndata*feature_dimensions values in data.
 */
void normalize_data(double *data, int feature_dimensions, int ndata);

#endif //HELPING_PROCEDURES_H
