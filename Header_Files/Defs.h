#ifndef MNIST_KDTREE_KMEANS_LSH_DEFS_H
#define MNIST_KDTREE_KMEANS_LSH_DEFS_H

#define DEBUG                         0

#define PI                            3.141592654
#define TINY_DELTA                    1.0e-05

/**
#define INFINITY                      987654321012.3
*/

#define MNIST_FEATURE_DIM             784
#define BIO_FEATURE_DIM               74
#define HIGGS_FEATURE_DIM             28

#define MNIST_NON_FEATURE_DIM         1
#define BIO_NON_FEATURE_DIM           3
#define HIGGS_NON_FEATURE_DIM         1

// Min and max values are for normalizing the data sets.
#define MNIST_FEATURE_MIN_VALUE       0
#define MNIST_FEATURE_MAX_VALUE       255
#define BIO_FEATURE_MIN_VALUE         -1096.0
#define BIO_FEATURE_MAX_VALUE         64129.4
#define HIGGS_FEATURE_MIN_VALUE       -2.969725
#define HIGGS_FEATURE_MAX_VALUE       20.937159

#define MNIST_TRAIN_SIZE              60000
#define MNIST_TEST_SIZE               10000
#define BIO_TRAIN_SIZE                145751
#define BIO_TEST_SIZE                 139658
#define HIGGS_TRAIN_SIZE              10500000
#define HIGGS_TEST_SIZE               10000

#define MNIST_m                       3
#define BIO_m                         33
#define HIGGS_m                       40

#define MNIST_w                       6.0
#define BIO_w                         50.0
#define HIGGS_w                       90.0

#define min(x,y)          ((x)>(y) ? (y) : (x))
#define max(x,y)          ((x)<(y) ? (y) : (x))

#endif //MNIST_KDTREE_KMEANS_LSH_DEFS_H
