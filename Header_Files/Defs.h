#ifndef MNIST_KDTREE_KMEANS_LSH_DEFS_H
#define MNIST_KDTREE_KMEANS_LSH_DEFS_H

#define DEBUG             0

#define PI                3.141592654
#define INFINITY          987654321012.3
#define TINY_DELTA        1.0e-05

#define FEATURE_DIM       784
#define LABEL_DIM         1
#define FEATURE_MIN_VALUE 0
#define FEATURE_MAX_VALUE 255
#define TRAIN_SIZE        60000
#define TEST_SIZE         10000

#define min(x,y)          ((x)>(y) ? (y) : (x))
#define max(x,y)          ((x)<(y) ? (y) : (x))

#endif //MNIST_KDTREE_KMEANS_LSH_DEFS_H
