/**
 *
 * The following program ...
 *
 */

#include "Header_Files/Headers.h"
#include "Header_Files/Defs.h"
#include "Header_Files/LSH.h"
#include "Header_Files/kdtree.h"
#include "Header_Files/kmeans.h"

/**************************************************************************************
 **************************************************************************************/

void execute_kdtree(int *train_labels, double *train_features, int *test_labels, double *test_features) {
  int ndata = TRAIN_SIZE, dim = FEATURE_DIM, k = 10, i, j, correct_labeling_count = 0;

  int *cluster_assign = malloc(ndata * sizeof(cluster_assign));
  int *cluster_size = malloc(k * sizeof(cluster_size));
  int *cluster_start = malloc(k * sizeof(cluster_start));

  // Initialize cluster assignments
  for(i = 0; i < ndata; i++) {
    cluster_assign[i] = -1;
  }

  // Initialize cluster start and cluster size
  for(i = 0; i < k; i++) {
    cluster_size[i] = 0;
    cluster_start[i] = 0;
  }

  double **cluster_bdry = malloc(k * sizeof(double*));
  double **cluster_centroid = malloc(k * sizeof(double*));
  for(i = 0; i < k; i++) {
    cluster_bdry[i] = malloc((dim*2) * sizeof(double));
    cluster_centroid[i] = malloc(dim * sizeof(double));
  }

  // Initialize cluster boundaries
  for(i = 0; i < k; i++) {
    for(j = 0; j < dim*2; j++) {
      if(j % 2 == 0) { cluster_bdry[i][j] = (double) INT_MAX; }
      else { cluster_bdry[i][j] = (double) INT_MIN; }
    }
  }

  // Initialize cluster centroids
  for(i = 0; i < k; i++) {
    for(j = 0; j < dim; j++) {
      cluster_centroid[i][j] = 0.0;
    }
  }

  printf("\nBuilding KDTree...\n");
  kdtree(dim, TRAIN_SIZE, train_features, train_labels, k,
         cluster_size, cluster_start, cluster_bdry,
         cluster_centroid, cluster_assign);

  printf("\nPerforming searches using test data...\n");

  double *query = malloc(dim * sizeof(double));

  int h = 0;
  for(i = 0; i < TEST_SIZE; i++) {
    for(j = i * FEATURE_DIM; j < i * FEATURE_DIM + FEATURE_DIM; j++) {
      query[h] = test_features[j];
      h++;
    }

    h = 0;

    search_kdtree(dim, TRAIN_SIZE, train_features, train_labels, test_labels,
                  k, i, cluster_size, cluster_start, cluster_bdry,
                  query, &correct_labeling_count);

    free(query);
    query = malloc(dim * sizeof(double));
  }

  printf("Accuracy of labeling = %.2f%%\n", ((double) correct_labeling_count / (double) TEST_SIZE) * 100);
}

void execute_kmeans(int *train_labels, double *train_features, int *test_labels, double *test_features) {
  int ndata = TRAIN_SIZE, dim = FEATURE_DIM, k = 10,
      i, j, correct_labeling_count = 0, thresh_hold = 2;

  int *cluster_size = malloc(k * sizeof(double));
  int *cluster_start = malloc(k * sizeof(double));
  int *cluster_assign = malloc(ndata * sizeof(double));

  // Initialize cluster assignments
  for (i = 0; i < ndata; i++) {
    cluster_assign[i] = -1;
  }

  // Initialize cluster start
  for (i = 0; i < k; i++) {
    cluster_start[i] = 0;
  }

  double *cluster_radius = malloc(k * sizeof(double));
  double **cluster_centroid = malloc(k * sizeof(double *));
  for (i = 0; i < k; i++) {
    cluster_centroid[i] = malloc(dim * sizeof(double));
    cluster_radius[i] = 0.0;
  }

  // Initialize cluster centroids
  for (i = 0; i < k; i++) {
    for (j = 0; j < dim; j++) {
      cluster_centroid[i][j] = 0.0;
    }
  }

  printf("\nForming %d clusters using K-means...\n", k);
  int numIterations = kmeans(dim, TRAIN_SIZE, train_features, train_labels,
                             k, cluster_size, cluster_start, cluster_radius,
                             cluster_centroid, cluster_assign, thresh_hold);
  printf("Number of iterations for K-means clustering = %d\n", numIterations);

  double *query = malloc(dim * sizeof(double));

  printf("\nPerforming searches using test data...\n");

  int h = 0;
  for(i = 0; i < TEST_SIZE; i++) {
    for(j = i * FEATURE_DIM; j < i * FEATURE_DIM + FEATURE_DIM; j++) {
      query[h] = test_features[j];
      h++;
    }

    h = 0;

    search_clusters(dim, ndata, train_features, train_labels, test_labels,
                    k, i, cluster_size, cluster_start, cluster_radius, cluster_centroid,
                    query, &correct_labeling_count);

    free(query);
    query = malloc(dim * sizeof(double));
  }

  printf("Accuracy of labeling = %.2f%%\n", ((double) correct_labeling_count / (double) TEST_SIZE) * 100);
}

void execute_LSH(int *train_labels, double *train_features, int *test_labels, double *test_features) {
  int dim = FEATURE_DIM, ndata = TRAIN_SIZE, m = 3, i, j;
  int *num_clusters = &i; int correct_labeling_count = 0;
  double w = 7;

  double *b = malloc(m * sizeof(double));
  for(i = 0; i < m; i++) {
    b[i] = 0.0;
  }

  double **r = malloc(m * sizeof(double *));
  for(i = 0; i < m; i++) {
    r[i] = malloc(dim * sizeof(double));
  }

  for(i = 0; i < m; i++) {
    for(j = 0; j < dim; j++) {
      r[i][j] = gauss_rand();
    }
  }

  int **H = malloc(ndata * sizeof(int *));
  for(i = 0; i < ndata; i++) {
    H[i] = malloc(m * sizeof(int));
  }

  printf("\nGenerating clusters via LSH...\n");
  cluster *clusters = LSH(dim, ndata, train_features, m, r, b, w, num_clusters, H);

  printf("\nTotal cluster count = %d\n\n", *num_clusters);

  printf("Performing searches using test data...\n");

  double *q_pt = malloc(dim * sizeof(double));
  int k = 0;
  for(i = 0; i < TEST_SIZE; i++) {
    for(j = i*FEATURE_DIM; j < i*FEATURE_DIM+FEATURE_DIM; j++) {
      q_pt[k] = test_features[j];
      k++;
    }
    k = 0;

    int *q_pt_hash = hash_q_pt(dim, q_pt, m, r, b, w);

    search_clusters_for_apprx_neighbors(dim, train_features, train_labels, test_labels, q_pt,
                                        q_pt_hash, i, clusters, m, &correct_labeling_count);

    free(q_pt);
    q_pt = malloc(dim * sizeof(double));
  }

  printf("Accuracy of labeling = %.2f%%\n", ((double) correct_labeling_count / (double) TEST_SIZE) * 100);
}

void print_time(clock_t begin, clock_t end)
{
  double total_sec = (double) (end-begin)/CLOCKS_PER_SEC, remaining_secs = 0.0, remaining_mins = 0.0;

  if(total_sec >= 60) {
    double mins = floor(total_sec / 60);

    if(mins >= 60) {
      double hours = floor(mins / 60);
      remaining_mins = mins - (60 * hours);
      printf("\nTotal execution time: %.0f hour", hours);
      if(hours > 1) { printf("s %.1f min(s)\n\n", remaining_mins); }
      else { printf(" %.1f min(s)\n\n", remaining_mins); }
    }
    else {
      remaining_secs = total_sec - (mins * 60);
      printf("\nTotal execution time: %d min", (int) mins);
      if(mins > 1) { printf("s %.1f sec(s)\n\n", remaining_secs); }
      else { printf(" %.1f sec(s)\n\n", remaining_secs); }
    }
  }
  else {
    printf("\nTotal execution time: %.1f sec(s)\n\n", (double) (end-begin)/CLOCKS_PER_SEC);
  }
}

int main(int argc, char **argv)
{
  if(argc < 2) {
    printf("\nExecute using the following parameters: ./MNIST_kdtree_kmeans_LSH 1 || ./MNIST_kdtree_kmeans_LSH 2 || "
               "./MNIST_kdtree_kmeans_LSH 3\n");
    printf("where the second argument specifies the algorithm used\n\n");
    printf("1: kdtree\n2: kmeans\n3: LSH\n");

    return 1;
  }

  clock_t begin = clock();
  srand(time(NULL));

  int *train_labels = malloc(sizeof(int) * TRAIN_SIZE);
  double *train_features = malloc(sizeof(double) * TRAIN_SIZE * FEATURE_DIM);

  int *test_labels = malloc(sizeof(int) * TEST_SIZE);
  double *test_features = malloc(sizeof(double) * TEST_SIZE * FEATURE_DIM);

  read_binary_dataset("train.bin", TRAIN_SIZE, train_labels, train_features);
  read_binary_dataset("test.bin", TEST_SIZE, test_labels, test_features);

  normalize_data(train_features, FEATURE_DIM, TRAIN_SIZE);
  normalize_data(test_features, FEATURE_DIM, TEST_SIZE);

  if(strtol(argv[1], NULL, 0) == 1) { execute_kdtree(train_labels, train_features, test_labels, test_features); }

  if (strtol(argv[1], NULL, 0) == 2) { execute_kmeans(train_labels, train_features, test_labels, test_features); }

  if(strtol(argv[1], NULL, 0) == 3) { execute_LSH(train_labels, train_features, test_labels, test_features); }

  clock_t end = clock();

  print_time(begin, end);

  return 0;
}
