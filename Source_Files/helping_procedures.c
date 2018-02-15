#include "../Header_Files/Headers.h"
#include "../Header_Files/Defs.h"
#include "../Header_Files/helping_procedures.h"
#include "../Header_Files/LSH.h"
#include "../Header_Files/kdtree.h"
#include "../Header_Files/bkmeans_j.h"
#include "../Header_Files/kdtree_median.h"
#include "../Header_Files/bkmeans_z.h"

double randMToN(double M, double N)
{
  return M + (rand() / ( RAND_MAX / (N-M) ) ) ;
}

void writeResults(int dim, int ndata, double* data, int* cluster_assign)
{
  int i;
  FILE *file;

  file = fopen("data.txt", "w");
  fprintf(file, "%d\n", dim);
  fprintf(file, "%d\n", ndata);

  for (i = 0; i < dim * ndata; i++) {
    fprintf(file, "%lf\n", data[i]);
  }

  for (i = 0; i < ndata; i++) {
    fprintf(file, "%d\n", cluster_assign[i]);
  }

  fclose(file);
}

void execute_LSH(int dim, int train_size, double *train_data, int data_set)
{
  // m is the hash size
  int m, i, j, correct_labeling_count = 0, cluster_count = 0;
  double w;

  if(data_set == 1) { m = 3; w = 6.0; } // MNIST data set
  if(data_set == 2) { m = 30; w = 60.0; } // BIO data set
  if(data_set == 3) { m = 40; w = 90.0; } // HIGGS data set

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

  printf("\nGenerating clusters via LSH...\n");
  Tree clusters = LSH(dim, train_size, train_data, m, r, b, w);

  cluster_count = get_cluster_count(clusters);

  if(DEBUG) {
    int *data_pts = malloc(train_size * sizeof(int));
    for(i = 0; i < train_size; i++) { data_pts[i] = -1; }
    verify_data_pts_clustered(clusters, data_pts, train_size);

    for(i = 0; i < train_size; i++) {
      if(data_pts[i] == -1) {
        printf("\n** Point %d not clustered **\n", i);
      }
    }
  }

  printf("\nTotal cluster count = %d\n\n", cluster_count);

  if(DEBUG) { write_LSH_clusters_info(clusters, dim, m, w, cluster_count); }

//  printf("Performing searches using test data...\n\n");
//
//  double *q_pt = malloc(dim * sizeof(double));
//  int k = 0, pts_searched = 0;
//  for(i = 0; i < TEST_SIZE; i++) {
//    for(j = i*FEATURE_DIM; j < i*FEATURE_DIM+FEATURE_DIM; j++) {
//      q_pt[k] = test_features[j];
//      k++;
//    }
//    k = 0;
//
//    int *q_pt_hash = hash_q_pt(dim, q_pt, m, r, b, w);
//
//    pts_searched += search_clusters_for_apprx_neighbors(dim, train_features, train_labels, test_labels, q_pt,
//                                                        q_pt_hash, i, clusters, m, &correct_labeling_count);
//
//    free(q_pt);
//    q_pt = malloc(dim * sizeof(double));
//  }
//
//  printf("Average amount of points searched per query = %.2f\n", (double) pts_searched / (double) TEST_SIZE);
//
//  printf("Accuracy of labeling = %.2f%%\n", ((double) correct_labeling_count / (double) TEST_SIZE) * 100);
}

void execute_kdtree(int dim, int k, int train_size, double *train_data, int data_set)
{
  int i, j, correct_labeling_count = 0;

  int *cluster_assign = malloc(train_size * sizeof(cluster_assign));
  int *cluster_size = malloc(k * sizeof(cluster_size));
  int *cluster_start = malloc(k * sizeof(cluster_start));

  // Initialize cluster assignments
  for(i = 0; i < train_size; i++) {
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

  printf("\nGenerating KDTree with %d clusters...\n", k);
  kdtree(dim, train_size, train_data, k,
         cluster_size, cluster_start, cluster_bdry,
         cluster_centroid, cluster_assign);

  printf("\nKDTree clustering complete.\n");

//  printf("\nPerforming searches using test data...\n");
//
//  double *query = malloc(dim * sizeof(double));
//
//  int h = 0;
//  for(i = 0; i < TEST_SIZE; i++) {
//    for(j = i * FEATURE_DIM; j < i * FEATURE_DIM + FEATURE_DIM; j++) {
//      query[h] = test_features[j];
//      h++;
//    }
//
//    h = 0;
//
//    search_kdtree(dim, ndata, train_features, train_labels, test_labels,
//                  k, i, cluster_size, cluster_start, cluster_bdry,
//                  query, &correct_labeling_count);
//
//    free(query);
//    query = malloc(dim * sizeof(double));
//  }
//
//  printf("Accuracy of labeling = %.2f%%\n", ((double) correct_labeling_count / (double) TEST_SIZE) * 100);
}

void execute_bkmeans_j(int dim, int k, int train_size, double *train_data, int data_set)
{
  int i, j, num_clusters = 0, correct_labeling_count = 0, num_iterations = 0;

  int *cluster_size = malloc(k * sizeof(double));
  int *cluster_start = malloc(k * sizeof(double));
  int *cluster_assign = malloc(train_size * sizeof(double));

  // Initialize cluster assignments
  for (i = 0; i < train_size; i++) {
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

  if(DEBUG) {
    int debug_size = 1000000, debug_dim = 2, debug_k = k;
    double *data = malloc(debug_size * debug_dim * sizeof(double));
    for(i = 0; i < debug_size * debug_dim; i++) {
      data[i] = randMToN(0, 1000);
    }
    printf("\nForming %d clusters via Bisecting K-means_j...\n", debug_k);
    num_iterations = bisecting_kmeans(debug_dim, debug_size, data, debug_k, cluster_size,
                                      cluster_start, cluster_radius, cluster_centroid,
                                      cluster_assign);

    writeResults(debug_dim, debug_size, data, cluster_assign);
  }
  else {
    printf("\nForming %d clusters via Bisecting K-means_j...\n", k);
    num_iterations = bisecting_kmeans(dim, train_size, train_data, k, cluster_size,
                                         cluster_start, cluster_radius, cluster_centroid,
                                         cluster_assign);
  }



  printf("\nNumber of iterations for bisecting K-means clustering = %d\n", num_iterations);

//  if (! DEBUG) {
//    printf("\nPerforming searches using test data...\n");
//
//    double *query = malloc(dim * sizeof(double));
//
//    int h = 0;
//    for(i = 0; i < TEST_SIZE; i++) {
//      for(j = i * FEATURE_DIM; j < i * FEATURE_DIM + FEATURE_DIM; j++) {
//        query[h] = test_features[j];
//        h++;
//      }
//
//      h = 0;
//
//      search_clusters_bkm(dim, ndata, train_features, train_labels, test_labels,
//                          k, i, cluster_size, cluster_start, cluster_radius, cluster_centroid,
//                          query, &correct_labeling_count);
//
//      free(query);
//      query = malloc(dim * sizeof(double));
//    }
//
//    printf("Accuracy of labeling = %.2f%%\n", ((double) correct_labeling_count / (double) TEST_SIZE) * 100);
//  }
}

void execute_kdtree_median(int *train_labels, double *train_features, int *test_labels, double *test_features)
{
  int ndata = 1000, dim = 2, kk = 16, i, j;

  double *data = malloc(ndata * dim * sizeof(double));
  for(i = 0; i < ndata * dim; i++) {
    data[i] = randMToN(0, 100);
  }

  int *cluster_assign = malloc(ndata * sizeof(cluster_assign));
  int *cluster_size = malloc(kk * sizeof(cluster_size));
  int *cluster_start = malloc(kk * sizeof(cluster_start));

  // Initialize cluster assignments
  for(i = 0; i < ndata; i++) {
    cluster_assign[i] = -1;
  }

  // Initialize cluster start and cluster size
  for(i = 0; i < kk; i++) {
    cluster_size[i] = 0;
    cluster_start[i] = 0;
  }

  double **cluster_bdry = malloc(kk * sizeof(double*));
  double **cluster_centroid = malloc(kk * sizeof(double*));
  for(i = 0; i < kk; i++) {
    cluster_bdry[i] = malloc((dim*2) * sizeof(double));
    cluster_centroid[i] = malloc(dim * sizeof(double));
  }

  // Initialize cluster boundaries
  for(i = 0; i < kk; i++) {
    for(j = 0; j < dim*2; j++) {
      if(j % 2 == 0) { cluster_bdry[i][j] = (double) INT_MAX; }
      else { cluster_bdry[i][j] = (double) INT_MIN; }
    }
  }

  // Initialize cluster centroids
  for(i = 0; i < kk; i++) {
    for(j = 0; j < dim; j++) {
      cluster_centroid[i][j] = 0.0;
    }
  }

  double *buf = malloc(ndata * sizeof(double));
  double *datum = malloc(dim * sizeof(double));

  printf("\nBuilding KDTree (median)...\n");
  kdtree_hybrid(dim, ndata, data, kk,
                cluster_start, cluster_size,
                cluster_bdry, cluster_centroid,
                cluster_assign, datum, buf);

  writeResults(dim, ndata, data, cluster_assign);

//  printf("\nPerforming searches...\n");
//
//  double *query = malloc(dim * sizeof(double));
//  double *result_pt = malloc(dim * sizeof(double));
//
//  int h = 0, pts_searched;
//  for(i = 0; i < TEST_SIZE; i++) {
//    pts_searched = 0;
//
//    for(j = i * FEATURE_DIM; j < i * FEATURE_DIM + FEATURE_DIM; j++) {
//      query[h] = test_features[j];
//      h++;
//    }
//
//    h = 0;
//
//    pts_searched = search_kdtree_hybrid(dim, ndata, train_features, kk,
//                                        cluster_start, cluster_size, cluster_bdry,
//                                        query, result_pt);
//
//    //if(i == 999 || i == 1999 || i == 2999 || i == 3999 || i == 4999) {
//    printf("%d.\tpoints searched = %d\n", i+1, pts_searched);
//    //}
//
//    free(query);
//    free(result_pt);
//    query = malloc(dim * sizeof(double));
//    result_pt = malloc(dim * sizeof(double));
//  }

  printf("\n");
}

void execute_bkmeans_z(int *train_labels, double *train_features, int *test_labels, double *test_features)
{
//  int ndata, dim, kk, i, j, num_clusters = 0, correct_labeling_count = 0;
//
//  double *data;
//  if(DEBUG) {
//    ndata = 100000; dim = 2; kk = 350;
//    data = malloc(ndata * dim * sizeof(double));
//    for(i = 0; i < ndata * dim; i++) {
//      data[i] = randMToN(0, 100);
//    }
//  }
//  else {
//    ndata = TRAIN_SIZE; dim = FEATURE_DIM; kk = 10;
//    data = train_features;
//  }
//
//  int *cluster_size = malloc(kk * sizeof(double));
//  int *cluster_start = malloc(kk * sizeof(double));
//  int *cluster_assign = malloc(ndata * sizeof(double));
//
//  // Initialize cluster assignments
//  for(i = 0; i < ndata; i++) {
//    cluster_assign[i] = -1;
//  }
//
//  // Initialize cluster start and cluster size
//  for(i = 0; i < kk; i++) {
//    cluster_start[i] = 0;
//    cluster_size[i] = 0;
//  }
//
//  double *cluster_radius = malloc(kk * sizeof(double));
//  double *cluster_center = malloc(kk * sizeof(double));
//  for(i = 0; i < kk; i++) {
//    cluster_center[i] = 0.0;
//    cluster_radius[i] = 0.0;
//  }
//
//  double *datum = malloc(dim * sizeof(double));
//
//  double *cluster_ssd = malloc(kk * sizeof(double));
//
//  printf("\nForming clusters via Bisecting K-means_z...\n\n");
//  num_clusters = bkmeans_z(10, kk, dim, 0, ndata, data,
//                           cluster_assign, datum,
//                           cluster_center, cluster_radius,
//                           cluster_start, cluster_size, cluster_ssd);
//
//  if(DEBUG) { writeResults(dim, ndata, data, cluster_assign); }
//
//  printf("Number of clusters = %d\n", num_clusters);
//
////  printf("\nPerforming searches...\n");
////
////  double *query = malloc(dim * sizeof(double));
////  double *result_pt = malloc(dim * sizeof(double));
////
////  int h = 0, pts_searched;
////  for(i = 0; i < TEST_SIZE; i++) {
////    pts_searched = 0;
////
////    for(j = i * FEATURE_DIM; j < i * FEATURE_DIM + FEATURE_DIM; j++) {
////      query[h] = test_features[j];
////      h++;
////    }
////
////    h = 0;
////
////    pts_searched = search_kdtree_hybrid(dim, ndata, train_features, kk,
////                                        cluster_start, cluster_size, cluster_bdry,
////                                        query, result_pt);
////
////    //if(i == 999 || i == 1999 || i == 2999 || i == 3999 || i == 4999) {
////    printf("%d.\tpoints searched = %d\n", i+1, pts_searched);
////    //}
////
////    free(query);
////    free(result_pt);
////    query = malloc(dim * sizeof(double));
////    result_pt = malloc(dim * sizeof(double));
////  }
////
////  printf("\n");
}

void read_MNIST_binary_dataset(char *file_path, int size, int *non_feature_data, double *feature_data)
{
  FILE *file = fopen(file_path, "rb");

  if(file == NULL) {
    perror("Error");
    exit(1);
  }

  double number_label = -1.0; // 0 - 9

  int i;
  for (i = 0; i < size; i++) {
    // Read the label
    fread(&number_label, sizeof(double), 1, file);
    non_feature_data[i] = (int) number_label;

    // Read the features
    fread(&feature_data[i * MNIST_FEATURE_DIM], sizeof(double), MNIST_FEATURE_DIM, file);
  }

  fclose(file);
}

void read_BIO_binary_dataset(char *file_path, int size, int *non_feature_data, double *feature_data)
{
  FILE *file = fopen(file_path, "rb");

  if(file == NULL) {
    perror("Error");
    exit(1);
  }

  double block_id = -1.0, example_id = -1.0, example_class = -1.0;

  int i, j = 0;
  for(i = 0; i < size; i++) {
    // Read the BLOCK ID
    fread(&block_id, sizeof(double), 1, file);
    non_feature_data[j] = (int) block_id;

    // Read the EXAMPLE ID
    fread(&example_id, sizeof(double), 1, file);
    non_feature_data[j + 1] = (int) example_id;

    // Read the class of the example (1 or 0)
    fread(&example_class, sizeof(double), 1, file);
    non_feature_data[j + 2] = (int) example_class;

    // Read the features
    fread(&feature_data[i * BIO_FEATURE_DIM], sizeof(double), BIO_FEATURE_DIM, file);

    j += BIO_NON_FEATURE_DIM;
  }

  fclose(file);
}

void read_HIGGS_binary_dataset(char *file_path, int size, int *non_feature_data, double *feature_data)
{
  FILE *file = fopen(file_path, "rb");

  if(file == NULL) {
    perror("Error");
    exit(1);
  }

  double class_label = -1.0; // 1 for signal or 0 for background

  int i, j = 0;

  for(i = 0; i < size; i++) {
    // Read the class label
    fread(&class_label, sizeof(double), 1, file);
    non_feature_data[j] = (int) class_label;

    // Read the features
    fread(&feature_data[i * HIGGS_FEATURE_DIM], sizeof(double), HIGGS_FEATURE_DIM, file);

    j += HIGGS_NON_FEATURE_DIM;
  }

  fclose(file);
}

void normalize_MNIST_data(double *data, int feature_dimensions, int ndata)
{
  int i;
  for(i = 0; i < ndata * feature_dimensions; i++) {
    data[i] = (data[i] - MNIST_FEATURE_MIN_VALUE) / (MNIST_FEATURE_MAX_VALUE - MNIST_FEATURE_MIN_VALUE);
  }
}

void normalize_BIO_data(double *data, int feature_dimensions, int ndata)
{
  int i;
  for(i = 0; i < ndata * feature_dimensions; i++) {
    data[i] = (data[i] - BIO_FEATURE_MIN_VALUE) / (BIO_FEATURE_MAX_VALUE - BIO_FEATURE_MIN_VALUE);
  }
}

void normalize_HIGGS_data(double *data, int feature_dimensions, int ndata)
{
  int i;
  for(i = 0; i < ndata * feature_dimensions; i++) {
    data[i] = (data[i] - HIGGS_FEATURE_MIN_VALUE) / (HIGGS_FEATURE_MAX_VALUE - HIGGS_FEATURE_MIN_VALUE);
  }
}
