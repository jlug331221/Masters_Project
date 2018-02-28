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

void execute_LSH(int data_set, int dim,
                 int train_size, double *train_data,
                 int test_size, double *test_data,
                 int *train_non_feature_data, int *test_non_feature_data)
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
  printf("\nLSH generated %d clusters.\n", cluster_count);

  if(DEBUG) {
    int *data_pts = malloc(train_size * sizeof(int));
    for (i = 0; i < train_size; i++) { data_pts[i] = -1; }
    verify_data_pts_clustered(clusters, data_pts, train_size);

    for (i = 0; i < train_size; i++) {
      if (data_pts[i] == -1) {
        printf("\n** Point %d not clustered **\n", i);
      }
    }

    write_LSH_clusters_info(clusters, dim, m, w, cluster_count);
  }
}

void execute_kdtree(int data_set, int dim, int k,
                    int train_size, double *train_data,
                    int test_size, double *test_data,
                    int *train_non_feature_data, int *test_non_feature_data)
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
}

void execute_bkmeans_j(int data_set, int dim, int k,
                       int train_size, double *train_data,
                       int test_size, double *test_data,
                       int *train_non_feature_data, int *test_non_feature_data)
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

void perform_queries(int clustering_algorithm, int data_set,
                     double *train_feature_data, double *test_feature_data,
                     int *train_non_feature_data, int *test_non_feature_data)
{
  char perform_queries = '\0';
  printf("\nWould you like to perform queries against the ");

  if(data_set == 1) { printf("MNIST "); }
  if(data_set == 2) { printf("BIO "); }
  if(data_set == 3) { printf("HIGGS "); }

  printf("data set?\nEnter 'y' or 'n': ");
  scanf("%s", &perform_queries);

  while(perform_queries != 'y' && perform_queries != 'n') {
    printf("\nPlease enter 'y' or 'n' if you do or do not want to perform queries: ");
    scanf("%s", &perform_queries);
  }

  if(perform_queries == 'n') { return; }

  search_clusters(clustering_algorithm, data_set,
                  train_feature_data, test_feature_data,
                  train_non_feature_data, test_non_feature_data);
}

void search_clusters(int clustering_algorithm, int data_set,
                     double *train_feature_data, double *test_feature_data,
                     int *train_non_feature_data, int *test_non_feature_data)
{
  printf("TODO: Search clusters\n");
}
