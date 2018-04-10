#include "../Header_Files/Headers.h"
#include "../Header_Files/Defs.h"
#include "../Header_Files/helping_procedures.h"
#include "../Header_Files/LSH.h"
#include "../Header_Files/kdtree.h"
#include "../Header_Files/bkmeans_j.h"
#include "../Header_Files/kdtree_median.h"
#include "../Header_Files/bkmeans_z.h"

int debug_size = 50000, debug_dim = 2;
double *debug_data = NULL;

Tree clusters = NULL; // used in LSH
double *b = NULL;     // used in LSH - initialized in execute_LSH()
double **r = NULL;    // used in LSH - initialized in execute_LSH()
int m = 0; // hash size, used in LSH - initialized in execute_LSH()
double w = 0.0; // used in LSH - initialized in execute_LSH()

int K = -1;

int *cluster_assign = NULL; // used in Kdtree & Bkmeans
int *cluster_size = NULL;   // used in Kdtree & Bkmeans
int *cluster_start = NULL;  // used in Kdtree & Bkmeans

double **cluster_bdry = NULL; // used in Kdtree - initialized in execute_kdtree()

double *cluster_radius = NULL;    // used in Bkmeans - initialized in execute_bkmeans()
double **cluster_centroid = NULL; // used in Kdtree & Bkmeans

void print_execution_error_message()
{
  printf("\nExecute using the following: ./main <clustering> <dataset>\n\n");

  printf("The second argument specifies the algorithm used to cluster the data and is given "
             "as an integer.\n\n");
  printf("\t1: LSH\n\t2: Kdtree\n\t3: BKmeans\n");

  printf("\nThe third argument specifies the dataset to cluster and is given as an integer.\n\n");
  printf("\t1: MNIST\n\t2: BIO\n\t3: HIGGS\n\n");
}

void print_execution_time(clock_t begin, clock_t end, char message[50])
{
  double total_sec = (double) (end-begin)/CLOCKS_PER_SEC,
      remaining_secs = 0.0, remaining_mins = 0.0;

  if(total_sec >= 60) {
    double mins = floor(total_sec / 60);

    if(mins >= 60) {
      double hours = floor(mins / 60);
      remaining_mins = mins - (60 * hours);
      printf("\n%s execution time: %.0f hour", message, hours);
      if(hours > 1) { printf("s %.1f min(s)\n", remaining_mins); }
      else { printf(" %.1f min(s)\n", remaining_mins); }
    }
    else {
      remaining_secs = total_sec - (mins * 60);
      printf("\n%s execution time: %d min", message, (int) mins);
      if(mins > 1) { printf("s %.1f sec(s)\n", remaining_secs); }
      else { printf(" %.1f sec(s)\n", remaining_secs); }
    }
  }
  else {
    printf("\n%s execution time: %.1f sec(s)\n", message,
           (double) (end-begin)/CLOCKS_PER_SEC);
  }
}

void fetch_datasets(int data_set, double **train_feature_data, double **test_feature_data,
                    int **train_non_feature_data, int **test_non_feature_data)
{
  if(data_set == 1) { // MNIST data set
    printf("\nReading in the MNIST training and testing data sets...\n");

    *train_feature_data = malloc(MNIST_TRAIN_SIZE * MNIST_FEATURE_DIM * sizeof(double));
    *train_non_feature_data = malloc(MNIST_TRAIN_SIZE * MNIST_NON_FEATURE_DIM * sizeof(int)); // train number labels

    *test_feature_data = malloc(MNIST_TEST_SIZE * MNIST_FEATURE_DIM * sizeof(double));
    *test_non_feature_data = malloc(MNIST_TEST_SIZE * MNIST_NON_FEATURE_DIM * sizeof(int)); // test number labels

    read_MNIST_binary_dataset("MNIST_train.bin", MNIST_TRAIN_SIZE, *train_non_feature_data, *train_feature_data);
    read_MNIST_binary_dataset("MNIST_test.bin", MNIST_TEST_SIZE, *test_non_feature_data, *test_feature_data);
  }

  if(data_set == 2) { // BIO data set
    printf("\nReading in the BIO training and testing data sets...\n");

    *train_feature_data = malloc(BIO_TRAIN_SIZE * BIO_FEATURE_DIM * sizeof(double));
    *train_non_feature_data = malloc(BIO_TRAIN_SIZE * BIO_NON_FEATURE_DIM * sizeof(int));

    *test_feature_data = malloc(BIO_TEST_SIZE * BIO_FEATURE_DIM * sizeof(double));
    *test_non_feature_data = malloc(BIO_TEST_SIZE * BIO_NON_FEATURE_DIM * sizeof(int));

    read_BIO_binary_dataset("BIO_train.bin", BIO_TRAIN_SIZE, *train_non_feature_data, *train_feature_data);
    read_BIO_binary_dataset("BIO_test.bin", BIO_TEST_SIZE, *test_non_feature_data, *test_feature_data);
  }

  if(data_set == 3) { // HIGGS data set
    printf("\nReading in the HIGGS training and testing data sets...\n");

    *train_feature_data = malloc(HIGGS_TRAIN_SIZE * HIGGS_FEATURE_DIM * sizeof(double));
    *train_non_feature_data = malloc(HIGGS_TRAIN_SIZE * HIGGS_NON_FEATURE_DIM * sizeof(int));

    *test_feature_data = malloc(HIGGS_TEST_SIZE * HIGGS_FEATURE_DIM * sizeof(double));
    *test_non_feature_data = malloc(HIGGS_TEST_SIZE * HIGGS_NON_FEATURE_DIM * sizeof(int));

    read_HIGGS_binary_dataset("HIGGS_train.bin", HIGGS_TRAIN_SIZE, *train_non_feature_data, *train_feature_data);
    read_HIGGS_binary_dataset("HIGGS_test.bin", HIGGS_TEST_SIZE, *test_non_feature_data, *test_feature_data);
  }
}

void cluster_and_search(int clustering_algorithm, int data_set, char normalize_data,
                        double *train_feature_data, double *test_feature_data,
                        int *train_non_feature_data, int *test_non_feature_data)
{
  if(normalize_data == 'y' && data_set == 1) {
    printf("Normalizing the MNIST data set...\n");
    normalize_data_values(train_feature_data, MNIST_FEATURE_DIM, MNIST_TRAIN_SIZE,
                          MNIST_FEATURE_MIN_VALUE, MNIST_FEATURE_MAX_VALUE);
    normalize_data_values(test_feature_data, MNIST_FEATURE_DIM, MNIST_TEST_SIZE,
                          MNIST_FEATURE_MIN_VALUE, MNIST_FEATURE_MAX_VALUE);
  }

  if(normalize_data == 'y' && data_set == 2) {
    printf("Normalizing the BIO data set...\n");
    normalize_data_values(train_feature_data, BIO_FEATURE_DIM, BIO_TRAIN_SIZE,
                          BIO_FEATURE_MIN_VALUE, BIO_FEATURE_MAX_VALUE);
    normalize_data_values(test_feature_data, BIO_FEATURE_DIM, BIO_TEST_SIZE,
                          BIO_FEATURE_MIN_VALUE, BIO_FEATURE_MAX_VALUE);
  }

  if(normalize_data == 'y' && data_set == 3) {
    printf("Normalizing the HIGGS data set...\n");
    normalize_data_values(train_feature_data, HIGGS_FEATURE_DIM, HIGGS_TRAIN_SIZE,
                          HIGGS_FEATURE_MIN_VALUE, HIGGS_FEATURE_MAX_VALUE);
    normalize_data_values(test_feature_data, HIGGS_FEATURE_DIM, HIGGS_TEST_SIZE,
                          HIGGS_FEATURE_MIN_VALUE, HIGGS_FEATURE_MAX_VALUE);
  }

  printf("\n==========================================================\n");

  clock_t clustering_start = clock();

  if(clustering_algorithm == 1) { // LSH clustering
    if(data_set == 1) { // MNIST data set
      clusters = execute_LSH(data_set, MNIST_FEATURE_DIM, MNIST_TRAIN_SIZE, train_feature_data);
    }

    if(data_set == 2) { // BIO data set
      clusters = execute_LSH(data_set, BIO_FEATURE_DIM, BIO_TRAIN_SIZE, train_feature_data);
    }

    if(data_set == 3) { // HIGGS data set
      clusters = execute_LSH(data_set, HIGGS_FEATURE_DIM, HIGGS_TRAIN_SIZE, train_feature_data);
    }
  }

  if(clustering_algorithm == 2) { // KDTree clustering
    if(DEBUG) { printf("\nDEBUG -> Enter the desired number of clusters: "); }
    else { printf("\nEnter the desired number of clusters: "); }

    scanf("%d", &K);

    if(data_set == 1) { // MNIST data set
      execute_kdtree(MNIST_FEATURE_DIM, K, MNIST_TRAIN_SIZE, train_feature_data);
    }

    if(data_set == 2) { // BIO data set
      execute_kdtree(BIO_FEATURE_DIM, K, BIO_TRAIN_SIZE, train_feature_data);
    }

    if(data_set == 3) { // HIGGS data set
      execute_kdtree(HIGGS_FEATURE_DIM, K, HIGGS_TRAIN_SIZE, train_feature_data);
    }
  }

  if(clustering_algorithm == 3) { // BKmeans clustering
    if(DEBUG) { printf("\nDEBUG -> Enter the desired number of clusters: "); }
    else { printf("\nEnter the desired number of clusters: "); }
    scanf("%d", &K);

    if(data_set == 1) { // MNIST data set
      execute_bkmeans_j(data_set, MNIST_FEATURE_DIM, K,
                        MNIST_TRAIN_SIZE, train_feature_data,
                        MNIST_TEST_SIZE, test_feature_data,
                        train_non_feature_data, test_non_feature_data);
    }

    if(data_set == 2) { // BIO data set
      execute_bkmeans_j(data_set, BIO_FEATURE_DIM, K,
                        BIO_TRAIN_SIZE, train_feature_data,
                        BIO_TEST_SIZE, test_feature_data,
                        train_non_feature_data, test_non_feature_data);
    }

    if(data_set == 3) { // HIGGS data set
      execute_bkmeans_j(data_set, HIGGS_FEATURE_DIM, K,
                        HIGGS_TRAIN_SIZE, train_feature_data,
                        HIGGS_TEST_SIZE, test_feature_data,
                        train_non_feature_data, test_non_feature_data);
    }
  }

  clock_t clustering_end = clock();

  print_execution_time(clustering_start, clustering_end, "Clustering");

  printf("\n==========================================================\n");

  if(! DEBUG) {
    perform_search_queries(clustering_algorithm, data_set,
                           train_feature_data, test_feature_data,
                           train_non_feature_data, test_non_feature_data);
  }
}

Tree execute_LSH(int data_set, int dim, int train_size, double *train_feature_data)
{
  int i, j, correct_labeling_count = 0, cluster_count = 0;

  if(data_set == 1) { m = MNIST_m; w = MNIST_w; }   // MNIST data set
  if(data_set == 2) { m = BIO_m; w = BIO_w; } // BIO data set
  if(data_set == 3) { m = HIGGS_m; w = HIGGS_w; } // HIGGS data set

  b = malloc(m * sizeof(double));
  for(i = 0; i < m; i++) {
    b[i] = 0.0;
  }

  r = malloc(m * sizeof(double *));
  for(i = 0; i < m; i++) {
    r[i] = malloc(dim * sizeof(double));
  }

  for(i = 0; i < m; i++) {
    for(j = 0; j < dim; j++) {
      r[i][j] = gauss_rand();
    }
  }

  if(DEBUG) {
    return execute_debug_LSH();
  }
  else {
    printf("\nGenerating clusters via LSH...\n");
    Tree clusters = LSH(dim, train_size, train_feature_data, m, r, b, w);

    cluster_count = get_cluster_count(clusters);
    printf("\nLSH generated %d clusters.\n", cluster_count);

    return clusters;
  }
}

void execute_kdtree(int dim, int k, int train_size, double *train_feature_data)
{
  int i, j;

  cluster_assign = malloc(train_size * sizeof(cluster_assign));
  cluster_size = malloc(k * sizeof(cluster_size));
  cluster_start = malloc(k * sizeof(cluster_start));

  // Initialize cluster assignments
  for(i = 0; i < train_size; i++) {
    cluster_assign[i] = -1;
  }

  // Initialize cluster start and cluster size
  for(i = 0; i < k; i++) {
    cluster_size[i] = 0;
    cluster_start[i] = 0;
  }

  cluster_bdry = malloc(k * sizeof(double *));
  cluster_centroid = malloc(k * sizeof(double *));
  for(i = 0; i < k; i++) {
    cluster_bdry[i] = malloc((dim * 2) * sizeof(double));
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

  if(DEBUG) {
    debug_data = read_debug_data();

    printf("\nDEBUG -> Generating KDTree with %d clusters...\n", k);
    kdtree(debug_dim, debug_size, debug_data, k,
           cluster_size, cluster_start, cluster_bdry,
           cluster_centroid, cluster_assign);

    write_results(debug_dim, debug_size, debug_data, cluster_assign);
  }
  else {
    printf("\nGenerating KDTree with %d clusters...\n", k);
    kdtree(dim, train_size, train_feature_data, k,
           cluster_size, cluster_start, cluster_bdry,
           cluster_centroid, cluster_assign);
  }

  printf("\nKDTree clustering complete.\n");
}

void execute_bkmeans_j(int data_set, int dim, int k,
                       int train_size, double *train_feature_data,
                       int test_size, double *test_data,
                       int *train_non_feature_data, int *test_non_feature_data)
{
  int i, j, num_iterations = 0;

  cluster_assign = malloc(train_size * sizeof(double));
  cluster_size = malloc(k * sizeof(double));
  cluster_start = malloc(k * sizeof(double));

  // Initialize cluster assignments
  for (i = 0; i < train_size; i++) {
    cluster_assign[i] = -1;
  }

  // Initialize cluster start
  for (i = 0; i < k; i++) {
    cluster_size[i] = 0;
    cluster_start[i] = 0;
  }

  cluster_radius = malloc(k * sizeof(double));
  cluster_centroid = malloc(k * sizeof(double *));
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
    debug_data = read_debug_data();
    printf("\nDEBUG -> Forming %d clusters via Bisecting K-means_j...\n", k);
    num_iterations = bisecting_kmeans(debug_dim, debug_size, debug_data, k, cluster_size,
                                      cluster_start, cluster_radius, cluster_centroid,
                                      cluster_assign);

    write_results(debug_dim, debug_size, debug_data, cluster_assign);
  }
  else {
    printf("\nForming %d clusters via Bisecting K-means_j...\n", k);
    num_iterations = bisecting_kmeans(dim, train_size, train_feature_data, k, cluster_size,
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

void normalize_data_values(double *data, int feature_dimensions, int ndata,
                           int feature_min_value, int feature_max_value)
{
  int i;
  for(i = 0; i < ndata * feature_dimensions; i++) {
    data[i] = (data[i] - feature_min_value) / (feature_max_value - feature_min_value);
  }
}

void perform_search_queries(int clustering_algorithm, int data_set,
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

  printf("\n**********************************************************\n\nPerforming queries using test data...\n");

  clock_t query_start = clock();

  search_clusters(clustering_algorithm, data_set,
                  train_feature_data, test_feature_data,
                  train_non_feature_data, test_non_feature_data);

  clock_t query_end = clock();

  print_execution_time(query_start, query_end, "Test data query");

  printf("\n**********************************************************\n");
}

void search_clusters(int clustering_algorithm, int data_set,
                     double *train_feature_data, double *test_feature_data,
                     int *train_non_feature_data, int *test_non_feature_data)
{
  if(clustering_algorithm == 1) {
    if(data_set == 1) {
      LSH_search_clusters_for_approx_neighbors(clusters, MNIST_FEATURE_DIM, MNIST_TEST_SIZE,
                                               m, w, b, r,
                                               train_feature_data, test_feature_data,
                                               train_non_feature_data, test_non_feature_data);
    }
    if(data_set == 2) {
      LSH_search_clusters_for_approx_neighbors(clusters, BIO_FEATURE_DIM, BIO_TEST_SIZE,
                                               m, w, b, r,
                                               train_feature_data, test_feature_data,
                                               train_non_feature_data, test_non_feature_data);
    }
    if(data_set == 3) {
      LSH_search_clusters_for_approx_neighbors(clusters, HIGGS_FEATURE_DIM, HIGGS_TEST_SIZE,
                                               m, w, b, r,
                                               train_feature_data, test_feature_data,
                                               train_non_feature_data, test_non_feature_data);
    }
  }
  if(clustering_algorithm == 2) {
    if(data_set == 1) {
      kdtree_search_clusters_for_approx_neighbors(MNIST_FEATURE_DIM, MNIST_TEST_SIZE, K,
                                                  train_feature_data, test_feature_data,
                                                  train_non_feature_data, test_non_feature_data,
                                                  cluster_size, cluster_start, cluster_bdry);
    }
    if(data_set == 2) {
      kdtree_search_clusters_for_approx_neighbors(BIO_FEATURE_DIM, BIO_TEST_SIZE, K,
                                                  train_feature_data, test_feature_data,
                                                  train_non_feature_data, test_non_feature_data,
                                                  cluster_size, cluster_start, cluster_bdry);
    }
    if(data_set == 3) {
      kdtree_search_clusters_for_approx_neighbors(HIGGS_FEATURE_DIM, HIGGS_TEST_SIZE, K,
                                                  train_feature_data, test_feature_data,
                                                  train_non_feature_data, test_non_feature_data,
                                                  cluster_size, cluster_start, cluster_bdry);
    }
  }
  if(clustering_algorithm == 3) {
    if(data_set == 1) {
      bkmeans_search_clsuters_for_approx_neighbors(MNIST_FEATURE_DIM, MNIST_TEST_SIZE, K,
                                                   train_feature_data, test_feature_data,
                                                   train_non_feature_data, test_non_feature_data,
                                                   cluster_size, cluster_start,
                                                   cluster_radius, cluster_centroid);
    }
    if(data_set == 2) {
      bkmeans_search_clsuters_for_approx_neighbors(BIO_FEATURE_DIM, BIO_TEST_SIZE, K,
                                                   train_feature_data, test_feature_data,
                                                   train_non_feature_data, test_non_feature_data,
                                                   cluster_size, cluster_start,
                                                   cluster_radius, cluster_centroid);
    }
    if(data_set == 3) {
      bkmeans_search_clsuters_for_approx_neighbors(HIGGS_FEATURE_DIM, HIGGS_TEST_SIZE, K,
                                                   train_feature_data, test_feature_data,
                                                   train_non_feature_data, test_non_feature_data,
                                                   cluster_size, cluster_start,
                                                   cluster_radius, cluster_centroid);
    }
  }
}

void print_search_results(int test_size, double total_closest_neighbor_distance, double total_pts_searched)
{
  printf("\nQuery testing size = %d\n", test_size);
  printf("\nAverage distance to the approximate neighbor = %.1e\n",
         total_closest_neighbor_distance / (double) test_size);
  printf("\nAverage points searched per query = %.1lf\n", total_pts_searched / (double) test_size);
}

double randMToN(double M, double N)
{
  return M + (rand() / ( RAND_MAX / (N-M) ) ) ;
}

void write_results(int dim, int ndata, double *data, int *cluster_assign)
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

Tree execute_debug_LSH()
{
  int i, j, debug_m = 6; double debug_w = 18.0;

  double *debug_b = malloc(debug_m * sizeof(double));
  for(i = 0; i < debug_m; i++) {
    debug_b[i] = 0.0;
  }

  double **debug_r = malloc(debug_m * sizeof(double *));
  for(i = 0; i < debug_m; i++) {
    debug_r[i] = malloc(debug_dim * sizeof(double));
  }

  for(i = 0; i < debug_m; i++) {
    for(j = 0; j < debug_dim; j++) {
      debug_r[i][j] = gauss_rand();
    }
  }

  debug_data = read_debug_data();

  int *debug_cluster_assign = malloc(debug_size * sizeof(debug_cluster_assign));
  for(i = 0; i < debug_size; i++) {
    debug_cluster_assign[i] = -1;
  }

  printf("\nDEBUG -> Generating clusters via LSH...\n");
  Tree debug_clusters = LSH(debug_dim, debug_size, debug_data, debug_m, debug_r, debug_b, debug_w);

  int debug_cluster_count = get_cluster_count(debug_clusters);
  write_LSH_clusters_info(debug_clusters, debug_dim, debug_m, debug_w, debug_cluster_count);

  debug_LSH_generate_cluster_assign(debug_cluster_assign, debug_cluster_count);

  write_results(debug_dim, debug_size, debug_data, debug_cluster_assign);

  return debug_clusters;
}

double* read_debug_data()
{
  FILE *file = fopen("../debug_data.txt", "r");

  if(file == NULL) {
    perror("Error");
    exit(1);
  }

  char *data_number = malloc(256 * sizeof(char));
  double d;
  double *debug_data = malloc(debug_size * debug_dim * sizeof(double));

  int i;
  for(i = 0; i < debug_size * debug_dim; i++) {
    fscanf(file, "%s", data_number);
    debug_data[i] = strtod(data_number, NULL);
  }

  free(data_number);
  fclose(file);

  return debug_data;
}