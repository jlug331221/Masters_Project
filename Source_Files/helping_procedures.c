#include "../Header_Files/Headers.h"
#include "../Header_Files/Defs.h"
#include "../Header_Files/helping_procedures.h"
#include "../Header_Files/kdtree.h"
#include "../Header_Files/kmeans.h"
#include "../Header_Files/LSH.h"
#include "../Header_Files/kdtree_median.h"
#include "../Header_Files/bkmeans.h"

double randMToN(double M, double N){
  return M + (rand() / ( RAND_MAX / (N-M) ) ) ;
}

void writeResults(int dim, int ndata, double* data, int* cluster_assign)
{
  int i;
  FILE* file;

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

void execute_kdtree(int *train_labels, double *train_features, int *test_labels, double *test_features)
{
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
  kdtree(dim, ndata, train_features, train_labels, k,
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

    search_kdtree(dim, ndata, train_features, train_labels, test_labels,
                  k, i, cluster_size, cluster_start, cluster_bdry,
                  query, &correct_labeling_count);

    free(query);
    query = malloc(dim * sizeof(double));
  }

  printf("Accuracy of labeling = %.2f%%\n", ((double) correct_labeling_count / (double) TEST_SIZE) * 100);
}

void execute_kmeans(int *train_labels, double *train_features, int *test_labels, double *test_features)
{
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

void execute_LSH(int *train_labels, double *train_features, int *test_labels, double *test_features)
{
  int dim = FEATURE_DIM, ndata = TRAIN_SIZE, m = 3, i, j, z = 0;
  int *num_clusters = &z; int correct_labeling_count = 0;
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

  printf("\nGenerating clusters via LSH...\n");
  cluster *clusters = LSH(dim, ndata, train_features, m, r, b, w, num_clusters);

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

void execute_kdtree_median(int *train_labels, double *train_features, int *test_labels, double *test_features)
{
  int ndata = 32, dim = 2, kk = 2, i, j;

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

  printf("\nBuilding KDTree Hybrid (median)...\n");
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

void execute_bkmeans(int *train_labels, double *train_features, int *test_labels, double *test_features)
{
  int ndata = 1000, dim = 2, kk = 32, i, num_clusters = 0;

  double *data = malloc(ndata * dim * sizeof(double));
  for(i = 0; i < ndata * dim; i++) {
    data[i] = randMToN(0, 100);
  }

  int *cluster_size = malloc(kk * sizeof(double));
  int *cluster_start = malloc(kk * sizeof(double));
  int *cluster_assign = malloc(ndata * sizeof(double));

  // Initialize cluster assignments
  for(i = 0; i < ndata; i++) {
    cluster_assign[i] = -1;
  }

  // Initialize cluster start and cluster size
  for(i = 0; i < kk; i++) {
    cluster_start[i] = 0;
    cluster_size[i] = 0;
  }

  double *cluster_radius = malloc(kk * sizeof(double));
  double *cluster_center = malloc(kk * sizeof(double));
  for(i = 0; i < kk; i++) {
    cluster_center[i] = 0.0;
    cluster_radius[i] = 0.0;
  }

  double *datum = malloc(dim * sizeof(double));

  double *cluster_ssd = malloc(kk * sizeof(double));

  printf("\nForming clusters...\n\n");
  num_clusters = bkmeans(5, kk, dim, ndata, 0, ndata, data,
                         cluster_assign, datum,
                         cluster_center, cluster_radius,
                         cluster_start, cluster_size, cluster_ssd);

  writeResults(dim, ndata, data, cluster_assign);

  printf("Number of clusters = %d\n", num_clusters);
}

void read_binary_dataset(char *path, int size, int *labels, double *features)
{
  FILE* file = fopen(path, "rb");

  if(file == NULL) {
    perror("Error");
    exit(1);
  }

  double label = -1.0;

  int i;
  for (i = 0; i < size; i++) {
    // Read the label
    fread(&label, sizeof(double), 1, file);
    labels[i] = (int)label;

    // Read the features
    fread(&features[i * FEATURE_DIM], sizeof(double), FEATURE_DIM, file);
  }

  fclose(file);
}

void normalize_data(double *data, int feature_dimensions, int ndata)
{
  int i;
  for(i = 0; i < ndata*feature_dimensions; i++) {
    data[i] = (data[i] - FEATURE_MIN_VALUE) / (FEATURE_MAX_VALUE - FEATURE_MIN_VALUE);
  }
}
