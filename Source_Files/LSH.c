#include "../Header_Files/Headers.h"
#include "../Header_Files/Defs.h"
#include "../Header_Files/LSH.h"

cluster* LSH(int dim, int ndata, double *data,
             int m, double **r, double *b, double w,
             int *num_clusters)
{
  int i, j, pt = -1, *pt_hash = malloc(m * sizeof(int));
  double *data_pt_vector = malloc(dim * sizeof(double));

  cluster *clusters = malloc(sizeof(cluster));
  clusters->cluster_hash = malloc(m * sizeof(int));
  clusters->data_pts = malloc(sizeof(data_pt));
  clusters->data_pts->d_pt = -1;
  clusters->data_pts->next = NULL;
  clusters->next = NULL;

  for(i = 0; i < ndata; i++) {
    pt = i;
    for(j = 0; j < dim; j++) {
      data_pt_vector[j] = data[i * dim + j];
    }

    hash_pt(dim, data_pt_vector, m, r, b, w, pt_hash);

    clusters = add_pt_to_cluster(clusters, pt, pt_hash, m, num_clusters);
  }

  free(data_pt_vector); free(pt_hash);

  return clusters;
}

void hash_pt(int dim, double *data_pt_vector, int m, double **r, double *b, double w, int *pt_hash)
{
  int i, j;
  double *vector = malloc(dim * sizeof(double));

  for(i = 0; i < m; i++) {
    for(j = 0; j < dim; j++) {
      vector[j] = r[i][j];
    }

    pt_hash[i] = (int) floor((dot_product(dim, data_pt_vector, vector) - b[i]) / w);
  }

  free(vector);
}

int* hash_q_pt(int dim, double *q_pt_vector, int m, double **r, double *b, double w)
{
  int i, j;
  double *vector = malloc(dim * sizeof(double));
  int *q_pt_hash = malloc(m * sizeof(double));

  for(i = 0; i < m; i++) {
    for(j = 0; j < dim; j++) {
      vector[j] = r[i][j];
    }

    q_pt_hash[i] = (int) floor((dot_product(dim, q_pt_vector, vector) - b[i]) / w);
  }

  free(vector);

  return q_pt_hash;
}

double dot_product(int dim, const double *vector_a, const double *vector_b)
{
  int i;
  double result = 0.0;
  for(i = 0; i < dim; i++) {
    result += vector_a[i] * vector_b[i];
  }

  return result;
}

cluster* add_pt_to_cluster(cluster* clusters, int pt, int *pt_hash, int m, int *num_clusters)
{
  int j, k;
  bool matching_hash = true;
  cluster* current_cluster = clusters;

  if(pt == 0) { // First point, create first cluster with pt_hash for pt and add pt to cluster data_pts.
    for(j = 0; j < m; j++) {
      clusters->cluster_hash[j] = pt_hash[j];
    }
    clusters->data_pts->d_pt = pt;

    return clusters;
  }

  while(current_cluster != NULL) {
    for(j = 0; j < m; j++) {
      if(pt_hash[j] != current_cluster->cluster_hash[j]) {
        matching_hash = false;
      }
    }

    if(matching_hash) {
      data_pt *new_data_pt = malloc(sizeof(data_pt));
      new_data_pt->d_pt = pt;
      new_data_pt->next = current_cluster->data_pts;

      current_cluster->data_pts = new_data_pt;

      return clusters;
    }
    else {
      if(current_cluster->next == NULL) { // Reached end of cluster list without a match; create new cluster
        cluster *last_cluster = current_cluster;

        cluster* new_cluster = malloc(sizeof(cluster));
        new_cluster->cluster_hash = malloc(m * sizeof(int));
        for(k = 0; k < m; k++) {
          new_cluster->cluster_hash[k] = pt_hash[k];
        }
        new_cluster->data_pts = malloc(sizeof(data_pt));
        new_cluster->data_pts->d_pt = pt;
        new_cluster->data_pts->next = NULL;
        new_cluster->next = NULL;

        last_cluster->next = new_cluster;

        *num_clusters = *num_clusters + 1;

        return clusters;
      }
      else { // Have not reached the end of the cluster list; keep looking for a match
        matching_hash = true;
        current_cluster = current_cluster->next;
      }
    }
  }

  return clusters;
}

void search_clusters_for_apprx_neighbors(int dim, double *train_features, int *train_labels,
                                         int *test_labels, double *q_pt, int *q_pt_hash, int query_index,
                                         cluster *clusters, int m, int *correct_labeling_count)
{
  int i, pts_searched = 0, correct_label = 0;
  cluster *current_cluster = clusters;
  bool matching_hash = true;
  double closest_neighbor_dist = (double) INT_MAX, distance = 0.0;

  while(current_cluster != NULL) {
    for(i = 0; i < m; i++) {
      if(q_pt_hash[i] != current_cluster->cluster_hash[i]) {
        matching_hash = false;
      }
    }

    if(matching_hash) {
      int neighbor_data_pt = -1, closest_neighbor_pt = -1;
      data_pt *neighbors = current_cluster->data_pts;
      while(neighbors != NULL) {
        neighbor_data_pt = neighbors->d_pt;

        distance = calc_dist_to_neighbor(dim, train_features, q_pt, neighbor_data_pt);
        if(distance < closest_neighbor_dist) {
          closest_neighbor_dist = distance;
          closest_neighbor_pt = neighbor_data_pt;
        }

        neighbors = neighbors->next;
        pts_searched++;
      }

      if(train_labels[closest_neighbor_pt] == test_labels[query_index]) {
        correct_label++; *correct_labeling_count += correct_label;
      }

      return;
    }
    else if(current_cluster->next == NULL) {
      return;
    }
    else { // Keep iterating through clusters list to check the next hash value.
      matching_hash = true;
      current_cluster = current_cluster->next;
    }
  }
}

double calc_dist_to_neighbor(int dim, double *data, double *q_pt, int neighbor_data_pt)
{
  int i;
  double distance = 0.0;
  for(i = 0; i < dim; i++) {
    distance += (q_pt[i] - data[neighbor_data_pt * dim + i]) * (q_pt[i] - data[neighbor_data_pt * dim + i]);
  }

  return sqrt(distance);
}

double gauss_rand()
{
  static double U, V;
  static int phase = 0;
  double Z;

  if(phase == 0) {
    U = (rand() + 1.) / (RAND_MAX + 2.);
    V = rand() / (RAND_MAX + 1.);
    Z = sqrt(-2 * log(U)) * sin(2 * PI * V);
  } else
    Z = sqrt(-2 * log(U)) * cos(2 * PI * V);

  phase = 1 - phase;

  return Z;
}

void print_clusters_info(int m, cluster *clusters)
{
  int i = 1, j, cluster_count = 0;
  cluster *current_cluster = clusters;
  data_pt *data_pts = NULL;

  while(current_cluster != NULL) {
    data_pts = current_cluster->data_pts;
    printf("\nCluster %d hash\t=\t", i);
    for(j = 0; j < m; j++) {
      printf("%d\t", current_cluster->cluster_hash[j]);
    }

    printf("\nCluster points\t=\t");

    while(data_pts != NULL) {
      printf("%d\t", data_pts->d_pt);
      data_pts = data_pts->next;
    }

    printf("\n");

    current_cluster = current_cluster->next;
    i++;
    cluster_count++;
  }

  printf("\nTotal cluster count = %d\n\n", cluster_count);
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


