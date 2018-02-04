#include "../Header_Files/Headers.h"
#include "../Header_Files/Defs.h"
#include "../Header_Files/LSH.h"

Tree LSH(int dim, int ndata, const double *data,
             int m, double **r, double *b, double w)
{
  int i, j, *pt_hash = malloc(m * sizeof(int));
  double *data_pt_vector = malloc(dim * sizeof(double));

  Tree clusters = NULL;

  for(i = 0; i < ndata; i++) {
    for(j = 0; j < dim; j++) {
      data_pt_vector[j] = data[i * dim + j];
    }

    hash_pt(dim, data_pt_vector, m, r, b, w, pt_hash);

    clusters = add_pt_to_cluster(clusters, i, pt_hash, m);
  }

  free(data_pt_vector); free(pt_hash);

  return clusters;
}

void hash_pt(int dim, double *data_pt_vector, int m, double **r, const double *b, double w, int *pt_hash)
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

Tree add_pt_to_cluster(Tree clusters, int pt, int *pt_hash, int m)
{
  clusters = insert(clusters, pt, m, pt_hash);

  return clusters;
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

int* hash_q_pt(int dim, double *q_pt_vector, int m, double **r, const double *b, double w)
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

int search_clusters_for_apprx_neighbors(int dim, double *train_features, int *train_labels,
                                         int *test_labels, double *q_pt, int *q_pt_hash,
                                         int test_query_index, Tree clusters, int m,
                                         int *correct_labeling_count)
{
  int i, pts_searched = 0, correct_label = 0, closest_neighbor_pt = -1;
  double closest_neighbor_dist = (double) INT_MAX, current_neighbor_pt_distance = 0.0;

  Data_pt *neighbors = find_neighbors(clusters, m, q_pt_hash);

  while(neighbors != NULL) {
    current_neighbor_pt_distance = calc_dist_to_neighbor(dim, train_features, q_pt, neighbors->d_pt);

    if(current_neighbor_pt_distance < closest_neighbor_dist) {
      closest_neighbor_dist = current_neighbor_pt_distance;
      closest_neighbor_pt = neighbors->d_pt;
    }

    neighbors = neighbors->next;
    pts_searched++;
  }

  if(train_labels[closest_neighbor_pt] == test_labels[test_query_index]) { *correct_labeling_count += 1; }

  return pts_searched;
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
