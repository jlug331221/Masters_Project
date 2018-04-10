#include "../Header_Files/Headers.h"
#include "../Header_Files/LSH_cluster_ADT.h"

int get_height(Tree T)
{
  if(T == NULL) {
    return 0;
  }
  else {
    return T->ht;
  }
}

int max_height(int x, int y)
{
  if(x > y) { return x; }

  return y;
}

int set_height(Tree T)
{
  T->ht = 1 + max_height(get_height(T->left), get_height(T->right));
}

Tree single_rotate_left(Tree T)
{
  Tree R = T->right;
  T->right = R->left;
  set_height(T);

  R->left = T;
  T = R;
  set_height(T);

  return T;
}

Tree single_rotate_right(Tree T)
{
  Tree L = T->left;
  T->left = L->right;
  set_height(T);

  L->right = T;
  T = L;
  set_height(T);

  return T;
}

Tree double_rotate_left(Tree T)
{
  T->right = single_rotate_right(T->right);
  T = single_rotate_left(T);

  return T;
}

Tree double_rotate_right(Tree T)
{
  T->left = single_rotate_left(T->left);
  T = single_rotate_right(T);

  return T;
}

Tree rotate_left(Tree T)
{
  Tree R = T->right;
  int zag = get_height(R->left);
  int zig = get_height(R->right);

  if(zig > zag) {
    T = single_rotate_left(T);
  }
  else {
    T = double_rotate_left(T);
  }

  return T;
}

Tree rotate_right(Tree T)
{
  Tree L = T->left;
  int zig = get_height(L->left);
  int zag = get_height(L->right);

  if(zig > zag) {
    T = single_rotate_right(T);
  }
  else {
    T = double_rotate_right(T);
  }

  return T;
}

Tree rebalance(Tree T)
{
  int h_left = get_height(T->left);
  int h_right = get_height(T->right);

  if(h_right > h_left + 1) {
    T = rotate_left(T);
  }
  else if(h_left > h_right + 1) {
    T = rotate_right(T);
  }
  else {
    set_height(T);
  }

  return T;
}

Tree insert(Tree T, int d_pt, int hash_size, int *d_pt_hash)
{
  if(T == NULL) {
    T = malloc(sizeof(Cluster));
    T->cluster_hash = malloc(hash_size * sizeof(int));

    int i;
    for(i = 0; i < hash_size; i++) {
      T->cluster_hash[i] = d_pt_hash[i];
    }

    T->data_pts = malloc(sizeof(Data_pt));
    T->data_pts->d_pt = d_pt;
    T->data_pts->next = NULL;
    T->left = NULL;
    T->right = NULL;
    T->ht = 1 + max_height(get_height(T->left), get_height(T->right));

    return T;
  }
  else if(compare_hash(hash_size, d_pt_hash, T->cluster_hash) < 0) {
    T->left = insert(T->left, d_pt, hash_size, d_pt_hash);
    T = rebalance(T);
    return T;
  }
  else if(compare_hash(hash_size, d_pt_hash, T->cluster_hash) > 0) {
    T->right = insert(T->right, d_pt, hash_size, d_pt_hash);
    T = rebalance(T);
    return T;
  }
  else { // equal hash values, add d_pt to cluster node
    if(T->data_pts == NULL) {
      T->data_pts = malloc(sizeof(Data_pt));
      T->data_pts->d_pt = d_pt;
      T->data_pts->next = NULL;
    }
    else { // Add d_pt to front of data points list
      Data_pt *new_data_pt = malloc(sizeof(Data_pt));
      new_data_pt->d_pt = d_pt;
      new_data_pt->next = T->data_pts;

      T->data_pts = new_data_pt;
    }

    return T;
  }
}

int compare_hash(int hash_size, const int *d_pt_hash, const int *cluster_node_hash)
{
  int i;
  for(i = 0; i < hash_size; i++) {
    if(d_pt_hash[i] < cluster_node_hash[i]) { return -1; }
    if(d_pt_hash[i] > cluster_node_hash[i]) { return 1; }
  }

  return 0;
}

Data_pt* find_neighbors(Tree T, int hash_size, int *q_pt_hash)
{
  if(T == NULL) { return NULL; }

  if(compare_hash(hash_size, q_pt_hash, T->cluster_hash) == 0) {
    return T->data_pts;
  }

  if(compare_hash(hash_size, q_pt_hash, T->cluster_hash) < 0) {
    return find_neighbors(T->left, hash_size, q_pt_hash);
  }

  if(compare_hash(hash_size, q_pt_hash, T->cluster_hash) > 0) {
    return find_neighbors(T->right, hash_size, q_pt_hash);
  }
}

int get_cluster_count(Tree T)
{
  if(T == NULL) { return 0; }
  else { return 1 + get_cluster_count(T->left) + get_cluster_count(T->right); }
}

void write_LSH_clusters_info(Tree T, int dim, int hash_size, double w, int cluster_count)
{
  int i;
  FILE *file1, *file2;

  file1 = fopen("../cmake-build-debug/LSH_clusters.dat", "w");
  fprintf(file1, "cluster count = %d\n", cluster_count);
  fprintf(file1, "hash size = %d\n", hash_size);
  fprintf(file1, "w = %.01lf\n\n", w);

  file2 = fopen("../cmake-build-debug/LSH_cluster_assign.dat", "w");

  write_cluster_node_info(file1, file2, T, hash_size);

  fclose(file1);
  fclose(file2);
}

void write_cluster_node_info(FILE *file1, FILE *file2, Tree T, int hash_size)
{
  if(T == NULL) { return; }

  write_cluster_node_info(file1, file2, T->left, hash_size);

  // Print cluster node info (hash value of cluster node and data points at cluster node
  int i;
  fprintf(file1, "Cluster node hash: ");
  for(i = 0; i < hash_size; i++) {
    fprintf(file1, "%4d", T->cluster_hash[i]);
  }
  fprintf(file1, "\n");

  fprintf(file1, "Data points: ");
  Data_pt *tmp = T->data_pts;
  while(tmp != NULL) {
    fprintf(file1, "%6d", tmp->d_pt);

    fprintf(file2, "%d ", tmp->d_pt);

    tmp = tmp->next;
  }
  fprintf(file1, "\n*******************************************************************\n\n");
  fprintf(file2, "-1\n");

  write_cluster_node_info(file1, file2, T->right, hash_size);
}

void debug_LSH_generate_cluster_assign(int *debug_cluster_assign, int debug_cluster_count)
{
  FILE *file = fopen("../cmake-build-debug/LSH_cluster_assign.dat", "r");

  if(file == NULL) {
    perror("Error");
    exit(1);
  }

  int cluster = 0;
  char *data_pt = malloc(7 * sizeof(char));

  while(cluster < debug_cluster_count) {
    fscanf(file, "%s", data_pt);
    while((int) strtod(data_pt, NULL) != -1) {
      debug_cluster_assign[(int) strtod(data_pt, NULL)] = cluster;

      fscanf(file, "%s", data_pt);
    }

    cluster++;
  }

  fclose(file);
}

void verify_data_pts_clustered(Tree T, int *data_pts, int ndata) {
  Data_pt *tmp = T->data_pts;
  int i;

  while(tmp != NULL) {
    data_pts[tmp->d_pt] = tmp->d_pt;

    tmp = tmp->next;
  }

  if(T->left != NULL) { verify_data_pts_clustered(T->left, data_pts, ndata); }
  if(T->right != NULL) { verify_data_pts_clustered(T->right, data_pts, ndata); }
}
