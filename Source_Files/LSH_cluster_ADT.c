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

void single_rotate_left(Tree T)
{
  Tree R = T->right;
  T->right = R->left;
  set_height(T);

  R->left = T;
  T = R;
  set_height(T);
}

void single_rotate_right(Tree T)
{
  Tree L = T->left;
  T->left = L->right;
  set_height(T);

  L->right = T;
  T = L;
  set_height(T);
}

void double_rotate_left(Tree T)
{
  single_rotate_right(T->right);
  single_rotate_left(T);
}

void double_rotate_right(Tree T)
{
  single_rotate_left(T->left);
  single_rotate_right(T);
}

void rotate_left(Tree T)
{
  Tree R = T->right;
  int zag = get_height(R->left);
  int zig = get_height(R->right);

  if(zig > zag) {
    single_rotate_left(T);
  }
  else {
    double_rotate_left(T);
  }
}

void rotate_right(Tree T)
{
  Tree L = T->left;
  int zig = get_height(L->left);
  int zag = get_height(L->right);

  if(zig > zag) {
    single_rotate_right(T);
  }
  else {
    double_rotate_right(T);
  }
}

void rebalance(Tree T)
{
  int h_left = get_height(T->left);
  int h_right = get_height(T->right);

  if(h_right > h_left + 1) {
    rotate_left(T);
  }
  else if(h_left > h_right + 1) {
    rotate_right(T);
  }
  else {
    set_height(T);
  }
}

Tree insert(Tree T, int d_pt, int hash_size, int *d_pt_hash)
{
  if(T == NULL) {
    T = malloc(sizeof(Cluster));
    T->cluster_hash = malloc(hash_size*sizeof(int));

    int i;
    for(i = 0; i < hash_size; i++) {
      T->cluster_hash[i] = d_pt_hash[i];
    }

    T->data_pts = malloc(sizeof(data_pt));
    T->data_pts->d_pt = d_pt;
    T->data_pts->next = NULL;
    T->left = NULL;
    T->right = NULL;
    T->ht = 1 + max_height(get_height(T->left), get_height(T->right));

    return T;
  }
  else if(compare_hash(hash_size, d_pt_hash, T->cluster_hash) < 0) {
    T->left = insert(T->left, d_pt, hash_size, d_pt_hash);
    rebalance(T);
  }
  else if(compare_hash(hash_size, d_pt_hash, T->cluster_hash) > 0) {
    T->right = insert(T->right, d_pt, hash_size, d_pt_hash);
    rebalance(T);
  }
  else { // equal hash values, add d_pt to cluster node
    if(T->data_pts == NULL) {
      T->data_pts = malloc(sizeof(data_pt));
      T->data_pts->d_pt = d_pt;
      T->data_pts->next = NULL;
    }
    else { // Add d_pt to front of data points list
      data_pt *new_data_pt = malloc(sizeof(data_pt));
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

int get_cluster_count(Tree T)
{
  if(T == NULL) { return 0; }
  else { return 1 + get_cluster_count(T->left) + get_cluster_count(T->right); }
}

void write_LSH_clusters_info(Tree T, int dim, int hash_size, double w, int cluster_count)
{
  int i;
  FILE *file;

  file = fopen("../LSH_Cluster_Data_Info/clusters.dat", "w");
  fprintf(file, "cluster count = %d\n", cluster_count);
  fprintf(file, "hash size = %d\n", hash_size);
  fprintf(file, "w = %.01lf\n\n", w);

  write_cluster_node_info(file, T, hash_size);

  fclose(file);
}

void write_cluster_node_info(FILE *file, Tree T, int hash_size)
{
  if(T == NULL) { return; }

  write_cluster_node_info(file, T->left, hash_size);

  // Print cluster node info (hash value of cluster node and data points at cluster node
  int i;
  fprintf(file, "Cluster node hash: ");
  for(i = 0; i < hash_size; i++) {
    fprintf(file, "%4d", T->cluster_hash[i]);
  }
  fprintf(file, "\n");

  fprintf(file, "Data points: ");
  data_pt *tmp = T->data_pts;
  while(tmp != NULL) {
    fprintf(file, "%6d", tmp->d_pt);

    tmp = tmp->next;
  }
  fprintf(file, "\n*******************************************************************\n\n");

  write_cluster_node_info(file, T->right, hash_size);
}
