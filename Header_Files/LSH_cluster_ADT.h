#ifndef MASTERS_PROJECT_LSH_CLUSTER_ADT_H
#define MASTERS_PROJECT_LSH_CLUSTER_ADT_H

/**
 * This implementation of a balanced binary search tree is taken from the lecture notes given by
 * Dr. Abrahamson, East Carolina University. All credit goes to Dr. Abrahamson.
 */

struct Cluster;
typedef struct Cluster *Tree;

/**
 * Return the height of cluster tree C.
 */
int get_height(Tree T);

int max_height(int x, int y);

/**
 * Structure used for the data points in a cluster.
 */
typedef struct data_pt {
  int d_pt;
  struct data_pt *next;
} data_pt;

/**
 * Structure of a cluster. Each cluster node has a hash and associated data points.
 *
 * Clusters are structured as a balanced binary tree.
 */
typedef struct Cluster {
  int ht; // Height of this cluster node
  int *cluster_hash;
  data_pt *data_pts;
  struct Cluster *left;
  struct Cluster *right;
} Cluster;

/**
 * Sets C->ht to the height of C.
 *
 * C must not be empty.
 */
int set_height(Tree T);

/**
 * Perform a single rotation from right to left at the root of T.
 */
void single_rotate_left(Tree T);

/**
 * Performs a single rotation from left to right at the root of T.
 */
void single_rotate_right(Tree T);

/**
 * Performs a double rotation from right to left at the root of T.
 */
void double_rotate_left(Tree T);

/**
 * Performs a double rotation from left to right at the root of T.
 */
void double_rotate_right(Tree T);

/**
 * Performs a rotation from right to left at the root of T, with either a single or double rotation.
 */
void rotate_left(Tree T);

/**
 * Performs a rotation from left to right at the root of T, with either a single or double rotation.
 */
void rotate_right(Tree T);

/**
 * Rebalance by performing a rotation at T if required. Set the height field of T correctly (even if no
 * rotation is performed).
 */
void rebalance(Tree T);

/**
 * Inserts data point d_pt with d_pt_hash[hash_size] into the binary search tree. If d_pt_hash[hash_size] already
 * exists in the tree, d_pt is added to the cluster. If not, a new cluster node is added to T.
 * Rebalancing is done after insertion if necessary.
 */
void insert(int d_pt, int hash_size, int *d_pt_hash, Tree T);

/**
 * Return -1 if cluster_hash < cluster_node_hash
 * Return 1 if cluster_hash > cluster_node_hash
 * Return 0 if cluster_hash == cluster_node_hash
 */
int compare_hash(int hash_size, const int *cluster_hash, const int *cluster_node_hash)

#endif //MASTERS_PROJECT_LSH_CLUSTER_ADT_H
