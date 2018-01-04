#ifndef KDTREE_MEDIAN_H
#define KDTREE_MEDIAN_H

double calc_median(int size, double *data, int *cluster_assign);

int bipartition_median(int dim, int i0, int im, double *data,
                       int cluster_start[2], int cluster_size[2],
                       double *cluster_bdry[2], double *cluster_centroid[2],
                       int *cluster_assign, double *datum, double *buf);

int kdtree_hybrid(int dim, int ndata, double *data, int kk,
                  int *cluster_start, int *cluster_size,
                  double **cluster_bdry, double **cluster_centroid,
                  int *cluster_assign, double *datum, double *buf);

int search_kdtree_hybrid(int dim, int ndata, double *data, int kk,
                         int *cluster_start, int *cluster_size, double **cluster_bdry,
                         double *query, double *result_pt);

#endif //KDTREE_MEDIAN_H
