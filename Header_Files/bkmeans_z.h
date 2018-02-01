#ifndef BKMEANS_H
#define BKMEANS_H

double calc_dist_square(int dim, double *datum1, double *datum2);

/*************************************************************************
 * array sizes:                                                          *
 *       cluster_assign[ndata], datum[dim], center[2*dim]                *
 * Input:                                                                *
 *       center[0] - centroid of the cluster                             *
 *       radius_pt[0] - farthest pt to the centroid                      *
 * Output:                                                               *
 *       radius_pt, center, start, size, ssd                             *
 * buffers: cluster_assign[], datum, cluster_center0[]                   *
 * radius_pt[2*dim]: - the radius pt of the cluster[k]                   *
 *************************************************************************/
int two_means(int iterat_limit, int dim, int i0, int im, double *data,   // input
              int *cluster_assign, double *datum, double *center0,       // buffers
              double *radius_pt, double *center, int start[2],
              int size[2], double ssd[2]); // output

int bkmeans_z(int iterat_limit, int kk, int dim, int i0_in, int im_in, double *data, // input
              int *cluster_assign, double *datum,                                    // buffers
              double *cluster_center, double *cluster_radius,                        // output
              int *cluster_start, int *cluster_size, double *cluster_ssd);           // output

/******************************************************************************
 kk :               number of clusters, i.e. the K in K-mean.
 cluster_center[kk*dim]: input  -- stores initial kk centers
                         output -- stores kk centers
 cluster_radius[kk]:output -- the radius of each output cluster
 cluster_start[kk]: output -- the index of the 1st in each cluster
 cluster_size[kk]:  output -- the num of datapoints in each cluster
 dataset_size and mem_capacity in unit of dim*sizeof(double),
                               i.e. in unit of data items
*******************************************************************************/
int kmeans_z(int iterat_limit, int kk, int dim, int i0, int im, double *data,   // input
             int *cluster_assign, double *datum, double *cluster_center0,       // buffers
             double *radius_pt, double *cluster_center, double *cluster_radius, // output
             int *cluster_start, int *cluster_size, double *cluster_ssd);       // output

#endif //BKMEANS_H
