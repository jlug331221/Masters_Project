#include "../Header_Files/Headers.h"
#include "../Header_Files/Defs.h"
#include "../Header_Files/kdtree_median.h"
#include "../Header_Files/bkmeans.h"

double calc_median(int size, double *data, int *cluster_assign)
{
  int i, j , iteration, i0, i1, i2,
      start0, start1, start2, end0, end1, end2, size0, size1, size2 ;
  double median, tmp, tmp1, tmp2, end_pt, end_pt0 ;

  median = data[0] ;
  for(i=1; i<size; i++) median += data[i];
  median *= 1.0/(double)size ; /* median is centroid now */

  size0=0; size1=0; size2=0;
  for(i=0; i<size; i++) {
    if( data[i] < median - TINY_DELTA) {cluster_assign[i] =0; size0++ ;}
    else if(data[i]>median+TINY_DELTA) {cluster_assign[i] =2; size2++ ;}
    else {cluster_assign[i] = 1 ; size1++ ;}
  }
  if((size0+size1+size2)!=size) {
    printf("size0+size1+size2)!=size\n"); return 0.0 ;
  }

  start0=0;       end0 = start0+size0;
  start1=end0;    end1 = size ;
  while(start0 < end0) {
    while((start0<end0) && (cluster_assign[start0]==0)) start0++ ;
    while((start1<end1) && (cluster_assign[start1]!=0)) start1++ ;
    if( (start0<end0) && (start1<end1) ) {
      tmp=data[start0]; data[start0]=data[start1]; data[start1]=tmp;
      start0++;   start1++;
    }
    else if(start0<end0) {
      printf("Error: start0<end0 && start1>=end1\n");   return 0.0;
    }
    else if( start1<end1 ) {
      printf("Error: start0>=end0 && start1<end1\n");   return 0.0;
    }
    else {;}
  }

  start1=size0;   end1 = start1+size1;
  start2=end1;    end2 = start2+size2;
  while(start1 < end1) {
    while((start1<end1) && (cluster_assign[start1]==1)) start1++ ;
    while((start2<end2) && (cluster_assign[start2]==2)) start2++ ;
    if( (start1<end1) && (start2<end2) ) {
      tmp=data[start1]; data[start1]=data[start2]; data[start2]=tmp;
      start1++;   start2++;
    }
    else if(start1<end1) {
      printf("Error: start1<end1 && start2>=end2\n");   return 0.0;
    }
    else if( start2<end2 ) {
      printf("Error: start1>=end1 && start2<end2\n");   return 0.0;
    }
    else {;}
  }

  if( size0+size1+1 < size2 ) {
    end_pt0 = 0.0;
    for(i=size0; i<size; i++) end_pt0 += data[i];
    end_pt0 *= 1.0/(size1+size2); /* end_pt0 is centroid of cluster2 */
    tmp = 1.0*((double)size2-(double)size0)/(double)size ;
    end_pt = median ;
    median = 0.75*end_pt + 0.25*end_pt0 ;
  }
  else if(size0 > size1+size2+1) {
    end_pt0 = 0.0;
    for(i=0; i<size0+size1; i++) end_pt0 += data[i];
    end_pt0 *= 1.0/(size0+size1); /* end_pt0 is centroid of cluster 0*/
    tmp = 1.0*((double)size0-(double)size2)/(double)size ;
    end_pt = median ;
    median = 0.75*end_pt + 0.25*end_pt0 ;
  }
  else {
    return median ;
  }

  iteration = 1 ;
  while( (size0+size1<size2) || (size0 > size1+size2) ) {
    iteration++ ;
    size0=0; size1=0; size2=0;

    for(i=0; i<size; i++) {
      if( data[i] < median - TINY_DELTA) {cluster_assign[i] =0; size0++ ;}
      else if(data[i]>median+TINY_DELTA) {cluster_assign[i] =2; size2++ ;}
      else {cluster_assign[i] = 1 ; size1++ ;}
    }
    if((size0+size1+size2)!=size) {
      printf("Error in sizes, %d, %d, %d, %d\n", size0, size1, size2, size);
      return 0 ;
    }

    if(size0+size1<size2) {
      end_pt0 = max(end_pt0, end_pt) ;
      end_pt = median ;
      median = 0.5*end_pt + 0.5*end_pt0 ;
    }
    else if(size0 > size1+size2) {
      end_pt0 = min(end_pt0, end_pt) ;
      end_pt = median ;
      median = 0.5*end_pt + 0.5*end_pt0 ;
    }
    else {
      return median ; }
  }
}

/*********************************************************************
 * cluster_bdry: cluster_bdry[2][2*dim], for each of the two clusters.
 *           AS output, for each cluster, each dim has
 *           a min and max bdry for each cluster.
 *   Sizes of arrays:
 *         buf[ndata], datum[dim], cluster_centroid[2][dim]
 *   Input:
 *         cluster_bdry[0] contains the bdry of the input dataset
 *         cluster_centroid[0] contains the centroid of the input set
 *   Output:
 *         cluster_start[], cluster_size[], cluster_bdry[], cluster_centroid[]
 *   Buffers:
 *         datum[], buf[], cluster_assign[]
 ***********************************************************************/
int bipartition_median(int dim, int i0, int im, double *data,
                       int cluster_start[2], int cluster_size[2],
                       double *cluster_bdry[2], double *cluster_centroid[2],
                       int *cluster_assign, double *datum, double *buf)

{
  int i, j, k, j_max, membership, size, start0, end0, start1, end1 ;
  double tmp, var_max, partition_pt, val_max, val_min ;

/* partition along dimension with largest variance */
  for(j=0; j<dim; j++) datum[j]= 0.0 ;
  for(i=i0; i<im; i++) {
    for(j=0; j<dim; j++) {
      tmp = data[dim*i+j] - cluster_centroid[0][j] ;
      datum[j] += tmp*tmp ;
    }
  }
  var_max=datum[0];  j_max=0;
  for(j=1; j<dim; j++)
    if(datum[j] > var_max) {var_max=datum[j]; j_max=j;}


/* partition at an approximate median along dimension j_max */
  for(i=i0; i<im; i+=8) buf[(i-i0)/8] = data[dim*i+j_max] ;
  size=(im-i0)/8 ;
  partition_pt = calc_median(size, buf, cluster_assign) ;

  cluster_size[0] = 0;   cluster_size[1] = 0;
  for(j=0; j<dim; j++) {
    cluster_centroid[0][j]=0.0;     cluster_centroid[1][j]=0.0;
  }
  for(i=i0; i<im; i++) {
    for(j=0; j<dim; j++) datum[j] = data[dim*i+j] ;
    membership = ( (datum[j_max] < partition_pt) ? 0 : 1 ) ;
    for(j=0; j<dim; j++) cluster_centroid[membership][j] += datum[j] ;
    cluster_assign[i] = membership ;
    cluster_size[membership]++ ;
  }
  for(k=0; k<2; k++) {
    if( cluster_size[k]>0 )
      for(j=0;j<dim;j++) cluster_centroid[k][j]/=(double)cluster_size[k] ;
  }
  memcpy(cluster_bdry[1], cluster_bdry[0], 2*dim*sizeof(double)) ;
  cluster_bdry[0][2*j_max+1]= partition_pt ;
  cluster_bdry[1][2*j_max]  = partition_pt ;

/*** Data re-ordering ***/
  cluster_start[0] = i0 ;    cluster_start[1] = i0+cluster_size[0];
  start0=cluster_start[0];   end0 = start0+cluster_size[0] ;
  start1=cluster_start[1];   end1 = start1+cluster_size[1] ;
  while(start0 < end0) {
    while((start0<end0) && (cluster_assign[start0]==0)) start0++ ;
    while((start1<end1) && (cluster_assign[start1]==1)) start1++ ;
    if( (start0<end0) && (start1<end1) ) {
      memcpy(datum, data+start0*dim, dim*sizeof(double)) ;
      memcpy(data+start0*dim, data+start1*dim, dim*sizeof(double)) ;
      memcpy(data+start1*dim, datum, dim*sizeof(double)) ;
      start0++;   start1++;
    }
    else if( (start0<end0)&&(start1==end1) ) {
      printf("Error: start0<end0 && start1==end1\n");   return 0;
    }
    else if( (start0==end0)&&(start1<end1) ) {
      printf("Error: start0==end0 && start1<end1\n");   return 0;
    }
    else {;}
  }
/*** End of data re-ordering ***/

  return 1;
} /****** End of bipartition_median ******/

int kdtree_hybrid(int dim, int ndata, double *data, int kk,
           int *cluster_start, int *cluster_size,
           double **cluster_bdry, double **cluster_centroid,
           int *cluster_assign, double *datum, double *buf)
/* buf size: buf[ndata]     */
{
  int i0, im, i, j, k, k_chosen, level, width, level_kid, width_kid,
      start[2], size[2] ;
  double tmp, dist_max, dist,  *bdry[2], *centroid[2] ;

  bdry[0] = (double *)calloc((2*dim), sizeof(double));
  bdry[1] = (double *)calloc((2*dim), sizeof(double));
  centroid[0] = (double *) calloc((dim),  sizeof(double));
  centroid[1] = (double *) calloc((dim),  sizeof(double));

  for(j=0; j<dim; j++) {
    centroid[0][j] = 0.0 ;
    bdry[0][2*j]  = -INFINITY ;
    bdry[0][2*j+1]=  INFINITY ;
  }
  for(i=0; i<ndata; i++)
    for(j=0; j<dim; j++) centroid[0][j] += data[dim*i+j] ;
  for(j=0; j<dim; j++) centroid[0][j] /= (double)ndata ;

  i0 = 0 ;
  im = ndata ;
  bipartition_median(dim, i0, im, data, start, size, bdry, centroid,
                     cluster_assign, datum, buf);
  cluster_start[0]= start[0];   cluster_start[kk/2]= start[1];
  cluster_size[0] = size[0] ;   cluster_size[kk/2] = size[1] ;
  memcpy(cluster_bdry[0],      bdry[0],   2*dim*sizeof(double)) ;
  memcpy(cluster_bdry[kk/2], bdry[1], 2*dim*sizeof(double)) ;
  memcpy(cluster_centroid[0],   centroid[0],  dim*sizeof(double)) ;
  memcpy(cluster_centroid[kk/2], centroid[1], dim*sizeof(double)) ;

  level = 2;         width=kk/level;
  level_kid=2*level; width_kid = kk/level_kid ;
  while( level < kk ) {
    for(k_chosen=0; k_chosen<kk; k_chosen+=width) {
      k=k_chosen;
      i0 = cluster_start[k_chosen];
      im = i0+cluster_size[k_chosen];
      for(j=0; j<dim+dim; j++) bdry[0][j] = cluster_bdry[k][j];
      for(j=0; j<dim; j++) centroid[0][j] = cluster_centroid[k][j];
      bipartition_median(dim, i0, im, data, start, size, bdry, centroid,
                  cluster_assign, datum, buf) ;
      cluster_start[k]= start[0];   cluster_start[k+width_kid]= start[1];
      cluster_size[k] = size[0] ;   cluster_size[k+width_kid] = size[1] ;
      memcpy(cluster_bdry[k],      bdry[0],   2*dim*sizeof(double)) ;
      memcpy(cluster_bdry[k+width_kid], bdry[1], 2*dim*sizeof(double)) ;
      memcpy(cluster_centroid[k],   centroid[0],  dim*sizeof(double)) ;
      memcpy(cluster_centroid[k+width_kid], centroid[1], dim*sizeof(double)) ;
      /* End of k=k_chosen */
    }
    level *= 2;        width=kk/level;
    level_kid=2*level; width_kid = kk/level_kid ;
  } /************ End of while() loop ************/


  for(k=0; k<level; k++) {
    if(cluster_size[k] == 0)
      printf("In kd-tree: cluster %d empty!\n");
  }

  free(bdry[0]);  free(bdry[1]);  free(centroid[0]);  free(centroid[1]);
  return 1 ;
} /****** End of kdtree_hybrid ******/

int search_kdtree_hybrid(int dim, int ndata, double *data, int kk,
                         int *cluster_start, int *cluster_size, double **cluster_bdry,
                         double *query, double *result_pt)
/* This function returns the number of data points checked */
{
  int    i, j, k, k_min, n_pts, istart, iend, i_min ;
  double tmp, tmp1, tmp2, temp, *datum, *d_cluster2query, d_c2q_min, dist_min ;
  datum = (double *) calloc( dim, sizeof(double) ) ;
  d_cluster2query = (double *) calloc( kk, sizeof(double) ) ;

  for(k=0; k<kk; k++) {
    d_cluster2query[k] = 0.0 ;
    for(j=0; j<dim; j++) {
      tmp = query[j] ;
      if( tmp<cluster_bdry[k][2*j] ) {
        temp = cluster_bdry[k][2*j] - tmp ;
        d_cluster2query[k] += temp*temp ;
      }
      else if( tmp>cluster_bdry[k][2*j+1] ) {
        temp = tmp - cluster_bdry[k][2*j+1] ;
        d_cluster2query[k] += temp*temp ;
      }
      else {;}
    }
  }

  d_c2q_min = d_cluster2query[0] ;   k_min = 0 ;
  for(k=1; k<kk; k++) /* Find the closest cluster */
    if(d_cluster2query[k] < d_c2q_min) {d_c2q_min=d_cluster2query[k]; k_min=k;}

  if(d_c2q_min==0.0) {
    temp = INFINITY ;
    for(j=0; j<dim; j++) {
      tmp1 = cluster_bdry[k_min][2*j] ;
      tmp2 = cluster_bdry[k_min][2*j+1] ;
      tmp = min((tmp2-query[j]),(query[j]-tmp1));
      if(tmp < temp) temp = tmp ; /* temp: min distance to a bdry */
    }
  }
  istart = cluster_start[k_min];  iend = istart+cluster_size[k_min];
  dist_min = INFINITY ;
  for(i=istart; i<iend; i++) { /* loop thru all data in cluster k_min */
    for(j=0; j<dim; j++) datum[j] = data[i*dim+j] ;
    tmp = calc_dist_square(dim, query, datum) ;
    if(tmp<dist_min) { dist_min=tmp; i_min=i; }
  }
  n_pts = cluster_size[k_min] ;

  if(sqrt(dist_min) <= temp) { /* closest point found */
    for(j=0; j<dim; j++) result_pt[j] = data[i_min*dim+j] ;
  }
  else {
    for(k=0; k<k_min; k++) {
      if(d_cluster2query[k]<dist_min) {
        istart = cluster_start[k];  iend = istart+cluster_size[k];
        for(i=istart; i<iend; i++) {
          for(j=0; j<dim; j++) datum[j] = data[i*dim+j] ;
          tmp = calc_dist_square(dim, query, datum) ;
          if(tmp<dist_min) { dist_min=tmp; i_min=i; }
        }
        n_pts += cluster_size[k] ;
      }
    }
    for(k=k_min+1; k<kk; k++) {
      if(d_cluster2query[k]<dist_min) {
        istart = cluster_start[k];  iend = istart+cluster_size[k];
        for(i=istart; i<iend; i++) {
          for(j=0; j<dim; j++) datum[j] = data[i*dim+j] ;
          tmp = calc_dist_square(dim, query, datum) ;
          if(tmp<dist_min) { dist_min=tmp; i_min=i; }
        }
        n_pts += cluster_size[k] ;
      }
    }
    for(j=0; j<dim; j++) result_pt[j] = data[i_min*dim+j] ;
  }

  free(d_cluster2query) ;
  return n_pts ;
}
