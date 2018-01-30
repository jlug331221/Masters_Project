#include "../Header_Files/Headers.h"
#include "../Header_Files/bkmeans_z.h"

double calc_dist_square(int dim, double *datum1, double *datum2)
{
  int    i;
  double tmp, dist = 0.0;

  for(i=0; i<dim; i++) { tmp = datum1[i]-datum2[i] ;   dist += tmp*tmp ; }
  return dist;   /* Square of Euclidean distance */
}  /****** end of calc_dist_square  ******/

int two_means(int iterat_limit, int dim, int i0, int im, double *data,    /* line of input   */
              int *cluster_assign, double *datum, double *center0,        /* line of buffers */
              double *radius_pt, double *center, int start[2],
              int size[2], double ssd[2])
{
  int    i, j, k, i_max, iterations, change, offset,
      start0, start1, end0, end1 ;
  double tmp, dist_max, dist, dist0, dist1 ;

  for(j=0;j<dim;j++) {/*** Choosing initial pair of centers ***/
    center[dim+j] = radius_pt[j] ;
    center[j] = 2.0*center[j] - radius_pt[j] ;
  } /*** End of choosing initial pair of centers ***/

  for(i=i0; i<im; i++) cluster_assign[i] = 0 ;
  change = 1 ;
  iterations = 0 ;
  while( (iterations<iterat_limit) && (change!=0) ) {
    iterations++ ;
    change = 0 ;
    size[0]=0 ;   size[1]=0 ;
    for(j=0;j<dim;j++) {
      center0[j]    = center[j] ;
      center0[dim+j]= center[dim+j];
      center[j]=0.0;  center[dim+j]=0.0;
    }/* center0 needed for calculating center at AAA */

    for(i=i0; i<im; i++) {
      memcpy(datum, data+i*dim, dim*sizeof(double)) ;
      dist0 = calc_dist_square(dim, datum, center0);
      dist1 = calc_dist_square(dim, datum, center0+dim);
      k = ((dist0 < dist1) ? 0 : 1 ) ;
      if( cluster_assign[i] != k ) change++ ;
      cluster_assign[i] = k ;
      size[k]++ ;
      for(j=0; j<dim; j++) center[k*dim+j] += datum[j]; // AAA
    } /****** Now, each datum has been assigned to a cluster ******/

    if(size[0] > 0)
      for(j=0; j<dim; j++) center[j] /= size[0] ;
    if(size[1] > 0)
      for(j=0; j<dim; j++) center[dim+j] /= size[1] ;

  }/****** End of while(iterations < iterat_limit) ******/


  /****** Data Re-ordering ******/
  start[0]=i0;       start[1] = i0+size[0];
  start0= start[0];  end0 = start0 + size[0] ;
  start1= start[1];  end1 = start1 + size[1] ;
  while(start0 < end0) {/*** end0 and end1 MUST NOT EXCEED im ***/
    while( (start0<end0) && (cluster_assign[start0]==0) ) start0++ ;
    while( (start1<end1) && (cluster_assign[start1]==1) ) start1++ ;
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
  }/*** End of Data Re-ordering ***/

  /****** Compute sse[k] and radius_pt[k] ******/
  ssd[0]= 0.0 ;    ssd[1] = 0.0 ;
  for(k=0; k<2; k++) {
    dist = 0.0 ; /*** dist holds radius of the cluster ***/
    for(i=start[k]; i<start[k]+size[k]; i++) {
      tmp = calc_dist_square(dim, data+i*dim, center+k*dim) ;
      ssd[k] += tmp;
      if(dist < tmp) {
        dist=tmp;    memcpy(radius_pt+k*dim, data+i*dim, dim*sizeof(double));
      }
    }
  }

  return 1 ;
} /****************** End of function two_means() ******************/

int bkmeans(int iterat_limit, int kk, int dim, int i0_in, int im_in, double *data, // line of input
            int *cluster_assign, double *datum,                                    // line of buffers
            double *cluster_center, double *cluster_radius,
            int *cluster_start, int *cluster_size, double *cluster_ssd)
{
  int    i0, im, i, j, k, k_max, nclusters, start[2], size[2] ;
  double tmp, dist_max, dist, ssd_initial, ssd[2],
      *center, *radius_pt, *cluster_center0, *cluster_radius_pt ;

  center  = (double *) calloc((2*dim), sizeof(double)) ;
  radius_pt=(double *) calloc((2*dim), sizeof(double)) ;
  cluster_center0  = (double *) calloc((kk*dim), sizeof(double)) ;
  cluster_radius_pt =(double *) calloc((kk*dim), sizeof(double)) ;

  /*** Compute centroid of the whole dataset ***/
  for(j=0; j<dim; j++) cluster_center[j] = 0.0 ;
  for(i=i0_in; i<im_in; i++)
    for(j=0; j<dim; j++) cluster_center[j] += data[i*dim+j];
  for(j=0; j<dim; j++) cluster_center[j] /= (double)(im_in-i0_in) ;
  /*** End of computing centroid of the whole dataset ***/

  /*** Compute farthest pt to centroid. Store it in radius_pt[0] ***/
  for(k=0; k<kk; k++) cluster_radius[k] = 0.0 ;
  ssd_initial = 0.0 ;
  for(i=i0_in; i<im_in; i++) {
    tmp = calc_dist_square(dim, data+i*dim, cluster_center) ;
    ssd_initial += tmp ;
    if(cluster_radius[0] < tmp) {
      cluster_radius[0] = tmp ;
      memcpy(cluster_radius_pt, data+i*dim, dim*sizeof(double)) ;
    }
  } /*** End of computing farthest pt to centroid ***/

  for(k=0; k<kk; k++) cluster_ssd[k] = 0.0 ;
  two_means(iterat_limit, dim, i0_in, im_in, data,
            cluster_assign, datum, cluster_center0,
            cluster_radius_pt,cluster_center,cluster_start,cluster_size,cluster_ssd) ;
  nclusters = 2 ;

  while( nclusters < kk ) {
    tmp = 0.0 ;
    for(k=0; k<nclusters; k++) {  /* Find cluster with largest ssd */
      if( cluster_ssd[k] > tmp ) { tmp=cluster_ssd[k]; k_max=k; }
    }

    /*** Split the cluster into 2 clusters ***/
    i0 = cluster_start[k_max];
    im = i0 + cluster_size[k_max];
    for(j=0;j<dim;j++) radius_pt[j]= cluster_radius_pt[k_max*dim+j] ;
    for(j=0;j<dim;j++) center[j]   = cluster_center[k_max*dim+j] ;
    two_means(iterat_limit, dim, i0, im, data,
              cluster_assign, datum, cluster_center0,
              radius_pt, center, start, size, ssd) ;
  #if 1
    for(k=nclusters-1; k>k_max; k--) {
      cluster_start[k+1]= cluster_start[k];
      cluster_size[k+1] = cluster_size[k] ;
      cluster_ssd[k+1]  = cluster_ssd[k] ;
      memcpy(cluster_radius_pt+(k+1)*dim,cluster_radius_pt+k*dim,dim*sizeof(double)) ;
      memcpy(cluster_center + (k+1)*dim, cluster_center + k*dim, dim*sizeof(double)) ;
    }
    k = k_max ;
    cluster_start[k]  = start[0] ;
    cluster_start[k+1]= start[1] ;
    cluster_size[k]  = size[0] ;
    cluster_size[k+1]= size[1] ;
    cluster_ssd[k]  = ssd[0] ;
    cluster_ssd[k+1]= ssd[1] ;
    memcpy(cluster_radius_pt+k*dim,    radius_pt,    dim*sizeof(double));
    memcpy(cluster_radius_pt+(k+1)*dim,radius_pt+dim,dim*sizeof(double));
    memcpy(cluster_center + k*dim,     center,     dim*sizeof(double)) ;
    memcpy(cluster_center + (k+1)*dim, center+dim, dim*sizeof(double)) ;
  #else
    cluster_start[k_max]  = start[0] ;
      cluster_start[nclusters]= start[1] ;
      cluster_size[k_max]  = size[0] ;
      cluster_size[nclusters]= size[1] ;
      cluster_ssd[k_max]  = ssd[0] ;
      cluster_ssd[nclusters]= ssd[1] ;
      memcpy(cluster_radius_pt+k_max*dim,    radius_pt,    dim*sizeof(double));
      memcpy(cluster_radius_pt+nclusters*dim,radius_pt+dim,dim*sizeof(double));
      memcpy(cluster_center + k_max*dim,     center,     dim*sizeof(double)) ;
      memcpy(cluster_center + nclusters*dim, center+dim, dim*sizeof(double)) ;
  #endif
    nclusters++;
  } /************ End of while() loop ************/

  for(k=0; k<nclusters; k++) {
    if(cluster_size[k] == 0) { printf("AAA: in bkmeans: cluster %d empty!\n"); }
  }

  if(kk > 2) {
    iterat_limit = 3;
    nclusters = kmeans_z(iterat_limit, kk, dim, i0_in, im_in, data,
                         cluster_assign, datum, cluster_center0, cluster_radius_pt,
                         cluster_center,cluster_radius,cluster_start,cluster_size,cluster_ssd);
  }

  free(center); free(radius_pt);  free(cluster_center0); free(cluster_radius_pt);
  return nclusters ;
} /****************** End of bkmeans() ******************/

/******************************************************************************
 kk : number of clusters, i.e. the K in K-mean.
 cluster_center[kk*dim]: input  -- stores initial kk centers
                         output -- stores kk centers
 cluster_radius[kk]:output -- the radius of each output cluster
 cluster_start[kk]: output -- the index of the 1st in each cluster
 cluster_size[kk]:  output -- the num of datapoints in each cluster
 cluster_ssd[k]:    output -- the sum of square distances to cluster centers
*******************************************************************************/
int kmeans_z(int iterat_limit, int kk, int dim, int i0, int im, double *data,
             int *cluster_assign, double *datum, double *cluster_center0,
             double *radius_pt, double *cluster_center, double *cluster_radius,
             int *cluster_start, int *cluster_size, double *cluster_ssd)
{
  int i, j, k, k_max, k_chosen, iterations, membership, start0, end0, start1, end1,
      position, change, *cluster_size0, *radius_index, k_max_used, nclusters ;
  double  tmp, dist_min ;

  nclusters = kk ;
  cluster_size0 = (int *) calloc(kk, sizeof(int)) ;
  radius_index  = (int *) calloc(kk, sizeof(int)) ;

  /****** Start of k-means iterations ******/
  change = 1 ;
  iterations = 0 ;
  while( (iterations < iterat_limit) && (change!=0) ) {
    iterations++ ;
    change = 0 ;
    for(k=0; k<kk; k++) {
      cluster_size0[k] = cluster_size[k] ;
      cluster_size[k]=0 ;
      cluster_radius[k]=0.0;
      // memcpy(cluster_center0+k*dim, cluster_center+k*dim, dim*sizeof(double)) ;
      for(j=0;j<dim;j++) cluster_center0[k*dim+j]=cluster_center[k*dim+j];
      for(j=0;j<dim;j++) cluster_center[k*dim+j] = 0.0 ;
    }/* cluster_center0 needed for calculating cluster_center at BBB */

    for(i=i0; i<im; i++) {
      memcpy(datum, data+i*dim, dim*sizeof(double));
      dist_min = 9876543210123.0 ; /* Find closest center to datum */
      for(k=0; k<kk; k++) {
        if(cluster_size0[k]>0) {
          tmp = calc_dist_square(dim, datum, cluster_center0+k*dim);
          if(tmp < dist_min) { membership=k ;    dist_min=tmp ; }
        }
      }
      if(cluster_assign[i] != membership) change++ ;
      cluster_assign[i] = membership ;
      cluster_size[membership]++ ;
      for(j=0;j<dim;j++) cluster_center[membership*dim+j] += datum[j];//BBB
      if(dist_min > cluster_radius[membership]) {
        cluster_radius[membership] = dist_min ;
        radius_index[membership] = i ;
        memcpy(radius_pt+dim*membership, datum, dim*sizeof(double));
      }
    }/*** Each data item has been assigned to a cluster ***/

    k_max=0;   tmp=cluster_radius[0];
    for(k=1; k<kk; k++) { /*** Find the cluster with the largest radius ***/
      if(tmp<cluster_radius[k]){ k_max=k; tmp=cluster_radius[k];}
    }
    k_max_used = 0 ;
    for(k=0; k<kk; k++) {
      if(cluster_size[k] > 0)
        for(j=0;j<dim;j++) cluster_center[k*dim+j] /= (double)cluster_size[k];
      #if 0
        else if(!k_max_used) {
            k_max_used = 1 ;
            i = radius_index[k_max] ;

            cluster_size[k_max]--;
            cluster_size[k]++ ;
            cluster_assign[i] = k ;
            memcpy(cluster_center+k*dim, radius_pt+k_max*dim, dim*sizeof(double));
         }
      #else
      else {
        k_chosen = (k_max_used + k_max) % kk;
        if(cluster_size[k_chosen] > 2) {
          i = radius_index[k_chosen];
          cluster_size[k]++;
          cluster_size[k_chosen]--;

          cluster_assign[i] = k;
          for(j = 0; j < dim; j++)
            cluster_center[k*dim+j] = radius_pt[k_chosen*dim+j];
        }
        k_max_used ++ ;
      }
      #endif
    }
  } /*** End of while( iterations < iterat_limit ) ***/

  /****** Data Re-ordering ******/
  position = 0 ;
  for(k=0;k<kk;k++) { cluster_start[k]=position;   position += cluster_size[k]; }
  for(k=0; k<kk-1; k++) {
    start0= cluster_start[k] ;
    end0  = start0 + cluster_size[k];
    start1= cluster_start[k+1];
    while(start0 < end0) {
      while((start0<end0)&&(cluster_assign[start0]==k)) start0++ ;
      if(start0 < end0) {
        while((start1<im)&&(cluster_assign[start1]!=k)) start1++ ;
        if(start1==im) { printf("\nError: start1==im.\n");  return 0; }
        memcpy(datum, data+start0*dim, dim*sizeof(double)) ;
        memcpy(data+start0*dim, data+start1*dim, dim*sizeof(double)) ;
        memcpy(data+start1*dim, datum, dim*sizeof(double)) ;
        cluster_assign[start1] = cluster_assign[start0] ;
        cluster_assign[start0] = k;
        start0++;   start1++;
      }
    }
  } /****** End of loop for(k=0; k<kk-1; k++). End of Data Re-ordering ******/

  int num_empty = 0 ;
  for(k = 0; k < kk; k++) {/*** Calculate sse, radius, delete empty clusters ***/
    if(cluster_size[k]==0) { num_empty++ ;    nclusters-- ; }
    else {
      cluster_ssd[k] = 0.0 ;        cluster_radius[k]=0.0;
      start0 = cluster_start[k] ;   end0 = start0 + cluster_size[k] ;
      for(i=start0; i<end0; i++) {
        tmp = calc_dist_square(dim, data+i*dim, cluster_center+k*dim) ;
        if(tmp > cluster_radius[k]) cluster_radius[k] = tmp ;
        cluster_ssd[k] += tmp ;
      }
      if(num_empty > 0) {
        cluster_start[k-num_empty] = cluster_start[k] ;
        cluster_size[k-num_empty]  = cluster_size[k] ;
        memcpy(cluster_center+(k-num_empty)*dim, cluster_center+k*dim, dim*sizeof(double)) ;
        cluster_radius[k-num_empty]= cluster_radius[k] ;
        cluster_ssd[k-num_empty]   = cluster_ssd[k] ;
      }
    }
  }
  for(k=0; k<nclusters; k++) cluster_radius[k] = sqrt(cluster_radius[k]) ;

  free(cluster_size0) ;   free(radius_index) ;
  return nclusters ;
} /****************** End of function kmeans() ******************/

/****************** End of File bkmeans.c ******************/