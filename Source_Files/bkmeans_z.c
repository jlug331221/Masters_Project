#include "../Header_Files/Headers.h"
#include "../Header_Files/bkmeans_z.h"

double calc_dist_square_z(int dim, double *datum1, double *datum2)
{
  int    i;
  double dist = 0.0;

  for(i=0; i<dim; i++) dist += (datum1[i]-datum2[i])*(datum1[i]-datum2[i]) ;
  return dist;   /* Square of Euclidean distance */
}  /****** end of calc_dist_square  ******/

int two_means(int iterat_limit, int dim, int i0, int im, double *data,    /* line of input   */
              int *cluster_assign, double *datum, double *center0,        /* line of buffers */
              double *sseDimwise, double *center, int start[2], int size[2])
{
  int    i, j, k, i_max, j_max, iterations, change, offset,
      start0, start1, end0, end1 ;
  double tmp, tmp0, tmp1, dist_max, dist, dist0, dist1, var_max ;

  /*** Compute dimesnion of largest variance ***/
  var_max = sseDimwise[0];     j_max = 0;
  for(j=1; j<dim; j++)
    if(sseDimwise[j] > var_max) {var_max=sseDimwise[j]; j_max=j;}
  tmp = center[j_max] ;
  /*** End of computing dimesnion of largest variance ***/

  /*** Choosing initial pair of centers ***/
  size[0] = 0;   size[1] = 0;
  for(j=0; j<dim; j++) { center[j]=0.0;     center[dim+j]=0.0; }
  for(i=i0; i<im; i++) {
    for(j=0; j<dim; j++) datum[j] = data[dim*i+j] ;
    k = ( (datum[j_max] < tmp) ? 0 : 1 ) ;
    for(j=0; j<dim; j++) center[k*dim+j] += datum[j] ;
    cluster_assign[i] = k ;
    size[k]++ ;
  }
  if(size[0] > 0)
    for(j=0; j<dim; j++) center[j] /= size[0] ;
  else { printf("AAA: size[0]=0.\n");  return 0; }
  if(size[1] > 0)
    for(j=0; j<dim; j++) center[dim+j] /= size[1] ;
  else { printf("AAA: size[1]=0.\n");  return 0; }
  /*** End of choosing initial pair of centers ***/


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
    }/* center0[] needed for calculating center[] at Line00000 */

    for(i=i0; i<im; i++) {
      for(j=0; j<dim; j++) datum[j] = data[i*dim+j] ;
      dist0 = 0.0 ;      dist1 = 0.0 ;
      for(j=0; j<dim; j++) {
        tmp0 = datum[j]-center0[j] ;      dist0 += tmp0*tmp0 ;
        tmp1 = datum[j]-center0[dim+j] ;  dist1 += tmp1*tmp1 ;
      }
      k = ((dist0 <= dist1) ? 0 : 1 ) ;
      if( cluster_assign[i] != k ) change++ ;
      cluster_assign[i] = k ;
      size[k]++ ;
      for(j=0; j<dim; j++) center[k*dim+j] += datum[j]; // Line00000
    } /****** Now, each datum has been assigned to a cluster ******/

    if(size[0] > 0)
      for(j=0; j<dim; j++) center[j] /= size[0] ;
    else { printf("BBB: size[0]=0. Iteration =%d\n", iterations);  return 0; }

    if(size[1] > 0)
      for(j=0; j<dim; j++) center[dim+j] /= size[1] ;
    else { printf("BBB: size[1]=0. Iteration =%d\n", iterations);  return 0; }
  }/****** End of while(iterations < iterat_limit) ******/


  /****** Data Re-ordering ******/
  start[0]=i0;       start[1] = i0+size[0];
  start0= start[0];  end0 = start0 + size[0] ;
  start1= start[1];  end1 = start1 + size[1] ;
  while(start0 < end0) {/*** end0 and end1 MUST NOT EXCEED im ***/
    while( (start0<end0) && (cluster_assign[start0]==0) ) start0++ ;
    while( (start1<end1) && (cluster_assign[start1]==1) ) start1++ ;
    if( (start0<end0) && (start1<end1) ) {
      for(j=0; j<dim; j++) datum[j]           = data[start0*dim+j] ;
      for(j=0; j<dim; j++) data[start0*dim+j] = data[start1*dim+j] ;
      for(j=0; j<dim; j++) data[start1*dim+j] = datum[j] ;
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

  /****** Compute sseDimwise[k] ******/
  for(k=0; k<2; k++) {
    for(j=0; j<dim; j++) sseDimwise[k*dim+j] = 0.0 ;
    for(i=start[k]; i<start[k]+size[k]; i++) {
      for(j=0; j<dim; j++) {
        tmp = data[i*dim+j]-center[k*dim+j] ;
        sseDimwise[k*dim+j] += tmp*tmp ;
      }
    }
  }

  return 1 ;
} /****************** End of function two_means() ******************/

int bkmeans_z(int iterat_limit, int kk, int dim, int i0_in, int im_in, double *data, // line of input
              int *cluster_assign, double *datum,                                    // line of buffers
              double *cluster_center, double *cluster_radius,
              int *cluster_start, int *cluster_size, double *cluster_sse)
{
  int    i0, im, i, i_max, j, k, k_max, nclusters, start[2], size[2] ;
  double tmp, dist_max, dist, sse_initial, *center, *sseDimwise,
      *cluster_sseDimwise, *cluster_center0, *cluster_radius_pt ;

  center    = (double *) calloc((2*dim), sizeof(double)) ;
  sseDimwise= (double *) calloc((2*dim), sizeof(double)) ;
  cluster_sseDimwise= (double *) calloc((kk*dim), sizeof(double)) ;
  cluster_center0   = (double *) calloc((kk*dim), sizeof(double)) ;
  cluster_radius_pt = (double *) calloc((kk*dim), sizeof(double)) ;


/*** Compute centroid and dimwise sse of the whole dataset ***/
  for(j=0; j<dim; j++) {cluster_center[j]=0.0 ; cluster_sseDimwise[j]=0.0;}
  for(i=i0_in; i<im_in; i++) {
    for(j=0; j<dim; j++) datum[j] = data[i*dim+j] ;
    for(j=0; j<dim; j++) cluster_center[j] += datum[j];
    for(j=0; j<dim; j++) cluster_sseDimwise[j] += datum[j]*datum[j] ;
  }
  tmp = (double)(im_in-i0_in) ;
  for(j=0; j<dim; j++) {
    cluster_center[j] /= tmp ;
    cluster_sseDimwise[j] -= tmp*cluster_center[j]*cluster_center[j] ;
  }
/*** End of computing centroid and dimwise sse of the whole dataset ***/


  for(k=0; k<kk; k++) cluster_sse[k] = 0.0 ;
  two_means(iterat_limit, dim, i0_in, im_in, data,
            cluster_assign, datum, cluster_center0,
            cluster_sseDimwise,cluster_center,cluster_start,cluster_size) ;
  cluster_sse[0] = 0.0 ;
  cluster_sse[1] = 0.0 ;
  for(j=0; j<dim; j++) cluster_sse[0] += cluster_sseDimwise[j] ;
  for(j=0; j<dim; j++) cluster_sse[1] += cluster_sseDimwise[dim+j] ;
  nclusters = 2 ;

  while( nclusters < kk ) {
    tmp = cluster_sse[0] ;    k_max = 0 ;
    for(k=1; k<nclusters; k++) {  /* Find cluster with largest sse */
      if( cluster_sse[k] > tmp ) { tmp=cluster_sse[k]; k_max=k; }
    }

    /*** Split the cluster into 2 clusters ***/
    i0 = cluster_start[k_max];
    im = i0 + cluster_size[k_max];
    for(j=0; j<dim; j++) sseDimwise[j]= cluster_sseDimwise[k_max*dim+j] ;
    for(j=0; j<dim; j++) center[j]    = cluster_center[k_max*dim+j] ;
    two_means(iterat_limit, dim, i0, im, data,
              cluster_assign, datum, cluster_center0,
              sseDimwise, center, start, size) ;
#if 1
    for(k=nclusters-1; k>k_max; k--) {
      cluster_start[k+1]= cluster_start[k];
      cluster_size[k+1] = cluster_size[k] ;
      for(j=0;j<dim;j++) cluster_sseDimwise[(k+1)*dim+j]=cluster_sseDimwise[k*dim+j] ;
      for(j=0;j<dim;j++) cluster_center[(k+1)*dim+j]    =cluster_center[k*dim+j] ;
    }
    k = k_max ;
    cluster_start[k]  = start[0] ;
    cluster_size[k]  = size[0] ;
    for(j=0;j<dim;j++) cluster_sseDimwise[k*dim+j] = sseDimwise[j] ;
    for(j=0;j<dim;j++) cluster_center[k*dim+j]     = center[j] ;
    cluster_sse[k] = 0.0 ;
    for(j=0; j<dim; j++) cluster_sse[k] += sseDimwise[j] ;

    k++ ;
    cluster_start[k]= start[1] ;
    cluster_size[k]= size[1] ;
    for(j=0;j<dim;j++) cluster_sseDimwise[k*dim+j]= sseDimwise[dim+j] ;
    for(j=0;j<dim;j++) cluster_center[k*dim+j]    = center[dim+j] ;
    cluster_sse[k]= 0.0 ;
    for(j=0; j<dim; j++) cluster_sse[k]+=sseDimwise[dim+j] ;
#else
    cluster_start[k_max]  = start[0] ;
      cluster_start[nclusters]= start[1] ;
      cluster_size[k_max]  = size[0] ;
      cluster_size[nclusters]= size[1] ;
      for(j=0;j<dim;j++) cluster_sseDimwise[k_max*dim+j]     = sseDimwise[j] ;
      for(j=0;j<dim;j++) cluster_sseDimwise[nclusters*dim+j] = sseDimwise[dim+j] ;
      for(j=0;j<dim;j++) cluster_center[k_max*dim+j]     = center[j] ;
      for(j=0;j<dim;j++) cluster_center[nclusters*dim+j] = center[dim+j] ;
      cluster_sse[k_max]  = 0.0 ;     cluster_sse[nclusters]= 0.0 ;
      for(j=0; j<dim; j++) {
         cluster_sse[k_max] += sseDimwise[j] ;
         cluster_sse[nclusters] += sseDimwise[dim+j] ;
      }
#endif
    nclusters++;
  } /************ End of while() loop ************/

  for(k=0; k<nclusters; k++) {
    if(cluster_size[k] == 0) { printf("CCC: in bkmeans: cluster %d empty!\n", k); }
  }

  if(kk > 2) {
    iterat_limit = 3;
    nclusters = kmeans_z(iterat_limit, kk, dim, i0_in, im_in, data,
                       cluster_assign, datum, cluster_center0, cluster_radius_pt,
                       cluster_center,cluster_radius,cluster_start,cluster_size,cluster_sse);
  }

  free(center); free(sseDimwise);  free(cluster_sseDimwise);
  free(cluster_center0); free(cluster_radius_pt);
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
           int *cluster_assign, double *datum, double *cluster_center0, double *radius_pt,
           double *cluster_center, double *cluster_radius,
           int *cluster_start, int *cluster_size, double *cluster_sse)
{
  int i, j, k, k_max, k_chosen, iterations, membership, start0, end0, start1, end1,
      position, change, *cluster_size0, *radius_index, k_max_used, nclusters ;
  double  tmp, dist, dist_min ;

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
      for(j=0;j<dim;j++) cluster_center0[k*dim+j]=cluster_center[k*dim+j];

      cluster_size[k]=0 ;
      cluster_radius[k]=0.0;
      for(j=0;j<dim;j++) cluster_center[k*dim+j] = 0.0 ;
    }/* cluster_center0 needed for calculating cluster_center at DDD */

    for(i=i0; i<im; i++) {
      for(j=0; j<dim; j++) datum[j] = data[i*dim+j] ;
      dist_min = 9876543210123.0 ; /* Find closest center to datum */
      for(k=0; k<kk; k++) {
        if(cluster_size0[k]>0) {
          dist = 0.0 ;
          for(j=0; j<dim; j++) {
            tmp = datum[j]-cluster_center0[k*dim+j] ;   dist += tmp*tmp;
          }
          if(dist < dist_min) { membership=k ;    dist_min=dist ; }
        }
      }
      if(cluster_assign[i] != membership) change++ ;
      cluster_assign[i] = membership ;
      cluster_size[membership]++ ;
      for(j=0;j<dim;j++) cluster_center[membership*dim+j] += datum[j];//DDD
      if(dist_min > cluster_radius[membership]) {
        cluster_radius[membership] = dist_min ;
        radius_index[membership] = i ;
        for(j=0; j<dim; j++) radius_pt[dim*membership+j] = datum[j] ;
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
            for(j=0; j<dim; j++) cluster_center[k*dim+j] = radius_pt[k_max*dim+j] ;
         }
#else
      else {
        k_chosen = (k_max_used + k_max) % kk;
        if(cluster_size[k_chosen] > 2) {
          i = radius_index[k_chosen];
          cluster_size[k]++;
          cluster_size[k_chosen]--;

          cluster_assign[i] = k;
          for(j=0; j<dim; j++)
            cluster_center[k*dim+j] = radius_pt[k_chosen*dim+j];
        }
        k_max_used ++ ;
      }
#endif
    }
  }/*** End of while( iterations < iterat_limit ) ***/


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
        for(j=0; j<dim; j++) datum[j]           = data[start0*dim+j] ;
        for(j=0; j<dim; j++) data[start0*dim+j] = data[start1*dim+j] ;
        for(j=0; j<dim; j++) data[start1*dim+j] = datum[j] ;
        cluster_assign[start1] = cluster_assign[start0] ;
        cluster_assign[start0] = k;
        start0++;   start1++;
      }
    }
  }/****** End of loop for(k=0; k<kk-1; k++). End of Data Re-ordering ******/

  int num_empty = 0 ;
  for(k = 0; k < kk; k++) {/*** Calculate sse, radius, delete empty clusters ***/
    if(cluster_size[k]==0) { num_empty++ ;    nclusters-- ; }
    else {
      cluster_sse[k] = 0.0 ;        cluster_radius[k]=0.0;
      start0 = cluster_start[k] ;   end0 = start0 + cluster_size[k] ;
      for(i=start0; i<end0; i++) {
        dist = 0.0 ;
        for(j=0; j<dim; j++) {
          tmp = data[i*dim+j]-cluster_center[k*dim+j] ;   dist += tmp*tmp;
        }
        if(dist > cluster_radius[k]) cluster_radius[k] = dist ;
        cluster_sse[k] += dist ;
      }
      if(num_empty > 0) {
        cluster_start[k-num_empty] = cluster_start[k] ;
        cluster_size[k-num_empty]  = cluster_size[k] ;
        for(j=0; j<dim; j++)
          cluster_center[(k-num_empty)*dim+j] = cluster_center[k*dim+j] ;
        cluster_radius[k-num_empty]= cluster_radius[k] ;
        cluster_sse[k-num_empty]   = cluster_sse[k] ;
      }
    }
  }
  for(k=0; k<nclusters; k++) cluster_radius[k] = sqrt(cluster_radius[k]) ;


  free(cluster_size0) ;   free(radius_index) ;
  return nclusters ;
} /****************** End of function kmeans() ******************/

/****************** End of File bkmeans.c ******************/