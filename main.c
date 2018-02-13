/**
 *
 * The following program ...
 *
 */

#include "Header_Files/Headers.h"
#include "Header_Files/Defs.h"
#include "Header_Files/helping_procedures.h"

void print_time(clock_t begin, clock_t end)
{
  double total_sec = (double) (end-begin)/CLOCKS_PER_SEC, remaining_secs = 0.0, remaining_mins = 0.0;

  if(total_sec >= 60) {
    double mins = floor(total_sec / 60);

    if(mins >= 60) {
      double hours = floor(mins / 60);
      remaining_mins = mins - (60 * hours);
      printf("\nTotal execution time: %.0f hour", hours);
      if(hours > 1) { printf("s %.1f min(s)\n\n", remaining_mins); }
      else { printf(" %.1f min(s)\n\n", remaining_mins); }
    }
    else {
      remaining_secs = total_sec - (mins * 60);
      printf("\nTotal execution time: %d min", (int) mins);
      if(mins > 1) { printf("s %.1f sec(s)\n\n", remaining_secs); }
      else { printf(" %.1f sec(s)\n\n", remaining_secs); }
    }
  }
  else {
    printf("\nTotal execution time: %.1f sec(s)\n\n", (double) (end-begin)/CLOCKS_PER_SEC);
  }
}

int main(int argc, char **argv)
{
  if(argc < 3 || strtol(argv[1], NULL, 0) > 6 || strtol(argv[1], NULL, 0) <= 0) {
    printf("\nExecute using the following: ./main <clustering> <dataset>\n\n");

    printf("The second argument specifies the algorithm used to cluster the data.\n\n");
    printf("\t1: LSH\n\t2: kdtree\n\t3: BKmeans_j\n\t4: kdtree_median\n\t5: BKmeans_z\n");

    printf("\nThe the third argument specifies the dataset to cluster.\n\n");
    printf("\t1: MNIST\n\t2: BIO\n\t3: HIGGS\n\n");

    return 1;
  }

  clock_t begin = clock();

  srand(time(NULL));

  double *train_feature_data = NULL, *test_feature_data = NULL;
  int *train_non_feature_data = NULL, *test_non_feature_data = NULL;

  if(strtol(argv[2], NULL, 0) == 1) { // MNIST data set
    printf("\nReading in the MNIST data set...\n");

    train_feature_data = malloc(MNIST_TRAIN_SIZE * MNIST_FEATURE_DIM * sizeof(double));
    train_non_feature_data = malloc(MNIST_TRAIN_SIZE * MNIST_NON_FEATURE_DIM * sizeof(int)); // train number labels

    //test_feature_data = malloc(MNIST_TEST_SIZE * MNIST_FEATURE_DIM * sizeof(double));
    //test_non_feature_data = malloc(MNIST_TEST_SIZE * MNIST_NON_FEATURE_DIM * sizeof(int)); // test number labels

    read_MNIST_binary_dataset("MNIST_train.bin", MNIST_TRAIN_SIZE, train_non_feature_data, train_feature_data);
    //read_MNIST_binary_dataset("MNIST_test.bin", MNIST_TEST_SIZE, test_non_feature_data, test_feature_data);
  }

  if(strtol(argv[2], NULL, 0) == 2) { // BIO data set
    printf("\nReading in the BIO data set...\n");

    train_feature_data = malloc(BIO_TRAIN_SIZE * BIO_FEATURE_DIM * sizeof(double));
    train_non_feature_data = malloc(BIO_TRAIN_SIZE * BIO_NON_FEATURE_DIM * sizeof(int));

    //test_feature_data = malloc(BIO_TEST_SIZE * BIO_FEATURE_DIM * sizeof(double));
    //test_non_feature_data = malloc(BIO_TEST_SIZE * BIO_NON_FEATURE_DIM * sizeof(int));

    read_BIO_binary_dataset("BIO_train.bin", BIO_TRAIN_SIZE, train_non_feature_data, train_feature_data);
    //read_BIO_binary_dataset("BIO_test.bin", BIO_TEST_SIZE, test_non_feature_data, test_feature_data);

    //double max = (double) INT_MIN, min = (double) INT_MAX;
    //for(int i = 0; i < BIO_TRAIN_SIZE * BIO_FEATURE_DIM; i++) {
    //  if(train_feature_data[i] > max) { max = train_feature_data[i]; }
    //  if(train_feature_data[i] < min) { min = train_feature_data[i]; }
    //}
    //printf("** min value of BIO data = %lf || max value of BIO data = %lf ** \n", min, max);
  }

  if(strtol(argv[2], NULL, 0) == 3) { // HIGGS data set
    printf("\nReading in the HIGGS data set...\n");

    train_feature_data = malloc(HIGGS_TRAIN_SIZE * HIGGS_FEATURE_DIM * sizeof(double));
    train_non_feature_data = malloc(HIGGS_TRAIN_SIZE * HIGGS_NON_FEATURE_DIM * sizeof(int));

    //test_feature_data = malloc(HIGGS_TEST_SIZE * HIGGS_FEATURE_DIM * sizeof(double));
    //test_non_feature_data = malloc(HIGGS_TEST_SIZE * HIGGS_NON_FEATURE_DIM * sizeof(int));

    read_HIGGS_binary_dataset("HIGGS_train.bin", HIGGS_TRAIN_SIZE, train_non_feature_data, train_feature_data);
    //read_HIGGS_binary_dataset("HIGGS_test.bin", HIGGS_TEST_SIZE, test_non_feature_data, test_feature_data);

    //double max = (double) INT_MIN, min = (double) INT_MAX;
    //for(int i = 0; i < BIO_TRAIN_SIZE * BIO_FEATURE_DIM; i++) {
    //  if(train_feature_data[i] > max) { max = train_feature_data[i]; }
    //  if(train_feature_data[i] < min) { min = train_feature_data[i]; }
    //}
    //printf("** min value of HIGGS data = %lf || max value of HIGGS data = %lf ** \n", min, max);
  }

  if(strtol(argv[1], NULL, 0) == 1) { // LSH clustering
    if(strtol(argv[2], NULL, 0) == 1) { // MNIST data set is normalized before LSH clustering
      normalize_MNIST_data(train_feature_data, MNIST_FEATURE_DIM, MNIST_TRAIN_SIZE);
      //normalize_MNIST_data(test_feature_data, MNIST_FEATURE_DIM, MNIST_TEST_SIZE);

      execute_LSH(MNIST_FEATURE_DIM, MNIST_TRAIN_SIZE, train_feature_data, (int) strtol(argv[2], NULL, 0));
    }

    if(strtol(argv[2], NULL, 0) == 2) { // BIO data set is normalized before LSH clustering
      normalize_BIO_data(train_feature_data, BIO_FEATURE_DIM, BIO_TRAIN_SIZE);
      //normalize_BIO_data(test_feature_data, BIO_FEATURE_DIM, BIO_TEST_SIZE);

      execute_LSH(BIO_FEATURE_DIM, BIO_TRAIN_SIZE, train_feature_data, (int) strtol(argv[2], NULL, 0));
    }

    if(strtol(argv[2], NULL, 0) == 3) { // HIGGS data set is normalized before LSH clustering
      normalize_HIGGS_data(train_feature_data, HIGGS_FEATURE_DIM, HIGGS_TRAIN_SIZE);
      normalize_HIGGS_data(train_feature_data, HIGGS_FEATURE_DIM, HIGGS_TEST_SIZE);

      execute_LSH(HIGGS_FEATURE_DIM, HIGGS_TRAIN_SIZE, train_feature_data, (int) strtol(argv[2], NULL, 0));
    }
  }

//  if(strtol(argv[1], NULL, 0) == 2) { execute_kdtree(train_labels, train_data, test_labels, test_data); }
//
//  if(strtol(argv[1], NULL, 0) == 3) { execute_bkmeans_j(train_labels, train_data, test_labels, test_data); }
//
//  if(strtol(argv[1], NULL, 0) == 4) { execute_kdtree_median(test_labels, train_data, test_labels, test_data); }
//
//  if(strtol(argv[1], NULL, 0) == 5) { execute_bkmeans_z(train_labels, train_data, test_labels, test_data); }

  clock_t end = clock();

  print_time(begin, end);

  return 0;
}
