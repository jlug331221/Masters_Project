/**
 *
 * The following program ...
 *
 */

#include "Header_Files/Headers.h"
#include "Header_Files/Defs.h"
#include "Header_Files/helping_procedures.h"

void print_execution_error_message()
{
  printf("\nExecute using the following: ./main <clustering> <dataset>\n\n");

  printf("The second argument specifies the algorithm used to cluster the data.\n\n");
  printf("\t1: LSH\n\t2: Kdtree\n\t3: BKmeans\n");

  printf("\nThe third argument specifies the dataset to cluster.\n\n");
  printf("\t1: MNIST\n\t2: BIO\n\t3: HIGGS\n\n");
}

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

/**
 * Fetch the training and testing data sets.
 */
void fetch_datasets(int data_set, double **train_feature_data, double **test_feature_data,
                    int **train_non_feature_data, int **test_non_feature_data)
{
  if(data_set == 1) { // MNIST data set
    printf("\nReading in the MNIST training and testing data sets...\n");

    *train_feature_data = malloc(MNIST_TRAIN_SIZE * MNIST_FEATURE_DIM * sizeof(double));
    *train_non_feature_data = malloc(MNIST_TRAIN_SIZE * MNIST_NON_FEATURE_DIM * sizeof(int)); // train number labels

    *test_feature_data = malloc(MNIST_TEST_SIZE * MNIST_FEATURE_DIM * sizeof(double));
    *test_non_feature_data = malloc(MNIST_TEST_SIZE * MNIST_NON_FEATURE_DIM * sizeof(int)); // test number labels

    read_MNIST_binary_dataset("MNIST_train.bin", MNIST_TRAIN_SIZE, *train_non_feature_data, *train_feature_data);
    read_MNIST_binary_dataset("MNIST_test.bin", MNIST_TEST_SIZE, *test_non_feature_data, *test_feature_data);
  }

  if(data_set == 2) { // BIO data set
    printf("\nReading in the BIO training and testing data sets...\n");

    *train_feature_data = malloc(BIO_TRAIN_SIZE * BIO_FEATURE_DIM * sizeof(double));
    *train_non_feature_data = malloc(BIO_TRAIN_SIZE * BIO_NON_FEATURE_DIM * sizeof(int));

    *test_feature_data = malloc(BIO_TEST_SIZE * BIO_FEATURE_DIM * sizeof(double));
    *test_non_feature_data = malloc(BIO_TEST_SIZE * BIO_NON_FEATURE_DIM * sizeof(int));

    read_BIO_binary_dataset("BIO_train.bin", BIO_TRAIN_SIZE, *train_non_feature_data, *train_feature_data);
    read_BIO_binary_dataset("BIO_test.bin", BIO_TEST_SIZE, *test_non_feature_data, *test_feature_data);
  }

  if(data_set == 3) { // HIGGS data set
    printf("\nReading in the HIGGS training and testing data sets...\n");

    *train_feature_data = malloc(HIGGS_TRAIN_SIZE * HIGGS_FEATURE_DIM * sizeof(double));
    *train_non_feature_data = malloc(HIGGS_TRAIN_SIZE * HIGGS_NON_FEATURE_DIM * sizeof(int));

    *test_feature_data = malloc(HIGGS_TEST_SIZE * HIGGS_FEATURE_DIM * sizeof(double));
    *test_non_feature_data = malloc(HIGGS_TEST_SIZE * HIGGS_NON_FEATURE_DIM * sizeof(int));

    read_HIGGS_binary_dataset("HIGGS_train.bin", HIGGS_TRAIN_SIZE, *train_non_feature_data, *train_feature_data);
    read_HIGGS_binary_dataset("HIGGS_test.bin", HIGGS_TEST_SIZE, *test_non_feature_data, *test_feature_data);
  }
}

/**
 * Perform clustering_algorithm on data_set and normalize the data set if normalize_data == 'y'.
 *
 * Then, perform search queries against the clustered data if desired by the user.
 */
void perform_clustering_and_search_queries(int clustering_algorithm, int data_set,
                                           char normalize_data,
                                           double *train_feature_data, double *test_feature_data,
                                           int *train_non_feature_data, int *test_non_feature_data)
{
  int k = -1;

  if(normalize_data == 'y' && data_set == 1) {
    printf("\nNormalizing the MNIST data set...\n");
    normalize_MNIST_data(train_feature_data, MNIST_FEATURE_DIM, MNIST_TRAIN_SIZE);
    normalize_MNIST_data(test_feature_data, MNIST_FEATURE_DIM, MNIST_TEST_SIZE);
  }

  if(normalize_data == 'y' && data_set == 2) {
    printf("\nNormalizing the BIO data set...\n");
    normalize_BIO_data(train_feature_data, BIO_FEATURE_DIM, BIO_TRAIN_SIZE);
    normalize_BIO_data(test_feature_data, BIO_FEATURE_DIM, BIO_TEST_SIZE);
  }

  if(normalize_data == 'y' && data_set == 3) {
    printf("\nNormalizing the HIGGS data set...\n");
    normalize_HIGGS_data(train_feature_data, HIGGS_FEATURE_DIM, HIGGS_TRAIN_SIZE);
    normalize_HIGGS_data(test_feature_data, HIGGS_FEATURE_DIM, HIGGS_TEST_SIZE);
  }

  if(clustering_algorithm == 1) { // LSH clustering
    if(data_set == 1) { // MNIST data set is normalized before LSH clustering
      execute_LSH(data_set, MNIST_FEATURE_DIM,
                  MNIST_TRAIN_SIZE, train_feature_data,
                  MNIST_TEST_SIZE, test_feature_data,
                  train_non_feature_data, test_non_feature_data);
    }

    if(data_set == 2) { // BIO data set is normalized before LSH clustering
      execute_LSH(data_set, BIO_FEATURE_DIM,
                  BIO_TRAIN_SIZE, train_feature_data,
                  BIO_TEST_SIZE, test_feature_data,
                  train_non_feature_data, test_non_feature_data);
    }

    if(data_set == 3) { // HIGGS data set is normalized before LSH clustering
      execute_LSH(data_set, HIGGS_FEATURE_DIM,
                  HIGGS_TRAIN_SIZE, train_feature_data,
                  HIGGS_TEST_SIZE, test_feature_data,
                  train_non_feature_data, test_non_feature_data);
    }
  }

  if(clustering_algorithm == 2) { // KDTree clustering
    printf("\nEnter the desired number of clusters: ");
    scanf("%d", &k);

    if(data_set == 1) { // MNIST data set
      execute_kdtree(data_set, MNIST_FEATURE_DIM, k,
                     MNIST_TRAIN_SIZE, train_feature_data,
                     MNIST_TEST_SIZE, test_feature_data,
                     train_non_feature_data, test_non_feature_data);
    }

    if(data_set == 2) { // BIO data set
      execute_kdtree(data_set, BIO_FEATURE_DIM, k,
                     BIO_TRAIN_SIZE, train_feature_data,
                     BIO_TEST_SIZE, test_feature_data,
                     train_non_feature_data, test_non_feature_data);
    }

    if(data_set == 3) { // HIGGS data set
      execute_kdtree(data_set, HIGGS_FEATURE_DIM, k,
                     HIGGS_TRAIN_SIZE, train_feature_data,
                     HIGGS_TEST_SIZE, test_feature_data,
                     train_non_feature_data, test_non_feature_data);
    }
  }

  if(clustering_algorithm == 3) { // BKmeans clustering
    printf("\nEnter the desired number of clusters: ");
    scanf("%d", &k);

    if(data_set == 1) { // MNIST data set
      execute_bkmeans_j(data_set, MNIST_FEATURE_DIM, k,
                        MNIST_TRAIN_SIZE, train_feature_data,
                        MNIST_TEST_SIZE, test_feature_data,
                        train_non_feature_data, test_non_feature_data);
    }

    if(data_set == 2) { // BIO data set
      execute_bkmeans_j(data_set, BIO_FEATURE_DIM, k,
                        BIO_TRAIN_SIZE, train_feature_data,
                        BIO_TEST_SIZE, test_feature_data,
                        train_non_feature_data, test_non_feature_data);
    }

    if(data_set == 3) { // HIGGS data set
      execute_bkmeans_j(data_set, HIGGS_FEATURE_DIM, k,
                        HIGGS_TRAIN_SIZE, train_feature_data,
                        HIGGS_TEST_SIZE, test_feature_data,
                        train_non_feature_data, test_non_feature_data);
    }
  }
}

int main(int argc, char **argv)
{
  if(argc < 3 || strtol(argv[1], NULL, 0) >= 4 || strtol(argv[1], NULL, 0) <= 0 ||
      strtol(argv[2], NULL, 0) >= 4 || strtol(argv[2], NULL, 0) <= 0) {
    print_execution_error_message();
    return 1;
  }

  clock_t begin = clock();

  srand(time(NULL));

  int data_set = (int) strtol(argv[2], NULL, 0), clustering_algorithm = (int) strtol(argv[1], NULL, 0);
  char normalize_data = '\0';

  printf("Would you like to normalize the data set?\nEnter 'y' or 'n': ");
  scanf("%s", &normalize_data);

  while(normalize_data != 'y' && normalize_data != 'n') {
    printf("\nPlease enter 'y' or 'n' if you do or do not want to normalize the dataset: ");
    scanf("%s", &normalize_data);
  }

  double *train_feature_data = NULL, *test_feature_data = NULL;
  int *train_non_feature_data = NULL, *test_non_feature_data = NULL;

  fetch_datasets(data_set, &train_feature_data, &test_feature_data,
                   &train_non_feature_data, &test_non_feature_data);

  perform_clustering_and_search_queries(clustering_algorithm, data_set,
                                        normalize_data,
                                        train_feature_data, test_feature_data,
                                        train_non_feature_data, test_non_feature_data);

  clock_t end = clock();

  print_time(begin, end);

  return 0;
}
