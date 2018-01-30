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
  if(argc < 2 || strtol(argv[1], NULL, 0) > 5 || strtol(argv[1], NULL, 0) <= 0) {
    printf("\nExecute using the following: ./main <clustering>\n");
    printf("The second argument specifies the algorithm used to cluster the data.\n\n");
    printf("\t1: kdtree\n\t2: kmeans\n\t3: LSH\n\t4: kdtree-median\n\t5: bisecting kmeans\n");

    return 1;
  }

  clock_t begin = clock();
  srand(time(NULL));

  int *train_labels = malloc(sizeof(int) * TRAIN_SIZE);
  double *train_features = malloc(sizeof(double) * TRAIN_SIZE * FEATURE_DIM);

  int *test_labels = malloc(sizeof(int) * TEST_SIZE);
  double *test_features = malloc(sizeof(double) * TEST_SIZE * FEATURE_DIM);

  read_binary_dataset("train.bin", TRAIN_SIZE, train_labels, train_features);
  read_binary_dataset("test.bin", TEST_SIZE, test_labels, test_features);

  normalize_data(train_features, FEATURE_DIM, TRAIN_SIZE);
  normalize_data(test_features, FEATURE_DIM, TEST_SIZE);

  if(strtol(argv[1], NULL, 0) == 1) { execute_kdtree(train_labels, train_features, test_labels, test_features); }

  if(strtol(argv[1], NULL, 0) == 2) { execute_kmeans(train_labels, train_features, test_labels, test_features); }

  if(strtol(argv[1], NULL, 0) == 3) { execute_LSH(train_labels, train_features, test_labels, test_features); }

  if(strtol(argv[1], NULL, 0) == 4) { execute_kdtree_median(test_labels, train_features, test_labels, test_features); }

  if(strtol(argv[1], NULL, 0) == 5) {
    // My implementation of bisecting-Kmeans
    execute_bisecting_kmeans(train_labels, train_features, test_labels, test_features);
  }
//  if(strtol(argv[1], NULL, 0) == 5) { execute_bkmeans(train_labels, train_features, test_labels, test_features); }

  clock_t end = clock();

  print_time(begin, end);

  return 0;
}
