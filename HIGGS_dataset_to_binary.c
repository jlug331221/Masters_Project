#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define TRAIN_SIZE  10500000
#define TEST_SIZE   500000
#define FEATURE_DIM 28

/**
 *
 * Data is defined as: Class Label (1 for signal or 0 for background) | FEATURES
 *
 */

double* read_train_dataset(char *file_path, int size)
{
  int i;
  FILE *f = fopen(file_path, "r");

  if(f == NULL) {
    perror("Error");
    exit(1);
  }

  double *data = malloc(size * (FEATURE_DIM + 1) * sizeof(double));
  char *value = malloc(256 * sizeof(char));

  for(i = 0; i < size; i++) {
    // Read in data to value
    fscanf(f, "%256[^,],", value);
    data[i] = strtod(value, NULL);
  }

  fclose(f);

  return data;
}

double* read_test_dataset(char *file_path, int size)
{
  int i;
  FILE *f = fopen(file_path, "r");

  if(f == NULL) {
    perror("Error");
    exit(1);
  }

  double *data = malloc(size * (FEATURE_DIM + 1) * sizeof(double));
  char *value = malloc(256 * sizeof(char));

  for(i = size * (FEATURE_DIM + 1); i < (TRAIN_SIZE + TEST_SIZE) * (FEATURE_DIM + 1); i++) {
    // Read in data to value
    fscanf(f, "%256[^,],", value);
    data[i] = strtod(value, NULL);
  }

  fclose(f);

  return data;
}

void write_binary_dataset(char *file_path, double *data, size_t size)
{
  FILE *bf = fopen(file_path, "wb");

  if(bf == NULL) {
    perror("Error");
    exit(1);
  }

  size_t elements_written = fwrite(data, sizeof(double), size * (FEATURE_DIM + 1), bf);

  printf("\n%d elements written to %s\n", (int) elements_written, file_path);

  fclose(bf);
}

read_binary_dataset(char *file_path, int size, double *data)
{
  int i;
  FILE *bf = fopen(file_path, "rb");

  if(bf == NULL) {
    perror("Error");
    exit(1);
  }

  for(i = 0; i < size; i++) {
    fread(&data[i * (FEATURE_DIM + 1)], sizeof(double), (FEATURE_DIM + 1), bf);
  }

  fclose(bf);
}

void generate_HIGGS_binary_datasets()
{
  printf("\nReading in training and testing datasets...\n");
  double *HIGGS_train_data = read_train_dataset("HIGGS.csv", TRAIN_SIZE);
  double *HIGGS_test_data = read_test_dataset("HIGGS.csv", TRAIN_SIZE);

  printf("Converting to binary...\n");
  write_binary_dataset("HIGGS_train.bin", HIGGS_train_data, TRAIN_SIZE);
  write_binary_dataset("HIGGS_test.bin", HIGGS_test_data, TEST_SIZE);

  free(HIGGS_train_data);
}

int main()
{
  clock_t begin = clock();

  generate_HIGGS_binary_datasets();

  double *HIGGS_train_data = malloc(TRAIN_SIZE * (FEATURE_DIM + 1) * sizeof(double));
  double *HIGGS_test_data = malloc(TEST_SIZE * (FEATURE_DIM + 1) * sizeof(double));

  printf("\nReading in data from binary files...\n");

  read_binary_dataset("HIGGS_train.bin", TRAIN_SIZE, HIGGS_train_data);
  read_binary_dataset("HIGGS_test.bin", TEST_SIZE, HIGGS_test_data);

  clock_t end = clock();

  printf("\nDone\n\n");

  printf("Total time to read in data and convert to binary: %.2f secs\n",
         (double) (end - begin) / CLOCKS_PER_SEC);
}
