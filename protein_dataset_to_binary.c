#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define TRAIN_SIZE  145751
#define TEST_SIZE   139658
#define FEATURE_DIM 74

/**
 *
 * Data is defined as: BLOCK_ID | EXAMPLE_ID | class(homologous vs non-homologous) | FEATURES
 *
 */

double* read_dataset(char *file_path, int size)
{
  int i, j, k = 0;
  FILE *f = fopen(file_path, "r");

  if(f == NULL) {
    perror("Error");
    exit(1);
  }


  double *data = malloc(size * (FEATURE_DIM + 3) * sizeof(double));
  char *value = malloc(256 * sizeof(char));

  for(i = 0; i < size; i++) {
    for(j = 0; j < FEATURE_DIM + 3; j++) {
      // Read in data to value
      fscanf(f, "%s", value);
      data[k] = strtod(value, NULL);
      k++;
    }
  }

  free(value);
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

  size_t elements_written = fwrite(data, sizeof(double), size * (FEATURE_DIM + 3), bf);

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
    fread(&data[i * (FEATURE_DIM + 3)], sizeof(double), (FEATURE_DIM + 3), bf);
  }

  fclose(bf);
}

void generate_protein_homology_binary_datasets()
{
  printf("\nReading in training and testing datasets...\n");
  double *bio_train_data = read_dataset("BIO_train.dat", TRAIN_SIZE);
  double *bio_test_data = read_dataset("BIO_test.dat", TEST_SIZE);

  printf("Converting to binary...\n");
  write_binary_dataset("BIO_train.bin", bio_train_data, TRAIN_SIZE);
  write_binary_dataset("BIO_test.bin", bio_test_data, TEST_SIZE);

  free(bio_train_data);
  free(bio_test_data);
}

int main()
{
  clock_t begin = clock();

  generate_protein_homology_binary_datasets();

  double *bio_train_data = malloc(TRAIN_SIZE * (FEATURE_DIM + 3) * sizeof(double));
  double *bio_test_data = malloc(TEST_SIZE * (FEATURE_DIM + 3) * sizeof(double));

  printf("\nTest reading in data from binary files...\n");

  read_binary_dataset("BIO_train.bin", TRAIN_SIZE, bio_train_data);
  read_binary_dataset("BIO_test.bin", TEST_SIZE, bio_test_data);

  clock_t end = clock();

  free(bio_train_data);
  free(bio_test_data);

  printf("\nDone\n\n");

  printf("Total time to read in data and convert to binary: %.2f secs\n",
         (double) (end - begin) / CLOCKS_PER_SEC);
}
