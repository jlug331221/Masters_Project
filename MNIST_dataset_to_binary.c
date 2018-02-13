#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define TRAIN_SIZE  60000
#define TEST_SIZE   10000
#define FEATURE_DIM 784
#define LABEL_DIM   1

double* read_dataset(char *file_path, int size)
{
  FILE *file = fopen(file_path, "r");

  if(file == NULL) {
    perror("Error");
    exit(1);
  }

  int i, j, k;
  double *data = malloc(size * (FEATURE_DIM + LABEL_DIM) * sizeof(double));
  char *current_pixel = malloc(256 * sizeof(char));
  char *current_label = malloc(256 * sizeof(char));

  k = 0;
  for (i = 0; i < size; i++) {
    // Read the label
    fscanf(file, "%s", current_label);
    data[k] = strtod(current_label, NULL);
    k++;

    // Read the features
    for (j = 0; j < FEATURE_DIM; j++) {
      fscanf(file, "%s", current_pixel);
      data[k] = strtod(current_pixel, NULL);
      k++;
    }
  }

  free(current_pixel);
  free(current_label);
  fclose(file);

  return data;
}

void write_binary_dataset(char *file_path, double *data, size_t size)
{
  FILE *file = fopen(file_path, "wb");

  if(file == NULL) {
    perror("Error");
    exit(1);
  }

  size_t elements_written = fwrite(data, sizeof(double), size * (FEATURE_DIM + LABEL_DIM), file);

  printf("\n%d elements written to %s\n", (int) elements_written, file_path);

  fclose(file);
}

void read_binary_dataset(char *path, int size, int *labels, double *features)
{
  int i;
  FILE *file = fopen(path, "rb");

  if(file == NULL) {
    perror("Error");
    exit(1);
  }

  double label = -1.0;
  for(i = 0; i < size; i++) {
    // Read the label
    fread(&label, sizeof(double), 1, file);
    labels[i] = (int) label;

    // Read the features
    fread(&features[i * FEATURE_DIM], sizeof(double), FEATURE_DIM, file);
  }

  fclose(file);
}

void generate_binary_datasets()
{
  printf("\nReading in training and testing datasets...\n");

  double *MNIST_train_data = read_dataset("MNIST_train.txt", TRAIN_SIZE);
  double *MNIST_test_data = read_dataset("MNIST_test.txt", TEST_SIZE);

  printf("Converting to binary...\n");
  write_binary_dataset("MNIST_train.bin", MNIST_train_data, TRAIN_SIZE);
  write_binary_dataset("MNIST_test.bin", MNIST_test_data, TEST_SIZE);

  free(MNIST_train_data);
  free(MNIST_test_data);
}

int main(int argc, char* argv[])
{
  clock_t begin = clock();

  generate_binary_datasets();

  int *train_labels = malloc(sizeof(int) * TRAIN_SIZE);
  double *train_features = malloc(sizeof(double) * TRAIN_SIZE * FEATURE_DIM);

  int *test_labels = malloc(sizeof(int) * TRAIN_SIZE);
  double *test_features = malloc(sizeof(double) * TRAIN_SIZE * FEATURE_DIM);

  printf("\nTest reading in data from binary files...\n");

  read_binary_dataset("MNIST_train.bin", TRAIN_SIZE, train_labels, train_features);
  read_binary_dataset("MNIST_test.bin", TEST_SIZE, test_labels, test_features);

  clock_t end = clock();

  free(train_features); free(train_labels);
  free(test_features); free(test_labels);

  printf("\nDone\n\n");

  printf("Total time to read in data, convert to binary and read in data from binary: %.2f secs\n",
         (double) (end - begin) / CLOCKS_PER_SEC);

  return 0;
}
