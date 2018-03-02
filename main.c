#include "Header_Files/Headers.h"
#include "Header_Files/Defs.h"
#include "Header_Files/helping_procedures.h"

int main(int argc, char **argv)
{
  if(argc < 3 || strtol(argv[1], NULL, 0) >= 4 || strtol(argv[1], NULL, 0) <= 0 ||
      strtol(argv[2], NULL, 0) >= 4 || strtol(argv[2], NULL, 0) <= 0) {
    print_execution_error_message();
    return 1;
  }

  srand(time(NULL));

  int data_set = (int) strtol(argv[2], NULL, 0);
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

  int clustering_algorithm = (int) strtol(argv[1], NULL, 0);
  cluster_and_search(clustering_algorithm, data_set,
                     normalize_data,
                     train_feature_data, test_feature_data,
                     train_non_feature_data, test_non_feature_data);

  return 0;
}
