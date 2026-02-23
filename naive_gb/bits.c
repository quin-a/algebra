#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Takes a string that represents an int "255" and returns the value integer
// value
int one_arg(char *str) {
  int len = 0;
  while (str[len] != '\0') {
    if (str[len] < '0' || str[len] > '9') {
      printf("Argv error: Please input a positive integer to see its base-2 "
             "rep!\n");
      exit(1);
    }
    len++;
  }
  int arg = 0;
  for (int i = 0; i < len; ++i) {
    arg += (int)(str[len - 1 - i] - '0') * (int)pow(10, i);
  }
  return arg;
}

// Return codes:
// -1 indicates error in the string
// 0 indicates no shift
// 1 indicates right shift
// 2 indicates left shift
int check_shift(char *str) {
  if (str[0] < '0' || str[0] > '9')
    return -1;

  int retflag = 0;
  int len = 0;
  while (str[len] != '\0') {
    if (str[len] < '0' || str[len] > '9') {
      if (str[len] == '<' && str[len + 1] == '<' && retflag == 0) {
        retflag = 2;
        len = len + 2;
        continue;
      } else if (str[len] == '>' && str[len + 1] == '>' && retflag == 0) {
        retflag = 1;
        len = len + 2;
        continue;
      } else {
        return -1;
      }
    }
    ++len;
  }
  return retflag;
}

int main(int argc, char **argv) {
  if (argc < 2) {
    printf("Argc error: Please enter an integer to see its base-2 rep!\n");
    exit(1);
  }
  // argv points to start of [ptr0, ptr1, ...]
  // argv + 1 dereferences to prt1
  // ptr1 points to start of ['2', '5', '5', ...]
  // ptr1 dereferences to char '2'

  int arg;
  if (argc == 2 && check_shift(argv[1]) == 0) {
    arg = one_arg(argv[1]);
  } else if (argc == 2 && check_shift(argv[1]) == 2) {
    // have 234<<345 or whatever,
    char *str1;
    int idx = 0;
    while (argv[1][idx] != '<') {
      str1[idx] = argv[1][idx];
      ++idx;
    }
    str1[idx] = '\0';
    char *str2;
    idx = idx + 2;
    int idy = 0;
    while (argv[1][idx] != '\0') {
      str2[idy] = argv[1][idx];
      ++idx;
      ++idy;
    }
    str2[idy] = '\0';
    arg = one_arg(str1) << one_arg(str2);
  } else if (argc == 2 && check_shift(argv[1]) == 1) {
    // have 234>>345 or whatever
    char *str1;
    int idx = 0;
    while (argv[1][idx] != '>') {
      str1[idx] = argv[1][idx];
      ++idx;
    }
    str1[idx] = '\0';
    char *str2;
    idx = idx + 2;
    int idy = 0;
    while (argv[1][idx] != '\0') {
      str2[idy] = argv[1][idx];
      ++idx;
      ++idy;
    }
    str2[idy] = '\0';
    arg = one_arg(str1) >> one_arg(str2);
  } else if (argc == 4) {
    // should be 234 >> 345 or 234 << 345
    int arg1 = one_arg(argv[1]);
    int arg2 = one_arg(argv[3]);
    if (argv[2][0] == '>' && argv[2][1] == '>' && argv[2][2] == '\0') {
      arg = arg1 >> arg2;
    } else if (argv[2][0] == '<' && argv[2][1] == '<' && argv[2][2] == '\0') {
      arg = arg1 << arg2;
    } else {
      printf("Error in phow to see dual pane gdbarsing multiple arguments!\n");
      exit(1);
    }
  } else {
    // error state
    printf("Error: exmple usage: \"./bits 10\" or \"./bits 10<<1\" or \"./bits "
           "10 << 1\"\n");
    exit(1);
  }

  int idx = 0;
  char base2[64];
  while (arg != 0) {
    base2[idx] = (char)(arg & 1) + '0';
    arg = arg >> 1;
    ++idx;
  }

  while (idx > 0) {
    printf("%c", base2[--idx]);
  }
  printf("\n");
  return 0;
}
