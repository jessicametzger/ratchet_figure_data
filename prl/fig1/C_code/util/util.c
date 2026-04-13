#include <stdbool.h>
#define _USE_MATH_DEFINES
#include <math.h>

#define M_PI_HALF (M_PI / 2.0)
#define M_TWO_PI (M_PI * 2.0)

void ERROR(char* msg){
  printf("%s",msg);
  exit(1);
}


void split_vec(int Nexpect, double** output, const char* input, char delimiter){
    int count = 1;
    for (const char *c = input; *c != '\0'; c++) if (*c == delimiter) count ++;
    if (count!=Nexpect) {
      char msg[200];
      sprintf(msg,"input list has length %d, not equal to Nexpect=%d\n", count, Nexpect);
      ERROR(msg);
    }
    
    char *temp = strdup(input); // Make a copy of the input string
    char *token = strtok(temp, ","); // Split by commas
    count = 0; // Initialize count of doubles
    while (token != NULL) {
        output[0][count]  = strtod(token, NULL); // Convert token to double
        token             = strtok(NULL, ","); // Get the next token
        count++;
    }
}
