#include <stdbool.h>
#define _USE_MATH_DEFINES
#include <math.h>

#define M_PI_HALF (M_PI / 2.0)
#define M_TWO_PI (M_PI * 2.0)

void ERROR(char* msg){
  printf("%s",msg);
  exit(1);
}

/// whether point x,y is inside triangle defined by given coordinates
bool inside_triangle(double x, double y, double x1, double y1, double x2, double y2, double x3, double y3){
  double a,b,c,d;
  d = (y2-y3)*(x1-x3) - (x2-x3)*(y1-y3);
  a = ((y2-y3)*(x-x3) - (x2-x3)*(y-y3))/d;
  b = ((x1-x3)*(y-y3) - (y1-y3)*(x-x3))/d;
  c = 1-a-b;
  return (a>=0) && (b>=0) && (c>=0);
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
