#define min(a,b) (((a)<(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define ceild(n,d) ceil(((double)(n))/((double)(d)))

int ** c;int ** ck; int zz = 2; int ** F; int N; char * RNA;  
#include "mem.h"
int paired(int i, int j) {return 0;}

int main(int argc, char *argv[]){
    int ll,p,q,l=0;
    int kind=1;
    N = 8;
    int check=1;

    if(kind==1 || check){
#pragma scop
      for (int c0 = 1; c0 < N - 1; c0 += 1)
        for (int c1 = -N + c0 + 1; c1 < 0; c1 += 1)
          for (int c2 = c0 - c1 + 1; c2 <= min(N, 2 * c0 - c1 + 1); c2 += 1)
          {
            if (2 * c0 >= c1 + c2)
              c[-c1][c2] += c[-c1][c2 - 1] + paired(-c0 + c2 - 1, c2) ? c[-c1][-c0 + c2 - 1 - 1] + c[-c0 + c2 - 1 + 1][c2 - 1] : 0;
            // S_0(-c1, c2, -c0 + c2 - 1);
            c[-c1][c2] += c[-c1][c2 - 1] + paired(c0 - c1, c2) ? c[-c1][c0 - c1 - 1] + c[c0 - c1 + 1][c2 - 1] : 0;
            // S_0(-c1, c2, c0 - c1);
          }
#pragma endscop
    }
    return 0;
}

