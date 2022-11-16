#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <omp.h>
#include <math.h>
#include <time.h>


#define min(a,b) (((a)<(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define ceild(n,d) ceil(((double)(n))/((double)(d)))

int ** c;
int ** ck;

int zz = 2;

int ** F;  //only ACGU

int N;
int DIM;

char * RNA;  //only ACGU


void rand_seq(char*a, int size){
  int i, tmp;
  srand(time(NULL));
  for(i=0; i<size; i++)
  {
      tmp = rand()%4;

      switch(tmp){
          case 0 : a[i] = 'A'; break;
          case 1 : a[i] = 'G'; break;
          case 2 : a[i] = 'C'; break;
          case 3 : a[i] = 'U'; break;
      }

  }

}

int **mem(int size)
{
  int i;
  int **S;
  S = (int **)malloc(size * sizeof(int *));
  for (i = 0; i < size; i++)
    S[i] = (int *)malloc(size * sizeof(int));

  return S;
}

int pared(int i, int j) {
   char nt1 = RNA[i];
   char nt2 = RNA[j];
         if ((nt1 == 'A' && nt2 == 'U') || (nt1 == 'U' && nt2 == 'A') ||
             (nt1 == 'G' && nt2 == 'C') || (nt1 == 'C' && nt2 == 'G') ||
             (nt1 == 'G' && nt2 == 'U') || (nt1 == 'U' && nt2 == 'G')){

            return 1;}
         else
            return 0;
}
int paired(int i, int j) {
    return pared(i, j);
}

void test_fun(int kind);
void do_one_step(int step_size, int kind)
{
    N = step_size;
    DIM = N+1;
    F =  mem(DIM);//add free to the end
    c = mem(DIM);
    ck = mem(DIM);
   for(int i=0; i<DIM; i++)
   {
    for(int j=0; j<DIM; j++){
     c[i][j] = 1;
     ck[i][j] = 1;
    }
   }
    RNA =  (char*) malloc(DIM * sizeof(char*)); 
    rand_seq(RNA, N);
    test_fun(kind);
    free(RNA);
}
void do_one_step_all_kinds(int step_size)
{
  do_one_step(step_size, 1);
  //do_one_step(step_size, 5);
  //do_one_step(step_size, 11);
  do_one_step(step_size, 12);
  do_one_step(step_size, 13);
}

int main(int argc, char *argv[])
{
  int num_proc = 1;
  int kind = 1;
  srand(time(NULL));
  for (num_proc = 1; num_proc <= 32; num_proc++)
  {
    printf("num_proc = %d\n", num_proc);
    omp_set_num_threads(num_proc);
    for (int i = 500; i <= 3000; i += 500)
    {
      do_one_step_all_kinds(i);
          printf("\n");
    }
  }
  return 0;
}

void test_fun(int kind)
{


    int i,j,k,ll,p,q,l=1;
    int c0, c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12;

    int t1, t2, t3, t4, t5, t6,t7,t8,t9,t10;
    int lb, ub, lbp, ubp, lb2, ub2;
    register int lbv, ubv;
    //printf(" -exp(Ebp/RT) = %5.3f\n", ERT);

    double start = omp_get_wtime();
    //  compute the partition functions Q and Qbp
    if(kind==1 ){
        #pragma scop
        for (int i = N - 2; i >= 1; i--) {
            for (int j = i + 2; j <= N; j++) {
                for (int k = i; k <= j - 1; k++) {
                    c[i][j] += paired(k, j) ? c[i][k - 1] + c[k + 1][j - 1] : 0;
                }
                c[i][j] = c[i][j] + c[i][j - 1];
            }
        }
       #pragma endscop
    }

    if (kind == 11)
    {
#pragma scop
      for (i = N - 2; i >= 1; i--)
      {
        for (j = i + 2; j <= N; j++)
        {
          for (k = i; k <= j - l; k++)
          {
            c[i][j] += pared(k-1, j-1) ? c[i][k - 1] + c[k + 1][j - 1] : 0;
          }
          c[i][j] = c[i][j] + c[i][j - 1];
        }
      }
#pragma endscop
    }

    if (kind == 12)
    {
      for (int w0 = floord(-N + 34, 160) - 1; w0 < floord(7 * N - 10, 80); w0 += 1)
      {
        int hMin = min((N - 2) / 16, w0 + floord(N - 80 * w0 + 46, 240) + 1);
#pragma omp parallel for
        for (int h0 = max(max(0, w0 - (N + 40) / 40 + 2), w0 + floord(-4 * w0 - 3, 9) + 1); 
        h0 <= hMin ; h0 += 1)
        {
          for (int h1 = max(max(max(5 * w0 - 9 * h0 - 3, -((N + 29) / 32)), w0 - h0 - (N + 40) / 40 + 1), -((N - 16 * h0 + 30) / 32)); h1 <= min(-1, 5 * w0 - 7 * h0 + 8); h1 += 1)
          {
            for (int i0 = max(max(1, 16 * h0), 20 * w0 - 20 * h0 - 4 * h1); i0 <= min(min(16 * h0 + 15, N + 32 * h1 + 30), 40 * w0 - 40 * h0 - 8 * h1 + 69); i0 += 1)
            {
              for (int i1 = max(max(32 * h1, -40 * w0 + 40 * h0 + 40 * h1 + i0 - 38), -N + i0 + 1); i1 <= min(32 * h1 + 31, -40 * w0 + 40 * h0 + 40 * h1 + 2 * i0 + 1); i1 += 1)
              {
                for (int i2 = max(40 * w0 - 40 * h0 - 40 * h1, i0 - i1 + 1); i2 <= min(min(N, 40 * w0 - 40 * h0 - 40 * h1 + 39), 2 * i0 - i1 + 1); i2 += 1)
                {
                  {
                    if (2 * i0 >= i1 + i2)
                    {
                      c[-i1][i2] += (pared((-i0 + i2 - 1), (i2)) ? (c[-i1][-i0 + i2 - 2] + c[-i0 + i2][i2 - 1]) : 0);
                    }
                    c[-i1][i2] += (pared((i0 - i1), (i2)) ? (c[-i1][i0 - i1 - 1] + c[i0 - i1 + 1][i2 - 1]) : 0);
                    if (i1 + i2 == i0 + 1)
                    {
                      c[-i1][i0 - i1 + 1] = (c[-i1][i0 - i1 + 1] + c[-i1][i0 - i1]);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    // c[i][j] +=  c[i][j-1] + paired(k,j) ?  c[i][k-1] + c[k+1][j-1] : 0;
    if (kind == 13)
    {
      for (c1 = 0; c1 < N + floord(N - 3, 128) - 2; c1 += 1)
#pragma omp parallel for
        for (c3 = max(0, -N + c1 + 3); c3 <= c1 / 129; c3 += 1)
          for (c4 = 0; c4 <= 1; c4 += 1)
          {
            if (c4 == 1)
            {
              for (c9 = N - c1 + 129 * c3; c9 <= min(N, N - c1 + 129 * c3 + 127); c9 += 1)
                for (c10 = max(0, -c1 + 64 * c3 - c9 + (N + c1 + c3 + c9 + 1) / 2 + 1); c10 <= 1; c10 += 1)
                {
                  if (c10 == 1)
                  {
                    c[(N - c1 + c3 - 2)][c9] = c[(N - c1 + c3 - 2)][c9] + c[(N - c1 + c3 - 2)][c9 - 1];
                  }
                  else
                  {
                    for (c11 = N - c1 + 129 * c3 + 1; c11 < c9; c11 += 1)
                      c[(N - c1 + c3 - 2)][c9] += paired(c11, c9) + c[(N - c1 + c3 - 2)][c11 - 1] + c[c11 + 1][c9 - 1] + 0;
                  }
                }
            }
            else
            {
              for (c5 = 0; c5 <= 8 * c3; c5 += 1)
                for (c9 = N - c1 + 129 * c3; c9 <= min(N, N - c1 + 129 * c3 + 127); c9 += 1)
                  for (c11 = N - c1 + c3 + 16 * c5 - 2; c11 <= min(min(N - c1 + 129 * c3, N - c1 + c3 + 16 * c5 + 13), c9 - 1); c11 += 1)
                    c[(N - c1 + c3 - 2)][c9] += paired(c11, c9) + c[(N - c1 + c3 - 2)][c11 - 1] + c[c11 + 1][c9 - 1] + 0;
            }
          }
    }
    if (kind == 5)
    {
      for (int c0 = 1; c0 < N - 1; c0 += 1)
        for (int c1 = -N + c0 + 1; c1 < 0; c1 += 1)
          for (int c2 = c0 - c1 + 1; c2 <= min(N, 2 * c0 - c1 + 1); c2 += 1)
          {
            if (2 * c0 >= c1 + c2)
              c[-c1][c2] += c[-c1][c2 - 1] + pared(-c0 + c2 - 1, c2) ? c[-c1][-c0 + c2 - 1 - 1] + c[-c0 + c2 - 1 + 1][c2 - 1] : 0;
            // S_0(-c1, c2, -c0 + c2 - 1);
            c[-c1][c2] += c[-c1][c2 - 1] + pared(c0 - c1, c2) ? c[-c1][c0 - c1 - 1] + c[c0 - c1 + 1][c2 - 1] : 0;
            // S_0(-c1, c2, c0 - c1);
          }
    }

    double stop = omp_get_wtime();
    printf("%.4f;",stop - start);

}

