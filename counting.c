#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <omp.h>
#include <math.h>
#include <time.h>

#define KIND_ORIGINAL 1
#define KIND_NON3D 11
#define KIND_TRACO 12
#define KIND_3D 13
#define KIND_3D_16 14
#define KIND_PLUTO 15
#define KIND_DAPT 16

#define min(a, b) (((a) < (b)) ? (a) : (b))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define max(a, b) (((a) > (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define floord(n, d) floor(((double)(n)) / ((double)(d)))
#define ceild(n, d) ceil(((double)(n)) / ((double)(d)))

int **c;
int **ck;

int zz = 2;

int **F; // only ACGU

int N;
int DIM;

int *RNA; // only ACGU

int getValue(const char c)
{
  /*
   *  'A' => 0 -> bitowo 0001 -> 1
   *  'G' => 1 -> bitowo 0010 -> 2
   *  'C' => 2 -> bitowo 0100 -> 4
   *  'U' => 3 -> bitowo 1000 -> 8
  */

  if(c=='A')    return 1;
  if(c=='G')    return 2;
  if(c=='C')    return 4;
  if(c=='U')    return 8;
  return 16;
}

void rand_seq(int *seq, int size)
{
  int i, tmp;
  int time_seed = 654322; //= time(NULL);
  srand(time_seed);
  for (i = 0; i < size; i++)
  {
    tmp = rand() % 4;

    switch (tmp)
    {
    case 0:
      seq[i] = getValue('A');
      break;
    case 1:
      seq[i] = getValue('G');
      break;
    case 2:
      seq[i] = getValue('C');
      break;
    case 3:
      seq[i] = getValue('U');
      break;
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

int match(const int e1, const int e2)
{
  /*
   *  'A' => 0 -> bitowo 0001 -> 1
   *  'G' => 1 -> bitowo 0010 -> 2
   *  'C' => 2 -> bitowo 0100 -> 4
   *  'U' => 3 -> bitowo 1000 -> 8
  */
  //const bool match =
  //  (e1 == 0 && e2 == 3) || (e1 == 3 && e2 == 0) ||
  //  (e1 == 1 && e2 == 2) || (e1 == 2 && e2 == 1) ||
  //  (e1 == 1 && e2 == 3) || (e1 == 3 && e2 == 1);
  //return match;
  const int match =
    (e1 + e2 == 9) ||
    (e1 + e2 == 6) ||
    (e1 + e2 == 10) ;
  return match;

  //(e1 == "A" && e2 == "U") ||
  //(e1 == "U" && e2 == "A") ||
  //(e1 == "G" && e2 == "C") ||
  //(e1 == "C" && e2 == "G") ||
  //(e1 == "G" && e2 == "U") ||
  //(e1 == "U" && e2 == "G");

}

#define pared(i, j) (match(RNA[i], RNA[j]))

//int pared(int i, int j)
// {
//   char nt1 = RNA[i];
//   char nt2 = RNA[j];
//   if ((nt1 == 'A' && nt2 == 'U') || (nt1 == 'U' && nt2 == 'A') ||
//       (nt1 == 'G' && nt2 == 'C') || (nt1 == 'C' && nt2 == 'G') ||
//       (nt1 == 'G' && nt2 == 'U') || (nt1 == 'U' && nt2 == 'G'))
//   {

//     return 1;
//   }
//   else
//     return 0;
// }
int paired(int i, int j)
{
  return pared(i, j);
}

void test_fun(int kind);
void free2dmem(int **S, int size)
{
  for (int i = 0; i < size; i++)
    free(S[i]);
  free(S);
}
void do_one_step(int step_size, int kind)
{
  N = step_size;
  DIM = N + 1;
  F = mem(DIM); // add free to the end
  c = mem(DIM);
  ck = mem(DIM);
  for (int i = 0; i < DIM; i++)
  {
    for (int j = 0; j < DIM; j++)
    {
      c[i][j] = 1;
      ck[i][j] = 1;
    }
  }
  RNA = (int *)malloc(DIM * sizeof(char *));
  rand_seq(RNA, N);
  test_fun(kind);
  free(RNA);
  free2dmem(F, DIM);
  free2dmem(c, DIM);
  free2dmem(ck, DIM);
}
void do_one_step_all_kinds(int step_size)
{
  printf("%4d;", step_size);
  do_one_step(step_size, KIND_ORIGINAL);
  do_one_step(step_size, KIND_NON3D);
  do_one_step(step_size, KIND_PLUTO);
  do_one_step(step_size, KIND_TRACO);
  do_one_step(step_size, KIND_3D);
  do_one_step(step_size, KIND_3D_16);
  do_one_step(step_size, KIND_DAPT);
}

int main(int argc, char *argv[])
{
  int num_proc = 1;
  int kind = 1;
  srand(time(NULL));
  //int step_size = 50;
  //omp_set_num_threads(16);
  //do_one_step(step_size, KIND_ORIGINAL);
  //do_one_step(step_size, KIND_NON3D);
  //do_one_step(step_size, KIND_PLUTO);
  //do_one_step(step_size, KIND_3D);
  //do_one_step(step_size, KIND_3D_16);

  for (num_proc = 16; num_proc <= 32; num_proc*=2)
  {
    printf("num_proc = %d\nN___;ORIGIN;NON3D_;PLUTO_;TRACO_;3D____;3D_16_;DAPT__;\n", num_proc);
  
    omp_set_num_threads(num_proc);
    for (int i = 500; i <= 12000; i += 500)
    {
      do_one_step_all_kinds(i);
      printf("\n");
    }
  }
  return 0;
}

void test_fun(int kind)
{

  int i, j, k, ll, p, q, l = 1;
  int c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12;
  int h0, h1, i0, i1, i2;
  int t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
  int lb, ub, lbp, ubp, lb2, ub2;
  register int hMin, lbv, ubv;
  // printf(" -exp(Ebp/RT) = %5.3f\n", ERT);

  double start = omp_get_wtime();
  //  compute the partition functions Q and Qbp
  if (kind == KIND_ORIGINAL)
  {
#pragma scop
    for (int i = N - 2; i >= 1; i--)
    {
      for (int j = i + 2; j <= N; j++)
      {
        for (int k = i; k <= j - 1; k++)
        {
          c[i][j] += paired(k, j) ? c[i][k - 1] + c[k + 1][j - 1] : 0;
        }
        c[i][j] = c[i][j] + c[i][j - 1];
      }
    }
#pragma endscop
  }

  if (kind == KIND_NON3D)
  { //
#pragma scop
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
#pragma endscop
  }
  if(kind == KIND_3D_16)
  {
    for (int w0 = -1; w0 < floord(N - 1, 8); w0 += 1) {
      hMin = min(floord(N - 2, 16), (w0 + 1) / 2);
  #pragma omp parallel for shared(hMin)
  for (h0 = max(w0 - (N + 16) / 16 + 2, (w0 + 1) / 3); h0 <= hMin; h0 += 1) {
    for (h1 = max(max(w0 - h0 - (N + 16) / 16 + 1, h0 - (N + 14) / 16), -((N + 13) / 16)); h1 < 0; h1 += 1) {
      for (i0 = max(max(1, 16 * h0), 8 * w0 - 8 * h0); i0 <= min(min(16 * h0 + 15, 16 * w0 - 16 * h0 + 29), N + 16 * h1 + 14); i0 += 1) {
        for (i1 = max(max(16 * h1, -16 * w0 + 16 * h0 + 16 * h1 + i0 - 14), -N + i0 + 1); i1 <= min(16 * h1 + 15, -16 * w0 + 16 * h0 + 16 * h1 + 2 * i0 + 1); i1 += 1) {
          for (i2 = max(16 * w0 - 16 * h0 - 16 * h1, i0 - i1 + 1); i2 <= min(min(N, 16 * w0 - 16 * h0 - 16 * h1 + 15), 2 * i0 - i1 + 1); i2 += 1) {
            {
              if (2 * i0 >= i1 + i2) {
                c[-i1][i2] += (pared((-i0 + i2 - 1), (i2)) ? (c[-i1][-i0 + i2 - 2] + c[-i0 + i2][i2 - 1]) : 0);
              }
              c[-i1][i2] += (pared((i0 - i1), (i2)) ? (c[-i1][i0 - i1 - 1] + c[i0 - i1 + 1][i2 - 1]) : 0);
              if (i1 + i2 == i0 + 1) {
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
  if (kind == KIND_3D)
  {
    for (int w0 = floord(-N + 34, 160) - 1; w0 < floord(7 * N - 10, 80); w0 += 1)
    {
      hMin = min((N - 2) / 16, w0 + floord(N - 80 * w0 + 46, 240) + 1);
#pragma omp parallel for shared(hMin)
      for (h0 = max(max(0, w0 - (N + 40) / 40 + 2), w0 + floord(-4 * w0 - 3, 9) + 1);
           h0 <= hMin; h0 += 1)
      {
        for (h1 = max(max(max(5 * w0 - 9 * h0 - 3, -((N + 29) / 32)), w0 - h0 - (N + 40) / 40 + 1), -((N - 16 * h0 + 30) / 32)); h1 <= min(-1, 5 * w0 - 7 * h0 + 8); h1 += 1)
        {
          for (i0 = max(max(1, 16 * h0), 20 * w0 - 20 * h0 - 4 * h1); i0 <= min(min(16 * h0 + 15, N + 32 * h1 + 30), 40 * w0 - 40 * h0 - 8 * h1 + 69); i0 += 1)
          {
            for (i1 = max(max(32 * h1, -40 * w0 + 40 * h0 + 40 * h1 + i0 - 38), -N + i0 + 1); i1 <= min(32 * h1 + 31, -40 * w0 + 40 * h0 + 40 * h1 + 2 * i0 + 1); i1 += 1)
            {
              for (i2 = max(40 * w0 - 40 * h0 - 40 * h1, i0 - i1 + 1); i2 <= min(min(N, 40 * w0 - 40 * h0 - 40 * h1 + 39), 2 * i0 - i1 + 1); i2 += 1)
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
  if (kind == KIND_TRACO)
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
                    c[(N - c1 + c3 - 2)][c9] += pared(c11, c9) ? c[(N - c1 + c3 - 2)][c11 - 1] + c[c11 + 1][c9 - 1] : 0;
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
  if (kind == KIND_DAPT)
  {
    for (int w0 = floord(-N - 14, 32); w0 < floord(N, 32); w0 += 1)
    {
#pragma omp parallel for
      for (int h0 = max(-((N + 13) / 16), w0 - (N + 32) / 32 + 1); h0 <= min(-1, 2 * w0 + 2); h0 += 1)
      {
        for (int i0 = max(max(-N + 2, -32 * w0 + 32 * h0 - 29), 16 * h0); i0 <= 16 * h0 + 15; i0 += 1)
        {
          for (int i1 = max(32 * w0 - 32 * h0, -i0 + 2); i1 <= min(N, 32 * w0 - 32 * h0 + 31); i1 += 1)
          {
            {
              for (int i3 = -i0; i3 < i1; i3 += 1)
              {
                c[-i0][i1] += (paired((i3), (i1)) ? (c[-i0][i3 - 1] + c[i3 + 1][i1 - 1]) : 0);
              }
              c[-i0][i1] = (c[-i0][i1] + c[-i0][i1 - 1]);
            }
          }
        }
      }
    }
  }
  if (kind == KIND_PLUTO)
  {
    if (N >= 3)
    {
      for (t1 = 3; t1 <= N; t1++)
      {
        lbp = 0;
        ubp = floord(t1 - 2, 32);
#pragma omp parallel for private(lbv, ubv, t3, t4, t5)
        for (t2 = lbp; t2 <= ubp; t2++)
        {
          for (t3 = t2; t3 <= floord(t1, 32); t3++)
          {
            if ((t1 >= 32 * t3 + 1) && (t1 <= 32 * t3 + 31))
            {
              for (t4 = max(1, 32 * t2); t4 <= min(t1 - 2, 32 * t2 + 31); t4++)
              {
                for (t5 = max(32 * t3, t4); t5 <= t1 - 1; t5++)
                {
                  c[t4][t1] += pared(t5, t1) ? c[t4][t5 - 1] + c[t5 + 1][t1 - 1] : 0;
                  ;
                }
                c[t4][t1] = c[t4][t1] + c[t4][t1 - 1];
                ;
              }
            }
            if (t1 >= 32 * t3 + 32)
            {
              for (t4 = max(1, 32 * t2); t4 <= min(t1 - 2, 32 * t2 + 31); t4++)
              {
                for (t5 = max(32 * t3, t4); t5 <= 32 * t3 + 31; t5++)
                {
                  c[t4][t1] += pared(t5, t1) ? c[t4][t5 - 1] + c[t5 + 1][t1 - 1] : 0;
                  ;
                }
              }
            }
            if (t1 == 32 * t3)
            {
              for (t4 = max(1, 32 * t2); t4 <= min(t1 - 2, 32 * t2 + 31); t4++)
              {
                if (t1 % 32 == 0)
                {
                  c[t4][t1] = c[t4][t1] + c[t4][t1 - 1];
                  ;
                }
              }
            }
          }
        }
      }
    }
  }

  double stop = omp_get_wtime();
  printf("%.4f;", stop - start);

  if(N<100)
  {
    printf("CK:\n");
    printf("N:%d, DIM:%d C:\n", N, DIM);
    for (int i = 0; i < N; i++)
    {
      for (int j = 0; j < N; j++)
      {
        printf("%d ", c[i][j]);
      }
      printf("\n");
    }
  }
  //write c to file with the timestamp as name
  // char filename[100];
  // char timestamp[50];
  // FILE *fp;
  // time_t t = time(NULL);
  // struct tm tm = *localtime(&t);
  // sprintf(timestamp, "%d_%d_%d_%d_%d_%d", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
  // sprintf(filename, "c_%d_%d_%s.txt", N, kind, timestamp);
  // fp = fopen(filename, "w");
  // for (int i = 0; i < N; i++)
  // {
  //   for (int j = 0; j < N; j++)
  //   {
  //     fprintf(fp, "%6d ", c[i][j]);
  //   }
  //   fprintf(fp, "\n");
  // }
  // fclose(fp);
  return;

  // for (i = 0; i < DIM; i++)
  //   for (j = 0; j < DIM; j++)
  //     if (c[i][j] != ck[i][j])
  //     {
  //       printf("err: %d %d %d %d\n", i, j, c[i][j], ck[i][j]);
  //       exit(0);
  //     }
}
