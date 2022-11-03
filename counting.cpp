// set check to 0 to measure time, set to 1 to check validity

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <omp.h>
#include <math.h>


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


#include "mem.h"

int paired(int i, int j) {
   char nt1 = RNA[i];
   char nt2 = RNA[j];
         if ((nt1 == 'A' && nt2 == 'U') || (nt1 == 'U' && nt2 == 'A') ||
             (nt1 == 'G' && nt2 == 'C') || (nt1 == 'C' && nt2 == 'G') ||
             (nt1 == 'G' && nt2 == 'U') || (nt1 == 'U' && nt2 == 'G')){

            return 1;}
         else
            return 0;
}



int main(int argc, char *argv[]){



    int num_proc=1;
    int i,j,k,ll,p,q,l=0;
    int c0, c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12;

    int t1, t2, t3, t4, t5, t6,t7,t8,t9,t10;
    int lb, ub, lbp, ubp, lb2, ub2;
    register int lbv, ubv;




    srand(time(NULL));



    if(argc > 1)
        num_proc = atoi(argv[1]);

    int kind=1;

    N = 8;
    DIM = 12;
    if(argc > 2)
        N = atoi(argv[2]);
    DIM = N+10;


    if(argc > 3)
        kind = atoi(argv[3]);



    omp_set_num_threads(num_proc);
    //printf(" -exp(Ebp/RT) = %5.3f\n", ERT);

    F =  mem();
    c = mem();
    ck = mem();

   for(i=0; i<DIM; i++)
    for(i=0; i<DIM; i++){
     c[i][j] = i+j;
     ck[i][j] = i+j;
    }

    RNA =  (char*) malloc(DIM * sizeof(char*));  //read from FASTA file
    rand_seq(RNA, N);

//for(i=0; i<DIM; i++)
 //printf("%c", RNA[i]);

printf("\n");
     int check=1;


    double start = omp_get_wtime();
    //  compute the partition functions Q and Qbp
    if(kind==1 || check){
        #pragma scop
        for (i = N-2;  i>=1; i--){
           for ( j=i+2; j<= N; j++){
             for ( k = i; k<=j-l; k++){
              ck[i][j] +=  ck[i][j-1] + paired(k,j) ?  ck[i][k-1] + ck[k+1][j-1] : 0;
           }
          }
        }
       #pragma endscop
    }
    if(kind==2) // pluto
    {
        /* Start of CLooG code */
/* Start of CLooG code */
if ((N >= 3) && (N >= l+1)) {
  for (t1=max(3,l+1);t1<=N;t1++) {
    lbp=0;
    ubp=min(floord(t1-2,16),floord(t1-l,16));
#pragma omp parallel for private(lbv,ubv,t3,t4,t5)
    for (t2=lbp;t2<=ubp;t2++) {
      for (t3=t2;t3<=floord(t1-l,16);t3++) {
        for (t4=max(1,16*t2);t4<=min(min(t1-2,t1-l),16*t2+15);t4++) {
          for (t5=max(16*t3,t4);t5<=min(t1-l,16*t3+15);t5++) {
            c[t4][t1] += c[t4][t1-1] + paired(t5,t1) ? c[t4][t5-1] + c[t5+1][t1-1] : 0;;
          }
        }
      }
    }
  }
}
/* End of CLooG code */


        /* End of CLooG code */
    }
    if(kind==3) // traco
    {

 

for( c1 = max(0, floord(l - 2, 8) - 1); c1 <= floord(N - 3, 8); c1 += 1)
   #pragma omp parallel for shared(c1) private(c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12) schedule(dynamic, 1)
  for( c3 = max(max(0, floord(l - 2, 16)), c1 - (N + 13) / 16 + 1); c3 <= c1 / 2; c3 += 1)
    for( c5 = 0; c5 <= min(c3 + floord(-l + 1, 16) + 1, floord(-l + N - 1, 16)); c5 += 1)
      for( c7 = max(-N + 16 * c1 - 16 * c3 + 2, l - N + 16 * c5); c7 <= min(-1, -N + 16 * c1 - 16 * c3 + 17); c7 += 1) {
        for( c9 = max(l + 16 * c5 - c7, 16 * c3 - c7 + 2); c9 <= min(min(N, l + 16 * c5 - c7 + 16), 16 * c3 - c7 + 17); c9 += 1) {
          if (c7 + c9 >= 16 * c3 + 3 && c7 + c9 >= l + 16 * c5 + 1)
            for( c11 = -c7; c11 < 16 * c5 - c7; c11 += 1)
              c[(-c7)][c9] += c[(-c7)][c9-1] + paired(c11,c9) ? c[(-c7)][c11-1] + c[c11+1][c9-1] : 0;
          for( c11 = 16 * c5 - c7; c11 <= min(16 * c5 - c7 + 15, -l + c9); c11 += 1)
            c[(-c7)][c9] += c[(-c7)][c9-1] + paired(c11,c9) ? c[(-c7)][c11-1] + c[c11+1][c9-1] : 0;
        }
        if (16 * c3 >= l + 16 * c5 + 15)
          for( c11 = 16 * c5 - c7; c11 <= 16 * c5 - c7 + 15; c11 += 1)
            c[(-c7)][(16*c3-c7+2)] += c[(-c7)][(16*c3-c7+2)-1] + paired(c11,(16*c3-c7+2)) ? c[(-c7)][c11-1] + c[c11+1][(16*c3-c7+2)-1] : 0;
      }




    }

    if(kind==4) // traco
    {
       for( c0 = max(0, floord(l - 2, 8) - 1); c0 <= floord(N - 3, 8); c0 += 1)
  #pragma omp parallel for  shared(c0) private(c1, c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12) schedule(dynamic, 1)
  for( c1 = (c0 + 1) / 2; c1 <= min(min(c0, c0 + (-l + 1/16) + 1), (N - 3) / 16); c1 += 1)
    for( c3 = max(l, 16 * c0 - 16 * c1 + 2); c3 <= min(N - 1, 16 * c0 - 16 * c1 + 17); c3 += 1)
      for( c4 = max(0, -c1 + (N - 1) / 16 - 1); c4 <= min((-l + N) / 16, -c1 + (-l + N + c3 - 2) / 16); c4 += 1)
        for( c6 = max(max(-N + 16 * c1 + 2, -N + c3), -16 * c4 - 15); c6 <= min(min(-1, -N + 16 * c1 + 17), -l + c3 - 16 * c4); c6 += 1)
          for( c10 = max(16 * c4, -c6); c10 <= min(16 * c4 + 15, -l + c3 - c6); c10 += 1)
            c[(-c6)][(c3-c6)] += c[(-c6)][(c3-c6)-1] + paired(c10,(c3-c6)) ? c[(-c6)][c10-1] + c[c10+1][(c3-c6)-1] : 0;
       
   if(1==0)    
    for( c0 = max(0, l + floord(l - 2, 16) - 2); c0 < N + floord(N - 3, 16) - 2; c0 += 1)
  #pragma omp parallel for shared(c0) private(c1, c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12) schedule(dynamic, 1)
  for( c1 = c0 - (c0 + 17) / 17 + 1; c1 <= min(min(N - 3, c0), c0 + (-l + 1)/16 + 1); c1 += 1)
    for( c3 = max(l, 16 * c0 - 16 * c1 + 2); c3 <= min(c1 + 2, 16 * c0 - 16 * c1 + 17); c3 += 1)
      for( c4 = (N - c1 - 2) / 16; c4 <= (-l + N - c1 + c3 - 2) / 16; c4 += 1)
        for( c10 = max(N - c1 - 2, 16 * c4); c10 <= min(-l + N - c1 + c3 - 2, 16 * c4 + 15); c10 += 1)
          c[(N-c1-2)][(N-c1+c3-2)] += c[(N-c1-2)][(N-c1+c3-2)-1] + paired(c10,(N-c1+c3-2)) ? c[(N-c1-2)][c10-1] + c[c10+1][(N-c1+c3-2)-1] : 0;
    }





    double stop = omp_get_wtime();
    printf("%.4f\n",stop - start);

    //printf("Q\n");
    //rna_array_print(Q);
    //printf("Qbp\n");
    //rna_array_print(Qbp);

    if(check)
    for(i=0; i<DIM; i++)
    for(j=0; j<DIM; j++)
     if(c[i][j] != ck[i][j]){
        printf("err: %d %d %d %d\n", i, j,c[i][j], ck[i][j]);
        exit(0);
     }


    return 0;

}


