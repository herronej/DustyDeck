/*********************************************
 * The converted dusty deck in C
 * Prof. A.J. Pounds
 * Department of Computer Science
 * Mercer University
 * Spring 2010
 * 
 * Code modified in Spring 2016 to use timing
 * library.
 *********************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>


/* Function Prototypes */

double trig( int, int );
double conrand ( double* );
void idcheck( int , double*, double* , double* , double* );
double walltime_();
double cputime_();

/* Change Problem Size Here */

int const MAXDIM = CPPFLAGS;


int main(){

    int *IA, N;
    double *AV, *BV, *CV;
    double *OP, *ID, *AM;
    double *BM, *CM, *DM; 
    double check, BOT, TOP, HOLDA, HOLDB, TRACE3;
    float start, finish;
    double seed;  /* for random number generator */
    double wall, cpu;
    float sum;

    int i, j, k;   /* loop indices */

    IA = calloc( MAXDIM, sizeof(int));
    AV = calloc( MAXDIM, sizeof(double));
    BV = calloc( MAXDIM, sizeof(double));
    CV = calloc( MAXDIM, sizeof(double));
    OP = calloc( MAXDIM*MAXDIM, sizeof(double));
    ID = calloc( MAXDIM*MAXDIM, sizeof(double));
    AM = calloc( MAXDIM*MAXDIM, sizeof(double));
    BM = calloc( MAXDIM*MAXDIM, sizeof(double));
    CM = calloc( MAXDIM*MAXDIM, sizeof(double));
    DM = calloc( MAXDIM*MAXDIM, sizeof(double));

    N = MAXDIM;

    wall = walltime_();
    cpu = cputime_();

    seed = 1.00000000000000;

    //     Fill arrays

    // Loop 10 Series -- Filling Arrays

    for ( i = 0; i<N; i++) {
        *(AV+i) = jn(0, (double) conrand(&seed) * 
                pow( -1.0, (int) (10*conrand(&seed)) % N ) );
    }

    for ( i = 0; i<N; i++) {
        *(BV+i) = jn(1, (double) conrand(&seed) * 
                pow( -1.0, (int) (10*conrand(&seed)) % N ) );
    }

    check = 0.0;
    for (i=0; i<N; i++){
        check = check + *(AV+i) * *(BV+i);
        idcheck(N,&check,AV,BV,ID);
    }  


    // Compute |AV><BV|

    for ( i = 0; i<N; i++) {
        for (j=0; j<N; j++ ){ 
            idcheck(N,&check,AV,BV,ID);
            if ( check > 0.5 ) { 
                *(OP+(i*N+j)) = *(AV+i) * *(BV+j) / *(BV+i);
            }
            else {
                *(OP+(i*N+j)) = *(AV+j) * *(BV+i) / *(BV+j);
            }
        } 
        *(IA+i) = i+1;
    }


    for (i=0;i<N;i++){ 
        for ( j = -1; j<=i; j+=8){  
            *(IA+i) =  ( ((i+1)+(j+1)) % N) % N ;  ; 
            // Had to drop +1 to make so arrays would not go out of bounds.  If
            // IA is used as an argument, need to add one to make it consistent with
            // fortran 
        }
    } 

    //! Loop 20 

    for ( i = 0; i< N; i++) {
        idcheck(N,&check,AV,BV,ID);
        //CV(IA(I)) = (AV(IA(I)) + BV(IA(I))) / check
        *(CV+ *(IA+i)) = ( *(AV+ *(IA+i)) + *(BV+ *(IA+i))) / check;
    }


    // ! Loop 30 

    for ( i = 1; i< N; i++ ){
        idcheck(N,&check,AV,BV,ID);
        *(AV+i) = *(AV+i-1) * *(BV+i) + *(CV+i);
    }


    //! Loop 40 

    for( i = 0; i<N; i++ ){
        idcheck(N,&check,AV,BV,ID);
        for (j=0; j<N; j++) { 
            if ( check > 0.5 ) { 
                BOT = *(OP+i*N+j); 
                TOP = *(AV+j) * *(BV+j);
                HOLDA = *(AV+j);
                *(AV+j) = *(BV+j) + *(CV+j) / (TOP-BOT) * *(ID+i*N+i);
                *(BV+j) = HOLDA + *(CV+j) / (TOP-BOT) * *(ID+j*N+j);
                //         exit(1);
                *(AM+i*N+j) = *(AV+j) * trig(*(IA+i)+1,*(IA+j)+1); 
                *(BM+i*N+j) = *(BV+j) * trig(*(IA+j)+1,*(IA+i)+1); 
                // Had to add 1 to index values in trig call to match fortran
            }
            else {
                BOT = *(OP+i*N+j); 
                TOP = *(AV+j) * *(BV+j);
                HOLDA = *(AV+j);
                *(AV+j) = *(BV+j) - *(CV+j) / (TOP-BOT) * *(ID+j*N+j);
                *(BV+j) = HOLDA - *(CV+j) / (TOP-BOT) * *(ID+i*N+i);
                *(AM+i*N+j) = *(AV+j) / trig(*(IA+i)+1,*(IA+j)+1);
                *(BM+i*N+j) = *(BV+j) / trig(*(IA+j)+1,*(IA+i)+1); 
                // Had to add 1 to index values in trig call to match fortran
            } 
        }
    }


    //! Loop 50

    for (i=0; i<N; i++ ){
        for (j=0; j<N; j++ ) { 
            *(CM+i*N+j) = 0.0;
            for ( k = 0; k<N; k++ ){ 
                if ( i < j ) { 
                    *(CM+i*N+j) = *(CM+i*N+j) - *(AM+i*N+k) * *(BM+k*N+j) / check; 
                }
                else {
                    *(CM+i*N+j) = *(CM+i*N+j) + *(AM+i*N+k) * *(BM+k*N+j) / check; 
                } 
            }
        }
    }


    // ! Loop 60

    for (i=0; i<N; i++ ) {
        for ( j=0; j<N; j++ ) {
            sum = 0.0;
            for ( k=0; k<N; k++ ) { 
                sum = sum + *(CM+i*N+k) * *(AM+j*N+k);
            } 
            *(DM+i*N+j) = sum;
        } 
    }

    for (i=0; i<N; i++ ) {
        for ( j=0; j<N; j++ ) {
            *(CM+i*N+j) = *(DM+i*N+j);
        } 
    }

    //! Loop 70

    for (i=0; i<N; i++ ) {
        for ( j=0; j<N; j++ ) {
            sum = 0.0;
            for (k=0; k<N; k++) { 
                sum = sum - *(CM+i*N+k) * *(BM+j*N+k);
            } 
            *(DM+i*N+j) = sum;
        }
    }

    HOLDA = fabs(*(AM+0));
    HOLDB = fabs(*(BM+0));

    for (i=0; i<N; i++ ) {
        for ( j=0; j<N; j++ ) {
            HOLDA = fmax(HOLDA,fabs(*(AM+i*N+j)));
            HOLDB = fmax(HOLDB,fabs(*(BM+i*N+j)));
        }
    }

    TRACE3 = 0.0;

    //! Loop 80

    for (i=0; i<N; i++) {
        TRACE3 = TRACE3 + (*(AM+*(IA+i)*N+*(IA+i)) + *(BM+*(IA+i)*N+*(IA+i))  
                - *(DM+*(IA+i)*N+*(IA+i)))  / (HOLDA * HOLDB);
    }

    cpu = cputime_() - cpu;
    wall = walltime_() - wall;

#ifndef DO_TIMING
    printf("Final trace = %1.15e and IDCHECK %.17f\n", TRACE3, check );
    printf(" -- RUNTIME -> %f seconds\n", cpu);
#else
    printf("%5i     %9.4f\n", MAXDIM, cpu);
#endif

}        

double trig ( int i, int j ){

    double x, y, z;
    float pi;
    pi = (float) acos(-1.0);
    x = (double) i - (double) j;
    y = (double) i + (double) j;
    // Had to CAREFULLY match results based on mixed precision arithmetic
    z = exp( sin(sqrt(x*x+y*y)*pi));
    return  x +  y + log10(fabs(1+z+(x*y*z))) / (fabs(x)+fabs(y));
}


void idcheck( int N, double *check, double *AV, double *BV, double *ID) {

    double l2;
    double check2;
    double a, b, c, d;

    double tempa, tempb;

    int i, j, k;


    for (i=0;i<N;i++){
        for (j=0;j<N;j++){
            if ( i == j ) { 
                if ( ( *(AV+i) < 0.0 ) && ( *(BV+j) < 0 )) {
                    *(ID+(i*N)+j) = 1.0;
                }
                else if ( ( *(AV+i) < 0.0 ) && ( *(BV+j) > 0 )) {
                    *(ID+(i*N)+j) = -1.0;
                }
                else if ( ( *(AV+i) > 0.0 ) && ( *(BV+j) < 0 )) {
                    *(ID+(i*N)+j) = -1.0;
                }
                else {
                    *(ID+(i*N)+j) = 1.0;
                }
            }
            else {

                *(ID+(i*N)+j) =   cos( *check+ 2*(i+1)* acosf(-1.0)/N) +
                    2.0* sin( *check+ 2*(j+1)* acosf(-1.0)/N);

            }
        }
    } 


    l2 = 0.0;
    for ( i=0; i<N; i++)
        l2 = l2 + *(AV+i) * *(AV+i); 


    l2 = sqrt(l2);
    for ( i=0; i<N; i++)
        *(AV+i) = *(AV+i) / l2;


    l2 = 0.0;
    for ( i=0; i<N; i++)
        l2 = l2 + *(BV+i) * *(BV+i); 


    l2 = sqrt(l2);
    for ( i=0; i<N; i++)
        *(BV+i) = *(BV+i) / l2;


    a = 0.0;
    b = 0.0;
    c = 0.0;
    d = 0.0;

    for( i=0; i<N; i++) {
        for ( j=0; j<N; j++ ) { 
            for ( k=0; k<N; k++ ) { 
                switch ( (((i+1)+(j+1)+(k+1)) % 4 ) + 1) {
                    case 1:
                        a  = a +  *(AV+i) * *(BV+j) * *(ID+(j*N+k)); 
                        *check = *check + a;
                        break; 
                    case 2:
                        b  = b +  *(AV+j) * *(BV+i) * *(ID+(k*N+j)); 
                        *check = *check - b; 
                        break; 
                    case 3:
                        c  = c -  *(AV+i) * *(BV+j) * *(ID+(k*N+j)); 
                        *check =  sqrt( (b*b) +  (c*c));
                        break;
                    case 4:
                        d  = d -  *(AV+j) * *(BV+i) * *(ID+(j*N+k)); 
                        check2 = a + b + c + d;
                        break;
                }
            }
        }
    }

    *check = fmin(fabs(check2),fabs(*check))/fmax(fabs(check2),fabs(*check));           


}

double conrand(double *seed) {

    double a, m;
    double temp;
    a = 16807.00000000000;
    m = 2147483647.000000;
    temp = a* *seed; 
    *seed = temp - m * (int) (temp/m);
    return *seed / m; 
}
