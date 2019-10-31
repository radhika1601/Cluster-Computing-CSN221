#include <math.h>
#include <string.h>
#include <stdio.h>
#ifdef _OPENACC
#include <openacc.h>
#endif /*_OPENACC*/
#include "common.h"

#include <mpi.h>

#define N 10000
#define M 10000

float A[N][M];
float Aref[N][M];
float newArray[N][M];

int main(int argc, char** argv)
{
    int iterator_max = 1000;
    
    const float pi  = 2.0 * asinf(1.0f);
    const float tol = 1.0e-5f;

    int rank = 0;
    int size = 1;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    memset(A, 0, N * M * sizeof(float));
    memset(Aref, 0, N * M * sizeof(float));
    
    // set boundary conditions
    for (int j = 0; j < N; j++)
    {
        float y0     = sinf( 2.0 * pi * j / (N-1));
        A[j][0]      = y0;
        A[j][M-1]    = y0;
        Aref[j][0]   = y0;
        Aref[j][M-1] = y0;
    }
    
#if _OPENACC
    int ngpus=acc_get_num_devices(acc_device);
    int devicenum=rank%ngpus;
    acc_set_device_num(devicenum,acc_device);

    // Call acc_init after acc_set_device_num to avoid multiple contexts on device 0 in multi GPU systems
    acc_init(acc_device);
#endif /*_OPENACC*/

    // Ensure correctness if N%size != 0
    int chunk_size = ceil( (1.0*N)/size );
    
    int jacobi_start = rank * chunk_size;
    int jacobi_end   = jacobi_start + chunk_size;
    
    // Do not process boundaries
    jacobi_start = max( jacobi_start, 1 );
    jacobi_end = min( jacobi_end, N - 1 );
    
    if ( rank == 0) printf("Jacobi relaxation Calculation: %d x %d mesh\n", N, M);

    if ( rank == 0) printf("Calculate reference solution and time serial execution.\n");
    StartTimer();
    laplace2d_serial( rank, iterator_max, tol );
    double runtime_serial = GetTimer();

    //Wait for all processes to ensure correct timing of the parallel version
    MPI_Barrier( MPI_COMM_WORLD );
    if ( rank == 0) printf("Parallel execution.\n");
    StartTimer();
    int iterator  = 0;
    float error = 1.0f;
    
    #pragma acc data copy(A) create(newArray)
    while ( error > tol && iterator < iterator_max )
    {
        error = 0.f;

        #pragma acc kernels
        for (int j = jacobi_start; j < jacobi_end; j++)
        {
            for( int i = 1; i < M-1; i++ )
            {
                newArray[j][i] = 0.25f * ( A[j][i+1] + A[j][i-1]
                                     + A[j-1][i] + A[j+1][i]);
                error = fmaxf( error, fabsf(newArray[j][i]-A[j][i]));
            }
        }
        float globalerror = 0.0f;
        MPI_Allreduce( &error, &globalerror, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD );
        error = globalerror;

        #pragma acc kernels
        for( int i = 1; i < M-1; i++ )
        {
                A[jacobi_start][i] = newArray[jacobi_start][i];
                A[jacobi_end-1][i] = newArray[jacobi_end-1][i];
        }
        
        #pragma acc kernels async
        for (int j = (jacobi_start+1); j < (jacobi_end-1); j++)
        {
            for( int i = 1; i < M-1; i++ )
            {
                A[j][i] = newArray[j][i];
            }
        }

        //Periodic boundary conditions
        int top    = (rank == 0) ? (size-1) : rank-1;
        int bottom = (rank == (size-1)) ? 0 : rank+1;

        #pragma acc host_data use_device( A )
        {
            //1. Sent row jacobi_start (first modified row) to top receive lower boundary (jacobi_end) from bottom
            MPI_Sendrecv( A[jacobi_start], M, MPI_FLOAT, top   , 0, A[jacobi_end], M, MPI_FLOAT, bottom, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );

            //2. Sent row (jacobi_end-1) (last modified row) to bottom receive upper boundary (jacobi_start-1) from top
            MPI_Sendrecv( A[(jacobi_end-1)], M, MPI_FLOAT, bottom, 0, A[(jacobi_start-1)], M, MPI_FLOAT, top   , 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
        }
        
        #pragma acc wait 
        if(rank == 0 && (iterator % 100) == 0) printf("%5d, %0.6f\n", iterator, error);
        
        iterator++;
    }
   // MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier( MPI_COMM_WORLD );
    int myid, numprocs;
    int namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];

    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Get_processor_name(processor_name, &namelen);

    fprintf(stdout, "Process %d of %d is on %s\n", myid, numprocs, processor_name);
    fflush(stdout);

    double runtime = GetTimer();
    
    if (check_results( rank, jacobi_start, jacobi_end, tol ) && rank == 0)
    {
        printf( "Num GPUs: %d\n", size );
        printf( "%dx%d: 1 GPU: %8.4f s, %d GPUs: %8.4f s, speedup: %8.2f, efficiency: %8.2f%\n", N,M, runtime_serial/ 1000.f, size, runtime/ 1000.f, runtime_serial/runtime, runtime_serial/(size*runtime)*100 );
    }
    MPI_Finalize();
    return 0;
}

#ifndef LAPLACE2D_SERIAL_H
#define LAPLACE2D_SERIAL_H

void laplace2d_serial( int rank, int iterator_max, float tol )
{
    int iterator  = 0;
    float error = 1.0f;
    #pragma acc data copy(Aref) create(newArray)
    while ( error > tol && iterator < iterator_max )
    {
        error = 0.f;

#pragma acc kernels
        for( int j = 1; j < N-1; j++)
        {
            for( int i = 1; i < M-1; i++ )
            {
                newArray[j][i] = 0.25f * ( Aref[j][i+1] + Aref[j][i-1]
                                     + Aref[j-1][i] + Aref[j+1][i]);
                error = fmaxf( error, fabsf(newArray[j][i]-Aref[j][i]));
            }
        }

#pragma acc kernels
        for( int j = 1; j < N-1; j++)
        {
            for( int i = 1; i < M-1; i++ )
            {
                Aref[j][i] = newArray[j][i];
            }
        }

        //Periodic boundary conditions
#pragma acc kernels
        for( int i = 1; i < M-1; i++ )
        {
                Aref[0][i]     = Aref[(N-2)][i];
                Aref[(N-1)][i] = Aref[1][i];
        }

        if(rank == 0 && (iterator % 100) == 0) printf("%5d, %0.6f\n", iterator, error);

        iterator++;
    }

}

int check_results( int rank, int jacobi_start, int jacobi_end, float tol )
{
    int result_correct = 1;
    for( int j = jacobi_start; j < jacobi_end && (result_correct == 1); j++)
    {
        for( int i = 1; i < M-1 && (result_correct == 1); i++ )
        {
            if ( fabs ( Aref[j][i] - A[j][i] ) >= tol )
            {
                printf("[MPI%d] ERROR: A[%d][%d] = %f does not match %f (reference)\n", rank, j,i, A[j][i], Aref[j][i]);
                result_correct = 0;
            }
        }
    }
#ifdef MPI_VERSION
    int global_result_correct = 0;
    MPI_Allreduce( &result_correct, &global_result_correct, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD );
    result_correct = global_result_correct;
#endif //MPI_VERSION
    return result_correct;
}

#endif // LAPLACE2D_SERIAL_H

