#include<math.h>
#include<iostream>
#include"Simulation.h"
#include"Helper.h"
include"mpi.h"
using namespace std;

Simulation::Simulation() {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
}

void Simulation::setParticles( double *x, double *y, double *z, int numberOfParticles ) {
    cout << "setParticles " << numberOfParticles << endl;
    this->x = x;
    this->y = y;
    this->z = z;
    this->numberOfParticles = numberOfParticles;
    active = new bool [ numberOfParticles ];
    for ( int i = 0; i < numberOfParticles; i++ )
        active[ i ] = true;
}

double Simulation::findMinSQ( int i, int *j ) {
    double dSQmin = 10000000.0;
    double dSQ;
    int kmin;
    for ( int k = 0; k < i; k++ ) {
        if ( ! active[ k ] ) continue;
        dSQ = Helper::getDistanceSQ(x,y,z,i,k);
        if ( dSQ < dSQmin ) {
            dSQmin = dSQ;
            kmin = k;
        }
    }
    for ( int k = i+1; k < numberOfParticles; k++ ) {
        if ( ! active[ k ] ) continue;
        dSQ = Helper::getDistanceSQ(x,y,z,i,k);
        if ( dSQ < dSQmin ) {
            dSQmin = dSQ;
            kmin = k;
        }
    }
    *j = kmin;

//    cout << "Closest to " << i << " is " << kmin << " -> " << dSQmin << endl;
    return dSQmin;
}

double Simulation::findMinSQ( int *i, int *j ) { 
    double dSQmin = 1000000.0;
    double dSQ;
    int iMin, jMin;
    int l;
    for ( int k = 0; k < numberOfParticles; k++ ) {
        if ( ! active[ k ] ) continue;
        dSQ = findMinSQ( k, &l );
        if ( dSQ < dSQmin ) {
            dSQmin = dSQ;
            iMin = k;
            jMin = l;
        }
    }
    *i = iMin;
    *j = jMin;

    cout << "Closest pair [" << iMin << ", " << jMin << "] " << dSQmin << endl;

    return dSQmin;
}

void Simulation::remove( int numberOfPairsToRemove ) {  // MPI this part
    int i, j;
    double middleX, middleY, middleZ;
    for ( int pairs = 0; pairs < numberOfPairsToRemove/size; pairs++ ) {
        findMinSQ( &i, &j );
        middleX = Helper::middle( x, i, j );
        middleY = Helper::middle( y, i, j );
        middleZ = Helper::middle( z, i, j );
        x[ i ] = middleX;
        y[ i ] = middleY;
        z[ i ] = middleZ;
        active[ j ] = false;
    }
}

void Simulation::calcAvgMinDistance( void ) { 
    double sum = 0.0;
    int effectiveParticles = 0;
    int avgMinDistancePerProc = 0;
    int l;
    double *avgMinDistanceBuf = new double [size];
    for ( int i = 0; i < numberOfParticles/size; i++ ) {
        if ( active[ i ] ) {
            sum += sqrt( findMinSQ( i, &l ) );
            effectiveParticles ++;
        }
    }
    // Return individual sums
    cout << "Suma " << sum << endl;


    avgMinDistancePerProc = sum / effectiveParticles;
    // Gather sum of sums
    MPI_Gather(&avgMinDistancePerProc,1,MPI_DOUBLE,avgMinDistanceBuf,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if (!rank)
        for (i=0 ; i<size ; i++)
            avgMinDistance += avgMinDistanceBuf[i];
}

double Simulation::getAvgMinDistance( void ) {
    return avgMinDistance;
}

void Simulation::shareData( void ) {
    MPI_Bcast(numberOfParticles,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Scatter(x,numberOfParticles/size,MPI_DOUBLE,x,numberOfParticles,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Scatter(y,numberOfParticles/size,MPI_DOUBLE,y,numberOfParticles,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Scatter(z,numberOfParticles/size,MPI_DOUBLE,z,numberOfParticles,MPI_DOUBLE,0,MPI_COMM_WORLD);        

    // Now all processess should have config
        this->x = x;
        this->y = y;
        this->z = z;
        this->numberOfParticles = numberOfParticles/size;
        active = new bool [ numberOfParticles/size ];
        for ( int i = 0; i < numberOfParticles/size; i++ )
            active[ i ] = true;
}