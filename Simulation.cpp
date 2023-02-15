#include<math.h>
#include<iostream>
#include"Simulation.h"
#include"Helper.h"
#include "mpi.h"
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

    //cout << "Closest to " << i << " is " << kmin << " -> " << dSQmin << endl;
    return dSQmin;
}

double Simulation::findMinSQ( int *i, int *j ) { 
    double dSQmin = 1000000.0;
    double dSQ;
    int iMin, jMin;
    int l;
    for ( int k = rank; k < numberOfParticles; k += size ) {
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

    //cout << "Closest pair [" << iMin << ", " << jMin << "] " << dSQmin << endl;

    return dSQmin;
}

void Simulation::remove( int numberOfPairsToRemove ) {  // MPI this part
    if (first){
        shareData();
        first = false;
    }
    int i, j;
    double middleX, middleY, middleZ;
    for ( int pairs = 0; pairs < numberOfPairsToRemove; pairs ++ ) {
        findMinSQ( &i, &j ); // This will be mpied too in this version, scattering the data at the start is a bad idea, the vector gets divided and the values get skewed
        middleX = Helper::middle( x, i, j );
        middleY = Helper::middle( y, i, j );
        middleZ = Helper::middle( z, i, j );
        x[ i ] = middleX;
        y[ i ] = middleY;
        z[ i ] = middleZ;
        active[ j ] = false;
    }
   if (!rank)
    calcAvgMinDistance();
}

void Simulation::calcAvgMinDistance( void ) { 
    double sum = 0.0;
    avgMinDistance = 0;
    int effectiveParticles = 0;
    int placeholderParticles = 0;
    int l;
    double *avgMinDistanceBuf = new double[size];
    int *effectiveParticlesBuf = new int[size];
    for ( int i = 0; i < numberOfParticles; i++ ) {
        if ( active[ i ] ) {
            sum += sqrt( findMinSQ( i, &l ) );
            effectiveParticles ++;
        }
    }
   

    if (!rank){
        avgMinDistance = sum/effectiveParticles;
    }

}

double Simulation::getAvgMinDistance( void ) {
    return avgMinDistance;
}

void Simulation::shareData( void ) {
    MPI_Bcast(&numberOfParticles,1,MPI_INT,0,MPI_COMM_WORLD);
    if (rank){
        x = new double [numberOfParticles];
        y = new double [numberOfParticles];
        z = new double [numberOfParticles];
    }        
    MPI_Bcast(x,numberOfParticles,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(y,numberOfParticles,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(z,numberOfParticles,MPI_DOUBLE,0,MPI_COMM_WORLD);        

    // Now all processess should have config

        //cout << "Set number of particles to " <<numberOfParticles << " in rank " << rank << endl;
        active = new bool [ numberOfParticles ];
        for ( int i = 0; i < numberOfParticles; i++ )
            active[ i ] = true;
}
