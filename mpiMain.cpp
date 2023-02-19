#include"Simulation.h"

#include <iostream>
#include <mpi.h>
#include <time.h>

using namespace std;

const int LEN = 40;
const int SIZE = LEN * LEN * LEN;
const int REPETITIONS = 5;
const int NUMBER_OF_PARTICLES_TO_REMOVE_ONCE = 2;
const double BOX_SIZE = 10.0;

int main(int ac, char **av) {

	MPI_Init(&ac, &av);

	int rank;

 	MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

	Simulation *sim = new Simulation();

	double *x = new double[ SIZE + REPETITIONS * NUMBER_OF_PARTICLES_TO_REMOVE_ONCE ];
	double *y = new double[ SIZE + REPETITIONS * NUMBER_OF_PARTICLES_TO_REMOVE_ONCE ];
	double *z = new double[ SIZE + REPETITIONS * NUMBER_OF_PARTICLES_TO_REMOVE_ONCE ];

	if ( ! rank ) {
	  int l;
  	  for ( int k = 0; k < LEN; k++ )
  	   for ( int j = 0; j < LEN; j++ )
  	     for ( int i = 0; i < LEN; i++ ) {
		l = i + ( LEN * k + j ) * LEN;
//		cout << " l[ " << i << ", " << j << ", " << k << " ] = " << l << endl;
		x[l] = i * BOX_SIZE;
		y[l] = j * BOX_SIZE;
		z[l] = k * BOX_SIZE;
  	     }
  	     
  	  double size = BOX_SIZE;
  	  double pos = size;
  	  for ( int i = SIZE; i < SIZE + REPETITIONS * NUMBER_OF_PARTICLES_TO_REMOVE_ONCE; i++ ) {
  	  	size *= 0.5;
  	  	pos -= size;
  	  	x[i] = 0.0;
  	  	y[i] = 0.0;
  	  	z[i] = pos;
  	  	cout << "z [ " << i << " ] = " << pos << endl;
  	  }  	  	     
	}
	if (!rank) {
	    sim->setParticles( x, y, z, SIZE + REPETITIONS * NUMBER_OF_PARTICLES_TO_REMOVE_ONCE );
	}
	for (int i = 0; i < REPETITIONS; i++) {
		sim->remove( NUMBER_OF_PARTICLES_TO_REMOVE_ONCE );

		if (!rank) {
		   cout << "Krok " << ( i + 1 ) << " -> " <<
		        sim->getAvgMinDistance() << endl;
		}

	} // REPETITIONS

	MPI_Finalize();
	return 0;
}