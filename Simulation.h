/*
 * Simulation.h
 *
 *  Created on: 5 maj 2015
 *      Author: oramus
 */

#ifndef SIMULATION_H_
#define SIMULATION_H_

class Simulation {
private:
	double *x, *y, *z;
	int numberOfParticles;
	bool *active;
	bool first = true;
	int rank, size;
	double avgMinDistance;
	double findMinSQ( int *i, int *j ); // odszukuje parę cząstek, które są najbliżej
	double findMinSQ( int i, int *j ); // odszukuje cząstkę najbliższą do i

public:
	Simulation();

	void setParticles( double *x, double *y, double *z, int numberOfParticles );
	
        void shareData(void);
	void remove( int numberOfPairsToRemove );
	void calcAvgMinDistance( void );
	
	double getAvgMinDistance( void );
};

#endif /* SIMULATION_H_ */
