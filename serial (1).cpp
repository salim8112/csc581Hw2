#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include <vector>
#include <iostream>

//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    int navg,nabsavg=0;
    double davg,dmin, absmin=1.0, absavg=0.0;

    if( find_option( argc, argv, "-h" ) >= 0 )
    {   
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 1000 );

    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );
    
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
	
    for( int step = 0; step < NSTEPS; step++ )
    {
	    navg = 0;
        davg = 0.0;
	    dmin = 1.0;
	    
	    double density = .0005;
        double size = sqrt( density * n );
        
        //----------
        // Create bins as matrix of Vectors
        int num = 18;
        std::vector<particle_t*> bins[num][num];
        
        for (int i = 0; i < num; i++)
        {
            for (int j = 0; j < num; j++)
            {
                bins[i][j] = std::vector<particle_t*>();
            }
        }
        double cell_size = size / num;
        
        // Fill the matrix of Vectors
        for (int i = 0; i < n; i++)
        {
            double x1 = particles[i].x;
            int bin_x = floor(x1 / cell_size);

            double y1 = particles[i].y;
            int bin_y = floor(y1 / cell_size);

            bins[bin_y][bin_x].push_back(&particles[i]);
        }
        
        // Compute forces

        for (int row  = 0; row < num; row++)
        {
            for (int col = 0; col < num; col++)
            {
                //Add current bin and neighboring bins to close_bins
                std::vector<std::vector<particle_t*> > close_bins;

                //Add neighbors
                for (int i = row - 1; i <= row + 1; i++) {
                    for (int j = col - 1; j <= col + 1; j++) {
                        if (i >= 0 && i < num && j >= 0 && j < num) {
                            close_bins.push_back(bins[i][j]);
                        }
                    }
                }

                //Find forces for each particle in bins[row][col]
                std::vector<particle_t*> currbin = bins[row][col];
               
                for (int p = 0; p < currbin.size(); p++ ) {
                    currbin[p]->ax = currbin[p]->ay = 0;

                    //Iterate through bins
                    for (int index = 0; index < close_bins.size(); index ++){
                        std::vector<particle_t*> neighbor = close_bins[index];
                        //iterate through neighbours
                        for (int t = 0; t < neighbor.size(); t++){
                            apply_force( *currbin[p], *neighbor[t], &dmin, &davg, &navg);
                        }
                    }
                }
            }
        }
        //
        //  move particles
        //
        for( int i = 0; i < n; i++ ) {
            move( particles[i] );		
        }

        if( find_option( argc, argv, "-no" ) == -1 )
        {
          //
          // Computing statistical data
          //
          if (navg) {
            absavg +=  davg/navg;
            nabsavg++;
          }
          if (dmin < absmin) absmin = dmin;
		
          //
          //  save if necessary
          //
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }
    }
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d, simulation time = %g seconds", n, simulation_time);

    if( find_option( argc, argv, "-no" ) == -1 )
    {
      if (nabsavg) absavg /= nabsavg;
    // 
    //  -the minimum distance absmin between 2 particles during the run of the simulation
    //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
    //  -A simulation were particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
    //
    //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
    //
    printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
    if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
    if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");     

    //
    // Printing summary data
    //
    if( fsum) 
        fprintf(fsum,"%d %g\n",n,simulation_time);
 
    //
    // Clearing space
    //
    if( fsum )
        fclose( fsum );    
    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
