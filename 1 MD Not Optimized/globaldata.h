
//  globaldata.h
//  Simulation Methods course, 2018
//  First Assignment: Molecular Dynamics (Brownian Dynamics) Simulation

#ifndef globaldata_h
#define globaldata_h

#include <stdio.h>

struct global_struct {
    double SX, SY;
    double halfSX, halfSY;

    int N_particles;

    double *particle_x;                        //particle x position
    double *particle_y;                        //particle y position
    double *particle_fx;                       //particle total force x direction
    double *particle_fy;                       //particle total force y direction
    int *particle_color;                    //particle color
    double *particle_direction;                //direction of external force

    double particle_driving_force;             //external driving force
    double particle_particle_screening_length;  //inter-partcile force screening length

    double dt;          //simulation step length

    int total_time;     //total running time 100,000
    int echo_time;      //echo to screen 1,000
    int movie_time;     //write to file 100
    int time;           //current time step 1.. 100,000

    FILE *moviefile;
};

extern struct global_struct global;

#endif /* globaldata_h */
