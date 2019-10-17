//  running.h
//  Simulation Methods course, 2018
//  First Assignment: Molecular Dynamics (Brownian Dynamics) Simulation

#include "running.h"
#include "globaldata.h"
#include "timing.h"
#include <stdlib.h>
#include <math.h>

#define PI 3.1415926535

//runs the simulation for the required steps
void run_simulation() {

    global.total_time = 100000;
    global.echo_time = 1000;
    global.movie_time = 100;

    get_starting_time();

    for (global.time = 0; global.time < global.total_time; global.time++) {

        calculate_external_forces_on_particles();
        calculate_pairwise_forces();

        move_particles();

        //echo time
        if (global.time % global.echo_time == 0) {
            printf("Timestep: %d / %d\n", global.time, global.total_time);
            fflush(stdout);
        }

        //movie write time
        if (global.time % global.movie_time == 0)
            write_cmovie_frame();

    }

    get_finishing_time();
    echo_running_time();
}


void calculate_external_forces_on_particles() {
    int i;

    for (i = 0; i < global.N_particles; i++) {
        global.particle_fx[i] += global.particle_direction[i] * global.particle_driving_force;
    }
}

void calculate_pairwise_forces() {
    int i, j;
    double r, r2, f;
    double dx, dy;


    for (i = 0; i < global.N_particles; i++)
        for (j = i + 1; j < global.N_particles; j++) {
            //calculate the distance between the particles
            distance_squared_folded_PBC(global.particle_x[i], global.particle_y[i],
                                        global.particle_x[j], global.particle_y[j], &r2, &dx, &dy);

            //the particles are close enough to interact
            if (r2 < 16.0) {

                r = sqrt(r2);

                if (r < 0.1) {
                    printf("WARNING:PARTICLES TOO CLOSE. LOWER CUTOFF FORCE USED\n");
                    f = 100.0;
                } else {
                    //calculate the force
                    f = 1 / r2 * exp(-r / global.particle_particle_screening_length);
                }

                //projection to the x,y axes
                f = f / r;

                global.particle_fx[i] = global.particle_fx[i] - f * dx;
                global.particle_fy[i] -= f * dy;

                global.particle_fx[j] += f * dx;
                global.particle_fy[j] += f * dy;
            }

        }
}

//calculates the shortest distance squared between 2 points in a PBC configuration
//this is squared because I want to save on sqrt with the lookup table
//also used by the Verlet rebuild flag check where I check how much a particle moved
void distance_squared_folded_PBC(double x0, double y0, double x1, double y1, double *r2_return, double *dx_return,
                                 double *dy_return) {
    double dr2;
    double dx, dy;

    dx = x1 - x0;
    dy = y1 - y0;

//PBC fold back
//if any distance is larger than half the box
//the copy in the neighboring box is closer
    if (dx > global.halfSX) dx -= global.SX;
    if (dx <= -global.halfSX) dx += global.SX;
    if (dy > global.halfSY) dy -= global.SY;
    if (dy <= -global.halfSY) dy += global.SY;

    dr2 = dx * dx + dy * dy;

    *r2_return = dr2;
    *dx_return = dx;
    *dy_return = dy;
}


//moves the particles one time step
void move_particles() {
    int i;
    double dx, dy;

    for (i = 0; i < global.N_particles; i++) {
        dx = global.particle_fx[i] * global.dt;
        dy = global.particle_fy[i] * global.dt;

        global.particle_x[i] += dx;
        global.particle_y[i] += dy;

        //if (i==300) printf("%lf\n",global.particle_x[i]);

        fold_particle_back_PBC(i);

        global.particle_fx[i] = 0.0;
        global.particle_fy[i] = 0.0;
    }
}

void fold_particle_back_PBC(int i) {

//fold back the particle into the PBC simulation box
//assumes it did not jump more than a box length
//if it did the simulation is already broken anyhow

    if (global.particle_x[i] < 0) global.particle_x[i] += global.SX;
    if (global.particle_y[i] < 0) global.particle_y[i] += global.SY;
    if (global.particle_x[i] >= global.SX) global.particle_x[i] -= global.SX;
    if (global.particle_y[i] >= global.SY) global.particle_y[i] -= global.SY;
}

void write_cmovie_frame() {
    int i;
    float floatholder;
    int intholder;

    intholder = global.N_particles;
    fwrite(&intholder, sizeof(int), 1, global.moviefile);

    intholder = global.time;
    fwrite(&intholder, sizeof(int), 1, global.moviefile);

    for (i = 0; i < global.N_particles; i++) {
        intholder = global.particle_color[i];
        fwrite(&intholder, sizeof(int), 1, global.moviefile);
        intholder = i;//ID
        fwrite(&intholder, sizeof(int), 1, global.moviefile);
        floatholder = (float) global.particle_x[i];
        fwrite(&floatholder, sizeof(float), 1, global.moviefile);
        floatholder = (float) global.particle_y[i];
        fwrite(&floatholder, sizeof(float), 1, global.moviefile);
        floatholder = 1.0;//cum_disp, cmovie format
        fwrite(&floatholder, sizeof(float), 1, global.moviefile);
    }

}
