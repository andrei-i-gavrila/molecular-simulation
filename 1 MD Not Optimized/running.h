//  running.h
//  Simulation Methods course, 2018
//  First Assignment: Molecular Dynamics (Brownian Dynamics) Simulation

#ifndef running_h
#define running_h

#include <stdio.h>

void run_simulation(void);

void calculate_external_forces_on_particles(void);

void calculate_pairwise_forces(void);

void move_particles(void);

void fold_particle_back_PBC(int i);

void distance_squared_folded_PBC(double x0, double y0, double x1, double y1,
                                 double *r2_return, double *dx_return, double *dy_return);

void write_cmovie_frame(void);

#endif /* running_h */
