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
void run_simulation()
{

global.total_time = 100000;
global.echo_time = 1000;
global.movie_time = 100;

rebuild_Verlet_list();//build for the first time

get_starting_time();

for(global.time=0;global.time<global.total_time;global.time++)
    {
    
    calculate_external_forces_on_particles();
    calculate_pairwise_forces();
   
    move_particles();
    check_Verlet_rebuild_condition_and_set_flag();
    //one particle moved enough to rebuild the Verlet list
    if (global.flag_to_rebuild_Verlet==1)
        rebuild_Verlet_list();
        
    //echo time
    if (global.time % global.echo_time == 0)
        {
        printf("Timestep: %d / %d\n",global.time,global.total_time);
        fflush(stdout);
        }
    
    //movie write time
    if (global.time % global.movie_time == 0)
        write_cmovie_frame();
        
    }
    
get_finishing_time();
echo_running_time();
}



void calculate_external_forces_on_particles()
{
int i;

for(i=0;i<global.N_particles;i++)
    {
    global.particle_fx[i] += global.particle_direction[i] * global.particle_driving_force;
    }
}

void calculate_pairwise_forces()
{
int i,j,ii;
double r,r2,f;
double dx,dy;

//old not optimized version:
/*
for(i=0;i<global.N_particles;i++)
    for(j=i+1;j<global.N_particles;j++)
    {
*/

for(ii=0;ii<global.N_Verlet;ii++)
    {
    //obtain the i,j from the Verlet list
    i = global.Verletlisti[ii];
    j = global.Verletlistj[ii];
    
    //calculate the distance between the particles
    distance_squared_folded_PBC(global.particle_x[i],global.particle_y[i],
            global.particle_x[j],global.particle_y[j],&r2,&dx,&dy);
    
    //the particles are close enough to interact
    if (r2<16.0)
        {
        //r2->find f/r in the table
        //fx = f/r * dx
        //fy = f/r * dy
        r = sqrt(r2);
        if (r<0.2)
            {
            printf("WARNING:PARTICLES TOO CLOSE. LOWER CUTOFF FORCE USED\n");
            f = 100.0;
            }
        else
            {
            //calculate the force
            f = 1/r2 * exp(-r / global.particle_particle_screening_length);
            }
            
        //projection to the x,y axes
        f = f/r;
        
        global.particle_fx[i] -= f*dx;
        global.particle_fy[i] -= f*dy;
    
        global.particle_fx[j] += f*dx;
        global.particle_fy[j] += f*dy;
        }
   
    }
}

//calculates the shortest distance squared between 2 points in a PBC configuration
//this is squared because I want to save on sqrt with the lookup table
//also used by the Verlet rebuild flag check where I check how much a particle moved
void distance_squared_folded_PBC(double x0,double y0,double x1,double y1, double *r2_return, double *dx_return,double *dy_return)
{
double dr2;
double dx,dy;

dx = x1 - x0;
dy = y1 - y0;

//PBC fold back
//if any distance is larger than half the box
//the copy in the neighboring box is closer
if (dx > global.halfSX) dx -= global.SX;
if (dx <= -global.halfSX) dx += global.SX;
if (dy > global.halfSY) dy -= global.SY;
if (dy <= -global.halfSY) dy += global.SY;

dr2 = dx*dx + dy*dy;

*r2_return = dr2;
*dx_return = dx;
*dy_return = dy;
}


//moves the particles one time step
void move_particles()
{
int i;
double dx,dy;

for(i=0;i<global.N_particles;i++)
    {
    dx = global.particle_fx[i] * global.dt;
    dy = global.particle_fy[i] * global.dt;
    
    global.particle_x[i] += dx;
    global.particle_y[i] += dy;
    
    global.particle_dx_so_far[i] += dx;
    global.particle_dy_so_far[i] += dy;
  
    fold_particle_back_PBC(i);
    
    global.particle_fx[i] = 0.0;
    global.particle_fy[i] = 0.0;
    }
}

void fold_particle_back_PBC(int i)
{

//fold back the particle into the pBC simulation box
//assumes it did not jump more thana  box length
//if it did the simulation is already broken anyhow

if (global.particle_x[i]<0) global.particle_x[i] += global.SX;
if (global.particle_y[i]<0) global.particle_y[i] += global.SY;
if (global.particle_x[i]>=global.SX) global.particle_x[i] -= global.SX;
if (global.particle_y[i]>=global.SY) global.particle_y[i] -= global.SY;
}

void write_cmovie_frame()
{
int i;
float floatholder;
int intholder;

intholder = global.N_particles;
fwrite(&intholder,sizeof(int),1,global.moviefile);

intholder = global.time;
fwrite(&intholder,sizeof(int),1,global.moviefile);

for (i=0;i<global.N_particles;i++)
    {
    intholder = global.particle_color[i];
    fwrite(&intholder,sizeof(int),1,global.moviefile);
    intholder = i;//ID
    fwrite(&intholder,sizeof(int),1,global.moviefile);
    floatholder = (float)global.particle_x[i];
    fwrite(&floatholder,sizeof(float),1, global.moviefile);
    floatholder = (float)global.particle_y[i];
    fwrite(&floatholder,sizeof(float),1, global.moviefile);
    floatholder = 1.0;//cum_disp, cmovie format
    fwrite(&floatholder,sizeof(float),1,global.moviefile);
    }

}

void check_Verlet_rebuild_condition_and_set_flag()
{
int i;
double dr2;

//check if any particle moved (potentially) enough to enter the inner Verlet shell
//the worst case scenario is when the center particle moved out and the external
//particle moved inward this is why we divide by 2.0 when we set the intershell up
global.flag_to_rebuild_Verlet = 0;

for(i=0;i<global.N_particles;i++)
        {
        dr2 = global.particle_dx_so_far[i] * global.particle_dx_so_far[i] +
                global.particle_dy_so_far[i] * global.particle_dy_so_far[i];
        
        if (dr2>=global.Verlet_intershell_squared)
            {
            global.flag_to_rebuild_Verlet = 1;
            break; //exit the for cycle, the Verlet will be rebuilt anyhow
            }
        }

}

//rebuilds the Verlet list
void rebuild_Verlet_list()
{
int i,j;
double dr2,dx,dy;
double estimation;

//initialize the Verlet list for the first time
if (global.N_Verlet_max==0)
    {
    //we are building the Verlet list for the first time in the simulation
    printf("Verlet list will be built for the first time\n");
     
    estimation = global.N_particles/(double)global.SX/(double)global.SY;
    printf("System density is %.3lf\n",estimation);
    
    estimation *= PI * global.Verlet_cutoff_distance * global.Verlet_cutoff_distance;
    printf("Particles in a R = %.2lf shell = %lf\n",
        global.Verlet_cutoff_distance,estimation);
    
    global.N_Verlet_max = (int)estimation * global.N_particles / 2;
    printf("Estimated N_Verlet_max = %d\n",global.N_Verlet_max);
    
    global.Verletlisti = (int *) malloc(global.N_Verlet_max * sizeof(int));
    global.Verletlistj = (int *) malloc(global.N_Verlet_max * sizeof(int));
    }

//build the Verlet list
global.N_Verlet = 0;

for(i=0;i<global.N_particles;i++)
    {
    for(j=i+1;j<global.N_particles;j++)
        {
        distance_squared_folded_PBC(global.particle_x[i],global.particle_y[i],
            global.particle_x[j],global.particle_y[j],&dr2,&dx,&dy);
            
        //if (global.time==0)
        //    if ((i==150)||(j==150)) printf("(%d %d %lf)\n",i,j,dr2);
            
        if (dr2<36.0)
            {
            global.Verletlisti[global.N_Verlet] = i;
            global.Verletlistj[global.N_Verlet] = j;
            
            //if (global.time==0)
            //if ((i==150)||(j==150)) printf("(%d %d)\n",i,j);
            
            global.N_Verlet++;
            if (global.N_Verlet>=global.N_Verlet_max)
                {
                printf("Verlet list reallocated from %d\n",global.N_Verlet_max);
                global.N_Verlet_max = (int)(1.1*global.N_Verlet);
                global.Verletlisti = (int *) realloc(global.Verletlisti ,global.N_Verlet_max * sizeof(int));
                global.Verletlistj = (int *) realloc(global.Verletlistj ,global.N_Verlet_max * sizeof(int));
                printf("New Verlet list max size = %d\n",global.N_Verlet_max);
                }
            }
        }
    }

printf("Counted Verlet list lenght = %d\n",global.N_Verlet);
fflush(stdout);

global.flag_to_rebuild_Verlet = 0;
for(i=0;i<global.N_particles;i++)
    {
    global.particle_dx_so_far[i] = 0.0;
    global.particle_dy_so_far[i] = 0.0;
    }

//recoloring

for(i=0;i<global.N_particles;i++)
 if (global.particle_direction[i]==1) global.particle_color[i] = 3;
 else global.particle_color[i] = 2;
 
for(i=0;i<global.N_Verlet;i++)
    {
    if (global.Verletlisti[i] == 30)
        global.particle_color[global.Verletlistj[i]] = 4;
    if (global.Verletlistj[i] == 30)
        global.particle_color[global.Verletlisti[i]] = 4;
    
    }

}
