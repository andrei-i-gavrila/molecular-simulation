//
// Created by andrei on 10/18/2019.
//

#include "tabulator.h"
#include "globaldata.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define EXP_BITS 3
#define MANTISSA_BITS 17
//http://www.sciencedirect.com.secure.sci-hub.tw/science/article/pii/S0010465599002179

int extract_hash(double d) {
    int *int_p;
    int index;
    int_p = (int *) &d;
    index = (*int_p & global.MASK_TABULATE);
    index = index >> 3;

    return index;
}

void initialize_table() {
    global.MASK_TABULATE = 1 << (EXP_BITS + 3);
    int N = 1 << 20;
    global.tabulate_forces = (double *) malloc(N * sizeof(double));

    double rMin = 0.1;
    double deltaR = (global.particle_particle_screening_length - rMin) / N;

    double r = 0;
    while (r < global.particle_particle_screening_length) {
        double rSq = r * r;
        int index = extract_hash(rSq);
        if (r < 0.2) {
            global.tabulate_forces[index] = 100.0 / r;
        } else {
            global.tabulate_forces[index] = 1 / rSq * exp(-r / global.particle_particle_screening_length) / r;
        }

        r += deltaR;
    }
}

double get_force(double r2) {
    return global.tabulate_forces[extract_hash(r2)];
}