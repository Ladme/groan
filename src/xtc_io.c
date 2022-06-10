// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#include "xtc_io.h"

void box_xtc2gro(const float box[3][3], box_t gro_box)
{
    gro_box[0] = box[0][0];
    gro_box[1] = box[1][1];
    gro_box[2] = box[2][2];
}

void reset_velocities(system_t *system)
{
    for (size_t i = 0; i < system->n_atoms; ++i) {
        memset(system->atoms[i].velocity, 0, 3 * sizeof(float));
    }
}

int read_xtc_step(XDRFILE *xtc, system_t *system)
{
    int step = 0;
    float time = 0.;
    float box[3][3] = {0};
    float prec = 0.;

    // allocate memory for an array of coordinates
    // this is a 2D array of [n_atoms][3] floats
    vec_t *coordinates = malloc(system->n_atoms * sizeof(vec_t));

    if (read_xtc(xtc, system->n_atoms, &step, &time, box, coordinates, &prec) != 0) {
        free(coordinates);
        return 1;
    }

    // update the system
    system->step = step;
    system->time = time;
    box_xtc2gro(box, system->box);
    // loop through all the atoms and update their coordinates
    for (size_t i = 0; i < system->n_atoms; ++i) {
        memcpy(system->atoms[i].position, coordinates[i], 3 * sizeof(float));
    }

    /*DEBUG
    for (size_t i = 0; i < system->n_atoms; ++i) {
        atom_t *atom = &(system->atoms[i]);
        printf("%f %f %f\n", atom->position[0], atom->position[1], atom->position[2]);
    }*/

    free(coordinates);

    return 0;
}


int validate_xtc(const char *filename, const int n_atoms)
{
    int xtc_atoms = n_atoms;
    read_xtc_natoms(filename, &xtc_atoms);
    return (xtc_atoms == n_atoms);
}


