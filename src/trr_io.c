// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#include "trr_io.h"

int read_trr_step(XDRFILE *trr, system_t *system)
{
    float box[3][3] = {0};

    // allocate memory for a arrays of coordinates, velocities and forces
    // these are 2D arrays of [n_atoms][3] floats
    vec_t *coordinates = calloc(system->n_atoms, sizeof(vec_t));
    vec_t *velocities  = calloc(system->n_atoms, sizeof(vec_t));
    vec_t *forces      = calloc(system->n_atoms, sizeof(vec_t));

    // read trr step using the xdrfile library
    if (read_trr(trr, system->n_atoms, &(system->step), &(system->time), &(system->lambda), box, coordinates, velocities, forces) != 0) {
        free(coordinates);
        free(velocities);
        free(forces);
        return 1;
    }

    // despite its name, this function also converts box dimensions from trr format to native groan format
    box_xtc2gro(box, system->box);

    // update the system
    for (size_t i = 0; i < system->n_atoms; ++i) {
        memcpy(system->atoms[i].position, coordinates[i], 3 * sizeof(float));
        memcpy(system->atoms[i].velocity, velocities[i],  3 * sizeof(float));
        memcpy(system->atoms[i].force,    forces[i],      3 * sizeof(float));
    }

    free(coordinates);
    free(velocities);
    free(forces);

    return 0;
}

int write_trr_step(XDRFILE *trr, const atom_selection_t *selection, int step, float time, box_t box, float lambda)
{
    float trr_box[3][3] = {{0.}};
    // convert box dimensions from groan format to native trr/xtc format
    box_gro2xtc(box, trr_box);

    // extract atom coordinates, velocities and forces from the system
    vec_t *coordinates = malloc(selection->n_atoms * sizeof(vec_t));
    vec_t *velocities  = malloc(selection->n_atoms * sizeof(vec_t));
    vec_t *forces      = malloc(selection->n_atoms * sizeof(vec_t));

    for (size_t i = 0; i < selection->n_atoms; ++i) {
        memcpy(coordinates[i], selection->atoms[i]->position, 3 * sizeof(float));
        memcpy(velocities[i],  selection->atoms[i]->velocity, 3 * sizeof(float));
        memcpy(forces[i],      selection->atoms[i]->force,    3 * sizeof(float));
    }

    // write trr step
    if (write_trr(trr, selection->n_atoms, step, time, lambda, trr_box, coordinates, velocities, forces) != 0) {
        free(coordinates);
        free(velocities);
        free(forces);
        return 1;
    }

    free(coordinates);
    free(velocities);
    free(forces);

    return 0;
}


int validate_trr(const char *filename, const int n_atoms)
{
    int trr_atoms = n_atoms;
    read_trr_natoms(filename, &trr_atoms);
    return (trr_atoms == n_atoms);
}