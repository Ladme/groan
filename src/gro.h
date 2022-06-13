// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#ifndef GRO_H
#define GRO_H

#include <stdio.h>
#include <stdlib.h>

/*
 * An array of three floats.
 */
typedef float vec_t[3];

/*
 * An array of nine floats.
 */
typedef float box_t[9];

/*
 * Specifies whether all atoms from system or 
 * only atoms in an index list should be used.
 */
typedef enum loop_mode {
    all = 0,       // all atoms
    selection = 1  // atoms from index list
} loop_mode_t;

/*
 * Structure containing all the available information about a specific atom.
 */
typedef struct atom {
    int residue_number;
    char residue_name[6];
    char atom_name[6];
    int atom_number;
    vec_t position;
    vec_t velocity;
} atom_t;

/*
 * Structure containing information about the system, or more specifically
 * about the simulation box, time-step of the simulation, and the atoms in the system.
 */
typedef struct system {
    box_t box;
    int step;
    float time;
    size_t n_atoms;
    atom_t atoms[];
} system_t;

/*
 * Structure containing an array of pointers to atoms and a number of atoms in this array.
 */
typedef struct atom_selection {
    size_t n_atoms;
    atom_t *atoms[];
} atom_selection_t;

#endif /* GRO_H */