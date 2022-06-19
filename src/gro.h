// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

/* This file contains structures used by the groan library */

#ifndef GRO_H
#define GRO_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "gen_structures/dyn_array.h"

/* Integer that fits any number that can be assigned to an atom or residue in a gro file. */
typedef uint32_t groint_t;

/* An array of three floats. */
typedef float vec_t[3];

/* An array of nine floats. */
typedef float box_t[9];

/* Structure containing all the available information about a specific atom. */
typedef struct atom {
    groint_t residue_number;
    char residue_name[6];
    char atom_name[6];
    groint_t atom_number;      /* this is an atom number as taken from gro file */
    size_t gmx_atom_number;    /* this is an atom number as gromacs uses it */
    vec_t position;
    vec_t velocity;
} atom_t;

/*
 * Structure containing information about the system, or more specifically
 * about the simulation box, time-step of the simulation, and the atoms in the system.
 */
typedef struct system {
    box_t box;           /* box dimensions */
    int step;            /* simulation step; we use int only because the xdrfile library uses int */
    float time;          /* simulation time in ps */
    size_t n_atoms;      /* number of atoms in the system */
    atom_t atoms[];      /* array of atoms in the system */
} system_t;

/*
 * Structure containing an array of pointers to atoms and a number of atoms in this array.
 */
typedef struct atom_selection {
    size_t n_atoms;
    atom_t *atoms[];
} atom_selection_t;

/* Macro shortcut for atom_selection_t.
 * To avoid confusion, this shortcut should never be used anywhere in the core groan library code.
 * It might however be useful in external or peripheral programs using groan library.
 */
#define select_t atom_selection_t

/*
 * Used to specify geometry for a geometric selection of atoms.
 * xcylinder = cylinder with its principal axis aligned to x axis
 * ycylinder = cylinder aligned to y axis
 * zcylinder = cylinder aligned to z axis
 * box       = rectangular box
 * sphere    = sphere [duh]
 */
typedef enum geometry {
        xcylinder,
        ycylinder,
        zcylinder,
        box,
        sphere
} geometry_t;

/*
 * Used to specify plane.
 */
typedef enum plane {
        xy,
        xz,
        yz
} plane_t;

#endif /* GRO_H */