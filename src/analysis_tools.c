// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#include "vector.h"
#include "analysis_tools.h"

float distance2D(const vec_t particle1, const vec_t particle2, const plane_t plane)
{   
    float dim1_d = 0;
    float dim2_d = 0;

    if (plane == xy) {
        dim1_d = particle1[0] - particle2[0];
        dim2_d = particle1[1] - particle2[1];
    } else if (plane == xz) {
        dim1_d = particle1[0] - particle2[0];
        dim2_d = particle1[2] - particle2[2];
    } else {
        dim1_d = particle1[1] - particle2[1];
        dim2_d = particle1[2] - particle2[2];
    }

    return sqrt( dim1_d*dim1_d + dim2_d*dim2_d );
}

float distance3D(const vec_t particle1, const vec_t particle2)
{
    float xd = particle1[0] - particle2[0];
    float yd = particle1[1] - particle2[1];
    float zd = particle1[2] - particle2[2];

    return sqrt( xd*xd + yd*yd + zd*zd );
}

int center_of_geometry(const atom_selection_t *selection, vec_t center)
{
    if (selection == NULL || selection->n_atoms == 0) return 1;

    for (size_t i = 0; i < selection->n_atoms; ++i) {
        vec_sum(center, selection->atoms[i]->position);
    }

    vec_div(center, selection->n_atoms);

    return 0;
}

void selection_translate(atom_selection_t *selection, box_t box, vec_t trans)
{
    for (size_t i = 0; i < selection->n_atoms; ++i) {
        atom_t *atom = selection->atoms[i];
        vec_sum(atom->position, trans);

        // check that atom is inside the box
        // if it's not, move it
        if (atom->position[0] > box[0]) atom->position[0] -= box[0];
        else if (atom->position[0] < 0) atom->position[0] += box[0];

        if (atom->position[1] > box[1]) atom->position[1] -= box[1];
        else if (atom->position[1] < 0) atom->position[1] += box[1];

        if (atom->position[2] > box[2]) atom->position[2] -= box[2];
        else if (atom->position[2] < 0) atom->position[2] += box[2];
    }
}