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

int center_of_geometry(const atom_selection_t *input_atoms, vec_t center)
{
    if (input_atoms == NULL || input_atoms->n_atoms == 0) return 1;

    for (size_t i = 0; i < input_atoms->n_atoms; ++i) {
        vec_sum(center, input_atoms->atoms[i]->position);
    }

    vec_div(center, input_atoms->n_atoms);

    return 0;
}