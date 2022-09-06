// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#include "analysis_tools.h"

#define M_PI 3.141592f
#define M_PI_X2 6.283184f

/* Simple function wrapping a coordinate into a simulation box. */
static inline void wrap_coordinate(float *x, const float dimension)
{
    while (*x > dimension) *x -= dimension;
    while (*x < 0) *x += dimension;
}

/* Simple function for modifying the distance based on the minimum image convention */
static inline void min_image(float *dx, const float halfbox, const float box)
{
    while (*dx > halfbox) *dx -= box;
    while (*dx < -halfbox) *dx += box;
}

float distance1D(const vec_t particle1, const vec_t particle2, const dimension_t dimension, const box_t box)
{
    float d = 0.0;
    if (dimension == x) {
        d = particle1[0] - particle2[0];
        min_image(&d, box[0] / 2, box[0]);
    } else if (dimension == y) {
        d = particle1[1] - particle2[1];
        min_image(&d, box[1] / 2, box[1]);
    } else {
        d = particle1[2] - particle2[2];
        min_image(&d, box[2] / 2, box[2]);
    }

    return d;
}

float distance2D(const vec_t particle1, const vec_t particle2, const plane_t plane, const box_t box)
{
    float dim1_d = 0.0;
    float dim2_d = 0.0;

    if (plane == xy) {
        dim1_d = particle1[0] - particle2[0];
        dim2_d = particle1[1] - particle2[1];
        min_image(&dim1_d, box[0] / 2, box[0]);
        min_image(&dim2_d, box[1] / 2, box[1]);
    } else if (plane == xz) {
        dim1_d = particle1[0] - particle2[0];
        dim2_d = particle1[2] - particle2[2];
        min_image(&dim1_d, box[0] / 2, box[0]);
        min_image(&dim2_d, box[2] / 2, box[2]);
    } else {
        dim1_d = particle1[1] - particle2[1];
        dim2_d = particle1[2] - particle2[2];
        min_image(&dim1_d, box[1] / 2, box[1]);
        min_image(&dim2_d, box[2] / 2, box[2]);
    }

    return sqrtf( dim1_d*dim1_d + dim2_d*dim2_d );
}

float distance2D_naive(const vec_t particle1, const vec_t particle2, const plane_t plane)
{   
    float dim1_d = 0.0;
    float dim2_d = 0.0;

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

    return sqrtf( dim1_d*dim1_d + dim2_d*dim2_d );
}

float distance3D(const vec_t particle1, const vec_t particle2, const box_t box)
{
    float xd = particle1[0] - particle2[0];
    float yd = particle1[1] - particle2[1];
    float zd = particle1[2] - particle2[2];

    min_image(&xd, box[0] / 2, box[0]);
    min_image(&yd, box[1] / 2, box[1]);
    min_image(&zd, box[2] / 2, box[2]);

    return sqrtf( xd*xd + yd*yd + zd*zd );
}

float distance3D_naive(const vec_t particle1, const vec_t particle2)
{
    float xd = particle1[0] - particle2[0];
    float yd = particle1[1] - particle2[1];
    float zd = particle1[2] - particle2[2];

    return sqrtf( xd*xd + yd*yd + zd*zd );
}

void calc_vector(vec_t result, const vec_t particle1, const vec_t particle2, const box_t box)
{
    register float boxx2 = box[0] / 2.;
    register float boxy2 = box[1] / 2.;
    register float boxz2 = box[2] / 2.;

    result[0] = pymod(particle2[0] - particle1[0] + boxx2, box[0]) - boxx2;
    result[1] = pymod(particle2[1] - particle1[1] + boxy2, box[1]) - boxy2;
    result[2] = pymod(particle2[2] - particle1[2] + boxz2, box[2]) - boxz2;
}

int center_of_geometry(const atom_selection_t *selection, vec_t center, box_t box)
{
    if (selection == NULL || selection->n_atoms == 0) return 1;

    // the following calculation approach is adapted from Bai, Linge; Breen, David (2008)
    // this should be able to calculate center of geometry for any distribution of atoms 
    // (except for completely homogeneous distribution)

    // reciprocal box sizes
    float rec_box_x = 1 / box[0];
    float rec_box_y = 1 / box[1];
    float rec_box_z = 1 / box[2];

    float sum_xi[3] = {0.0f};
    float sum_zeta[3] = {0.0f};

    for (size_t i = 0; i < selection->n_atoms; ++i) {
        atom_t *atom = selection->atoms[i];

        // make sure that each coordinate is inside the box
        float real_x = atom->position[0];
        float real_y = atom->position[1];
        float real_z = atom->position[2];
        wrap_coordinate(&real_x, box[0]);
        wrap_coordinate(&real_y, box[1]);
        wrap_coordinate(&real_z, box[2]);

        // then calculate magic angles
        float theta_x = real_x * rec_box_x * M_PI_X2;
        float theta_y = real_y * rec_box_y * M_PI_X2;
        float theta_z = real_z * rec_box_z * M_PI_X2;
        
        sum_xi[0]     += cosf(theta_x);
        sum_xi[1]     += cosf(theta_y);
        sum_xi[2]     += cosf(theta_z);
        sum_zeta[0]   += sinf(theta_x);
        sum_zeta[1]   += sinf(theta_y);
        sum_zeta[2]   += sinf(theta_z); 
        
    }

    // transform magic angles into real coordinates
    float final_theta_x = atan2f(-sum_zeta[0], -sum_xi[0]) + M_PI;
    float final_theta_y = atan2f(-sum_zeta[1], -sum_xi[1]) + M_PI;
    float final_theta_z = atan2f(-sum_zeta[2], -sum_xi[2]) + M_PI;

    center[0] = box[0] * (final_theta_x / M_PI_X2);
    center[1] = box[1] * (final_theta_y / M_PI_X2);
    center[2] = box[2] * (final_theta_z / M_PI_X2);

    //printf("%f %f %f\n", center[0], center[1], center[2]);

    return 0;
}

int center_of_geometry_naive(const atom_selection_t *selection, vec_t center)
{
    if (selection == NULL || selection->n_atoms == 0) return 1;

    for (size_t i = 0; i < selection->n_atoms; ++i) {
        vec_sum(center, selection->atoms[i]->position);
    }

    vec_div(center, selection->n_atoms);

    return 0;
}

void selection_translate(atom_selection_t *selection, vec_t trans, box_t box)
{
    for (size_t i = 0; i < selection->n_atoms; ++i) {
        atom_t *atom = selection->atoms[i];
        vec_sum(atom->position, trans);

        // check that atom is inside the box
        wrap_coordinate(&atom->position[0], box[0]);
        wrap_coordinate(&atom->position[1], box[1]);
        wrap_coordinate(&atom->position[2], box[2]);
    }
}