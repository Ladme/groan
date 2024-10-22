// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#ifndef ANALYSIS_TOOLS_H
#define ANALYSIS_TOOLS_H

#include <math.h>
#include <string.h>
#include "gro.h"
#include "selection.h"

#define M_PI 3.141592f
#define M_PI_X2 6.283184f

/*! @brief Python-like modulo function. */
inline float pymod(float n, float M)
{
    return fmodf( fmodf(n, M) + M, M);
}


/*! @brief Converts radians to degrees. */
inline float rad2deg(float rad)
{
    return rad * (180.0f / M_PI);
}


/*! @brief Returns oriented line distance between two points in space. Handles rectangular PBC.
 * 
 * @param particle1     pointer to an array of floats specifying the position of particle 1
 * @param particle2     pointer to an array of floats specifying the position of particle 2
 * @param dimension     specifies in which dimension the distance should be calculated
 * @param box           simulation box dimensions
 * 
 * @return Oriented line distance between particle1 and particle2
 */
float distance1D(const vec_t particle1, const vec_t particle2, const dimension_t dimension, const box_t box);


/*! @brief Returns plane distance between two points in space. Handles rectangular PBC.
 * 
 * @param particle1     pointer to an array of floats specifying the position of particle 1
 * @param particle2     pointer to an array of floats specifying the position of particle 2
 * @param plane         specifies in which plane the 2D distance should be calculated
 * @param box           simulation box dimensions
 * 
 * @return Plane distance between particle1 and particle2
 */
float distance2D(const vec_t particle1, const vec_t particle2, const plane_t plane, const box_t box);


/*! @brief Returns plane distance between two points in space. Ignores PBC and minimum image convention!
 * 
 * @param particle1     pointer to an array of floats specifying the position of particle 1
 * @param particle2     pointer to an array of floats specifying the position of particle 2
 * @param plane         specifies in which plane the 2D distance should be calculated
 * 
 * @return Plane distance between particle1 and particle2
 */
float distance2D_naive(const vec_t particle1, const vec_t particle2, const plane_t plane);


/*! @brief Returns distance between two points in space. Handles rectangular PBC.
 * 
 * @param particle1     pointer to a vector specifying the position of particle1
 * @param particle2     pointer to a vector specifying the position of particle2
 * @param box           simulation box dimensions
 * 
 * @return Distance between particle1 and particle2
 */
float distance3D(const vec_t particle1, const vec_t particle2, const box_t box);


/*! @brief Returns distance between two points in space. Ignores PBC and minimum image convention!
 * 
 * @param particle1     pointer to a vector specifying the position of particle1
 * @param particle2     pointer to a vector specifying the position of particle2
 * 
 * @return Distance between particle1 and particle2
 */
float distance3D_naive(const vec_t particle1, const vec_t particle2);

/*! @brief Calculates vector from particle1 to particle2. Handles rectangular PBC.
 * 
 * @paragraph Details
 * If the particle and its image are equidistant from the other particle,
 * the behavior of this function is undefined. In other words,
 * the function will calculate either a vector between the real particle and
 * the other particle or a vector between the particle image and the other particle.
 * This function does not guarantee that a specific vector from these two will be used.
 *
 * @param result        vector to which the result shall be saved
 * @param particle1     pointer to a vector specifying the position of particle1
 * @param particle2     pointer to a vector specifying the position of particle2
 * @param box           simulation box dimensions
 * 
 */
void calc_vector(vec_t result, const vec_t particle1, const vec_t particle2, const box_t box);

/*! @brief Calculates center of geometry for selected atoms. Handles rectangular PBC.
 *
 * 
 * @param selection             selection of atoms
 * @param center                pointer to an array for saving center of geometry
 * @param box                   current size of the simulation box
 * 
 * @return Zero, if successful; else non-zero.
 */
int center_of_geometry(const atom_selection_t *selection, vec_t center, box_t box);


/*! @brief Calculates center of geometry for selected atoms DISREGARDING PBC!
 * 
 * @param selection             selection of atoms
 * @param center                pointer to an array for saving center of geometry
 * 
 * @return Zero, if successful; else non-zero.
 */
int center_of_geometry_naive(const atom_selection_t *selection, vec_t center);


/*! @brief Calculates center of geometry for atoms selected using string query. Handles rectangular PBC.
 * 
 * @paragraph Details
 * This function creates an atom selection using the provided string query, then calculates
 * center of geometry of this selection and then deallocates the atom selection.
 * 
 * @param input_selection       atoms from which the analyzed atoms should be selected
 * @param query                 string specifying the selection of atoms
 * @param ndx_groups            dictionary containing ndx groups and their atoms
 * @param center                pointer to an array for saving center of geometry
 * @param box                   current size of the simulation box
 * 
 * @return Zero, if successful; else non-zero.
 */
int smart_center_of_geometry(
        const atom_selection_t *input_selection, 
        const char *query,
        const dict_t *ndx_groups, 
        vec_t center, 
        box_t box);


/*! @brief Translates all atoms of selection by trans. Handles rectangular PBC.
 *
 * @param selection             selection of atoms
 * @param trans                 translation vector
 * @param box                   box size
 * 
 */
void selection_translate(atom_selection_t *selection,  vec_t trans, box_t box);


/*! @brief Rotates a point around specified axis and origin counterclockwise.
 *
 * 
 * @param point                 position of the point in the xyz coordinate system
 * @param origin                position of the origin
 * @param axis                  axis to rotate the point around (x/y/z)
 * @param theta                 rotation around the axis (in degrees)
 *
 */
void rotate_point(vec_t point, const vec_t origin, const float theta, const dimension_t axis);


/*! @brief Rotates all atoms of selection around specified origin counterclockwise. Handles rectangular PBC.
 *
 * @param selection             selection of atoms to be rotated
 * @param origin                position of the origin
 * @param axis                  axis to rotate the point around (x/y/z)
 * @param theta                 rotation around the axis (in degrees)
 * @param box                   simulation box size
 *
 */
void selection_rotate(atom_selection_t *selection, const vec_t origin, const float theta, const dimension_t axis, box_t box);


/*! @brief Rotates all atoms of selection around specified origin counterclockwise. Does NOT wrap coordinates into the simulation box.
 *
 * @param selection             selection of atoms to be rotated
 * @param origin                position of the origin
 * @param axis                  axis to rotate the point around (x/y/z)
 * @param theta                 rotation around the axis
 *
 */
void selection_rotate_naive(atom_selection_t *selection, const vec_t origin, const float theta, const dimension_t axis);


/*! @brief Calculates angle between two vectors in degrees.
 *
 * @param vecA                  first vector
 * @param vecB                  second vector
 * 
 * @return Angle between the vectors in degrees.
 * 
 */
float calc_angle(const vec_t vecA, const vec_t vecB);


/*! @brief Returns distance between two points depending on the dimensionality of the problem.
 *
 * @param point1                coordinates of the first point
 * @param point2                coordinates of the second point
 * @param dim                   dimensionality of the problem
 * @param box                   simulation box size
 * @param oriented              should oriented distance be returned for 1D problems?
 * 
 * @paragraph Oriented distance
 * If oriented is set to 1 and dimensionality is x, y or z, 
 * _oriented_ distance is returned for the two points, i.e. point1.coordinate - point2.coordinate.
 * If oriented is set to 0 and dimensionality is x, y or z,
 * normal distance is returned for the two points, i.e. abs(point1.coordinate - point2.coordinate).
 * 
 * For higher dimensionality than 1, the value of oriented plays no role.
 * 
 * @paragraph Unknown dimensionality
 * In case unknown dimensionality is encountered, the behavior of this function is undefined.
 * 
 * @return Distance between the points.
 * 
 */
float calc_distance_dim(
        const vec_t point1, 
        const vec_t point2, 
        const dimensionality_t dim, 
        const box_t box, 
        const int oriented);


/*! @brief Sorts atoms in selection by their absolute distance from a reference point. Handles PBC.
 *
 * @paragraph Note on dimensionality.
 * `dim` is an enum with options dimensionality_xyz, dimensionality_xy, dimensionality_xz, dimensionality_yz,
 * dimensionality_x, dimensionality_y, dimensionality_z. 
 * 
 * This option specifies the dimensions in which the distance between the atoms and the reference point should be calculated. 
 * 
 * Even in case of 1D distance, an absolute (always positive) distance is used.
 *
 * @param selection             selection of atoms to sort
 * @param reference             reference point
 * @param dim                   what dimensions should be included in the calculation
 * @param box                   simulation box size
 * 
 */
void selection_sort_by_dist(atom_selection_t *selection, const vec_t reference, const dimensionality_t dim, box_t box);

#endif /* ANALYSIS_TOOLS_H */