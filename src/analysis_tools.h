// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#ifndef ANALYSIS_TOOLS_H
#define ANALYSIS_TOOLS_H

#include <math.h>
#include <string.h>
#include "gro.h"

/*! \brief Returns oriented line distance between two points in space. Handles rectangular PBC.
 * 
 * \param particle1     pointer to an array of floats specifying the position of particle 1
 * \param particle2     pointer to an array of floats specifying the position of particle 2
 * \param dimension     specifies in which dimension the distance should be calculated
 * \param box           simulation box dimensions
 * 
 * \return Oriented line distance between particle1 and particle2
 */
float distance1D(const vec_t particle1, const vec_t particle2, const dimension_t dimension, const box_t box);


/*! \brief Returns plane distance between two points in space. Handles rectangular PBC.
 * 
 * \param particle1     pointer to an array of floats specifying the position of particle 1
 * \param particle2     pointer to an array of floats specifying the position of particle 2
 * \param plane         specifies in which plane the 2D distance should be calculated
 * \param box           simulation box dimensions
 * 
 * \return Plane distance between particle1 and particle2
 */
float distance2D(const vec_t particle1, const vec_t particle2, const plane_t plane, const box_t box);


/*! \brief Returns plane distance between two points in space. Ignores PBC and minimum image convention!
 * 
 * \param particle1     pointer to an array of floats specifying the position of particle 1
 * \param particle2     pointer to an array of floats specifying the position of particle 2
 * \param plane         specifies in which plane the 2D distance should be calculated
 * 
 * \return Plane distance between particle1 and particle2
 */
float distance2D_naive(const vec_t particle1, const vec_t particle2, const plane_t plane);


/*! \brief Returns distance between two points in space. Handles rectangular PBC.
 * 
 * \param particle1     pointer to a vector specifying the position of particle1
 * \param particle2     pointer to a vector specifying the position of particle2
 * \param box           simulation box dimensions
 * 
 * \return Distance between particle1 and particle2
 */
float distance3D(const vec_t particle1, const vec_t particle2, const box_t box);


/*! \brief Returns distance between two points in space. Ignores PBC and minimum image convention!
 * 
 * \param particle1     pointer to a vector specifying the position of particle1
 * \param particle2     pointer to a vector specifying the position of particle2
 * 
 * \return Distance between particle1 and particle2
 */
float distance3D_naive(const vec_t particle1, const vec_t particle2);


/*! \brief Calculates center of geometry for selected atoms. Handles rectangular PBC.
 *
 * \param selection             selection of atoms
 * \param center                pointer to an array for saving center of geometry
 * \param box                   current size of the simulation box
 * 
 * \return Zero, if successful; else non-zero.
 */
int center_of_geometry(const atom_selection_t *selection, vec_t center, box_t box);


/*! \brief Calculates center of geometry for selected atoms DISREGARDING PBC!
 * 
 * \param selection             selection of atoms
 * \param center                pointer to an array for saving center of geometry
 * 
 * \return Zero, if successful; else non-zero.
 */
int center_of_geometry_naive(const atom_selection_t *selection, vec_t center);


/*! \brief Translates all atoms of selection by trans. Handles rectangular PBC.
 *
 * \param selection             selection of atoms
 * \param trans                 translation vector
 * \param box                   box size
 * 
 */
void selection_translate(atom_selection_t *selection,  vec_t trans, box_t box);

#endif /* ANALYSIS_TOOLS_H */