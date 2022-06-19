// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#ifndef ANALYSIS_TOOLS_H
#define ANALYSIS_TOOLS_H

#include <math.h>
#include <string.h>
#include "gro.h"

/* ! \brief Returns plane distance between two points in space.
 * 
 * \param particle1     pointer to an array of floats specifying the position of particle 1
 * \param particle2     pointer to an array of floats specifying the position of particle 2
 * \param plane         specifies in which plane the 2D distance should be calculated
 * 
 * \return Plane distance between particle1 and particle2
 */
float distance2D(const vec_t particle1, const vec_t particle2, const plane_t plane);

/* ! \brief Returns distance between two points in space.
 * 
 * \param particle1     pointer to a vector specifying the position of particle1
 * \param particle2     pointer to a vector specifying the position of particle2
 * 
 * \return Distance between particle1 and particle2
 */
float distance3D(const vec_t particle1, const vec_t particle2);

/* ! \brief Calculates center of geometry for selected atoms.
 *
 * \param selection             selection of atoms
 * \param center                pointer to an array for saving center of geometry
 * 
 * \return Zero, if successful; else non-zero.
 */
int center_of_geometry(const atom_selection_t *selection, vec_t center);

/* ! \brief Translates all atoms of selection by trans. Handles basic PBC.
 *
 * \param selection             selection of atoms
 * \param box                   box size
 * \param trans                 translation vector
 * 
 */
void selection_translate(atom_selection_t *selection, box_t box, vec_t trans);

#endif /* ANALYSIS_TOOLS_H */