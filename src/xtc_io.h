// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#ifndef XTC_IO_H
#define XTC_IO_H

#include <string.h>
#include "gro.h"
#include "xdrfile/xdrfile_xtc.h"

/* ! \brief Converts box dimensions from the xtc format into gro format.
 * 
 * Currently only supports rectangular boxes.
 * 
 * \param box           box dimensions in the xtc format
 * \param gro_box       box dimensions in the gro format
 */
void box_xtc2gro(float box[3][3], box_t gro_box);

/* ! \brief Sets velocities of all particles in a system to zero.
 *
 * \param system        pointer to a structure containing information about the system
 */
void reset_velocities(system_t *system);

/* ! \brief Reads a single trajectory step from an open xtc file and updates the system.
 * 
 * \param xtc           open XDRFILE structure corresponding to target xtc file
 * \param system        pointer to a structure containing information about the system
 * 
 * \return One if reading was successful, else non-zero.
 * Non-zero return code indicates that the file has been fully read.
 */
int read_xtc_step(XDRFILE *xtc, system_t *system);

/* ! \brief Checks that the number of atoms in xtc file matches the provided number.
 * 
 * \param filename      path to the xtc file
 * \param n_atoms       the expected number of atoms
 * 
 * \return Non-zero if the number of atoms in the xtc file matches the expected number.
 * Else returns zero. 
 */
int validate_xtc(const char *filename, const int n_atoms);

#endif /* XTC_IO_H */