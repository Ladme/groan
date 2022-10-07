// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#ifndef TRR_IO_H
#define TRR_IO_H

#include "gro.h"
#include "xtc_io.h"
#include "xdrfile/xdrfile_trr.h"

/*! @brief Reads a single step from an open trr file and updates the system.
 * 
 * @param xtc           open XDRFILE structure corresponding to target trr file
 * @param system        pointer to a structure containing information about the system
 * 
 * @return Zero if reading was successful, else non-zero.
 * Non-zero return code indicates that the file has been fully read.
 */
int read_trr_step(XDRFILE *xtc, system_t *system);

/*! @brief Writes the current positions, velocities, and forces of the selected atoms to trr file.
 * 
 * @param xtc           open XDRFILE structure corresponding to target trr file
 * @param selection     selection of atoms to write into trr file
 * @param step          current step of the simulation
 * @param time          time of the simulation frame
 * @param box           box size in gro format
 * @param lambda        gromacs lambda
 *  
 * @return Zero if writing has been successful, else non-zero.
 */
int write_trr_step(XDRFILE *trr, const atom_selection_t *selection, int step, float time, box_t box, float lambda);


/*! @brief Checks that the number of atoms in trr file matches the provided number.
 * 
 * @param filename      path to the trr file
 * @param n_atoms       the expected number of atoms
 * 
 * @return Non-zero if the number of atoms in the trr file matches the expected number.
 * Else returns zero. 
 */
int validate_trr(const char *filename, const int n_atoms);


#endif /* TRR_IO_H */