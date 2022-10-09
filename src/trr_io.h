// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#ifndef TRR_IO_H
#define TRR_IO_H

#include <string.h>
#include "gro.h"
#include "xtc_io.h"
#include "xdrfile/xdrfile_trr.h"

/*! @brief Reads a single step from an open trr file and updates the system.
 * 
 * @paragraph Note on missing information in the trr file
 * The trr file may not be complete, i.e. each frame may not contain coordinates AND velocities AND
 * forces acting on atoms. For instance, coordinates may be supplied in a specific frame but the frame may
 * be missing velocities of particles. Or, instead, velocities and forces may be supplied but not the coordinates.
 * 
 * There are two possible approaches to take when dealing with this situation.
 * Either the missing information may just not be updated in the system_t structure (e.g. if velocities
 * are missing, the atoms will just keep the velocities from the previous frame) OR the missing
 * information may be set to zero (e.g. if forces acting on atoms are missing, they will be set to 0 for all atoms).
 * 
 * For better or worse, groan library uses the second approach, i.e. the missing values are set to zero.
 * The reason for that is that we want the user to realize that something is wrong with their data.
 * If a zero velocity or force suddenly appears during the analysis, it is much more likely that the user
 * will notice the issue and solve it, either by skipping some frames in the analysis or, if they wish,
 * copy the missing information from the previous frame.
 *
 * @param trr           open XDRFILE structure corresponding to target trr file
 * @param system        pointer to a structure containing information about the system
 * 
 * @return Zero if reading was successful, else non-zero.
 * Non-zero return code indicates that the file has been fully read.
 */
int read_trr_step(XDRFILE *trr, system_t *system);


/*! @brief Writes the current positions, velocities, and forces of the selected atoms to trr file.
 * 
 * @param trr           open XDRFILE structure corresponding to target trr file
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