// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#ifndef ANALYSIS_TOOLS_H
#define ANALYSIS_TOOLS_H

#include <math.h>
#include <string.h>
#include "gro.h"

/*
 * Maximal length of the match string for the center_of_geometry function.
 */ 
#define MAX_MATCH_STRING_LEN 100

/*
 * Initial number of atom indices in selection array.
 * See function select_atoms().
 */
#define INITIAL_SELECTION_SIZE 64

/* ! \brief Splits string by delimiter and saves the substrings into an array. 
 * 
 * \param string        string to be split
 * \param array         pointer to an array into which the substring should be saved
 * \param delim         delimiter for spliting
 * 
 * \return The number of substrings resulting from spliting. A non-positive number in case of an error.
 */
size_t strsplit(char *string, char ***array, const char *delim);

/* ! \brief Compares string and the residue name of target atom.
 * 
 * \param atom          pointer to target \a atom_t structure
 * \param string        string to compare
 * 
 * \return Returns non-zero if string matches the residue name of atom or if the string is NULL. Else returns zero.
 */
int match_residue_name(const atom_t *atom, const char *string);

/* ! \brief Compares string and the atom name of target atom.
 * 
 * \param atom          pointer to target \a atom_t structure
 * \param string        string to compare
 * 
 * \return Returns non-zero if string matches the atom name of atom or if the string is NULL. Else returns zero.
 */
int match_atom_name(const atom_t *atom, const char *string);

/* ! \brief Selects atoms from the system.
 * 
 * Loops through the atoms in input_atom_ids, matching them against
 * the individual elements in the match_string.
 * Currently only support single match_function for all the match elements.
 * Note that the maximal length of the match_string is 99 characters.
 * You can change this by modifying MAX_MATCH_STRING_LEN.
 * 
 * The function also allocates memory for an array of atom indices (atom_ids) and saves
 * indices of the selected atoms into this array.
 * 
 * If input_atom_ids is NULL or n_input_atoms is zero, the function will loop through all atoms in the system.
 * 
 * \param system                system_t structure containing information about the system
 * \param input_atom_ids        an array of input atom indices
 * \param n_input_atoms         number of input atom indices
 * \param output_atom_ids       pointer to an array where the selected atom indices will be saved
 * \param match_string          string of elements separated by spaces
 * \param match_function        pointer to a function used for matching
 * 
 * \return The number of selected atoms.
 */
size_t select_atoms(
        const system_t *system,
        const size_t *input_atom_ids,
        const size_t n_input_atoms,
        size_t **output_atom_ids,
        const char *match_string,
        int (*match_function)(const atom_t *, const char *));


/* ! \brief Joins two atom selections.
 *
 * Allocates enough memory for target. 
 * Does NOT destroy or deallocate any of the input selections.
 * 
 * \param target                pointer to an array into which the concatenated selection should be saved
 * \param selection1            array of atom indices
 * \param n_selection1          number of elements in selection1
 * \param selection2            array of atom indices
 * \param n_selection2          number of elements in selection2
 * 
 * \return Length of the concatenated selection.
 */
size_t join_selections(
        size_t **target,
        const size_t *selection1,
        const size_t n_selection1,
        const size_t *selection2,
        const size_t n_selection2);

/* ! \brief Returns xy plane distance between two points in space.
 * 
 * \param particle1     pointer to an array of floats specifying the position of particle 1
 * \param particle2     pointer to an array of floats specifying the position of particle 2
 * 
 * \return xy plane distance between particle1 and particle2
 */
float distance2D(float *particle1, float *particle2);

/* ! \brief Returns distance between two points in space.
 * 
 * \param particle1     pointer to a vector specifying the position of particle1
 * \param particle2     pointer to a vector specifying the position of particle2
 * 
 * \return distance between particle1 and particle2
 */
float distance3D(vec_t particle1, vec_t particle2);

/* ! \brief Calculates center of geometry for selected atoms from specified coordinates.
 *
 * Loops through atom indices in atom_ids and calculates the center of geometry
 * from their coordinates. 
 * If coordinates are not supplied, information from system will be used.
 * System can be NULL, if you provide atom_ids and coordinates.
 * TODO: expand to cases when atom_ids is not supplied (loop through all atoms)
 * 
 * \param system                system_t structure containing information about the system
 * \param atom_ids              an array of atom indices to calculate center of geometry for 
 * \param n_selected            number of selected atoms
 * \param center                pointer to a vector into which the center of geometry will be saved
 * \param coordinates           x,y,z coordinates for all atoms in the system
 * 
 * \return Returns zero in case the calculation finished successfuly. Else returns non-zero.
 */
int center_of_geometry(
        const system_t *system,
        const size_t *atom_ids,
        const size_t n_selected,
        vec_t center,
        const vec_t *coordinates);



#endif /* ANALYSIS_TOOLS_H */