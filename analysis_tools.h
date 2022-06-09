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
 * Loops through the atoms in the system, matching them against
 * the individual elements in the match_string.
 * Currently only support single match_function for all the match elements.
 * Note that the maximal length of the match_string is 99 characters.
 * You can change this by modifying MAX_MATCH_STRING_LEN.
 * 
 * The function also allocates memory for an array of atom indices (atom_ids) and saves
 * indices of the selected atoms into this array.
 * 
 * \param system                system_t structure containing information about the system
 * \param atom_ids              pointer to an array where the selected atom indices will be saved
 * \param match_string          string of elements separated by spaces
 * \param match_function        pointer to a function used for matching
 * 
 * \return The number of selected atoms.
 */
size_t select_atoms(
        const system_t *system,
        size_t **atom_ids,
        const char *match_string,
        int (*match_function)(const atom_t *, const char *));



int center_of_geometry(
        const system_t *system,
        const size_t n_selected,
        const size_t *atom_ids,
        vec_t center,
        const vec_t *coordinates);

/* ! \brief Returns xy plane distance between two points in space.
 * 
 * \param particle1     pointer to an array of floats specifying the position of particle 1
 * \param particle2     pointer to an array of floats specifying the position of particle 2
 * 
 * \return xy plane distance between particle1 and particle2
 */
float distance2D(float *particle1, float *particle2);

#endif /* ANALYSIS_TOOLS_H */