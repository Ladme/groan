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

/*
 * Used to specify geometry for a geometric selection of atoms.
 * xcylinder = cylinder with its principal axis aligned to x axis
 * ycylinder = cylinder aligned to y axis
 * zcylinder = cylinder aligned to z axis
 * rectangular = rectangular box
 * sphere    = sphere [duh]
 */
typedef enum geometry {
        xcylinder,
        ycylinder,
        zcylinder,
        box,
        sphere
} geometry_t;

/*
 * Used to specify plane.
 */
typedef enum plane {
        xy,
        xz,
        yz
} plane_t;

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

/* ! \brief Adds target atom index to a list of atom indices.
 *
 * Handles reallocation.
 * 
 * If output_atoms is NULL, does not do anything.
 * 
 * \param output_atoms          pointer to an array of indices
 * \param allocated_atoms       number of elements for which the memory has been allocated in output_atoms
 * \param array_index           index of the output_atoms array to which the atom_id should be saved
 * \param atom_id               index of the atom to be added
 *
 */
void add_selected_atom(
        size_t **output_atoms,
        size_t *allocated_atoms,
        const size_t array_index,
        const size_t atom_id);

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


/* ! \brief Selects atoms based on specified geometric property.
 *
 * Loops through atoms of system or (if supplied) atoms of atom_ids.
 * For each atom checks whether it is inside a specified area
 * and if it is, adds it to output_atom_ids.
 * 
 * Atoms are selected relative to center. Use {0, 0, 0} as center for absolute reference.
 * 
 * If coordinates are NULL, coordinates from system are used.
 * 
 * The function handles memory allocation for output_atom_ids. If you do not want to save
 * indices of selected atoms, supply NULL as **output_atom_ids.
 * 
 * Available geometries:
 *      cylinder (xcylinder,ycylinder,zcylinder)
 *              > selects atoms inside a specified cylinder
 *              > geometry definition is float[3] = {radius, bottom of the cylinder, top of the cylinder}
 *      box
 *              > selects atoms inside a box
 *              > geometry definition is float[9] = {min_x, max_x, min_y, max_y, min_z, max_z}
 *      sphere
 *              > selects atoms inside a sphere
 *              > geometry definition is *float = radius
 * 
 * \param system                system_t structure containing information about the system
 * \param atom_ids              an array of input atom indices
 * \param n_selected            number of input atom indices
 * \param output_atom_ids       pointer to an array where the selected atom indices will be saved
 * \param coordinates           coordinates of the atoms from an xtc/trr file
 * \param center                reference coordinates
 * \param geometry              geometry type (see above)
 * \param geometry_definition   geometric description of the selection area (see above)
 * 
 * \return The number of selected atoms.
 * 
 */
size_t select_geometry(
        const system_t *system,
        const size_t *atom_ids,
        const size_t n_selected,
        size_t **output_atom_ids,
        const vec_t *coordinates,
        const vec_t center,
        const geometry_t geometry,
        const void *geometry_definition);

/* ! \brief Returns plane distance between two points in space.
 * 
 * \param particle1     pointer to an array of floats specifying the position of particle 1
 * \param particle2     pointer to an array of floats specifying the position of particle 2
 * \param plane         specifies in which plane the 2D distance should be calculated
 * 
 * \return Plane distance between particle1 and particle2
 */
float distance2D(const float *particle1, const float *particle2, const plane_t plane);

/* ! \brief Returns distance between two points in space.
 * 
 * \param particle1     pointer to a vector specifying the position of particle1
 * \param particle2     pointer to a vector specifying the position of particle2
 * 
 * \return Distance between particle1 and particle2
 */
float distance3D(const vec_t particle1, const vec_t particle2);

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