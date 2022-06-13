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

/* ! \brief Adds pointer to a target atom to target atom selection.
 *
 * Handles reallocation and increments the n_atoms counter in the selection.
 * 
 * \param output_atoms          selection of atoms
 * \param allocated_atoms       number of elements for which the memory has been allocated in output_atoms
 * \param atom                  pointer to an atom to be added to the selection
 *
 */
void selection_add_atom(
        atom_selection_t **output_atoms,
        size_t *allocated_atoms,
        atom_t *atom);

/* ! \brief Selects specified atoms from atom selection.
 * 
 * Loops through the atoms in input_atoms, matching them against
 * the individual elements in the match_string.
 * Currently only support single match_function for all the match elements.
 * Note that the maximal length of the match_string is 99 characters.
 * You can change this by modifying MAX_MATCH_STRING_LEN.
 * 
 * The function allocates memory for a new selection and returns a pointer to this selection.
 * 
 * If input_atom_ids is NULL or n_input_atoms is zero, the function will loop through all atoms in the system.
 * 
 * \param input_atoms           selection of atoms to choose from
 * \param match_string          string of elements separated by spaces
 * \param match_function        pointer to a function used for matching
 * 
 * \return Pointer to new atom selection.
 */
atom_selection_t *select_atoms(
        const atom_selection_t *input_atoms,
        const char *match_string,
        int (*match_function)(const atom_t *, const char *));

/* ! \brief Selects ALL atoms from system.
 *
 * Creates atom_selection structure for all atoms in the system.
 *
 * \param system                system_t structure containing information about the system
 * 
 * \return Pointer to atom selection.        
 */
atom_selection_t *select_system(system_t *system);

/* ! \brief Concatenates two atom selections.
 *
 * Allocates enough memory for target. 
 * Does NOT destroy or deallocate any of the input selections.
 * Does NOT care about duplicites in the output selection.
 * If you want to remove duplicate atoms, use selection_cat_unique().
 * 
 * \param selection1            atom selection n1
 * \param selection2            atom selection n2
 * 
 * \return Pointer to the concatenated atom selection.
 */
atom_selection_t *selection_cat(const atom_selection_t *selection1, const atom_selection_t *selection2);

/* ! \brief Concatenates two atom selections removing duplicates.
 *
 * Allocates enough memory for target. 
 * Does NOT destroy or deallocate any of the input selections.
 * Makes sure that there are no duplicate atoms in the concatenated selection.
 * 
 * Note that selection_cat_unique() is slower than selection_cat(). If you are sure that your
 * input selections do not contain duplicates, use selection_cat().
 * 
 * If you have to use selection_cat_unique(), it is a good idea (performance-wise) to use
 * the larger selection as the first argument and the smaller selection as the second argument.
 * 
 * \param selection1            atom selection n1
 * \param selection2            atom selection n2
 * 
 * \return Pointer to the concatenated atom selection.
 */
atom_selection_t *selection_cat_unique(const atom_selection_t *selection1, const atom_selection_t *selection2);

/* ! \brief Selects atoms based on specified geometric property.
 *
 * Selects atoms from input_atoms located inside a specified area and
 * add them to a new atom selection.
 * 
 * Atoms are selected relative to center. Use {0, 0, 0} as center for absolute reference.
 * 
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
 *              > geometry definition is float = radius
 * 
 * \param input_atoms           selection of atoms to choose from
 * \param center                reference coordinates
 * \param geometry              geometry type (see above)
 * \param geometry_definition   geometric description of the selection area (see above)
 * 
 * \return Pointer to new atom selection.
 * 
 */
atom_selection_t *select_geometry(
        const atom_selection_t *input_atoms,
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
 * \param input_atoms           selection of atoms
 * \param center                pointer to an array for saving center of geometry
 * 
 * \return Zero, if successful; else non-zero.
 */
int center_of_geometry(const atom_selection_t *input_atoms, vec_t center);

#endif /* ANALYSIS_TOOLS_H */