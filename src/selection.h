// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#ifndef SELECTION_H
#define SELECTION_H

#include <ctype.h>
#include <string.h>
#include "gro.h"

/*! @brief Splits string by delimiter and saves the substrings into an array. 
 * 
 * @paragraph Details
 * Modifies the original string! Note that the items in array are mere pointers
 * to the segments of the original line. Thus array can be deallocated by simple free.
 * 
 * @param string        string to be split
 * @param array         pointer to an array into which the substring should be saved
 * @param delim         delimiter for spliting
 * 
 * @return The number of substrings resulting from spliting. A non-positive number in case of an error.
 */
int strsplit(char *string, char ***array, const char *delim);


/*! @brief Strips leading and trailing whitespaces from string.
 * 
 * @paragraph NULL or no 0-terminator
 * The behavior of this function is UNDEFINED if string is NULL or if it is not terminated by \0.
 * 
 * @paragraph Acknowledgments
 * This function is a slightly modified version of jkramer's trim() from 
 * https://stackoverflow.com/questions/122616/how-do-i-trim-leading-trailing-whitespace-in-a-standard-way
 * 
 * @param string        string to strip
 * 
 */
void strstrip(char *string);

/*! @brief Removes all white spaces from a string.
 * 
 * @paragraph Acknowledgments
 * This function is a slightly modified version of a function from
 * https://stackoverflow.com/questions/1726302/remove-spaces-from-a-string-in-c
 * 
 * @param string        string to modify
 * 
 */
void strremwhite(char *string);

/*! @brief Compares string and the residue name of target atom.
 * 
 * @paragraph NULL atom pointer
 * The behavior of this function is UNDEFINED if the atom_t pointer points to NULL.
 * 
 * @param atom          pointer to target atom_t structure
 * @param string        string to compare
 * 
 * @return Returns non-zero if string matches the residue name of atom or if the string is NULL. Else returns zero.
 */
int match_residue_name(const atom_t *atom, const char *string);


/*! @brief Compares string and the residue number of target atom.
 * 
 * @paragraph NULL atom pointer
 * The behavior of this function is UNDEFINED if the atom_t pointer points to NULL.
 * 
 * @param atom          pointer to target atom_t structure
 * @param string        string to compare
 * 
 * @return Returns non-zero if string matches the residue number of atom or if the string is NULL. Else returns zero.
 */
int match_residue_num(const atom_t *atom, const char *string);


/*! @brief Compares string and the atom name of target atom.
 * 
 * @paragraph NULL atom pointer
 * The behavior of this function is UNDEFINED if the atom_t pointer points to NULL.
 *
 * @param atom          pointer to target atom_t structure
 * @param string        string to compare
 * 
 * @return Returns non-zero if string matches the atom name of atom or if the string is NULL. Else returns zero.
 */
int match_atom_name(const atom_t *atom, const char *string);


/*! @brief Compares the provided string containing a number with GMX_ATOM_NUMER. 
 *
 * @paragraph NULL atom pointer
 * The behavior of this function is UNDEFINED if the atom_t pointer points to NULL.
 * 
 * @param atom          pointer to target atom_t structure
 * @param string        string containing an integer
 *
 * @return Returns non-zero if the number in the string matches the gmx_atom_number or if the string is NULL. Else returns zero.
 */
int match_atom_num(const atom_t *atom, const char *string);


/*! @brief Creates an empty atom_selection structure and allocates memory for it.
 * 
 * @paragraph Details
 * Empty atom_selection means that the n_atoms variable is set to zero.
 * Enough memory to hold 'items' atoms is allocated, though the memory can be later reallocated,
 * e.g. by using selection_add_atom().
 *
 * @param items         number of atoms for which memory should be allocated
 * 
 * @return Pointer to empty atom selection.
 */
atom_selection_t *selection_create(size_t items);


/*! @brief Copies an atom selection to a new selection.
 * 
 * @paragraph Details
 * This indeed returns DEEP copy of the provided selection.
 * (The atom pointers inside the selection are still pointers, of course.)
 * 
 * @param selection     pointer to a selection to be copied
 * 
 * @return Pointer to new copy of the selection.
 */
atom_selection_t *selection_copy(const atom_selection_t *selection);


/*! @brief Same as selection_copy() but the input selection is deallocated. */
atom_selection_t *selection_copy_d(atom_selection_t *selection);


/*! @brief Empties an existing atom selection.
 *
 * @paragraph Details
 * Memory is not deallocated and is left as it was.
 * 
 * Actually, the only thing this function does is setting the n_atoms variable to 0.
 * The atom pointers are still present in the structure, but should not be
 * available to any proper groan function.
 * 
 * @param selection     pointer to a selection to be emptied
 *
 */
void selection_empty(atom_selection_t *selection);


/*! @brief Adds pointer to a target atom to target atom selection.
 *
 * @paragraph Details
 * Handles reallocation and increments the n_atoms counter in the selection.
 * 
 * @param output_atoms          selection of atoms
 * @param allocated_atoms       number of elements for which the memory has been allocated in output_atoms
 * @param atom                  pointer to an atom to be added to the selection
 *
 */
void selection_add_atom(
        atom_selection_t **output_atoms,
        size_t *allocated_atoms,
        atom_t *atom);

/*! @brief Adds atoms from target selection to another selection.
 *
 * @paragraph Details
 * Calls selection_add_atom() on every atom of the 'atoms_to_add' selection. 
 * Handles reallocation. Does not remove duplicate atoms.
 * 
 * @param output_atoms          selection of atoms
 * @param allocated_atoms       number of elements for which the memory has been allocated in output_atoms
 * @param atoms_to_add          selection of atoms that should be added to the output_atoms
 *
 */
void selection_add(
        atom_selection_t **output_atoms, 
        size_t *allocated_atoms, 
        atom_selection_t *atoms_to_add);

/*! @brief Selects specified atoms from atom selection.
 * 
 * @paragraph Details
 * Loops through the atoms in input_atoms, matching them against
 * the individual elements in the match_string.
 * Only supports single match_function for all the match elements.
 * Use smart_select() for more advanced queries.
 * 
 * The function allocates memory for a new selection and returns a pointer to this selection.
 * 
 * @param input_atoms           selection of atoms to choose from
 * @param match_string          string of elements separated by spaces
 * @param match_function        pointer to a function used for matching
 * 
 * @return Pointer to new atom selection.
 */
atom_selection_t *select_atoms(
        const atom_selection_t *input_atoms,
        const char *match_string,
        int (*match_function)(const atom_t *, const char *));


/*! @brief Same as select_atoms() but the input selection is deallocated. */
atom_selection_t *select_atoms_d(
        atom_selection_t *input_atoms,
        const char *match_string,
        int (*match_function)(const atom_t *, const char *));


/*! @brief Selects ALL atoms from system.
 *
 * Creates atom_selection structure for all atoms in the system.
 *
 * @param system                system_t structure containing information about the system
 * 
 * @return Pointer to atom selection.        
 */
atom_selection_t *select_system(system_t *system);


/*! @brief Concatenates two atom selections.
 *
 * @paragraph Details
 * Allocates enough memory for the output selection.
 * Does NOT destroy or deallocate any of the input selections.
 * Does NOT care about duplicates in the output selection.
 * If you want to remove duplicate atoms, use selection_cat_unique().
 * 
 * @param selection1            atom selection n1
 * @param selection2            atom selection n2
 * 
 * @return Pointer to the concatenated atom selection.
 */
atom_selection_t *selection_cat(const atom_selection_t *selection1, const atom_selection_t *selection2);


/* ! @brief Same as selection_cat() but the input selections are both deallocated. */
atom_selection_t *selection_cat_d(atom_selection_t *selection1, atom_selection_t *selection2);


/*! @brief Concatenates two atom selections while removing duplicate atoms.
 *
 * @paragraph Details
 * Allocates enough memory for the output selection.
 * Does NOT destroy or deallocate any of the input selections.
 * Makes sure that there are no duplicate atoms in the concatenated selection.
 * 
 * @paragraph Note on speed
 * Note that selection_cat_unique() is slower than selection_cat(). If you are sure that your
 * input selections do not contain duplicates, use selection_cat().
 * 
 * If you have to use selection_cat_unique(), it is a good idea (performance-wise) to use
 * the larger selection as the first argument and the smaller selection as the second argument.
 * 
 * You can achieve the same behavior as selection_cat_unique() with subsequent calls of
 * selection_cat() and selection_unique(). Selection_cat_unique() is however an order of magnitude faster.
 * 
 * @param selection1            atom selection n1
 * @param selection2            atom selection n2
 * 
 * @return Pointer to the concatenated atom selection.
 */
atom_selection_t *selection_cat_unique(const atom_selection_t *selection1, const atom_selection_t *selection2);


/*! @brief Same as selection_cat_unique() but the input selections are both deallocated. */
atom_selection_t *selection_cat_unique_d(atom_selection_t *selection1, atom_selection_t *selection2);


/*! @brief Creates an atom selection which is an intersection of the provided selections.
 *
 * @paragraph Details
 * Allocates enough memory for the output selection.
 * Does NOT destroy or deallocate any of the input selections.
 *
 * @param selection1            atom selection n1
 * @param selection2            atom selection n2
 *
 * @return Pointer to the intersected atom selection.
 */
atom_selection_t *selection_intersect(const atom_selection_t *selection1, const atom_selection_t *selection2);


/*! @brief Same as selection_intersect() but the input selections are both deallocated. */
atom_selection_t *selection_intersect_d(atom_selection_t *selection1, atom_selection_t *selection2);


/*! @brief Removes an atom from atom selection.
 *
 * @paragraph Details
 * If there are multiple identical atoms (meaning if there are multiple pointers to the same atom in system),
 * ALL of them will be removed.
 * 
 * The amount of memory allocated to selection is not changed by this function.
 * 
 * @param selection             atom selection to modify
 * @param remove                pointer to an atom to remove
 * 
 * @return Number of deleted atoms.
 */ 
size_t selection_remove_atom(atom_selection_t *selection, atom_t *remove);


/*! @brief Removes atoms which are part of selection_sub from selection_result. Legacy version: DO NOT USE.
 *
 * @paragraph Details
 * The amount of memory allocated to any selection is not changed by this function.
 * This is a pretty dumb and inefficient function as it simply consists of calling...
 * ...selection_intersect() on the two intersections and then looping through the atoms of the intersection...
 * ...removing them from the result using selection_remove_atom().
 * 
 * @param selection_result             atom selection to modify
 * @param selection_sub                atom selection to be subtracted
 * 
 * @return Number of deleted atoms.
 */ 
size_t selection_remove_legacy(atom_selection_t *selection_result, const atom_selection_t *selection_sub);


/*! @brief Removes atoms which are part of selection_sub from selection_result.
 *
 * @paragraph Details
 * The amount of memory allocated to any selection is not changed by this function.
 * 
 * @param selection_result             atom selection to modify
 * @param selection_sub                atom selection to be subtracted
 * 
 * @return Number of deleted atoms.
 */ 
size_t selection_remove(atom_selection_t *selection_result, const atom_selection_t *selection_sub);


/*! @brief Same as selection_remove() but the selection_sub is deallocated. */
size_t selection_remove_d(atom_selection_t *selection_result, atom_selection_t *selection_sub);


/*! @brief Removes all duplicit entries in the selection.
 *
 * @paragraph Details
 * Keeps the other atoms in the same order. Does not deallocate memory.
 * 
 * @param selection                     atom selection to modify
 * 
 * @return Number of deleted atoms.
 */
size_t selection_unique(atom_selection_t *selection);


/*! @brief Compares two atom selections.
 *
 * @paragraph Details
 * The function proclaims two selections as identical
 * if they contain the same atoms. 
 * The order of the atoms is not relevant.
 * 
 * The behavior of this function is UNDEFINED for non-unique atom selections.
 * 
 * @param selection1            atom selection n1
 * @param selection2            atom_selection n2
 * 
 * @return One, if the atom selections are identical. Else zero.
 */
int selection_compare(
        const atom_selection_t *selection1, 
        const atom_selection_t *selection2);


/*! @brief Compares two atom selections STRICTLY.
 *
 * @paragraph Details
 * The function proclaims two selections as identical
 * if they contain the same atoms IN THE SAME ORDER.
 * 
 * Non-unique selections can also be compared using this function.
 * 
 * @param selection1            atom selection n1
 * @param selection2            atom_selection n2
 * 
 * @return One, if the atom selections are strictly identical. Else zero.
 */
int selection_compare_strict(
        const atom_selection_t *selection1, 
        const atom_selection_t *selection2);


/*! @brief Renumbers atoms of a selection.
 *
 * @paragraph Details
 * Atoms and the corresponsing residues in a selection are renumbered based
 * on their position in the selection.
 * Be careful when using this function. It might be better to use selection_to_system().
 * 
 * Notes:
 * 1) The renumbering propagates to the system_t structure and
 * consequently to all the other selections associated with the system!
 * 
 * 2) The function maintains correct residue numbering even if the atoms
 * of the same residue are not located next to each other.
 * Gromacs will however have problems with such split residues.
 * Therefore, it might be useful to fix the selection beforehand by using selection_fixres().
 * 
 * 3) The behavior of this function is UNDEFINED when applied to 
 * selection containing duplicate atoms.
 * 
 * 4) For the above reasons, using this function is not recommended. 
 * If you want to properly renumber atoms in a selection, use selection_to_system().
 * 
 * @param selection             selection of atoms to renumber
 *  
 */ 
void selection_renumber(atom_selection_t *selection);


/*! @brief Sorts the atoms in a selection based on the atom number.
 *
 * @paragraph Limitations
 * The behavior of this function is UNDEFINED when applied to 
 * selection containing duplicate atoms.
 * 
 * @param selection             selection of atoms to sort
 * 
 */ 
void selection_sort(atom_selection_t *selection);


/*! @brief Sorts the atoms in a selection based on the gmx atom number.
 *
 * @paragraph Limitations
 * The behavior of this function is UNDEFINED when applied to 
 * selection containing duplicate atoms.
 * 
 * @param selection             selection of atoms to sort
 * 
 */ 
void selection_sort_gmx(atom_selection_t *selection);


/*! @brief Reverses the order of atoms in the atom selection.
 *
 * @param selection             selection of atoms to reverse
 * 
 */ 
void selection_reverse(atom_selection_t *selection);

/*! @brief Takes a slice of atoms from selection.
 *
 * @paragraph Indices
 * Slice start is the index of the first atom in the selection to be
 * included in the slice.
 * Slice end is the index right after the last atom in the selection
 * that is included in the slice.
 * For instance, if the selection contains atoms [A, B, C, D, E, F, G, H],
 * and slice start is 2 and slice end is 5, the output slice will be
 * [C, D, E].
 * 
 * selection_slice() does not care about atom numbers or residues,
 * it only cares about the positions of atoms in the selection itself.
 * 
 * In other words, this replicates the pythonic `list[a:b]`.
 * 
 * @paragraph Zero indices
 * If slice_start is zero, the slicing will start at the start of the selection.
 * If slice_end is zero, the slicing will end at the END of the selection 
 * (the last item will be included in the slice).
 * 
 * @paragraph Negative indices
 * Using negative index for slice_start or slice_end leads to the same
 * behavior as in python (-1 is the last index, -2 is the second to last index).
 * 
 * @paragraph Input selection
 * Input selection will not be modified.
 * 
 * @param selection            selection to slice
 * @param slice_start          index of the first item
 * @param slice_end            index of the next item after the last item 
 * 
 * @return Atom selection containing slice of the input selection. NULL in case of an error.
 */
atom_selection_t *selection_slice(const atom_selection_t *selection, int slice_start, int slice_end);

/*! @brief Fixes split residues in a selection.
 *
 * @paragraph Details
 * Changes the order of atoms in a selection in such a way
 * that all atoms of every residue are located right after each other.
 * Also sorts the atoms of all residues according to their gmx atom number.
 * (Yes, "gmx atom number", not "atom number".)
 * 
 * Does not change the order of residues themselves.
 * 
 * @paragraph Limitations
 * The behavior of this function is UNDEFINED when applied to
 * selection containing duplicate atoms or when it is applied
 * to selection with more than 99,999 residues.
 * 
 * @param selection             selection of atoms to fixres
 * 
 */ 
void selection_fixres(atom_selection_t *selection);


/*! @brief Checks whether the atom is a part of selection.
 * 
 * @param selection             selection of atoms to search in
 * @param atom                  pointer to atom to search for
 * 
 * @return One, if atom IS part of the selection. Else zero.
 */ 
int selection_isin(atom_selection_t *selection, atom_t *atom);


/*! @brief Calculates the number of unique residues in selection.
 * 
 * @paragraph Limitations
 * Be careful about the total number of residues in the system.
 * This function uses residue_number to indentify residues and
 * the residue_number is capped at 99,999.
 * 
 * This function will not behave correctly if you have more than
 * 99,999 residues in your system.
 *
 * @param selection             selection of atoms to search in
 * 
 * @return The number of unique residues in selection.
 */ 
size_t selection_getnres(const atom_selection_t *selection);


/*! @brief Collects unique residue names from a selection.
 *
 * @param selection             selection of atoms to search in
 * 
 * @return Pointer to a list_t structure.
 */ 
list_t *selection_getresnames(const atom_selection_t *selection);


/*! @brief Splits selection into an array of atom selections of the same residue.
 *
 * @paragraph Memory
 * The function allocates memory for the split array as well as for the newly constructed
 * atom selections in the array. To deallocate memory for this array, you have to
 * loop through the array, deallocating memory for each atom_selection_t structure using free().
 * Then you can deallocate the split array using free().
 * 
 * @paragraph Note on ordering
 * The selections in the split array are NOT sorted by the residue number.
 * Instead, they are sorted based on the position of the first atom of each residue in
 * the original selection.
 * 
 * @paragraph Limitations
 * The behavior of this function is undefined if the number of residues is higher than 99,999.
 *
 * @param selection             selection of atoms to split
 * @param split                 pointer to an array containing resulting selections
 * 
 * @return Number of unique residues in the selection.
 */ 
size_t selection_splitbyres(const atom_selection_t *selection, atom_selection_t ***split);


/*! @brief Creates a new system_t structure from provided atom selection.
 * 
 * @paragraph Details
 * Atoms pointed at by selection are DEEP copied to new system_t structure.
 * 
 * Duplicate atoms are removed from the new system.
 * Note that duplicate atoms are NOT atoms with the same atom number, 
 * residue number, position or whatever. Duplicate atoms are atoms
 * located in the same part of memory. (It would be more correct to call
 * these "duplicate atom pointers".) 
 * You may still have atoms with the same atom number, residue number, position, velocity etc. 
 * and they will be treated as separate atoms IF they are saved in different parts
 * of memory. Such atoms will NOT be removed as duplicates and they will be
 * renumbered based on their position in the selection.
 * 
 * Any split residues are fixed by sorting the residue atoms by gmx atom numbers.
 * In case of other residues, the function maintains the order of the atoms in the selection.
 * 
 * New atoms are assigned new gmx atom numbers based on their position in the new system.
 * 
 * The original selection is kept intact! The original selection thus
 * still points to the atoms of the original system. The original system
 * is also not changed in any way.
 * 
 * If you want to select atoms of the new system, use select_system().
 * 
 * @paragraph Limitations
 * The behavior of this function is UNDEFINED when applied to
 * selection with more than 99,999 residues.
 * 
 * @param selection             selection of atoms
 * @param box                   boxsize of the new system
 * @param step                  current simulation step
 * @param time                  current simulation time
 * 
 * @return Pointer to the new system_t structure.
 */
system_t *selection_to_system(
        const atom_selection_t *selection, 
        const box_t box, 
        const int step, 
        const float time);


/*! @brief Same as selection_to_system() but the input selection is deallocated. */
system_t *selection_to_system_d(
        atom_selection_t *selection, 
        const box_t box, 
        const int step, 
        const float time);

/*! @brief Creates a new system_t structure from provided atom selection without fixing residues.
 * 
 * @param selection             selection of atoms
 * @param box                   boxsize of the new system
 * @param step                  current simulation step
 * @param time                  current simulation time
 * 
 * @return Pointer to the new system_t structure.
 */
system_t *selection_to_system_nofixres(
        const atom_selection_t *selection, 
        const box_t box, 
        const int step, 
        const float time);

/*! @brief Selects atoms based on specified geometric property. Handles rectangular PBC.
 *
 * @paragraph Details
 * Selects atoms from input_atoms located inside a specified area and
 * add them to a new atom selection.
 * 
 * Atoms are selected relative to center. Use {0, 0, 0} as center for absolute reference.
 * 
 * @paragraph Supported geometries
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
 * @param input_atoms           selection of atoms to choose from
 * @param center                reference coordinates
 * @param geometry              geometry type (see above)
 * @param geometry_definition   geometric description of the selection area (see above)
 * @param system_box            simulation box dimensions
 * 
 * @return Pointer to new atom selection.
 * 
 */
atom_selection_t *select_geometry(
        const atom_selection_t *input_atoms,
        const vec_t center,
        const geometry_t geometry,
        const void *geometry_definition,
        const box_t system_box);


/*! @brief Same as select_geometry() but the input selection is deallocated. */
atom_selection_t *select_geometry_d(
        atom_selection_t *input_atoms,
        const vec_t center,
        const geometry_t geometry,
        const void *geometry_definition,
        const box_t system_box);


/*! @brief Selects atoms based on provided string query. 
 *
 * @paragraph Groan selection language
 * See https://github.com/Ladme/groan#groan-selection-language for more information about what queries are supported by smart_select().
 * 
 * @paragraph Note on memory deallocation
 * The selection returned by this function must ALWAYS be freed when it is not needed anymore.
 * This holds true even for selections originating from the ndx_groups dictionary (which you usually do not free separately).
 * 
 * @paragraph 'ndx_groups' can be NULL
 * If 'ndx_groups' is NULL, the ndx groups will not be used for selecting atoms but the parsing can still be successful
 * (if the ndx groups are not needed to understand the query).
 * 
 * @param selection             selection of atoms to choose from
 * @param query                 query to be parsed
 * @param ndx_groups            dictionary containing definitions of the ndx groups (use read_ndx() to obtain it)
 * 
 * @return Pointer to the atom selection. NULL in case the parsing fails. If query is NULL, returns copy of the input selection.
 * 
 */
atom_selection_t *smart_select(const atom_selection_t *selection, const char *query, const dict_t *ndx_groups);


/*! @brief Select atoms based on the provided geometry query.
 *
 * @paragraph Groan selection language
 * Use groan selection language for 'selection_query' and 'reference_query' (see https://github.com/Ladme/groan#groan-selection-language).
 * 
 * @paragraph Geometry query
 * The geometry query has the following format: GEOMETRY_TYPE OPTIONS...
 * The number and character of options depends on the GEOMETRY_TYPE. 
 * > xcylinder/ycylinder/zcylinder 
 *      >>> options: RADIUS MIN-MAX
 *      >>> example: zcylinder 1.5 -3-3
 * > sphere
 *      >>> options: RADIUS
 *      >>> example: sphere 3
 * > box
 *      >>> options: MINX-MAXX MINY-MAXY MINZ-MAXZ
 *      >>> example: box -1-1 -2--1 0-4
 * 
 * @paragraph Specifying center
 * The cylinder, sphere, or box used for geometry selection is placed to the center of geometry of reference atoms specified by 'reference_query'.
 * In case the 'reference_query' is NULL, a reference point in box origin {0, 0, 0} is used.
 * The 'reference_query' can also start with 'point' in which case the supplied reference point will be used.
 * For example, using "point 3.5 4.2 1.1" will place the reference point to x = 3.5 nm, y = 4.2 nm, z = 1.1 nm.
 * 
 * @paragraph NULL parameters
 * In case the 'input selection' is NULL, the function returns NULL.
 * In case the 'selection_query' is NULL, the entire input_selection is used for geometry selection.
 * In case the 'reference_query' is NULL, the box origin is used as a reference point for the geometry selection.
 * In case the 'geometry query' is NULL, the entire selection specified by 'selection_query' is returned.
 * In case the 'ndx_groups' is NULL, ndx groups will not be used for selecting atoms (see smart_select()).
 * In case the 'system_box' is NULL, the function returns NULL.
 * 
 * @param input_selection       selection of atoms to be used by the function
 * @param selection_query       query to specify atoms to use for geometry selection
 * @param reference_query       query to specify reference atoms
 * @param geometry_query        query to specify geometry to use for sleection
 * @param ndx_groups            dictionary containing definitions of the ndx groups (use read_ndx() to obtain it)
 * @param system_box            dimensions of the simulation box (get as system->box)
 * 
 * @return Pointer to atom selection specified by selection_query and geometry_query. NULL in case the parsing fails.
 * 
 */
atom_selection_t *smart_geometry(
        const atom_selection_t *input_selection,
        const char *selection_query, 
        const char *reference_query,
        const char *geometry_query,
        const dict_t *ndx_groups,
        box_t system_box);


/*! @brief Reads an ndx file creating atom selection for each index group.
 *
 * @paragraph Structure of the output dictionary
 * Index groups and their atoms are saved into dictionary, pointer to
 * which is returned by this function. Keys are group names as defined
 * in the index file and values are atom_selection_t structures. Even
 * groups with no atoms are included in this dictionary.
 * 
 * @paragraph Selecting group from the dictionary
 * Use function dict_get() and typecast the returning value as (atom_selection_t *).
 * 
 * @paragraph Deallocating the dictionary
 * The dictionary returned by this function must be deallocated when it is no longer
 * needed by using dict_destroy(). The atom selections from this dictionary
 * should NOT be freed separately.
 * 
 * @paragraph On checking sanity of the input
 * Note that the function does NOT check whether a group contains multiple identical atoms.
 * The function however DOES raise error (returns NULL), if it encounters an atom with an unknown atom number.
 * You can remove any duplicate atoms afterwards by applying selection_unique() to the individual atom selections.
 * 
 * @paragraph On atom sorting
 * This function SHOULD work also for systems in which atoms are not sorted
 * by their gmx_atom_number, this behavior is however untested!
 * (Unless you are doing really shady stuff, the above note does not concern you at all.)
 * 
 * @param filename              path to the ndx file
 * @param system                pointer to a system_t structure
 * 
 * @return Pointer to dictionary containing group_name->atom_selection pairs. NULL in case the reading fails.
 * 
 */
dict_t *read_ndx(const char *filename, system_t *system);

#endif /* SELECTION_H */