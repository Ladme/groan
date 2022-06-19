// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#include "selection.h"
#include "analysis_tools.h"

/* Simple function for qsort comparison of atom numbers */
static int compare_atomnum(const void *x, const void *y)
{
    const atom_t *atom1 = *((const atom_t **) x);
    const atom_t *atom2 = *((const atom_t **) y);

    return ( atom1->atom_number - atom2->atom_number );
}

/* Simple function for qsort comparison of gmx atom numbers */
static int compare_gmxatomnum(const void *x, const void *y)
{
    const atom_t *atom1 = *((const atom_t **) x);
    const atom_t *atom2 = *((const atom_t **) y);

    return ( atom1->gmx_atom_number - atom2->gmx_atom_number );
}

size_t strsplit(char *string, char ***array, const char *delim)
{
    // allocate memory for array
    size_t allocated_items = 2;
    *array = malloc(allocated_items * sizeof(char *));
    if (*array == NULL) return 0;

    // initialize strtok
    char *word;
    if ((word = strtok(string, delim)) == NULL) {
        free(*array);
        return 0;
    }

    (*array)[0] = word;

    // repeat strtok until the string is split completely
    size_t items = 1;
    while ((word = strtok(NULL, delim)) != NULL) {
        // reallocate, if needed
        if (items >= allocated_items) {
            allocated_items *= 2;
            char **new_array = realloc(*array, allocated_items * sizeof(char *));
            if (new_array == NULL) {
                free(*array);
                return 0;
            }

            (*array) = new_array; 
        }
        (*array)[items] = word;
        ++items;
    }

    return items;
}

int match_residue_name(const atom_t *atom, const char *string)
{
    if (string == NULL) return 1;
    return !strcmp(atom->residue_name, string);
}

int match_atom_name(const atom_t *atom, const char *string) 
{
    if (string == NULL) return 1;
    return !strcmp(atom->atom_name, string);
}

atom_selection_t *selection_create(size_t items)
{
    atom_selection_t *selection = malloc(sizeof(atom_selection_t) + items * sizeof(atom_t *));
    selection->n_atoms = 0;

    return selection;
}

atom_selection_t *selection_copy(const atom_selection_t *selection)
{
    atom_selection_t *output = selection_create(selection->n_atoms);
    memcpy(output->atoms, selection->atoms, selection->n_atoms * sizeof(atom_t *));
    output->n_atoms = selection->n_atoms;

    return output;
}

atom_selection_t *selection_copy_d(atom_selection_t *selection)
{
    atom_selection_t *output = selection_copy(selection);
    free(selection);
    return output;
}

void selection_empty(atom_selection_t *selection)
{
    // yes, this is all it does
    selection->n_atoms = 0;
}

void selection_add_atom(
        atom_selection_t **output_atoms,
        size_t *allocated_atoms,
        atom_t *atom)
{
    // reallocate array, if needed
    if ((*output_atoms)->n_atoms + 1 >= *allocated_atoms) {
        *allocated_atoms *= 2;
        *output_atoms = realloc(*output_atoms, sizeof(atom_selection_t) + *allocated_atoms * sizeof(atom_t *));
    }

    // add the atom_id to the selection
    (*output_atoms)->atoms[(*output_atoms)->n_atoms] = atom;
    // increase the number of atoms in selection
    (*output_atoms)->n_atoms++;
}

atom_selection_t *select_atoms(
        const atom_selection_t *input_atoms,
        const char *match_string,
        int (*match_function)(const atom_t *, const char *))
{
    // allocate memory for output_atoms
    size_t alloc_ids = INITIAL_SELECTION_SIZE;
    atom_selection_t *output_atoms = selection_create(alloc_ids);

    if (input_atoms == NULL || input_atoms->n_atoms == 0) return output_atoms;

    // split match_string into individual elements
    // copy the match_string into new string that can be modified
    char *to_match = malloc(MAX_MATCH_STRING_LEN);
    strncpy(to_match, match_string, MAX_MATCH_STRING_LEN - 1);
    to_match[MAX_MATCH_STRING_LEN - 1] = '\0';
    
    // next do the actual splitting
    char **elements = NULL;
    size_t n_elements = strsplit(to_match, &elements, " ");
    if (n_elements <= 0) return 0;

    // loop through atoms
    for (size_t i = 0; i < input_atoms->n_atoms; ++i) {
        // and for each atom try matching the match function to individual match elements
        for (size_t j = 0; j < n_elements; ++j) {
            if (match_function(input_atoms->atoms[i], elements[j])) {

                selection_add_atom(&output_atoms, &alloc_ids, input_atoms->atoms[i]);
                break;
            }
        }
    }

    // we have to free memory allocated for match elements
    free(elements);
    free(to_match);

    // return the pointer to output atoms
    return output_atoms;
}

atom_selection_t *select_atoms_d(
        atom_selection_t *input_atoms,
        const char *match_string,
        int (*match_function)(const atom_t *, const char *))
{
    atom_selection_t *output = select_atoms(input_atoms, match_string, match_function);
    free(input_atoms);
    return output;
}

atom_selection_t *select_system(system_t *system)
{
    atom_selection_t *output_atoms = selection_create(system->n_atoms);

    for (size_t i = 0; i < system->n_atoms; ++i) {
        output_atoms->atoms[i] = &(system->atoms[i]);
    }

    output_atoms->n_atoms = system->n_atoms;

    return output_atoms;
}


atom_selection_t *selection_cat(const atom_selection_t *selection1, const atom_selection_t *selection2)
{
    // create output atoms
    atom_selection_t *output_atoms = selection_create(selection1->n_atoms + selection2->n_atoms);
    output_atoms->n_atoms = selection1->n_atoms + selection2->n_atoms;

    // cat input selections
    memcpy(output_atoms->atoms, selection1->atoms, selection1->n_atoms * sizeof(atom_t *));
    memcpy(output_atoms->atoms + selection1->n_atoms, selection2->atoms, selection2->n_atoms * sizeof(atom_t *));

    return output_atoms;
}

atom_selection_t *selection_cat_d(atom_selection_t *selection1, atom_selection_t *selection2)
{
    atom_selection_t *output = selection_cat(selection1, selection2);
    free(selection1);
    if (selection1 != selection2) free(selection2);
    return output;
}


atom_selection_t *selection_cat_unique(const atom_selection_t *selection1, const atom_selection_t *selection2)
{
    // allocate enough memory to hold both selections in full
    // this may be VERY memory inefficient, if the selections significantly overlap, 
    // but it is also faster than reallocation
    atom_selection_t *output_atoms = selection_create(selection1->n_atoms + selection2->n_atoms);
    output_atoms->n_atoms = selection1->n_atoms;

    // copy first selection into the output selection
    memcpy(output_atoms->atoms, selection1->atoms, selection1->n_atoms * sizeof(atom_t *));

    // then loop through atoms of the selection2
    for (size_t i = 0; i < selection2->n_atoms; ++i) {
        int duplicate = 0;
        // and for each atom, check that it is NOT inside selection1
        for (size_t j = 0; j < selection1->n_atoms; ++j) {
            if (selection2->atoms[i] == selection1->atoms[j]) {
                duplicate = 1;
                break;
            }
        }

        // if the atom is not duplicate, add it to output selection
        if (!duplicate) {
            output_atoms->atoms[output_atoms->n_atoms] = selection2->atoms[i];
            // increase the counter of output atoms
            output_atoms->n_atoms++;
        }
    }

    return output_atoms;
}

atom_selection_t *selection_cat_unique_d(atom_selection_t *selection1, atom_selection_t *selection2)
{
    atom_selection_t *output = selection_cat_unique(selection1, selection2);
    free(selection1);
    if (selection1 != selection2) free(selection2);
    return output;
}

atom_selection_t *selection_intersect(const atom_selection_t *selection1, const atom_selection_t *selection2)
{
    // if both selections are the same, then just copy selection1 and skip all this
    if (selection1 == selection2) {
        return selection_copy(selection1);
    }

    size_t alloc_ids = INITIAL_SELECTION_SIZE;
    atom_selection_t *output_atoms = selection_create(alloc_ids);

    for (size_t i = 0; i < selection1->n_atoms; ++i) {
        atom_t *atom = selection1->atoms[i];
        for (size_t j = 0; j < selection2->n_atoms; ++j) {
            if (atom == selection2->atoms[j]) {
                selection_add_atom(&output_atoms, &alloc_ids, atom);
            }
        }
    }

    return output_atoms;
}

atom_selection_t *selection_intersect_d(atom_selection_t *selection1, atom_selection_t *selection2)
{
    atom_selection_t *output = selection_intersect(selection1, selection2);
    free(selection1);
    if (selection1 != selection2) free(selection2);
    return output;
}

size_t selection_remove_atom(atom_selection_t *selection, atom_t *remove)
{
    size_t deleted_atoms = 0;
    // loop through atoms searching for the atom to remove
    for (size_t i = 0; i < selection->n_atoms; ++i) {
        if (remove == selection->atoms[i]) {
            // if the atom is found, move all following atoms one position forward
            for (size_t j = i + 1; j < selection->n_atoms; ++j) {
                selection->atoms[j - 1] = selection->atoms[j];
            }

            // decrease the atom counter in the selection
            selection->n_atoms--;
            // increase the deleted atoms counter
            deleted_atoms++;
            // decrease the index, because we have to read this position again in the following search
            i--;
        }
    }

    return deleted_atoms;
}

size_t selection_remove(atom_selection_t *selection_result, const atom_selection_t *selection_sub)
{
    // if the selections are identical, just empty selection_result and return its original length
    if (selection_result == selection_sub) {
        size_t original_length = selection_result->n_atoms;
        selection_empty(selection_result);
        return original_length;
    }

    // find atoms to be deleted (intersection of the two atom selections)
    atom_selection_t *to_delete = selection_intersect(selection_result, selection_sub);

    // remove the selected atoms
    for (size_t i = 0; i < to_delete->n_atoms; ++i) {
        selection_remove_atom(selection_result, to_delete->atoms[i]);
    }

    size_t delete_n_atoms = to_delete->n_atoms;
    free(to_delete);
    return delete_n_atoms;
}

size_t selection_remove_d(atom_selection_t *selection_result, atom_selection_t *selection_sub)
{
    size_t result = selection_remove(selection_result, selection_sub);
    free(selection_sub);
    return result;
}

size_t selection_unique(atom_selection_t *selection)
{
    size_t deleted_atoms = 0;
    for (size_t i = 0; i < selection->n_atoms; ++i) {
        for (size_t j = i + 1; j < selection->n_atoms; ++j) {
            if (selection->atoms[i] == selection->atoms[j]) {
                for (size_t k = j + 1; k < selection->n_atoms; ++k) {
                    selection->atoms[k - 1] = selection->atoms[k];
                }
                deleted_atoms++;
                selection->n_atoms--;
                j--;
            }
        }
    }

    return deleted_atoms;
}

int selection_compare(
        const atom_selection_t *selection1, 
        const atom_selection_t *selection2)
{
    if (selection1 == selection2) return 1;
    if (selection1->n_atoms != selection2->n_atoms) return 0;

    for (size_t i = 0; i < selection1->n_atoms; ++i) {
        atom_t *atom = selection1->atoms[i];
        int found = 0;
        for (size_t j = 0; j < selection2->n_atoms; ++j) {
            if (atom == selection2->atoms[j]) {
                found = 1;
                break;
            }
        }

        if (!found) return 0;
    }

    return 1;
}

int selection_compare_strict(
        const atom_selection_t *selection1, 
        const atom_selection_t *selection2)
{
    if (selection1 == selection2) return 1;
    if (selection1->n_atoms != selection2->n_atoms) return 0;

    return !(memcmp(selection1->atoms, selection2->atoms, selection1->n_atoms * sizeof(atom_t *)));
}

void selection_renumber(atom_selection_t *selection) 
{
    // residue_numbers is an array for converting old residue numbers to new residue numbers
    // we allocate enough memory to hold n_atoms residue numbers...
    // ...because every atom in selection can be in its own residue
    size_t n_residues = 0;
    groint_t *residue_numbers = malloc(selection->n_atoms * sizeof(int));
    memset(residue_numbers, 0, selection->n_atoms * sizeof(int));

    for (size_t i = 0; i < selection->n_atoms; ++i) {
        atom_t *atom = selection->atoms[i];

        // changing atom number is quite simple
        // note that gro files don't support numbers larger than 99,999 and we must follow this
        // thus, we wrap the atom number to 1, if it reaches 100,000
        atom->atom_number = (i % 99999) + 1;

        // to change residue number, we must loop through the residue numbers array...
        // ...checking whether the residue number was already used or not
        // in this way, atoms corresponding to the same residue will have the same residue number...
        // ...after renumbering no matter where these atoms are positioned in the selection
        int found = 0;
        for (size_t j = 0; j < n_residues; ++j) {
            // if the residue number was already detected...
            if (atom->residue_number == residue_numbers[j]) {
                // ...new residue number will be the index of the old residue number + 1
                // note that we are not making sure that the residue number is lower than 100,000...
                // ...because there is no way that j (or n_residues below) is larger than 99,999...
                // ...if the residue numbers are read correctly from the gro file
                // (due to the nature of residue renumbering)
                atom->residue_number = j + 1;
                
                found = 1;
                break;
            }
        }

        // if the residue number was not found, we add it to the list of residue numbers...
        // ...and set the new residue number of the atom to the highest index plus 1
        if (!found) {
            residue_numbers[n_residues] = atom->residue_number;
            atom->residue_number = n_residues + 1;
            n_residues++;
        }
    }

    free(residue_numbers);
}

void selection_sort(atom_selection_t *selection)
{
    qsort(selection->atoms, selection->n_atoms, sizeof(atom_t *), &compare_atomnum);
}

void selection_sort_gmx(atom_selection_t *selection)
{
    qsort(selection->atoms, selection->n_atoms, sizeof(atom_t *), &compare_gmxatomnum);
}

void selection_fixres(atom_selection_t *selection)
{
    // this is a complicated function

    // first, we allocate memory for an array of atom selections
    // each atom selection will contain atoms corresponding to one residue
    // we allocate enough memory to fit n_atoms atom selections into the array because
    // all residues could be a single atom large
    atom_selection_t **residue_selections = malloc(selection->n_atoms * sizeof(atom_selection_t *));

    // we also have to create an array of allocations 
    // (how much memory was allocated for each atom selection)
    size_t *allocations = malloc(selection->n_atoms * sizeof(size_t));
    memset(allocations, 0, selection->n_atoms * sizeof(size_t));

    // we also create an array for remembering what residue numbers correspond to what residue indices
    groint_t *residue_numbers = malloc(selection->n_atoms * sizeof(size_t));

    // then we loop through the atoms, assigning them to the corresponding residues
    size_t n_residues = 0;
    for (size_t i = 0; i < selection->n_atoms; ++i) {
        atom_t *atom = selection->atoms[i];

        // we loop through residue numbers array, similarly to the selection_renumber function
        int found = 0;
        for (size_t j = 0; j < n_residues; ++j) {
            if (atom->residue_number == residue_numbers[j]) {
                selection_add_atom(&residue_selections[j], &allocations[j], atom);
                found = 1;
            }
        }

        // if no memory has been allocated for the atom selection corresponding to this residue, allocate default value
        if (!found) {
            allocations[n_residues] = INITIAL_SELECTION_SIZE;
            residue_selections[n_residues] = selection_create(allocations[n_residues]);
            // add the atom
            residue_selections[n_residues]->atoms[0] = atom;
            residue_selections[n_residues]->n_atoms++;
            // add the residue to the residue_numbers array
            residue_numbers[n_residues] = atom->residue_number;
            // increase the number of detected residues
            n_residues++;
        }
    }

    // after assigning all atoms to residues, sort all atoms of the residues
    // and assign all residue selections to the original selection
    int iterator = 0;
    for (size_t i = 0; i < n_residues; ++i) {
        // we sort the atoms based on the gmx atom number because that cannot be (easily) changed by the user...
        // ...and we are trying order the atoms in such a way that the gro file is actually usable by gromacs
        selection_sort_gmx(residue_selections[i]);
        memcpy(selection->atoms + iterator, residue_selections[i]->atoms, residue_selections[i]->n_atoms * sizeof(atom_t *));
        // move the selection->atoms pointer
        iterator += residue_selections[i]->n_atoms;
        free(residue_selections[i]);
    }

    // deallocate memory
    free(residue_selections);
    free(allocations);
    free(residue_numbers);
}

system_t *selection_to_system(
        const atom_selection_t *selection, 
        const box_t box, 
        const int step, 
        const float time)
{
    // create new atom selection
    atom_selection_t *new_selection = selection_copy(selection);

    // remove duplicate atoms
    selection_unique(new_selection);

    // fix residues
    selection_fixres(new_selection);

    // create new system_t structure
    system_t *new_system = malloc(sizeof(system_t) + new_selection->n_atoms * sizeof(atom_t));
    memset(new_system, 0, sizeof(system_t) + new_selection->n_atoms * sizeof(atom_t));
    new_system->n_atoms = new_selection->n_atoms;
    memcpy(new_system->box, box, sizeof(box_t));
    new_system->step = step;
    new_system->time = time;

    // copy atoms from the new_selection to new_system
    for (size_t i = 0; i < new_selection->n_atoms; ++i) {
        memcpy(&new_system->atoms[i], new_selection->atoms[i], sizeof(atom_t));
        // change gmx_atom_number
        new_system->atoms[i].gmx_atom_number = i + 1;
    }

    // reassign new_selection to be associated with the new system
    free(new_selection);
    new_selection = select_system(new_system);

    // renumber atoms in the new_system (this is little convoluted)
    selection_renumber(new_selection);

    // deallocate new_selection
    free(new_selection);

    // return pointer to the new system
    return new_system;
}

system_t *selection_to_system_d(
        atom_selection_t *selection, 
        const box_t box, 
        const int step, 
        const float time)
{
    system_t *system = selection_to_system(selection, box, step, time);
    free(selection);
    return system;
}

atom_selection_t *select_geometry(
        const atom_selection_t *input_atoms,
        const vec_t center,
        const geometry_t geometry,
        const void *geometry_definition)
{
    // allocate memory for output_atoms
    size_t alloc_ids = INITIAL_SELECTION_SIZE;
    atom_selection_t *output_atoms = selection_create(alloc_ids);

    if (input_atoms == NULL || input_atoms->n_atoms == 0) return output_atoms;

    // loop through atoms
    for (size_t i = 0; i < input_atoms->n_atoms; ++i) {
        atom_t *atom = input_atoms->atoms[i];

        // decide whether the atom should be included in the selection
        if (geometry == xcylinder || geometry == ycylinder || geometry == zcylinder) {
            const float *cyl_def = (float *) geometry_definition;
            
            // settings for zcylinder
            short index = 2;
            plane_t plane = xy;

            // settings for xcylinder
            if (geometry == xcylinder) {
                index = 0;
                plane = yz;
            } 
            // settings for ycylinder
            else if (geometry == ycylinder) {
                index = 1;
                plane = xz;
            }

            if (atom->position[index] - center[index] > cyl_def[1] && 
                atom->position[index] - center[index] < cyl_def[2] &&
                distance2D(atom->position, center, plane) < cyl_def[0]) {
                    selection_add_atom(&output_atoms, &alloc_ids, atom);
                }

        } else if (geometry == box) {
            const float *box_def = (float *) geometry_definition;

            // compare all three dimensions of the box (this is INTENTIONALLY unrolled loop)
            if (atom->position[0] - center[0] > box_def[0] && 
                atom->position[0] - center[0] < box_def[1] && 
                atom->position[1] - center[1] > box_def[2] &&
                atom->position[1] - center[1] < box_def[3] &&
                atom->position[2] - center[2] > box_def[4] &&
                atom->position[2] - center[2] < box_def[5]) {
                    selection_add_atom(&output_atoms, &alloc_ids, atom);
                }

        } else if (geometry == sphere) {
            const float *sphere_def = (float *) geometry_definition;

            if (distance3D(atom->position, center) < *sphere_def) {
                selection_add_atom(&output_atoms, &alloc_ids, atom);
            }
        }

    }

    return output_atoms;
}

atom_selection_t *select_geometry_d(
        atom_selection_t *input_atoms,
        const vec_t center,
        const geometry_t geometry,
        const void *geometry_definition)
{
    atom_selection_t *output = select_geometry(input_atoms, center, geometry, geometry_definition);
    free(input_atoms);
    return output;
}