// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#include "selection.h"
#include "analysis_tools.h"

/*! @brief Maximal number of query segments for smart_select(). These are two query segments: >resname POPC< && >name PO4< */
static const size_t MAX_QUERY_SEGMENTS = 50;

/*! @brief Initial number of atom indices in selection array. See function select_atoms(). */
static const size_t INITIAL_SELECTION_SIZE = 64;

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

int strsplit(char *string, char ***array, const char *delim)
{
    // allocate memory for array
    size_t allocated_items = 2;
    *array = malloc(allocated_items * sizeof(char *));
    if (*array == NULL) return 0;

    // initialize strtok
    char *word = NULL;
    if ((word = strtok(string, delim)) == NULL) {
        free(*array);
        *array = NULL;
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
                *array = NULL;
                return 0;
            }

            (*array) = new_array; 
        }
        (*array)[items] = word;
        ++items;
    }

    return items;
}

void strstrip(char *string)
{
    char *p = string;
    int l = strlen(p);

    if (l == 0) return;

    while(isspace(p[l - 1])) p[--l] = 0;
    while(* p && isspace(* p)) ++p, --l;

    memmove(string, p, l + 1);    
}

void strremwhite(char *string)
{
    char *d = string;

    do {
        while (isspace(*d)) {
            ++d;
        }
    *string = *d;
    string++;
    d++;
    } while (*string);
}

int match_residue_name(const atom_t *atom, const char *string)
{
    if (string == NULL) return 1;
    return !strcmp(atom->residue_name, string);
}

int match_residue_num(const atom_t *atom, const char *string)
{
    if (string == NULL) return 1;
    groint_t residue_number = 0;
    if (sscanf(string, "%u", &residue_number) != 1) {
        return 0;
    }

    return (atom->residue_number == residue_number);
}

int match_atom_name(const atom_t *atom, const char *string) 
{
    if (string == NULL) return 1;
    return !strcmp(atom->atom_name, string);
}

int match_atom_num(const atom_t *atom, const char *string)
{
    if (string == NULL) return 1;
    size_t atom_num = 0;
    if (sscanf(string, "%zu", &atom_num) != 1) {
        return 0;
    }
    return (atom->gmx_atom_number == atom_num);
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

void selection_add(
        atom_selection_t **output_atoms, 
        size_t *allocated_atoms, 
        atom_selection_t *atoms_to_add)
{
    for (size_t i = 0; i < atoms_to_add->n_atoms; ++i) {
        selection_add_atom(output_atoms, allocated_atoms, atoms_to_add->atoms[i]);
    }
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
    char *to_match = calloc(strlen(match_string) + 1, 1);
    strcpy(to_match, match_string);
    
    // next do the actual splitting
    char **elements = NULL;
    int n_elements = strsplit(to_match, &elements, " ");
    if (n_elements <= 0) {
        free(elements);
        free(to_match);
        return output_atoms;
    }

    // loop through atoms
    for (size_t i = 0; i < input_atoms->n_atoms; ++i) {
        // and for each atom try matching the match function to individual match elements
        for (int j = 0; j < n_elements; ++j) {
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

size_t selection_remove_legacy(atom_selection_t *selection_result, const atom_selection_t *selection_sub)
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

size_t selection_remove(atom_selection_t *selection_result, const atom_selection_t *selection_sub)
{
    // if the selections are identical, just empty selection_result and return its original length
    if (selection_result == selection_sub) {
        size_t original_length = selection_result->n_atoms;
        selection_empty(selection_result);
        return original_length;
    }

    size_t removed_atoms = 0;
    // loop through the selection_sub, removing matching atoms from selection_result
    for (size_t i = 0; i < selection_sub->n_atoms; ++i) {
        removed_atoms += selection_remove_atom(selection_result, selection_sub->atoms[i]);
    }

    return removed_atoms;
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

int selection_isin(atom_selection_t *selection, atom_t *atom)
{
    for (size_t i = 0; i < selection->n_atoms; ++i) {
        if (atom == selection->atoms[i]) return 1;
    }

    return 0;
}

size_t selection_getnres(atom_selection_t *selection)
{
    size_t total_residues = 0;
    // we allocate enough memory to hold N numbers where N is the number of atoms in selection
    groint_t *residue_numbers = calloc(selection->n_atoms, sizeof(groint_t));

    // loop through the selection and get residue number for each atom
    for (size_t i = 0; i < selection->n_atoms; ++i) {
        int found = 0;
        groint_t resid = selection->atoms[i]->residue_number;
        // search for the residue number in the array of residue numbers
        for (size_t j = 0; j < total_residues; ++j) {
            if (resid == residue_numbers[j]) {
                found = 1;
                break;
            }
        }

        // if the residue is unique, add the residue number into the array and increase the counter
        if (!found) {
            residue_numbers[total_residues] = resid;
            total_residues++;
        }
    }

    free(residue_numbers);

    return total_residues;
}


list_t *selection_getresnames(atom_selection_t *selection)
{
    list_t *resnames = list_create();
    if (resnames == NULL) return NULL;

    for (size_t i = 0; i < selection->n_atoms; ++i) {
        char *resname = selection->atoms[i]->residue_name;
        // if the residue name does not exist in the list
        if (list_index(resnames, resname) < 0) {
            list_append(&resnames, resname);
        }
    }
    
    return resnames;
}

size_t selection_splitbyres(atom_selection_t *selection, atom_selection_t ***split)
{
    // if selection does not exist or is empty, do not allocate any memory for split and return 0
    if (selection == NULL || selection->n_atoms <= 0) return 0;

    // allocate initial memory for the split array
    size_t alloc_residues = 8;
    size_t added_residues = 0;
    *split = calloc(alloc_residues, sizeof(atom_selection_t *));

    short *positions_assigned = calloc(selection->n_atoms, sizeof(short));

    for (size_t i = 0; i < selection->n_atoms; ++i) {
        // if this atom has already been assigned, skip it
        if (positions_assigned[i]) continue;
        atom_t *atom1 = selection->atoms[i];

        // reallocate memory for the split array, if needed
        if (added_residues >= alloc_residues) {
            alloc_residues *= 2;
            *split = realloc(*split, alloc_residues * sizeof(atom_selection_t *));
        }

        // allocate memory for the new atom_selection_t structure in the split array
        size_t alloc_atoms = 64;
        (*split)[added_residues] = selection_create(alloc_atoms);

        // add atom1 to the list of assigned positions
        positions_assigned[i] = 1;
        // and to the atom_selection
        selection_add_atom(&((*split)[added_residues]), &alloc_atoms, atom1);
        
        for (size_t j = i + 1; j < selection->n_atoms; ++j) {
            // if this atom has already been assigned, skip it
            if (positions_assigned[j]) continue;
            atom_t *atom2 = selection->atoms[j];

            if (atom1->residue_number == atom2->residue_number) {
                // add atom2 to the list of assigned positions
                positions_assigned[j] = 1;
                // and to the atom_selection
                selection_add_atom(&((*split)[added_residues]), &alloc_atoms, atom2);
            }
        }

        // increment the residue counter
        ++added_residues;
    }

    free(positions_assigned);

    return added_residues; 
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
        const void *geometry_definition,
        const box_t system_box)
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
            plane_t plane = xy;
            dimension_t line = z;

            // settings for xcylinder
            if (geometry == xcylinder) {
                plane = yz;
                line = x;
            } 
            // settings for ycylinder
            else if (geometry == ycylinder) {
                plane = xz;
                line = y;
            }

            float dist1d = distance1D(atom->position, center, line, system_box);

            if (dist1d > cyl_def[1] && dist1d < cyl_def[2] &&
                distance2D(atom->position, center, plane, system_box) < cyl_def[0]) {
                    selection_add_atom(&output_atoms, &alloc_ids, atom);
                }

        } else if (geometry == box) {
            const float *box_def = (float *) geometry_definition;

            float dist_x = distance1D(atom->position, center, x, system_box);
            float dist_y = distance1D(atom->position, center, y, system_box);
            float dist_z = distance1D(atom->position, center, z, system_box);

            // compare all three dimensions of the box
            if (dist_x > box_def[0] && dist_x < box_def[1] && 
                dist_y > box_def[2] && dist_y < box_def[3] &&
                dist_z > box_def[4] && dist_z < box_def[5]) {
                    selection_add_atom(&output_atoms, &alloc_ids, atom);
                }

        } else if (geometry == sphere) {
            const float *sphere_def = (float *) geometry_definition;

            if (distance3D(atom->position, center, system_box) < *sphere_def) {
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
        const void *geometry_definition,
        const box_t system_box)
{
    atom_selection_t *output = select_geometry(input_atoms, center, geometry, geometry_definition, system_box);
    free(input_atoms);
    return output;
}

/*! @brief Returns selection of atoms from 'all' which are NOT part of 'selection'. */
static atom_selection_t *selection_invert(const atom_selection_t *all, const atom_selection_t *selection)
{
    atom_selection_t *result = selection_copy(all);

    selection_remove(result, selection);

    return result;
}

/*! @brief Parses lexeme translating it to atom selection structure pointer to which is returned. */
static atom_selection_t *parse_lexeme(atom_selection_t *selection, char *lexeme, dict_t *ndx_groups)
{
    // check whether the lexeme contains 'not' or '!'
    int not = 0;
    size_t skip = 0;
    if (strlen(lexeme) >= 1 && memcmp(lexeme, "!", 1) == 0) {
        not = 1;
        skip = 2;
    } else if (strlen(lexeme) >= 3 && memcmp(lexeme, "not", 3) == 0) {
        not = 1;
        skip = 4;
    }
    
    atom_selection_t *result = NULL;
    // select all atoms from the selection
    if (strlen(lexeme + skip) >= 3 && memcmp(lexeme + skip, "all", 3) == 0) {
        result = selection_copy(selection);
    // select atoms based on residue names
    } else if (strlen(lexeme + skip) >= 8 && memcmp(lexeme + skip, "resname", 7) == 0) {
        result = select_atoms(selection, lexeme + skip + 7, &match_residue_name);
    // select atoms based on residue numbers
    } else if (strlen(lexeme + skip) >= 6 && memcmp(lexeme + skip, "resid", 5) == 0) {
        result = select_atoms(selection, lexeme + skip + 5, &match_residue_num);
    // select atoms based on atom names
    } else if (strlen(lexeme + skip) >= 5 && memcmp(lexeme + skip, "name", 4) == 0) {
        result = select_atoms(selection, lexeme + skip + 4, &match_atom_name);
    // select atoms based on atom numbers
    } else if (strlen(lexeme + skip) >= 7 && memcmp(lexeme + skip, "serial", 6) == 0) {
        result = select_atoms(selection, lexeme + skip + 6, &match_atom_num);
    // select atoms based on ndx groups
    } else if (ndx_groups != NULL) {
        // we have to replace the trailing space that is added during lexeme formation
        lexeme[strlen(lexeme) - 1] = 0;
        atom_selection_t *original = (atom_selection_t *) dict_get(ndx_groups, lexeme + skip);
        if (original == NULL) return NULL;

        // we have to copy the selection so we can then free it without disrupting the dictionary
        result = selection_copy(original);
    }

    // invert the selection, if 'not' or '!' is in the lexeme
    if (not && result != NULL) {
        atom_selection_t *result_inverted = selection_invert(selection, result);
        free(result);
        return result_inverted;
    }
    
    return result;
}

/*! @brief Expands 'a to b' or 'a - b' macro into a sequence that can be understood by the parser. 
 * 
 * @paragraph Memory Allocation
 * Allocates enough memory to hold the expanded sequence.
 * 
 * @return Zero in case of successful replacement. Else non-zero.
 */
static int expand_to(char **new_string, const char *original_string)
{
    // check whether there is any '-' or 'to' in the original string
    if (strchr(original_string, '-') == NULL && strstr(original_string, "to") == NULL) {
        *new_string = calloc(strlen(original_string) + 1, 1);
        strcpy(*new_string, original_string);
        return 0;
    }

    // split the original string
    char **split = NULL;
    char *string_to_split = calloc(strlen(original_string) + 1, 1);
    strcpy(string_to_split, original_string);
    size_t n_words = strsplit(string_to_split, &split, " \n\t");

    size_t string_allocated = 64;
    *new_string = calloc(string_allocated, 1);
    size_t string_len = 0;

    // loop through the individual words and detect keywords
    for (size_t i = 0; i < n_words; ++i) {
        // sanity check: -/to cannot be at the start or the end of the query
        if (strcmp(split[i], "-") == 0 || strcmp(split[i], "to") == 0) {
            
            if (i == 0 || i == n_words - 1) {
                free(string_to_split);
                free(split);
                return 1;
            }

            // try reading the starting and ending value of the loop
            int start = 0;
            int end = 0;
            if (sscanf(split[i - 1], "%d", &start) != 1 || sscanf(split[i + 1], "%d", &end) != 1) {
                free(string_to_split);
                free(split);
                return 1;
            }

            // loop must start at the lower value
            if (start > end) {
                free(string_to_split);
                free(split);
                return 1;
            }

            // we start from 'start + 1' and we end at 'end - 1' because the first and last number are included by the block in 'else'
            for (int j = start + 1; j < end; ++j) {
                char number[50] = "";
                sprintf(number, "%d", j);
                if (string_len + strlen(number) + 2 >= string_allocated) {
                    while (string_len + strlen(number) + 2 >= string_allocated) string_allocated *= 2;
                    *new_string = realloc(*new_string, string_allocated);
                }

                // add new number to the string
                strcpy( (*new_string) + string_len, number);
                string_len += strlen(number);
                // we have to add space at the end of the string
                (*new_string)[string_len] = ' ';
                (*new_string)[string_len + 1] = 0;
                ++string_len;
            }

        } else {
            if (string_len + strlen(split[i]) + 2 >= string_allocated) {
                while (string_len + strlen(split[i]) + 2 >= string_allocated) string_allocated *= 2;
                *new_string = realloc(*new_string, string_allocated);
            }
            
            strcpy( (*new_string) + string_len, split[i]);
            string_len += strlen(split[i]);
            // we have to add space to the end of the string
            (*new_string)[string_len] = ' ';
            (*new_string)[string_len + 1] = 0;
            ++string_len;

        }
    }

    free(string_to_split);
    free(split);

    //printf("TO expanded: %s, length expected: %d, length real: %d, allocated: %d\n", *new_string, string_len, strlen(*new_string), string_allocated);
    return 0;
}

static atom_selection_t *parse_query(atom_selection_t *selection, char *query, dict_t *ndx_groups)
{
    size_t query_len = strlen(query);
    // split the expanded query into individual lexemes
    char **split = NULL;
    size_t n_words = strsplit(query, &split, " \n\t");
    //printf("Split query into %d words\n", n_words);

    // prepare an array with token definitions
    atom_selection_t *tokens[MAX_QUERY_SEGMENTS];

    // loop through the lexemes, translating them to atom selections
    char *lexeme = calloc(2 * query_len + 1, 1);
    char *operators[MAX_QUERY_SEGMENTS];
    size_t counter = 0;
    size_t n_tokens = 0;
    size_t n_operators = 0;
    for (size_t i = 0; i < n_words; ++i) {
        // if parenthesis is detected
        if (strchr(split[i], '(')) {
            //printf("Detected an opening parenthesis.\n");
            char *block = calloc(2 * query_len + 1, 1);
            size_t block_len = 0;
            size_t par = 0;
            size_t j = i;
            // loop through all the words until we find a matching ')'
            for (; j < n_words; ++j) {
                //printf("Adding %s to block.\n", split[j]);

                // count the number of parenthesis in the word
                for (size_t k = 0; split[j][k] != 0; ++k) {
                    if (split[j][k] == '(') ++par;
                    else if (split[j][k] == ')') --par;
                }
                
                // copy the words into a new query 'block'
                size_t add_block_len = 0;
                if (j == i) {
                    if (par == 0) {
                        // remove the opening and the closing parenthesis
                        split[j][strlen(split[j])] = 0;
                        strcpy(block, split[j] + 1);
                        add_block_len = strlen(split[j]) - 2;
                    } else {
                        // remove the initial parenthesis from the first word
                        strcpy(block + block_len, split[j] + 1);
                        add_block_len = strlen(split[j]) - 1;
                    }
                } else if (par == 0) {
                    if (strlen(split[j]) != 1) {
                        strcpy(block + block_len, split[j]);
                        // remove the ending parenthesis from the last word
                        add_block_len = strlen(split[j]) - 1;
                    }
                } else {
                    // copy the entire word
                    strcpy(block + block_len, split[j]);
                    add_block_len = strlen(split[j]);
                }
                
                block_len += add_block_len;
                if (add_block_len > 0) {
                    block[block_len] = ' ';
                    block[block_len + 1] = 0;
                    ++block_len;
                }
                
                //printf("Current block: %s\n", block);

                if (par == 0) break;
            }

            // parse the block as a new query and save the output into a list of tokens
            atom_selection_t *parsed_block = parse_query(selection, block, ndx_groups);
            if (parsed_block == NULL) {
                //fprintf(stderr, "Could not parse block %s\n", block);
                free(split);
                free(block);
                free(lexeme);
                for (size_t i = 0; i < n_tokens; ++i) free(tokens[i]);
                return NULL;
            }

            // check whether there is 'not' operator in front of the block
            if (!strcmp(lexeme, "! ") || !strcmp(lexeme, "not ")) {
                //printf("Detected NOT operator.\n");
                tokens[n_tokens] = selection_invert(selection, parsed_block);
                free(parsed_block);
                memset(lexeme, 0, strlen(lexeme));
                counter = 0;
            } else {
                tokens[n_tokens] = parsed_block;
            }

            // if there are any characters in front of a parenthesis (other than '!' or 'not'), raise a syntax error
            if (strlen(lexeme) != 0) {
                ++n_tokens;
                free(split);
                free(block);
                free(lexeme);
                for (size_t i = 0; i < n_tokens; ++i) free(tokens[i]);
                return NULL;
            }
            
            ++n_tokens;
            // continue parsing the input at the end of the block
            i = j;

            free(block);
        } 

        else if ( ( strlen(split[i]) == 2 && (strcmp(split[i], "&&") == 0 || strcmp(split[i], "||") == 0 || strcmp(split[i], "or") == 0)) ||
             ( strlen(split[i]) == 3 && strcmp(split[i], "and") == 0) ) {

            // add the operator to the list of operators
            //printf("Adding operator %s\n", split[i]);
            operators[n_operators] = split[i];
            ++n_operators;

            if (strlen(lexeme) == 0) continue;

            //printf("Parsing lexeme %s\n", lexeme);

            // parse the lexeme
            atom_selection_t *parsed = parse_lexeme(selection, lexeme, ndx_groups);
            if (parsed == NULL) {
                //fprintf(stderr, "Could not parse lexeme %s\n", lexeme);
                free(split);
                free(lexeme);
                for (size_t i = 0; i < n_tokens; ++i) free(tokens[i]);
                return NULL;
            }

            tokens[n_tokens] = parsed;
            ++n_tokens;
            memset(lexeme, 0, strlen(lexeme));
            counter = 0;
        } else {
            //printf("Adding %s to lexeme\n", split[i]);
            strcpy(lexeme + counter, split[i]);
            counter += strlen(split[i]);
            lexeme[counter] = ' ';
            lexeme[counter + 1] = 0;
            counter += 1;
            //printf("Added %s to lexeme\n", split[i]);

            // if this is the last word, parse the current lexeme
            if (i == n_words - 1) {

                //printf("Parsing lexeme %s\n", lexeme);
                
                atom_selection_t *parsed = parse_lexeme(selection, lexeme, ndx_groups);
                if (parsed == NULL) {
                    //fprintf(stderr, "Could not parse lexeme %s\n", lexeme);
                    free(split);
                    free(lexeme);
                    for (size_t i = 0; i < n_tokens; ++i) free(tokens[i]);
                    return NULL;
                }

                tokens[n_tokens] = parsed;
                ++n_tokens;
            }
        }
    }

    //printf("Detected %lu query segments that were translated to tokens.\n", n_tokens);

    // check that the number of operators corresponds to the number of tokens
    if (n_tokens == 0 || n_operators + 1 != n_tokens) {
        free(lexeme);
        free(split);
        for (size_t i = 0; i < n_tokens; ++i) free(tokens[i]);
        return NULL;
    }

    // if there is only one token and no operators, return the first token
    if (n_tokens == 1) {
        free(lexeme);
        free(split);
        return tokens[0];
    }

    // combine atom selections
    atom_selection_t *final = tokens[0];
    for (size_t i = 1; i < n_tokens; ++i) {
        // AND operation
        if (strcmp(operators[i - 1], "&&") == 0 || strcmp(operators[i - 1], "and") == 0) {
            //printf("Performing AND operation\n");
            final = selection_intersect_d(final, tokens[i]);
        // OR operation
        } else if (strcmp(operators[i - 1], "||") == 0 || strcmp(operators[i - 1], "or") == 0) {
            //printf("Performing OR operation\n");
            // must be unique cat as we do not want duplicate atoms in the final selection!
            final = selection_cat_unique_d(final, tokens[i]);
        }
    }

    free(lexeme);
    free(split);
    

    return final;


}

atom_selection_t *smart_select(atom_selection_t *selection, const char *query, dict_t *ndx_groups)
{
    // check that the query is valid
    if (query == NULL) {
        atom_selection_t *copy = selection_copy(selection);
        return copy;
    };

    // check that the query ends with '\0'
    char *pos = strchr(query, 0);
    if (pos == NULL) return NULL;

    // check that the selection is valid
    if (selection == NULL) return NULL;

    // check that the number of '(' and ')' match each other
    int par_open = 0;
    int par_close = 0;
    for (size_t i = 0; query[i] != 0; ++i) {
        if (query[i] == '(') ++par_open;
        else if (query[i] == ')') ++par_close;
    }
    if (par_open != par_close) return NULL;

    // expand the 'to' macro
    char *query_expanded = NULL;
    if (expand_to(&query_expanded, query) != 0) {
        free(query_expanded);
        return NULL;
    }

    atom_selection_t *final = parse_query(selection, query_expanded, ndx_groups);
    free(query_expanded);

    return final;
}


dict_t *read_ndx(const char *filename, system_t *system)
{
    FILE *ndx = fopen(filename, "r");
    if (ndx == NULL) {
        return NULL;
    }

    dict_t *ndx_selections = dict_create();
    if (ndx_selections == NULL) {
        fclose(ndx);
        return NULL;
    }

    char line[1024] = "";
    char **split = NULL;
    char current_group[100] = "";
    atom_selection_t *current_selection = NULL;
    size_t alloc_atoms = 0;
    while (fgets(line, 1024, ndx) != NULL) {
        // remove newline character
        line[strcspn(line, "\n")] = 0;
        // split the line
        int n_items = strsplit(line, &split, " ");

        // negative n_items means that an error has occured
        if (n_items < 0) {
            dict_destroy(ndx_selections);
            free(split);
            fclose(ndx);
            return NULL;
        }

        // zero means that the line is completely empty
        else if (n_items == 0) {
            free(split);
            split = NULL;
            continue;
        }

        // one may also mean that the line is empty
        else if (n_items == 1 && strcmp(split[0], "\n") == 0) {
            free(split);
            split = NULL;
            continue;
        }

        // if three or more items are detected, try getting the group name
        else if (n_items >= 3 && strcmp(split[0], "[") == 0) {
            if (strcmp(split[n_items - 1], "]") == 0) {
                // if group name is detected, add the previous selection to the dictionary
                if (current_selection != NULL) {
                    dict_set(ndx_selections, current_group, current_selection, sizeof(atom_selection_t) + alloc_atoms * sizeof(atom_t *));
                    free(current_selection);
                }

                int offset = 0;
                for (int i = 1; i < n_items - 1; ++i) {
                    // add white space
                    if (i != 1) {
                        current_group[offset] = ' ';
                        ++offset;
                    }

                    strncpy(current_group + offset, split[i], 99 - offset);
                    offset += strlen(split[i]);
                }
                
                alloc_atoms = INITIAL_SELECTION_SIZE;
                current_selection = selection_create(alloc_atoms);
                free(split);
                split = NULL;
                continue;

            // if there is initial [, but no closing ], raise an error
            } else {
                free(current_selection);
                dict_destroy(ndx_selections);
                free(split);
                fclose(ndx);
                return NULL;
            }
        }

        // load all atoms to current selection
        for (int i = 0; i < n_items; ++i) {
            size_t atom_n = 0;

            // load atom number
            if (sscanf(split[i], "%lu", &atom_n) != 1) {
                free(current_selection);
                dict_destroy(ndx_selections);
                free(split);
                fclose(ndx);
                return NULL;   
            }

            // we test whether the atom actually has the correct gmx_atom_number
            // if not, we loop through the entire array searching for the correct atom
            // this may seem unnecessary, but the atoms in the system do not have to be sorted based on gmx_atom_number 
            if ((atom_n - 1) < system->n_atoms && system->atoms[atom_n - 1].gmx_atom_number == atom_n) {
                selection_add_atom(&current_selection, &alloc_atoms, &system->atoms[atom_n - 1]);
                continue;
            }

            // for broken systems, we have to use this
            for (size_t j = 0; j < system->n_atoms; ++j) {
                if (system->atoms[j].gmx_atom_number == atom_n) {
                    selection_add_atom(&current_selection, &alloc_atoms, &system->atoms[j]);
                    break;
                }
            }

            // if no suitable atom has been found
            free(current_selection);
            dict_destroy(ndx_selections);
            free(split);
            fclose(ndx);
            return NULL;
        }
        
        // free the array of splitted line
        free(split);
        split = NULL;

    }

    // add the last ndx group
    if (current_selection != NULL) {
        dict_set(ndx_selections, current_group, current_selection, sizeof(atom_selection_t) + alloc_atoms * sizeof(atom_t *));
        free(current_selection);
    }

    free(split);
    fclose(ndx);

    return ndx_selections;
}