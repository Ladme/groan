// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#include "analysis_tools.h"

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

void add_selected_atom(
        size_t **output_atoms,
        size_t *allocated_atoms,
        const size_t array_index,
        const size_t atom_id)
{
    // if output_atoms is NULL, don't do anything
    if (output_atoms == NULL) return;

    // reallocate array, if needed
    if (array_index >= *allocated_atoms) {
        *allocated_atoms *= 2;
        *output_atoms = realloc(*output_atoms, *allocated_atoms * sizeof(size_t));
    }

    // add the atom_id to the selection
    (*output_atoms)[array_index] = atom_id;
}

size_t select_atoms(
        const system_t *system,
        const size_t *input_atom_ids,
        const size_t n_input_atoms,
        size_t **output_atom_ids,
        const char *match_string,
        int (*match_function)(const atom_t *, const char *))
{
    if (system == NULL) return 0;

    size_t alloc_ids = INITIAL_SELECTION_SIZE;
    *output_atom_ids = malloc(alloc_ids * sizeof(size_t));
    

    // if string is NULL, return list of indices of all atoms
    if (match_string == NULL) {
        for (size_t i = 0; i < system->n_atoms; ++i) {
            (*output_atom_ids)[i] = i;
        }
        return system->n_atoms;
    }

    // split match_string into individual elements
    // copy the match_string into new string that can be modified
    char *to_match = malloc(MAX_MATCH_STRING_LEN);
    strncpy(to_match, match_string, MAX_MATCH_STRING_LEN - 1);
    to_match[MAX_MATCH_STRING_LEN - 1] = '\0';
    
    // next do the actual splitting
    char **elements = NULL;
    size_t n_elements = strsplit(to_match, &elements, " ");
    if (n_elements <= 0) return 0;

    // if the array of input atoms is not provided, use all atoms from the system
    size_t n_atoms = 0;
    loop_mode_t mode = all;
    if (n_input_atoms == 0 || input_atom_ids == NULL) {
        n_atoms = system->n_atoms;
    } else {
        n_atoms = n_input_atoms;
        mode = selection;
    }

    // loop through atoms
    size_t selection_iterator = 0;
    for (size_t i = 0; i < n_atoms; ++i) {
        
        // select the corresponding atom
        const atom_t *atom = NULL;
        if (mode == all) {
            atom = &(system->atoms[i]);
        } else {
            atom = &(system->atoms[input_atom_ids[i]]);
        }
        
        // and for each atom try matching the match function to individual match elements
        for (size_t j = 0; j < n_elements; ++j) {
            if (match_function(atom, elements[j])) {

                // if a match is found, add the atom index to selection
                if (mode == all) {
                    add_selected_atom(output_atom_ids, &alloc_ids, selection_iterator, i);
                } else {
                    add_selected_atom(output_atom_ids, &alloc_ids, selection_iterator, input_atom_ids[i]);
                }

                ++selection_iterator;
                break;
            }
        }
    }

    // we have to free memory allocated for match elements
    free(elements);
    free(to_match);

    // we return the number of selected atoms
    return selection_iterator;
}

size_t join_selections(
        size_t **target,
        const size_t *selection1,
        const size_t n_selection1,
        const size_t *selection2,
        const size_t n_selection2)
{
    *target = malloc(n_selection1 * sizeof(size_t) + n_selection2 * sizeof(size_t));

    memcpy(*target, selection1, n_selection1 * sizeof(size_t));
    memcpy(*target + n_selection1, selection2, n_selection2 * sizeof(size_t));

    return n_selection1 + n_selection2;
}

size_t select_geometry(
        const system_t *system,
        const size_t *atom_ids,
        const size_t n_selected,
        size_t **output_atom_ids,
        const vec_t *coordinates,
        const vec_t center,
        const geometry_t geometry,
        const void *geometry_definition)
{
    if (system == NULL && (atom_ids == NULL || n_selected == 0 || coordinates == NULL)) {
        return 1;
    }

    size_t alloc_ids = INITIAL_SELECTION_SIZE;
    if (output_atom_ids != NULL) *output_atom_ids = malloc(alloc_ids * sizeof(size_t));

    // if the array of atom_ids is not provided, use all atoms from the system
    size_t n_atoms = 0;
    loop_mode_t mode = all;
    if (n_selected == 0 || atom_ids == NULL) {
        n_atoms = system->n_atoms;
    } else {
        n_atoms = n_selected;
        mode = selection;
    }

    size_t selection_iterator = 0;
    for (size_t i = 0; i < n_atoms; ++i) {

        // select the corresponding atom index
        size_t atom_index = 0;
        if (mode == all) {
            atom_index = i;
        } else {
            atom_index = atom_ids[i];
        }

        // select the corresponding coordinates
        const float *coor = NULL;

        if (coordinates == NULL) {
            coor = system->atoms[atom_index].position;
        } else {
            coor = coordinates[atom_index];
        }

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

            if (coor[index] - center[index] > cyl_def[1] && 
                coor[index] - center[index] < cyl_def[2] &&
                distance2D(coor, center, plane) < cyl_def[0]) {
                    add_selected_atom(output_atom_ids, &alloc_ids, selection_iterator, atom_index);
                    selection_iterator++;
                }

        } else if (geometry == box) {
            const float *box_def = (float *) geometry_definition;
            short selected = 1;

            // loop through all three dimensions of the box
            for (short j = 0; j < 3; ++j) {
                if (coor[j] - center[j] < box_def[j*2] ||
                    coor[j] - center[j] > box_def[j*2 + 1]) {
                        selected = 0;
                        break;
                }
            }

            if (selected) {
                add_selected_atom(output_atom_ids, &alloc_ids, selection_iterator, atom_index);
                selection_iterator++;
            }
        } else if (geometry == sphere) {
            const float *sphere_def = (float *) geometry_definition;

            if (distance3D(coor, center) < *sphere_def) {
                add_selected_atom(output_atom_ids, &alloc_ids, selection_iterator, atom_index);
                selection_iterator++;
            }
        }

    }

    return selection_iterator;
}

float distance2D(const float *particle1, const float *particle2, const plane_t plane)
{   
    float dim1_d = 0;
    float dim2_d = 0;

    if (plane == xy) {
        dim1_d = particle1[0] - particle2[0];
        dim2_d = particle1[1] - particle2[1];
    } else if (plane == xz) {
        dim1_d = particle1[0] - particle2[0];
        dim2_d = particle1[2] - particle2[2];
    } else {
        dim1_d = particle1[1] - particle2[1];
        dim2_d = particle1[2] - particle2[2];
    }

    return sqrt( dim1_d*dim1_d + dim2_d*dim2_d );
}

float distance3D(const vec_t particle1, const vec_t particle2)
{
    float xd = particle1[0] - particle2[0];
    float yd = particle1[1] - particle2[1];
    float zd = particle1[2] - particle2[2];

    return sqrt( xd*xd + yd*yd + zd*zd );
}

int center_of_geometry(
        const system_t *system,
        const size_t *atom_ids,
        const size_t n_selected,
        vec_t center,
        const vec_t *coordinates)
{
    if (system == NULL && (atom_ids == NULL || n_selected == 0 || coordinates == NULL)) {
        return 1;
    }

    // if no coordinates are provided, use coordinates from system
    if (coordinates == NULL) {
        for (size_t i = 0; i < n_selected; ++i) {
            center[0] += system->atoms[atom_ids[i]].position[0];
            center[1] += system->atoms[atom_ids[i]].position[1];
            center[2] += system->atoms[atom_ids[i]].position[2];
        }
    }

    // else, use the supplied coordinates
    else {
        for (size_t i = 0; i < n_selected; i++) {
            center[0] += coordinates[atom_ids[i]][0];
            center[1] += coordinates[atom_ids[i]][1];
            center[2] += coordinates[atom_ids[i]][2];
        }
    }

    center[0] /= n_selected;
    center[1] /= n_selected;
    center[2] /= n_selected;

    return 0;
}