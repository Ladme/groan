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
    atom_selection_t *output_atoms = malloc(sizeof(atom_selection_t) + alloc_ids * sizeof(atom_t *));
    output_atoms->n_atoms = 0;

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
    size_t selection_iterator = 0;
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

atom_selection_t *select_system(system_t *system)
{
    atom_selection_t *output_atoms = malloc(sizeof(atom_selection_t) + system->n_atoms * sizeof(atom_t *));

    for (size_t i = 0; i < system->n_atoms; ++i) {
        output_atoms->atoms[i] = &(system->atoms[i]);
    }

    output_atoms->n_atoms = system->n_atoms;

    return output_atoms;
}


atom_selection_t *selection_cat(const atom_selection_t *selection1, const atom_selection_t *selection2)
{
    // allocate memory for output_atoms
    atom_selection_t *output_atoms = malloc(sizeof(atom_selection_t) + (selection1->n_atoms + selection2->n_atoms) * sizeof(atom_t *));
    output_atoms->n_atoms = selection1->n_atoms + selection2->n_atoms;

    // cat input selections
    memcpy(output_atoms->atoms, selection1->atoms, selection1->n_atoms * sizeof(atom_t *));
    memcpy(output_atoms->atoms + selection1->n_atoms, selection2->atoms, selection2->n_atoms * sizeof(atom_t *));

    return output_atoms;
}


atom_selection_t *select_geometry(
        const atom_selection_t *input_atoms,
        const vec_t center,
        const geometry_t geometry,
        const void *geometry_definition)
{
    // allocate memory for output_atoms
    size_t alloc_ids = INITIAL_SELECTION_SIZE;
    atom_selection_t *output_atoms = malloc(sizeof(atom_selection_t) + alloc_ids * sizeof(atom_t *));
    output_atoms->n_atoms = 0;

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

            // compare all three dimensions of the box
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

float distance2D(const vec_t particle1, const vec_t particle2, const plane_t plane)
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

int center_of_geometry(const atom_selection_t *input_atoms, vec_t center)
{
    if (input_atoms == NULL || input_atoms->n_atoms == 0) return 1;

    for (size_t i = 0; i < input_atoms->n_atoms; ++i) {
        center[0] += input_atoms->atoms[i]->position[0];
        center[1] += input_atoms->atoms[i]->position[1];
        center[2] += input_atoms->atoms[i]->position[2];
    }

    center[0] /= input_atoms->n_atoms;
    center[1] /= input_atoms->n_atoms;
    center[2] /= input_atoms->n_atoms;

    return 0;
}