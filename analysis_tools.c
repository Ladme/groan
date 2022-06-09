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

size_t select_atoms(
        const system_t *system,
        size_t **atom_ids,
        const char *match_string,
        int (*match_function)(const atom_t *, const char *))
{
    if (system == NULL) return 0;

    size_t alloc_ids = INITIAL_SELECTION_SIZE;
    *atom_ids = malloc(alloc_ids * sizeof(size_t));
    

    // if string is NULL, return list of indices of all atoms
    if (match_string == NULL) {
        for (size_t i = 0; i < system->n_atoms; ++i) {
            (*atom_ids)[i] = i;
        }
        return system->n_atoms;
    }

    // else, split match_string into individual elements
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
    for (size_t i = 0; i < system->n_atoms; ++i) {
        // and for each atom try matching the match function to individual match elements
        for (size_t j = 0; j < n_elements; ++j) {
            if (match_function(&(system->atoms[i]), elements[j])) {
                // reallocate array if needed
                if (selection_iterator >= alloc_ids) {
                    alloc_ids *= 2;
                    *atom_ids = realloc(*atom_ids, alloc_ids * sizeof(size_t));
                }

                (*atom_ids)[selection_iterator] = i;
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

int center_of_geometry(
        const system_t *system,
        const size_t n_selected,
        const size_t *atom_ids,
        vec_t center,
        const vec_t *coordinates)
{
    if (system == NULL || atom_ids == NULL || n_selected <= 0) return 1;

    // if no coordinates were provided, use coordinates from system
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

float distance2D(float *particle1, float *particle2)
{   
    float xd = particle1[0] - particle2[0];
    float yd = particle1[1] - particle2[1];
    return sqrt( xd*xd + yd*yd );
}