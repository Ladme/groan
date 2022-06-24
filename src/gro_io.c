// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#include "gro_io.h"

int isdecimal(const char *string) 
{
    short plusminus = 0;
    for (int i = 0; string[i] != '\0'; ++i) {
        if (!isdigit(string[i]) && !isblank(string[i])) {
            if (plusminus == 0 && (string[i] == '+' || string[i] == '-') ) ++plusminus;
            else return 0;
        }
    }

    return 1;
}

int isdecimalf(const char *string) 
{
    short point = 0;
    short plusminus = 0;
    for (int i = 0; string[i] != '\0'; ++i) {
        if (!isdigit(string[i]) && !isblank(string[i])) {
            if (point == 0 && string[i] == '.') ++point;
            else if (plusminus == 0 && (string[i] == '+' || string[i] == '-') ) ++plusminus;
            else return 0;
        }
    }

    return 1;
}

void get_fragment(const char *src, char *dest, const size_t start, const size_t len) 
{
    memcpy(dest, src + start, len);
    dest[len] = '\0';
}

int parse_int(const char *line, const size_t start, const size_t len, groint_t *element)
{
    char *segment = malloc(len + 1);
    get_fragment(line, segment, start, len);
    // check that only digits were loaded
    if (!isdecimal(segment)) {
        free(segment);
        return 1;
    }

    // TODO: check that the segment is not negative
    *element = (groint_t) atoi(segment);

    free(segment);
    return 0;
}

int parse_string(const char *line, const size_t start, const size_t len, char *element)
{
    char *segment = malloc(len + 1);
    get_fragment(line, segment, start, len);
    // trim string; this will only work if there is only one word
    sscanf(segment, "%s", element);

    free(segment);
    return 0;

}

int parse_float(const char *line, const size_t start, const size_t len, float *element)
{
    char *segment = malloc(len + 1);
    get_fragment(line, segment, start, len);
    // check that only digits or '.' were loaded
    if (!isdecimalf(segment)) {
        free(segment);
        return 1;
    }
    *element = atof(segment);
    
    free(segment);
    return 0;
}

int parse_gro_line(const char *line, atom_t *atom)
{
    // parse residue number
    if (parse_int(line, 0, 5, &(atom->residue_number)) != 0) return 1;

    // parse residue name
    if (parse_string(line, 5, 5, atom->residue_name) != 0) return 1;

    // parse atom name
    if (parse_string(line, 10, 5, atom->atom_name) != 0) return 1;

    // parse atom number
    if (parse_int(line, 15, 5, &(atom->atom_number)) != 0) return 1;

    // parse position
    for (short i = 0; i < 3; ++i) {
        if (parse_float(line, 20 + i * 8, 8, &(atom->position[i])) != 0) return 1;
    }

    // check whether velocity information is present
    short vel_present = 1;
    for (short i = 44; i < 68; ++i) {
        if (line[i] == '\0') vel_present = 0; 
    }

    // parse velocity
    if (vel_present) {
        for (short i = 0; i < 3; ++i) {
            if (parse_float(line, 44 + i * 8, 8, &(atom->velocity[i])) != 0) return 1;
        }
    }

    return 0;

}

system_t *load_gro(const char *filename)
{
    FILE *gro_file = fopen(filename, "r");
    if (gro_file == NULL) {
        fprintf(stderr, "File %s not found.\n", filename);
        return NULL;
    }

    char line[1024] = "";
    // read the first line of gro file and ignore it
    // TODO read optional time information from the first line
    if (fgets(line, 1024, gro_file) == NULL) {
        fclose(gro_file);
        fprintf(stderr, "Error. Gro file is empty.\n");
        return NULL;
    }

    // read the number of atoms
    if (fgets(line, 1024, gro_file) == NULL) {
        fclose(gro_file);
        fprintf(stderr, "Error. Could not read the number of atoms.\n");
        return NULL;
    }

    size_t n_atoms = 0;
    if (sscanf(line, "%lu", &n_atoms) != 1) {
        fclose(gro_file);
        fprintf(stderr, "Error. Could not read the number of atoms.\n");
        return NULL;
    }

    // allocate memory for the 'system' and set all bytes to zero
    system_t *system = malloc(sizeof(system_t) + n_atoms * sizeof(atom_t));
    if (system == NULL) {
        fclose(gro_file);
        fprintf(stderr, "Error. Could not allocate memory.\n");
        return NULL;
    }
    memset(system, 0, sizeof(system_t) + n_atoms * sizeof(atom_t));
    system->n_atoms = n_atoms;
    system->precision = 1000.;

    // read lines in gro file and load atoms
    for (size_t i = 0; i < n_atoms; ++i) {
        if (fgets(line, 1024, gro_file) == NULL) {
            free(system);
            fclose(gro_file);
            fprintf(stderr, "Error. Gro file ended unexpectedly.\n");
            return NULL;
        }

        if (parse_gro_line(line, &(system->atoms[i])) != 0) {
            free(system);
            fclose(gro_file);
            fprintf(stderr, "Error. Could not understand line: %s", line);
            return NULL;
        }

        // assign real gromacs atom number
        system->atoms[i].gmx_atom_number = i + 1;
    }

    // read information about the simulation box
    if (fgets(line, 1024, gro_file) == NULL) {
        free(system);
        fclose(gro_file);
        fprintf(stderr, "Error. Gro file is missing box information.\n");
        return NULL;
    }

    // TODO validate the input
    int loaded_values = sscanf(line, "%f %f %f %f %f %f %f %f %f",
    &(system->box[0]), &(system->box[1]), &(system->box[2]),
    &(system->box[3]), &(system->box[4]), &(system->box[5]), 
    &(system->box[6]), &(system->box[7]), &(system->box[8]));

    if (loaded_values != 9 && loaded_values != 3) {
        free(system);
        fclose(gro_file);
        fprintf(stderr, "Error. Could not obtain box size from: %s", line);
        return NULL;
    }

    fclose(gro_file);
    return system;

}

int write_gro(
        FILE *stream, 
        const atom_selection_t *atoms,
        const box_t boxsize,
        const write_mode_t write_mode,
        const char *comment)
{
    if (atoms == NULL || stream == NULL) return 1;

    fprintf(stream, "%s\n", comment);
    fprintf(stream, "%ld\n", atoms->n_atoms);

    // loop through atoms printing information about them
    for (size_t i = 0; i < atoms->n_atoms; ++i) {

        atom_t *atom = atoms->atoms[i];

        fprintf(stream, "%5u%-5s%5s%5u%8.3f%8.3f%8.3f",
        atom->residue_number, atom->residue_name, atom->atom_name, atom->atom_number,
        atom->position[0], atom->position[1], atom->position[2]);

        // print velocities if not forbidden
        if (write_mode != no_velocities) {
            fprintf(stream, "%8.4f%8.4f%8.4f", 
            atom->velocity[0], atom->velocity[1], atom->velocity[2]);
        }

        fprintf(stream, "\n");
    }

    // print box information
    for (int i = 0; i < 9; ++i) {
        fprintf(stream, " %9.5f", boxsize[i]);
    }
    fprintf(stream, "\n");

    return 0;
}