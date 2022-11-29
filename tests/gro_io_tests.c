// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#include "tests.h"

static void test_load_gro(void)
{
    printf("%-40s", "load_gro ");
    fflush(stdout);

    system_t *system = load_gro(INPUT_GRO_FILE);

    // check number of atoms
    assert(system->n_atoms == 48284);

    // check first atom in the list
    atom_t *atom = &system->atoms[0];

    assert(atom->residue_number == 1);
    assert(!strcmp(atom->residue_name, "LEU"));
    assert(!strcmp(atom->atom_name, "N"));
    assert(atom->atom_number == 1);
    assert(atom->gmx_atom_number == 1);
    assert(closef(atom->position[0], 5.028, 0.001));
    assert(closef(atom->position[1], 3.864, 0.001));
    assert(closef(atom->position[2], 6.219, 0.001));
    assert(closef(atom->velocity[0], -0.2376, 0.0001));
    assert(closef(atom->velocity[1],  0.2655, 0.0001));
    assert(closef(atom->velocity[2], -0.5587, 0.0001));

    // check randomly selected atom somewhere from the middle
    atom = &system->atoms[2134];
    assert(atom->residue_number == 36);
    assert(!strcmp(atom->residue_name, "POPE"));
    assert(!strcmp(atom->atom_name, "H8S"));
    assert(atom->atom_number == 2135);
    assert(atom->gmx_atom_number == 2135);
    assert(closef(atom->position[0], 2.367, 0.001));
    assert(closef(atom->position[1], 6.335, 0.001));
    assert(closef(atom->position[2], 5.366, 0.001));
    assert(closef(atom->velocity[0],  0.8083, 0.0001));
    assert(closef(atom->velocity[1], -2.0621, 0.0001));
    assert(closef(atom->velocity[2],  0.0627, 0.0001));
    
    // check the last atom
    atom = &system->atoms[48283];
    assert(atom->residue_number == 9207);
    assert(!strcmp(atom->residue_name, "NA"));
    assert(!strcmp(atom->atom_name, "NA"));
    assert(atom->atom_number == 48284);
    assert(atom->gmx_atom_number == 48284);
    assert(closef(atom->position[0], 1.593, 0.001));
    assert(closef(atom->position[1], 5.569, 0.001));
    assert(closef(atom->position[2], 6.361, 0.001));
    assert(closef(atom->velocity[0],  0.1499, 0.0001));
    assert(closef(atom->velocity[1], -0.0234, 0.0001));
    assert(closef(atom->velocity[2], -0.3537, 0.0001));

    free(system);
    printf("OK\n");
}

static void test_write_gro()
{
    printf("%-40s", "write_gro ");
    fflush(stdout);

    // full system
    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    FILE *output = fopen("temporary.gro", "w");
    write_gro(output, all, system->box, velocities, "Temporary gro file.");
    fclose(output);

    system_t *new_system = load_gro("temporary.gro");

    assert(system->n_atoms == new_system->n_atoms);
    assert(memcmp(system, new_system, sizeof(system_t) + system->n_atoms * sizeof(atom_t)) == 0);

    free(new_system);

    // part of the system
    select_t *part = select_atoms(all, "POPE", &match_residue_name);

    output = fopen("temporary.gro", "w");
    write_gro(output, part, system->box, velocities, "Temporary gro file.");
    fclose(output);

    new_system = load_gro("temporary.gro");

    assert(new_system->n_atoms == part->n_atoms);
    assert(new_system->box[0] == system->box[0]);
    assert(new_system->box[1] == system->box[1]);
    assert(new_system->box[2] == system->box[2]);

    for (size_t i = 0; i < part->n_atoms; ++i) {
        // temporarily change gmx atom number for `part` atoms so they can be compared with the new system
        size_t orig_gmx_atom_number = part->atoms[i]->gmx_atom_number;
        part->atoms[i]->gmx_atom_number = new_system->atoms[i].gmx_atom_number;
        assert(memcmp(part->atoms[i], &(new_system->atoms[i]), sizeof(atom_t)) == 0);
        part->atoms[i]->gmx_atom_number = orig_gmx_atom_number;
    }

    free(part);
    free(new_system);

    // empty system
    select_t *nothing = select_atoms(all, "XXXX", &match_residue_name);

    output = fopen("temporary.gro", "w");
    write_gro(output, nothing, system->box, velocities, "Temporary gro file.");
    fclose(output);

    new_system = load_gro("temporary.gro");
    
    assert(new_system->n_atoms == 0);
    assert(new_system->box[0] == system->box[0]);
    assert(new_system->box[1] == system->box[1]);
    assert(new_system->box[2] == system->box[2]);

    free(nothing);
    free(new_system);

    free(all);
    free(system);
    remove("temporary.gro");

    printf("OK\n");

}

void test_gro_io(void) 
{
    test_load_gro();
    test_write_gro();
}