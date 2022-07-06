// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#include "tests.h"

static void test_strsplit_space(void)
{
    printf("%-40s", "strsplit (space) ");

    char string[1024] = "This is an example string that should be split.";
    char **array = NULL;

    int n_split = strsplit(string, &array, " ");
    assert(n_split == 9);
    assert(!strcmp(array[0], "This"));
    assert(!strcmp(array[4], "string"));
    assert(!strcmp(array[8], "split."));

    free(array);
    printf("OK\n");
}

static void test_strsplit_spacetab(void)
{   
    printf("%-40s", "strsplit (spacetab) ");

    char string[1024] = "This is    an example string   that should be  split.";
    char **array = NULL;
    int n_split = strsplit(string, &array, " \t");
    assert(n_split == 9);
    assert(!strcmp(array[0], "This"));
    assert(!strcmp(array[4], "string"));
    assert(!strcmp(array[8], "split."));

    free(array);
    printf("OK\n");
}

static void test_load_gro(void)
{
    printf("%-40s", "load_gro ");

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

static void test_match_residue_name(void)
{
    printf("%-40s", "match_residue_name ");

    system_t *system = load_gro(INPUT_GRO_FILE);

    int atom_ids[5] = {12, 5061, 11349, 32542, 48191};
    char residue_names[5][5] = {"LEU", "POPE", "POPE", "SOL", "NA"};

    for (unsigned short i = 0; i < 5; ++i) {
        assert(match_residue_name(&system->atoms[atom_ids[i]], residue_names[i]));
        assert(match_residue_name(&system->atoms[atom_ids[i]], NULL));
    }

    free(system);
    printf("OK\n");
}

static void test_match_atom_name(void)
{
    printf("%-40s", "match_atom_name ");

    system_t *system = load_gro(INPUT_GRO_FILE);

    int atom_ids[5] = {12, 5061, 11349, 32542, 48191};
    char atom_names[5][5] = {"HD11", "H10X", "HB", "HW2", "NA"};

    for (unsigned short i = 0; i < 5; ++i) {
        assert(match_atom_name(&system->atoms[atom_ids[i]], atom_names[i]));
        assert(match_atom_name(&system->atoms[atom_ids[i]], NULL));
    }

    free(system);
    printf("OK\n");
}

static void test_match_atom_num(void)
{
    printf("%-40s", "match_atom_num ");

    system_t *system = load_gro(INPUT_GRO_FILE);

    int atom_ids[5] = {12, 5061, 11349, 32542, 48191};
    char gmx_numbers[5][6] = {"13", "5062", "11350", "32543", "48192"};

    for (unsigned short i = 0; i < 5; ++i) {
        assert(match_atom_num(&system->atoms[atom_ids[i]], gmx_numbers[i]));
        assert(match_atom_num(&system->atoms[atom_ids[i]], NULL));
    }

    free(system);
    printf("OK\n");
}

static void test_selection_create(void)
{
    printf("%-40s", "selection_create ");

    // we allocate enough memory to hold 10 atom pointers
    select_t *selection = selection_create(10);
    assert(selection->n_atoms == 0);

    system_t *system = load_gro(INPUT_GRO_FILE);

    // we assign a random atom pointer to the last allocated position in the selection
    // to test the allocation
    selection->atoms[9] = &system->atoms[10432];
    assert(selection->atoms[9]->atom_number == 10433);

    free(selection);
    free(system);
    printf("OK\n");
}

static void test_select_system(void)
{
    printf("%-40s", "select_system ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *selection = select_system(system);

    assert(selection->n_atoms == 48284);
    assert(selection->atoms[0] == &system->atoms[0]);
    assert(selection->atoms[48283] == &system->atoms[48283]);

    free(selection);
    free(system);
    printf("OK\n");
}

static void test_selection_copy(void)
{
    printf("%-40s", "selection_copy ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *selection = select_system(system);

    select_t *selection2 = selection_copy(selection);

    assert(selection->n_atoms == selection2->n_atoms);
    assert(selection->atoms[0] == selection2->atoms[0]);
    assert(selection->atoms[48283] == selection2->atoms[48283]);

    free(selection);
    free(selection2);
    free(system);
    printf("OK\n");
}

static void test_selection_empty(void)
{
    printf("%-40s", "selection_empty ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *selection = select_system(system);

    selection_empty(selection);

    assert(selection->n_atoms == 0);

    free(selection);
    free(system);
    printf("OK\n");
}

static void test_selection_add_atom(void)
{
    printf("%-40s", "selection_add_atom ");

    system_t *system = load_gro(INPUT_GRO_FILE);

    size_t allocated = 1;
    select_t *selection = selection_create(allocated);

    selection_add_atom(&selection, &allocated, &system->atoms[2134]);
    assert(!strcmp(selection->atoms[0]->atom_name, "H8S"));

    selection_add_atom(&selection, &allocated, &system->atoms[0]);
    assert(!strcmp(selection->atoms[1]->residue_name, "LEU"));

    assert(selection->n_atoms == 2);
    assert(allocated > 1);

    free(selection);
    free(system);
    printf("OK\n");
}

static void test_select_atoms_atomname(void)
{
    printf("%-40s", "select_atoms (atomname) ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *selection = select_atoms(all, "CA", &match_atom_name);

    assert(selection->n_atoms == 21);
    assert(selection->atoms[0] == &system->atoms[4]);
    assert(selection->atoms[20] == &system->atoms[312]);

    free(all);
    free(selection);
    free(system);
    printf("OK\n");
}

static void test_select_atoms_atomname_multiple(void)
{
    printf("%-40s", "select_atoms (atomname, multiple) ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *selection = select_atoms(all, "CA H1 HD21 HD22", &match_atom_name);

    assert(selection->n_atoms == 47);
    assert(selection->atoms[0] == &system->atoms[1]);
    assert(selection->atoms[46] == &system->atoms[330]);

    free(all);
    free(selection);
    free(system);
    printf("OK\n");
}

static void test_select_atoms_atomname_empty(void)
{
    printf("%-40s", "select_atoms (atomname, empty) ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *selection = select_atoms(all, "   ", &match_atom_name);

    assert(selection->n_atoms == 0);

    free(all);
    free(selection);
    free(system);
    printf("OK\n");
}

static void test_select_atoms_atomname_nomatch(void)
{
    printf("%-40s", "select_atoms (atomname, nomatch) ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *selection = select_atoms(all, "XGR", &match_atom_name);

    assert(selection->n_atoms == 0);

    free(all);
    free(selection);
    free(system);
    printf("OK\n");
}

static void test_select_atoms_resname(void)
{
    printf("%-40s", "select_atoms (resname) ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *selection = select_atoms(all, "NA", &match_residue_name);

    assert(selection->n_atoms == 67);
    assert(selection->atoms[0] == &system->atoms[48191]);
    assert(selection->atoms[66] == &system->atoms[48283]);

    free(all);
    free(selection);
    free(system);
    printf("OK\n");
}

static void test_select_atoms_resname_multiple(void)
{
    printf("%-40s", "select_atoms (resname, multiple) ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *selection = select_atoms(all, "NA POPG", &match_residue_name);

    assert(selection->n_atoms == 5401);
    assert(selection->atoms[0] == &system->atoms[16082]);
    assert(selection->atoms[5400] == &system->atoms[48283]);

    free(all);
    free(selection);
    free(system);
    printf("OK\n");
}

static void test_select_atoms_resname_empty(void)
{
    printf("%-40s", "select_atoms (resname, empty) ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *selection = select_atoms(all, " ", &match_residue_name);

    assert(selection->n_atoms == 0);

    free(all);
    free(selection);
    free(system);
    printf("OK\n");
}

static void test_select_atoms_resname_nomatch(void)
{
    printf("%-40s", "select_atoms (resname, nomatch) ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *selection = select_atoms(all, "GLU VAL ASN", &match_residue_name);

    assert(selection->n_atoms == 0);

    free(all);
    free(selection);
    free(system);
    printf("OK\n");
}

static void test_select_atoms_atomnum(void)
{
    printf("%-40s", "select_atoms (atomnum) ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *selection = select_atoms(all, "1 7 9 1465 43256", &match_atom_num);

    assert(selection->n_atoms == 5);
    assert(selection->atoms[0] == &system->atoms[0]);
    assert(selection->atoms[4] == &system->atoms[43255]);

    free(all);
    free(selection);
    free(system);
    printf("OK\n");
}

static void test_select_atoms_atomnum_nomatch(void)
{
    printf("%-40s", "select_atoms (atomnum, nomatch) ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *selection = select_atoms(all, "-864 9897674", &match_atom_num);

    assert(selection->n_atoms == 0);

    free(all);
    free(selection);
    free(system);
    printf("OK\n");
}

static void test_selection_cat(void)
{
    printf("%-40s", "selection_cat ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *selection1 = select_atoms(all, "LEU", &match_residue_name);

    select_t *selection2 = select_atoms_d(all, "CA", &match_atom_name);

    select_t *cat = selection_cat(selection1, selection2);

    assert(selection1->n_atoms + selection2->n_atoms == cat->n_atoms);
    for (size_t i = 0; i < selection1->n_atoms; ++i) {
        assert(selection1->atoms[i] == cat->atoms[i]);
    }
    for (size_t i = 0; i < selection2->n_atoms; ++i) {
        assert(selection2->atoms[i] == cat->atoms[selection1->n_atoms + i]);
    }

    free(selection1);
    free(selection2);
    free(cat);
    free(system);
    printf("OK\n");
}

static void test_selection_cat_unique(void)
{
    printf("%-40s", "selection_cat_unique ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *selection1 = select_atoms(all, "LEU", &match_residue_name);

    select_t *selection2 = select_atoms_d(all, "CA", &match_atom_name);

    select_t *cat = selection_cat_unique(selection1, selection2);

    assert(selection1->n_atoms + selection2->n_atoms - 12 == cat->n_atoms);
    for (size_t i = 0; i < selection1->n_atoms; ++i) {
        assert(selection1->atoms[i] == cat->atoms[i]);
    }
    
    for (size_t i = 0; i < cat->n_atoms; ++i) {
        for (size_t j = i + 1; j < cat->n_atoms; ++j) {
            assert(cat->atoms[i] != cat->atoms[j]);
        }
    }

    free(selection1);
    free(selection2);
    free(cat);
    free(system);
    printf("OK\n");
}

static void test_selection_intersect(void)
{
    printf("%-40s", "selection_intersect ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *selection1 = select_atoms(all, "LEU", &match_residue_name);

    select_t *selection2 = select_atoms_d(all, "CA", &match_atom_name);

    select_t *intersect = selection_intersect_d(selection1, selection2);

    assert(intersect->n_atoms == 12);
    for (size_t i = 0; i < intersect->n_atoms; ++i) {
        assert(match_atom_name(intersect->atoms[i], "CA"));
    }

    free(intersect);
    free(system);
    printf("OK\n");
}

static void test_selection_intersect_none(void)
{
    printf("%-40s", "selection_intersect (none) ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *selection1 = select_atoms(all, "LEU", &match_residue_name);

    select_t *selection2 = select_atoms_d(all, "HW1 HW2", &match_atom_name);

    select_t *intersect = selection_intersect_d(selection1, selection2);

    assert(intersect->n_atoms == 0);

    free(intersect);
    free(system);
    printf("OK\n");
}

static void test_selection_remove_atom(void)
{
    printf("%-40s", "selection_remove_atom ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *selection = select_atoms_d(all, "LEU", &match_residue_name);

    size_t orig_length = selection->n_atoms;

    size_t removed = selection_remove_atom(selection, &system->atoms[43]);

    assert(removed == 1);
    assert(selection->n_atoms == orig_length - removed);
    assert(selection->atoms[21]->gmx_atom_number == 45);

    free(selection);
    free(system);
    printf("OK\n");
}

static void test_selection_remove_atom_none(void)
{
    printf("%-40s", "selection_remove_atom (none) ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *selection = select_atoms_d(all, "LEU", &match_residue_name);

    size_t orig_length = selection->n_atoms;

    size_t removed = selection_remove_atom(selection, &system->atoms[25]);

    assert(removed == 0);
    assert(selection->n_atoms == orig_length);

    free(selection);
    free(system);
    printf("OK\n");
}

static void test_selection_remove_atom_all(void)
{
    printf("%-40s", "selection_remove_atom (all) ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *selection = select_atoms_d(all, "17564", &match_atom_num);

    size_t orig_length = selection->n_atoms;

    size_t removed = selection_remove_atom(selection, &system->atoms[17563]);

    assert(removed == orig_length);
    assert(selection->n_atoms == 0);

    free(selection);
    free(system);
    printf("OK\n");
}

static void test_selection_remove_atom_duplicates(void)
{
    printf("%-40s", "selection_remove_atom (duplicates) ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *selection = select_atoms_d(all, "LEU", &match_residue_name);

    select_t *selection_double = selection_cat_d(selection, selection);
    select_t *selection_quad = selection_cat_d(selection_double, selection_double);

    size_t orig_length = selection_quad->n_atoms;

    size_t removed = selection_remove_atom(selection_quad, &system->atoms[43]);

    assert(removed == 4);
    assert(selection_quad->n_atoms == orig_length - removed);

    free(selection_quad);
    free(system);
    printf("OK\n");
}

static void test_selection_remove(void)
{
    printf("%-40s", "selection_remove ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *selection1 = select_atoms(all, "LEU", &match_residue_name);

    select_t *selection2 = select_atoms_d(all, "CA", &match_atom_name);

    size_t orig_length = selection1->n_atoms;

    size_t removed = selection_remove(selection1, selection2);

    assert(removed == 12);
    assert(selection1->n_atoms == orig_length - removed);

    free(selection1);
    free(selection2);
    free(system);
    printf("OK\n");
}

static void test_selection_remove_none(void)
{
    printf("%-40s", "selection_remove (none) ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *selection1 = select_atoms(all, "LEU", &match_residue_name);

    select_t *selection2 = select_atoms_d(all, "OW", &match_atom_name);

    size_t orig_length = selection1->n_atoms;

    size_t removed = selection_remove(selection1, selection2);

    assert(removed == 0);
    assert(selection1->n_atoms == orig_length);

    free(selection1);
    free(selection2);
    free(system);
    printf("OK\n");
}

static void test_selection_remove_all(void)
{
    printf("%-40s", "selection_remove (all) ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *selection1 = select_atoms(all, "LEU", &match_residue_name);

    size_t orig_length = selection1->n_atoms;

    size_t removed = selection_remove(selection1, selection1);

    assert(removed == orig_length);
    assert(selection1->n_atoms == 0);

    free(selection1);
    free(all);
    free(system);
    printf("OK\n");
}

static void test_selection_remove_all2(void)
{
    printf("%-40s", "selection_remove (all2) ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *selection1 = select_atoms(all, "LEU", &match_residue_name);
    select_t *selection2 = select_atoms(all, "LEU", &match_residue_name);

    size_t orig_length = selection1->n_atoms;

    size_t removed = selection_remove(selection1, selection2);

    assert(removed == orig_length);
    assert(selection1->n_atoms == 0);

    free(selection1);
    free(selection2);
    free(all);
    free(system);
    printf("OK\n");
}

static void test_selection_compare(void)
{
    printf("%-40s", "selection_compare");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *selection1 = select_atoms(all, "CA", &match_atom_name);
    select_t *selection2 = select_atoms(all, "N", &match_atom_name);
    select_t *selection = selection_cat(selection1, selection2);

    select_t *selection3 = select_atoms(all, "N CA", &match_atom_name);

    assert(selection_compare(selection, selection3));
    assert(selection_compare(selection3, selection));
    assert(!selection_compare(selection1, selection2));
    assert(!selection_compare(selection1, selection3));
    assert(!selection_compare(selection1, selection));
    assert(!selection_compare(selection2, selection3));
    assert(!selection_compare(selection2, selection));

    free(selection);
    free(selection1);
    free(selection2);
    free(selection3);
    free(all);
    free(system);
    printf("OK\n");
}

static void test_selection_compare_empty(void)
{
    printf("%-40s", "selection_compare (empty)");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *selection1 = select_atoms(all, "CXJH", &match_atom_name);
    select_t *selection2 = select_atoms(all, "N", &match_atom_name);
    select_t *selection3 = select_atoms(all, "IFB", &match_atom_name);

    assert(!selection_compare(selection1, selection2));
    assert(selection_compare(selection1, selection3));

    free(selection1);
    free(selection2);
    free(selection3);
    free(all);
    free(system);
    printf("OK\n");
}

static void test_selection_compare_strict(void)
{
    printf("%-40s", "selection_compare_strict");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *selection1 = select_atoms(all, "CA", &match_atom_name);
    select_t *selection2 = select_atoms(all, "N", &match_atom_name);
    select_t *selection = selection_cat(selection1, selection2);

    select_t *selection3 = select_atoms(all, "N CA", &match_atom_name);

    select_t *selection_alt = selection_cat(selection2, selection1);

    select_t *selection4 = select_atoms(all, "CA N", &match_atom_name);


    assert(!selection_compare_strict(selection, selection3));
    assert(!selection_compare_strict(selection1, selection2));
    assert(!selection_compare_strict(selection1, selection3));
    assert(!selection_compare_strict(selection1, selection));
    assert(!selection_compare_strict(selection2, selection3));
    assert(!selection_compare_strict(selection2, selection));
    assert(!selection_compare_strict(selection_alt, selection3));
    assert(selection_compare_strict(selection3, selection4));
    assert(selection_compare_strict(selection4, selection3));

    free(selection);
    free(selection1);
    free(selection2);
    free(selection3);
    free(selection_alt);
    free(selection4);
    free(all);
    free(system);
    printf("OK\n");
}

static void test_selection_compare_strict_empty(void)
{
    printf("%-40s", "selection_compare_strict (empty)");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *selection1 = select_atoms(all, "CXJH", &match_atom_name);
    select_t *selection2 = select_atoms(all, "N", &match_atom_name);
    select_t *selection3 = select_atoms(all, "IFB", &match_atom_name);

    assert(!selection_compare_strict(selection1, selection2));
    assert(selection_compare_strict(selection1, selection3));

    free(selection1);
    free(selection2);
    free(selection3);
    free(all);
    free(system);
    printf("OK\n");
}

static void test_selection_unique(void)
{
    printf("%-40s", "selection_unique ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *selection = select_atoms_d(all, "LEU", &match_residue_name);
    size_t orig_length = selection->n_atoms;

    select_t *selection_double = selection_cat_d(selection, selection);
    select_t *selection_quad = selection_cat_d(selection_double, selection_double);

    assert(selection_unique(selection_quad) == orig_length * 3);
    assert(selection_quad->n_atoms == orig_length);

    free(selection_quad);
    free(system);
    printf("OK\n");
}

static void test_selection_renumber(void)
{
    printf("%-40s", "selection_renumber ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *selection1 = select_atoms(all, "SER", &match_residue_name);
    select_t *selection2 = select_atoms_d(all, "LEU", &match_residue_name);

    select_t *selection = selection_cat_d(selection1, selection2);
    selection_renumber(selection);

    assert(selection->atoms[0]->atom_number == 1);
    assert(selection->atoms[0]->residue_number == 1);
    assert(selection->atoms[0]->gmx_atom_number == 22);
    assert(system->atoms[21].atom_number == 1);
    assert(system->atoms[21].residue_number == 1);
    assert(system->atoms[21].gmx_atom_number == 22);

    assert(selection->atoms[21]->residue_number == 2);

    free(selection);
    free(system);
    printf("OK\n");
}

static void test_selection_sort(void)
{
    printf("%-40s", "selection_sort ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *selection1 = select_atoms(all, "SER", &match_residue_name);
    select_t *selection2 = select_atoms_d(all, "LEU", &match_residue_name);

    select_t *selection = selection_cat_d(selection1, selection2);
    selection_sort(selection);

    assert(selection->atoms[0]->atom_number == 1);
    assert(selection->atoms[0]->residue_number == 1);
    assert(selection->atoms[0]->gmx_atom_number == 1);
    assert(!strcmp(selection->atoms[0]->residue_name, "LEU"));
    

    free(selection);
    free(system);
    printf("OK\n");
}

static void test_selection_sort_renumber(void)
{
    printf("%-40s", "selection_sort (renumbered) ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *selection1 = select_atoms(all, "SER", &match_residue_name);
    select_t *selection2 = select_atoms(all, "LEU", &match_residue_name);

    select_t *selection = selection_cat_d(selection1, selection2);
    selection_renumber(selection);
    free(selection);
    selection = NULL;

    select_t *selection3 = select_atoms(all, "LEU", &match_residue_name);
    select_t *selection4 = select_atoms_d(all, "SER", &match_residue_name);

    select_t *selection_to_sort = selection_cat_d(selection3, selection4);
    selection_sort(selection_to_sort);

    assert(selection_to_sort->atoms[0]->atom_number == 1);
    assert(selection_to_sort->atoms[0]->residue_number == 1);
    assert(selection_to_sort->atoms[0]->gmx_atom_number == 22);
    assert(!strcmp(selection_to_sort->atoms[0]->residue_name, "SER"));

    assert(selection_to_sort->atoms[21]->residue_number == 2);

    free(selection_to_sort);
    free(system);
    printf("OK\n");
}

static void test_selection_sort_gmx(void)
{
    printf("%-40s", "selection_sort_gmx ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *selection1 = select_atoms(all, "SER", &match_residue_name);
    select_t *selection2 = select_atoms_d(all, "LEU", &match_residue_name);

    select_t *selection = selection_cat_d(selection1, selection2);
    selection_sort_gmx(selection);

    assert(selection->atoms[0]->atom_number == 1);
    assert(selection->atoms[0]->residue_number == 1);
    assert(selection->atoms[0]->gmx_atom_number == 1);
    assert(!strcmp(selection->atoms[0]->residue_name, "LEU"));
    

    free(selection);
    free(system);
    printf("OK\n");
}

static void test_selection_sort_gmx_renumber(void)
{
    printf("%-40s", "selection_sort_gmx (renumbered) ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *selection1 = select_atoms(all, "SER", &match_residue_name);
    select_t *selection2 = select_atoms(all, "LEU", &match_residue_name);

    select_t *selection = selection_cat_d(selection1, selection2);
    selection_renumber(selection);
    free(selection);
    selection = NULL;

    select_t *selection3 = select_atoms(all, "LEU", &match_residue_name);
    select_t *selection4 = select_atoms_d(all, "SER", &match_residue_name);

    select_t *selection_to_sort = selection_cat_d(selection3, selection4);
    selection_sort_gmx(selection_to_sort);
    
    assert(selection_to_sort->atoms[0]->atom_number == 100);
    assert(selection_to_sort->atoms[0]->residue_number == 10);
    assert(selection_to_sort->atoms[0]->gmx_atom_number == 1);
    assert(!strcmp(selection_to_sort->atoms[0]->residue_name, "LEU"));

    assert(selection_to_sort->atoms[21]->residue_number == 1);

    free(selection_to_sort);
    free(system);
    printf("OK\n");
}

static void test_selection_fixres(void)
{
    printf("%-40s", "selection_fixres ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *selection1 = select_atoms(all, "HD21 HD22 HD23", &match_atom_name);
    select_t *selection2 = select_atoms(all, "CD2 C", &match_atom_name);
    select_t *selection3 = select_atoms_d(all, "N", &match_atom_name);

    select_t *selection12 = selection_cat_d(selection1, selection2);
    select_t *selection123 = selection_cat_d(selection12, selection3);

    selection_fixres(selection123);

    assert(selection123->atoms[1]->atom_number == 16);
    assert(selection123->atoms[1]->residue_number == 1);
    assert(!strcmp(selection123->atoms[1]->atom_name, "CD2"));

    assert(selection123->atoms[73]->atom_number == 22);
    assert(selection123->atoms[73]->residue_number == 2);
    assert(!strcmp(selection123->atoms[73]->atom_name, "N"));
    assert(selection123->atoms[selection123->n_atoms - 1] == &system->atoms[15957]);
    

    free(selection123);
    free(system);
    printf("OK\n");
}

static void test_selection_isin(void)
{
    printf("%-40s", "selection_isin ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    for (size_t i = 0; i < system->n_atoms; ++i) {
        assert(selection_isin(all, &system->atoms[i]));
    }

    select_t *selectionLEU = select_atoms(all, "LEU", &match_residue_name);
    select_t *selection2 = select_atoms(all, "HD21 HD22 HD23", &match_atom_name);

    for (size_t i = 0; i < selection2->n_atoms; ++i) {
        assert(selection_isin(selectionLEU, selection2->atoms[i]));
    }

    select_t *selectionSER = select_atoms(all, "SER", &match_residue_name);
    
    for (size_t i = 0; i < selection2->n_atoms; ++i) {
        assert(!selection_isin(selectionSER, selection2->atoms[i]));
    }

    free(all);
    free(selectionLEU);
    free(selection2);
    free(selectionSER);
    free(system);
    printf("OK\n");
}

static void test_selection_to_system(void)
{
    printf("%-40s", "selection_to_system ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *selection1 = select_atoms(all, "HD21 HD22 HD23", &match_atom_name);
    select_t *selection2 = select_atoms(all, "CD2 C", &match_atom_name);
    select_t *selection3 = select_atoms_d(all, "N", &match_atom_name);

    select_t *selection12 = selection_cat_d(selection1, selection2);
    select_t *selection123 = selection_cat_d(selection12, selection3);

    system_t *new_system = selection_to_system(selection123, system->box, 0, 10.0f);

    // asserting the properties of the new system
    assert(new_system->n_atoms == selection123->n_atoms);
    assert(new_system->box[0] == system->box[0]);
    assert(new_system->box[1] == system->box[1]);
    assert(new_system->box[2] == system->box[2]);
    assert(new_system->step == 0);
    assert(closef(new_system->time, 10.0, 0.00001));

    assert(new_system->atoms[1].atom_number == 2);
    assert(new_system->atoms[1].residue_number == 1);
    assert(!strcmp(new_system->atoms[1].atom_name, "CD2"));

    assert(new_system->atoms[73].atom_number == 74);
    assert(new_system->atoms[73].residue_number == 13);
    assert(!strcmp(new_system->atoms[73].atom_name, "N"));

    // asserting the properties of the old system
    assert(selection123->atoms[1]->atom_number == 18);
    assert(selection123->atoms[1]->residue_number == 1);
    assert(!strcmp(selection123->atoms[1]->atom_name, "HD22"));

    assert(selection123->atoms[38]->atom_number == 31);
    assert(selection123->atoms[38]->residue_number == 2);
    assert(!strcmp(selection123->atoms[38]->atom_name, "C"));
    assert(selection123->atoms[selection123->n_atoms - 1] == &system->atoms[15957]);

    free(selection123);
    free(new_system);
    free(system);
    printf("OK\n");
}

static void test_select_geometry_sphere(void)
{
    printf("%-40s", "select_geometry (sphere) ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    vec_t center = {0.f};
    float g_definition = 4.0f;
    select_t *selection = select_geometry(all, center, sphere, &g_definition, system->box);

    assert(selection->n_atoms == 26253);
    assert(selection->atoms[0] == &system->atoms[582]);
    assert(selection->atoms[26252] == &system->atoms[48283]);

    free(all);
    free(selection);
    free(system);
    printf("OK\n");
}

static void test_select_geometry_box(void)
{
    printf("%-40s", "select_geometry (box) ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    vec_t center = {0.f};
    float g_definition[6] = {-2.5, 1.0, 0.0, 4.5, -0.5, 3.3};
    select_t *selection = select_geometry(all, center, box, g_definition, system->box);

    assert(selection->n_atoms == 4878);
    assert(selection->atoms[0] == &system->atoms[8361]);
    assert(selection->atoms[4877] == &system->atoms[48273]);

    free(all);
    free(selection);
    free(system);
    printf("OK\n");
}

static void test_select_geometry_zcylinder(void)
{
    printf("%-40s", "select_geometry (zcylinder) ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    vec_t center = {0.f};
    float g_definition[3] = {3.3, -2.1, 1.3};
    select_t *selection = select_geometry(all, center, zcylinder, g_definition, system->box);

    assert(selection->n_atoms == 11384);
    assert(selection->atoms[0] == &system->atoms[584]);
    assert(selection->atoms[11383] == &system->atoms[48282]);

    free(all);
    free(selection);
    free(system);
    printf("OK\n");
}

static void test_select_geometry_ycylinder(void)
{
    printf("%-40s", "select_geometry (ycylinder) ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    vec_t center = {0.f};
    float g_definition[3] = {3.3, -2.1, 1.3};
    select_t *selection = select_geometry(all, center, ycylinder, g_definition, system->box);

    assert(selection->n_atoms == 11704);
    assert(selection->atoms[0] == &system->atoms[2582]);
    assert(selection->atoms[11703] == &system->atoms[48283]);

    free(all);
    free(selection);
    free(system);
    printf("OK\n");
}

static void test_select_geometry_xcylinder(void)
{
    printf("%-40s", "select_geometry (xcylinder) ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    vec_t center = {0.f};
    float g_definition[3] = {3.3, -2.1, 1.3};
    select_t *selection = select_geometry(all, center, xcylinder, g_definition, system->box);

    assert(selection->n_atoms == 11622);
    assert(selection->atoms[0] == &system->atoms[584]);
    assert(selection->atoms[11621] == &system->atoms[48281]);

    free(all);
    free(selection);
    free(system);
    printf("OK\n");
}

static void test_read_ndx(void)
{
    printf("%-40s", "read_ndx (basic) ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    dict_t *ndx_groups = read_ndx(NDX_FILE, system);

    char groups[25][50] = {"System", "Protein", "Protein-H", "C-alpha",
                           "Backbone", "MainChain", "MainChain+Cb", "MainChain+H",
                           "SideChain", "SideChain-H", "Prot-Masses", "non-Protein",
                           "Other", "POPE", "POPG", "NA", "CL", "Water", "SOL",
                           "non-Water", "Membrane", "ION", "W_ION", "Protein_Membrane", "Empty"};
    
    size_t numbers[25] = {48284, 332, 151, 21,
                          64, 85, 106, 110,
                          222, 66, 332, 47952,
                          21084, 15750, 5334, 67, 26, 26775, 26775,
                          21509, 21084, 93, 26868, 21416, 0};


    for (int i = 0; i < 25; ++i) {
        select_t *selection = (atom_selection_t *) dict_get(ndx_groups, groups[i]);
        assert(selection != NULL);
        //printf("Selection %s: %lu\n", groups[i], selection->n_atoms);
        assert(selection->n_atoms == numbers[i]);
    }

    select_t *selection = (atom_selection_t *) dict_get(ndx_groups, "NonExistent");
    assert(selection == NULL);

    free(system);
    dict_destroy(ndx_groups);
    printf("OK\n");
}

static void test_read_ndx_advanced(void)
{
    printf("%-40s", "read_ndx (advanced) ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);
    dict_t *ndx_groups = read_ndx(NDX_FILE, system);

    // check protein
    select_t *protein = select_atoms(all, "LEU SER NHE", &match_residue_name);
    select_t *ndx_protein = (atom_selection_t *) dict_get(ndx_groups, "Protein");
    assert(ndx_protein != NULL);
    assert(selection_compare(protein, ndx_protein));
    
    // check protein without hydrogen atoms
    select_t *hydrogens = select_atoms(all, "H1 H2 H3 HA HB1 HB2 HG HD11 HD12 HD13 HD21 HD22 HD23 H", &match_atom_name);
    selection_remove(protein, hydrogens);
    select_t *ndx_protein_minus_h = (atom_selection_t *) dict_get(ndx_groups, "Protein-H");
    assert(ndx_protein_minus_h != NULL);
    assert(selection_compare(protein, ndx_protein_minus_h));
    free(protein);
    free(hydrogens);

    // check membrane
    select_t *pope = select_atoms(all, "POPE", &match_residue_name);
    select_t *popg = select_atoms(all, "POPG", &match_residue_name);
    select_t *ndx_pope = (atom_selection_t *) dict_get(ndx_groups, "POPE");
    select_t *ndx_popg = (atom_selection_t *) dict_get(ndx_groups, "POPG");
    assert(ndx_pope != NULL);
    assert(ndx_popg != NULL);
    assert(selection_compare(pope, ndx_pope));
    assert(selection_compare(popg, ndx_popg));
    select_t *membrane = selection_cat_d(pope, popg);
    select_t *ndx_membrane = (atom_selection_t *) dict_get(ndx_groups, "Membrane");
    assert(ndx_membrane != NULL);
    assert(selection_compare(membrane, ndx_membrane));
    free(membrane);

    // check water and non-water
    select_t *water = select_atoms(all, "SOL", &match_residue_name);
    select_t *ndx_water1 = (atom_selection_t *) dict_get(ndx_groups, "Water");
    select_t *ndx_water2 = (atom_selection_t *) dict_get(ndx_groups, "SOL");
    assert(ndx_water1 != NULL);
    assert(ndx_water2 != NULL);
    assert(selection_compare(water, ndx_water1));
    assert(selection_compare(water, ndx_water2));
    select_t *nonwater = selection_copy(all);
    selection_remove(nonwater, water);
    select_t *ndx_nonwater = (atom_selection_t *) dict_get(ndx_groups, "non-Water");
    assert(ndx_nonwater != NULL);
    assert(selection_compare(nonwater, ndx_nonwater));
    free(water);
    free(nonwater);

    // check ions
    select_t *na = select_atoms(all, "NA", &match_atom_name);
    select_t *cl = select_atoms(all, "CL", &match_residue_name);
    select_t *ndx_na = (atom_selection_t *) dict_get(ndx_groups, "NA");
    select_t *ndx_cl = (atom_selection_t *) dict_get(ndx_groups, "CL");
    assert(ndx_na != NULL);
    assert(ndx_cl != NULL);
    assert(selection_compare(na, ndx_na));
    assert(selection_compare(cl, ndx_cl));
    select_t *ion = selection_cat_d(na, cl);
    select_t *ndx_ion = (atom_selection_t *) dict_get(ndx_groups, "ION");
    assert(ndx_ion != NULL);
    assert(selection_compare(ion, ndx_ion));
    free(ion);

    free(system);
    free(all);
    dict_destroy(ndx_groups);
    printf("OK\n");
}

static void test_read_ndx_empty(void)
{
    printf("%-40s", "read_ndx (empty) ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    dict_t *ndx_groups = read_ndx(EMPTY_NDX_FILE, system);

    char groups[25][50] = {"System", "Protein", "Protein-H", "C-alpha",
                           "Backbone", "MainChain", "MainChain+Cb", "MainChain+H",
                           "SideChain", "SideChain-H", "Prot-Masses", "non-Protein",
                           "Other", "POPE", "POPG", "NA", "CL", "Water", "SOL",
                           "non-Water", "Membrane", "ION", "W_ION", "Protein_Membrane", "Empty"};
    
    for (int i = 0; i < 25; ++i) {
        select_t *selection = (atom_selection_t *) dict_get(ndx_groups, groups[i]);
        assert(selection == NULL);
    }

    free(system);
    dict_destroy(ndx_groups);
    printf("OK\n");
}

static void test_read_ndx_nonexistent(void)
{
    printf("%-40s", "read_ndx (non-existent) ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    dict_t *ndx_groups = read_ndx("index.ndx", system);

    assert(ndx_groups == NULL);

    free(system);
    printf("OK\n");
}

static void test_smart_select(void)
{
    printf("%-40s", "smart_select ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);
    dict_t *ndx_groups = read_ndx(NDX_FILE, system);

    assert(ndx_groups != NULL);

    // residue name selection
    select_t *manual_pope = select_atoms(all, "POPE", &match_residue_name);
    select_t *pope1 = smart_select(all, "resname POPE", ndx_groups);
    select_t *pope2 = smart_select(all, "resnamePOPE", NULL);
    select_t *pope3 = smart_select(all, "resname      POPE", ndx_groups);
    assert(pope1 != NULL);
    assert(pope2 != NULL);
    assert(pope3 != NULL);
    assert(selection_compare_strict(manual_pope, pope1));
    assert(selection_compare_strict(manual_pope, pope2));
    assert(selection_compare_strict(manual_pope, pope3));
    free(pope1);
    free(pope2);
    free(pope3);
    free(manual_pope);

    // atom name selection
    select_t *manual_hydrogens = select_atoms(all, "H1 H2 H3 HA HB1 HB2 HG HD11 HD12 HD13 HD21 HD22 HD23 H", &match_atom_name);
    select_t *hydrogens = smart_select(all, "name H1 H2 H3 HA HB1 HB2 HG HD11 HD12 HD13 HD21 HD22 HD23 H", NULL);
    assert(hydrogens != NULL);
    assert(selection_compare_strict(manual_hydrogens, hydrogens));
    free(manual_hydrogens);
    free(hydrogens);

    select_t *manual_o = select_atoms(all, "O", &match_atom_name);
    select_t *o = smart_select(all, "nameO", ndx_groups);
    assert(o != NULL);
    assert(selection_compare_strict(manual_o, o));
    free(manual_o);
    free(o);

    // atom number selection
    select_t *manual_serial = select_atoms(all, "73 6542 9875 23463 2 42653", &match_atom_num);
    select_t *serial = smart_select(all, "serial 73 6542 9875 23463 2 42653", ndx_groups);
    select_t *serial2 = smart_select(all, "serial 73 2 6542 9875 42653 23463", ndx_groups);
    select_t *serial3 = smart_select(all, "serial 73 73 73 73 73 2 6542 9875 42653 23463", NULL);
    select_t *serial4 = smart_select(all, "serial 76532 73 6542 9875 23463 2 987654 42653", ndx_groups);
    assert(serial != NULL);
    assert(serial2 != NULL);
    assert(serial3 != NULL);
    assert(serial4 != NULL);
    assert(selection_compare_strict(manual_serial, serial));
    assert(selection_compare_strict(manual_serial, serial2));
    assert(selection_compare_strict(manual_serial, serial3));
    assert(selection_compare_strict(manual_serial, serial4));
    free(manual_serial);
    free(serial);
    free(serial2);
    free(serial3);
    free(serial4);

    // index selection
    select_t *manual_sidechain = (atom_selection_t *) dict_get(ndx_groups, "SideChain");
    select_t *sidechain = smart_select(all, "SideChain", ndx_groups);
    assert(sidechain != NULL);
    assert(selection_compare_strict(manual_sidechain, sidechain));
    free(sidechain);

    select_t *manual_wion = (atom_selection_t *) dict_get(ndx_groups, "W_ION");
    select_t *wion = smart_select(all, "W_ION", ndx_groups);
    assert(wion != NULL);
    assert(selection_compare_strict(manual_wion, wion));
    free(wion);

    // select all
    select_t *smart_all = smart_select(all, NULL, NULL);
    assert(smart_all != NULL);
    assert(selection_compare_strict(all, smart_all));
    free(smart_all);

    free(system);
    free(all);
    dict_destroy(ndx_groups);
    printf("OK\n");
}

static void test_smart_select_fails(void)
{
    printf("%-40s", "smart_select (fails) ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);
    dict_t *ndx_groups = read_ndx(NDX_FILE, system);

    assert(ndx_groups != NULL);

    // query too long
    select_t *selection = smart_select(all, "name H1 H2 H3 HA HB1 HB2 HG HD11 HD12 HD13 HD21 HD22 HD23 H N CA C CB CG HN1 HN2 HN3 C12 H12A H12B C11 H11A H11B P O13 O14 O11 O12 C1 HS C2 O21 O22 C21 C22 H2R H2S C3 HX HY", NULL);
    assert(selection == NULL);

    // invalid selection
    select_t *selection2 = smart_select(selection, "resname POPE", ndx_groups);
    assert(selection2 == NULL);

    // nonexistent index group
    select_t *nonexistent = smart_select(all, "Nonexistent", ndx_groups);
    assert(nonexistent == NULL);

    // ndx groups not provided
    select_t *selection3 = smart_select(all, "Protein", NULL);
    assert(selection3 == NULL);

    // unknown query identifier
    select_t *selection4 = smart_select(all, "res POPE", ndx_groups);
    assert(selection4 == NULL);

    free(system);
    free(all);
    dict_destroy(ndx_groups);
    printf("OK\n");
}

/*  FILE *output = fopen("test.gro", "w");
    write_gro(output, selection, system->box, no_velocities, "Something");
    fclose(output);*/

void test_selection(void)
{
    test_strsplit_space();
    test_strsplit_spacetab();
    test_load_gro();

    test_match_residue_name();
    test_match_atom_name();
    test_match_atom_num();

    test_selection_create();
    test_select_system();
    test_selection_copy();
    test_selection_empty();
    test_selection_add_atom();

    test_select_atoms_atomname();
    test_select_atoms_atomname_multiple();
    test_select_atoms_atomname_empty();
    test_select_atoms_atomname_nomatch();

    test_select_atoms_resname();
    test_select_atoms_resname_multiple();
    test_select_atoms_resname_empty();
    test_select_atoms_resname_nomatch();

    test_select_atoms_atomnum();
    test_select_atoms_atomnum_nomatch();

    test_selection_cat();
    test_selection_cat_unique();

    test_selection_intersect();
    test_selection_intersect_none();

    test_selection_remove_atom();
    test_selection_remove_atom_none();
    test_selection_remove_atom_all();
    test_selection_remove_atom_duplicates();

    test_selection_remove();
    test_selection_remove_none();
    test_selection_remove_all();
    test_selection_remove_all2();

    test_selection_compare();
    test_selection_compare_empty();
    test_selection_compare_strict();
    test_selection_compare_strict_empty();
    
    test_selection_unique();
    test_selection_renumber();

    test_selection_sort();
    test_selection_sort_renumber();
    test_selection_sort_gmx();
    test_selection_sort_gmx_renumber();

    test_selection_fixres();
    test_selection_isin();

    test_selection_to_system();

    test_select_geometry_sphere();
    test_select_geometry_box();
    test_select_geometry_zcylinder();
    test_select_geometry_ycylinder();
    test_select_geometry_xcylinder();

    test_read_ndx();
    test_read_ndx_advanced();
    test_read_ndx_empty();
    test_read_ndx_nonexistent();

    test_smart_select();
    test_smart_select_fails();
}