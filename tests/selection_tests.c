// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#include "tests.h"

static void test_strsplit_space(void)
{
    printf("%-40s", "strsplit (space) ");
    fflush(stdout);

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
    fflush(stdout);

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

static void test_strstrip(void)
{
    printf("%-40s", "strstrip ");
    fflush(stdout);

    char string[1024] = "               Str\ning to st\nrip.         \n";
    strstrip(string);
    assert(!strcmp(string, "Str\ning to st\nrip."));

    printf("OK\n");
}

static void test_strstrip_zero(void)
{
    printf("%-40s", "strstrip (zero) ");
    fflush(stdout);

    char string[1024] = "";
    strstrip(string);
    assert(!strcmp(string, ""));

    printf("OK\n");
}

static void test_strstrip_white(void)
{
    printf("%-40s", "strstrip (whiteonly) ");
    fflush(stdout);

    char string[1024] = "    \n             \n  \n\n\n          ";
    strstrip(string);
    assert(!strcmp(string, ""));

    printf("OK\n");
}

static void test_strremwhite(void)
{
    printf("%-40s", "strremwhite ");
    fflush(stdout);

    char string[1024] = "bla - g . \t haha\n xx";
    strremwhite(string);
    assert(!strcmp(string, "bla-g.hahaxx"));

    printf("OK\n");
}

static void test_strremwhite_white(void)
{
    printf("%-40s", "strremwhite (whiteonly)");
    fflush(stdout);

    char string[1024] = "   \t\n\n ";
    strremwhite(string);
    assert(!strcmp(string, ""));

    printf("OK\n");
}

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

static void test_match_residue_name(void)
{
    printf("%-40s", "match_residue_name ");
    fflush(stdout);

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

static void test_match_residue_num(void)
{
    printf("%-40s", "match_residue_num ");
    fflush(stdout);

    system_t *system = load_gro(INPUT_GRO_FILE);

    int atom_ids[5] = {12, 5061, 11349, 32542, 48191};
    char residue_numbers[5][5] = {"1", "59", "110", "3898", "9115"};

    for (unsigned short i = 0; i < 5; ++i) {
        assert(match_residue_num(&system->atoms[atom_ids[i]], residue_numbers[i]));
        assert(match_residue_num(&system->atoms[atom_ids[i]], NULL));
    }

    free(system);
    printf("OK\n");
}

static void test_match_atom_name(void)
{
    printf("%-40s", "match_atom_name ");
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

    // we allocate enough memory to hold 10 atom pointers
    select_t *selection = selection_create(10);
    assert(selection->n_atoms == 0);

    system_t *system = load_gro(INPUT_GRO_FILE);

    // we assign an arbitrary atom pointer to the last allocated position in the selection
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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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

static void test_selection_add(void)
{
    printf("%-40s", "selection_add ");
    fflush(stdout);

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    size_t allocated = 1;
    select_t *selection1 = selection_create(allocated);

    select_t *selection2 = select_atoms(all, "NA", &match_residue_name);

    selection_add(&selection1, &allocated, selection2);
    selection_add(&selection1, &allocated, selection2);

    assert(selection1->n_atoms == 134);
    for (size_t i = 0; i < selection1->n_atoms; ++i) {
        assert(!strcmp(selection1->atoms[i]->residue_name, "NA"));
    }

    free(all);
    free(selection1);
    free(selection2);
    free(system);
    printf("OK\n");
}

static void test_select_atoms_atomname(void)
{
    printf("%-40s", "select_atoms (atomname) ");
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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

static void test_selection_reverse(void)
{
    printf("%-40s", "selection_reverse ");
    fflush(stdout);

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *all_copy = selection_copy(all);

    selection_reverse(all_copy);

    assert(selection_compare(all, all_copy));
    assert(!selection_compare_strict(all, all_copy));
    assert(all_copy->atoms[0]->gmx_atom_number == system->n_atoms);
    assert(all_copy->atoms[all_copy->n_atoms - 1]->gmx_atom_number == 1);

    selection_reverse(all_copy);

    assert(selection_compare_strict(all, all_copy));

    free(all_copy);

    select_t *all_minus_one = selection_copy(all);
    selection_remove_atom(all_minus_one, all_minus_one->atoms[0]);

    selection_reverse(all_minus_one);
    assert(all_minus_one->atoms[0]->gmx_atom_number == system->n_atoms);
    assert(all_minus_one->atoms[all_minus_one->n_atoms - 1]->gmx_atom_number == 2);

    free(all_minus_one);

    free(system);
    free(all);
    printf("OK\n");
}

static void test_selection_slice(void)
{
    printf("%-40s", "selection_slice ");
    fflush(stdout);

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    // simple slice
    select_t *slice1 = selection_slice(all, 4, 9);
    assert(slice1->n_atoms == 5);
    assert(slice1->atoms[0]->gmx_atom_number == 5);
    assert(slice1->atoms[4]->gmx_atom_number == 9);

    free(slice1);

    // full slice with zero
    select_t *slice2 = selection_slice(all, 0, 0);
    assert(slice2->n_atoms == all->n_atoms);
    assert(selection_compare_strict(slice2, all));

    free(slice2);

    // full slice with large end index
    select_t *slice3 = selection_slice(all, 0, 965432);
    assert(slice3->n_atoms == all->n_atoms);
    assert(selection_compare_strict(slice3, all));

    free(slice3);

    // slice with negative start
    select_t *slice4 = selection_slice(all, -9, 0);
    assert(slice4->n_atoms == 9);
    assert(slice4->atoms[0]->gmx_atom_number == 48276);
    assert(slice4->atoms[8]->gmx_atom_number == 48284);

    free(slice4);

    // slice with negative start and positive end
    select_t *slice5 = selection_slice(all, -9, 48282);
    assert(slice5->n_atoms == 7);
    assert(slice5->atoms[0]->gmx_atom_number == 48276);
    assert(slice5->atoms[6]->gmx_atom_number == 48282);

    free(slice5);

    // slice with negative start and negative end
    select_t *slice6 = selection_slice(all, -9, -2);
    assert(slice6->n_atoms == 7);
    assert(slice6->atoms[0]->gmx_atom_number == 48276);
    assert(slice6->atoms[6]->gmx_atom_number == 48282);

    free(slice6);

    // slice with positive start and negative end
    select_t *slice7 = selection_slice(all, 48275, -2);
    assert(slice7->n_atoms == 7);
    assert(slice7->atoms[0]->gmx_atom_number == 48276);
    assert(slice7->atoms[6]->gmx_atom_number == 48282);

    free(slice7);

    // slice is very negative
    select_t *slice8 = selection_slice(all, -65432, 0);
    assert(slice8->n_atoms == all->n_atoms);
    assert(selection_compare_strict(slice8, all));

    free(slice8);

    // failed slice: start is higher than end
    select_t *failed_slice1 = selection_slice(all, 36, 18);
    assert(failed_slice1 == NULL);

    // failed slice: negative start is higher than positive end
    select_t *failed_slice2 = selection_slice(all, -652, 18);
    assert(failed_slice2 == NULL);

    // failed slice: negative start is higher than negative end
    select_t *failed_slice3 = selection_slice(all, -18, -652);
    assert(failed_slice3 == NULL);

    // failed slice: end is way too negative
    select_t *failed_slice4 = selection_slice(all, 12, -64320);
    assert(failed_slice4 == NULL);


    free(system);
    free(all);
    printf("OK\n");
}

static void test_selection_fixres(void)
{
    printf("%-40s", "selection_fixres ");
    fflush(stdout);

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
    fflush(stdout);

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

static void test_selection_getnres(void)
{
    printf("%-40s", "selection_getnres ");
    fflush(stdout);

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    assert(selection_getnres(all) == 9207);

    select_t *membrane = select_atoms(all, "POPE POPG", &match_residue_name);
    assert(selection_getnres(membrane) == 168);

    free(all);
    free(membrane);
    free(system);
    printf("OK\n");
}

static void test_selection_getresnames(void)
{
    printf("%-40s", "selection_getresnames ");
    fflush(stdout);

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    list_t *all_resnames = selection_getresnames(all);

    assert(all_resnames->n_items == 8);
    assert(list_index(all_resnames, "LEU") >= 0);
    assert(list_index(all_resnames, "SER") >= 0);
    assert(list_index(all_resnames, "NHE") >= 0);
    assert(list_index(all_resnames, "POPE") >= 0);
    assert(list_index(all_resnames, "POPG") >= 0);
    assert(list_index(all_resnames, "SOL") >= 0);
    assert(list_index(all_resnames, "NA") >= 0);
    assert(list_index(all_resnames, "CL") >= 0);
    assert(list_index(all_resnames, "NAH") < 0);
    assert(list_index(all_resnames, "NONEXISTENT") < 0);

    list_destroy(all_resnames);
    free(all);
    free(system);
    printf("OK\n");
}

static void test_selection_splitbyres(void)
{
    printf("%-40s", "selection_splitbyres ");
    fflush(stdout);

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t **array = NULL;
    size_t n_residues = selection_splitbyres(all, &array);

    assert(n_residues == selection_getnres(all));

    char serine_atoms[11][4] = {"N", "H", "CA", "HA", "CB", "HB1", "HB2", "OG", "HG", "C", "O"};
    assert(array[1]->n_atoms == 11);
    for (size_t i = 0; i < array[1]->n_atoms; ++i) {
        assert(strcmp(serine_atoms[i], array[1]->atoms[i]->atom_name) == 0);
    }

    assert(strcmp(array[2]->atoms[0]->atom_name, "N") == 0);

    assert(strcmp(array[67]->atoms[0]->residue_name, "POPE") == 0);
    assert(array[67]->n_atoms == 125);

    assert(array[9100]->n_atoms == 3);
    assert(array[9206]->n_atoms == 1);
    assert(array[9206]->atoms[0]->atom_number == 48284);
    

    for (size_t i = 0; i < n_residues; ++i) {
        free(array[i]);
    }

    free(array);
    free(all);
    free(system);
    printf("OK\n");
}

static void test_selection_splitbyres_broken(void)
{
    printf("%-40s", "selection_splitbyres (broken) ");
    fflush(stdout);

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *selection1 = select_atoms(all, "HD21 HD22 HD23", &match_atom_name);
    select_t *selection2 = select_atoms(all, "CD2 C", &match_atom_name);
    select_t *selection3 = select_atoms_d(all, "N", &match_atom_name);

    select_t *selection12 = selection_cat_d(selection1, selection2);
    select_t *selection123 = selection_cat_d(selection12, selection3);

    select_t **array = NULL;
    size_t n_residues = selection_splitbyres(selection123, &array);

    assert(n_residues == selection_getnres(selection123));

    assert(array[0]->n_atoms == 6);
    assert(array[11]->n_atoms == 7);

    assert(strcmp(array[12]->atoms[0]->atom_name, "C") == 0);
    assert(strcmp(array[12]->atoms[1]->atom_name, "N") == 0);

    /*for (size_t i = 0; i < n_residues; ++i) {
        write_gro(stdout, array[i], system->box, no_velocities, "something");
    }*/
    
    for (size_t i = 0; i < n_residues; ++i) {
        free(array[i]);
    }

    free(array);
    free(selection123);
    free(system);
    printf("OK\n");
}

static void test_selection_splitbyres_empty(void)
{
    printf("%-40s", "selection_splitbyres (empty) ");
    fflush(stdout);

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *selection = select_atoms_d(all, "PO4", &match_atom_name);

    select_t **array = NULL;
    
    size_t n_residues = selection_splitbyres(selection, &array);

    assert(array == NULL);
    assert(n_residues == 0);

    free(selection);
    free(system);
    printf("OK\n");
}

static void test_selection_to_system(void)
{
    printf("%-40s", "selection_to_system ");
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

    system_t *system = load_gro(INPUT_GRO_FILE);
    dict_t *ndx_groups = read_ndx("index.ndx", system);

    assert(ndx_groups == NULL);

    free(system);
    printf("OK\n");
}

static void test_smart_select(void)
{
    printf("%-40s", "smart_select ");
    fflush(stdout);

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);
    dict_t *ndx_groups = read_ndx(NDX_FILE, system);

    assert(ndx_groups != NULL);

    // residue name selection
    select_t *manual_pope = select_atoms(all, "POPE", &match_residue_name);
    select_t *pope1 = smart_select(all, "resname POPE", ndx_groups);
    select_t *pope2 = smart_select(all, "resnamePOPE", NULL);
    select_t *pope3 = smart_select(all, "resname      POPE", ndx_groups);
    select_t *not_pope = smart_select(all, "not resname POPE", ndx_groups);
    assert(pope1 != NULL);
    assert(pope2 != NULL);
    assert(pope3 != NULL);
    assert(not_pope != NULL);
    assert(selection_compare_strict(manual_pope, pope1));
    assert(selection_compare_strict(manual_pope, pope2));
    assert(selection_compare_strict(manual_pope, pope3));
    assert(not_pope->n_atoms == all->n_atoms - pope1->n_atoms);
    for (size_t i = 0; i < not_pope->n_atoms; ++i) {
        assert(strcmp(not_pope->atoms[i]->residue_name, "POPE"));
    }

    free(pope1);
    free(pope2);
    free(pope3);
    free(manual_pope);
    free(not_pope);

    // residue number selection
    select_t *manual_residues = select_atoms(all, "8874 7734 4 5 6 9207 1", &match_residue_num);
    select_t *smart1 = smart_select(all, "resid 8874 7734 4 5 6 9207 1", NULL);
    select_t *smart2 = smart_select(all, "resid8874 7734 4 5 6 9207 1", NULL);
    select_t *not_resid = smart_select(all, "! resid 8874 7734 4 5 6 9207 1", NULL);
    assert(smart1 != NULL);
    assert(smart2 != NULL);
    assert(not_resid != NULL);
    assert(selection_compare_strict(manual_residues, smart1));
    assert(selection_compare_strict(manual_residues, smart2));
    assert(not_resid->n_atoms == all->n_atoms - smart1->n_atoms);
    for (size_t i = 0; i < not_resid->n_atoms; ++i) {
        assert(not_resid->atoms[i]->residue_number != 8874);
        assert(not_resid->atoms[i]->residue_number != 7734);
        assert(not_resid->atoms[i]->residue_number != 4);
        assert(not_resid->atoms[i]->residue_number != 5);
        assert(not_resid->atoms[i]->residue_number != 9207);
        assert(not_resid->atoms[i]->residue_number != 1);
    }
    free(manual_residues);
    free(smart1);
    free(smart2);
    free(not_resid);

    // atom name selection
    select_t *manual_hydrogens = select_atoms(all, "H1 H2 H3 HA HB1 HB2 HG HD11 HD12 HD13 HD21 HD22 HD23 H", &match_atom_name);
    select_t *hydrogens = smart_select(all, "name H1 H2 H3 HA HB1 HB2 HG HD11 HD12 HD13 HD21 HD22 HD23 H", NULL);
    select_t *not_hydrogens = smart_select(all, "not name H1 H2 H3 HA HB1 HB2 HG HD11 HD12 HD13 HD21 HD22 HD23 H", ndx_groups);
    select_t *not_hydrogens_from_hydrogens = smart_select(hydrogens, "not name H1 H2 H3 HA HB1 HB2 HG HD11 HD12 HD13 HD21 HD22 HD23 H", NULL);
    assert(hydrogens != NULL);
    assert(not_hydrogens != NULL);
    assert(not_hydrogens_from_hydrogens != NULL);
    assert(selection_compare_strict(manual_hydrogens, hydrogens));
    assert(not_hydrogens->n_atoms == all->n_atoms - hydrogens->n_atoms);
    for (size_t i = 0; i < not_hydrogens->n_atoms; ++i) {
        assert(strcmp(not_hydrogens->atoms[i]->atom_name, "H1"));
        assert(strcmp(not_hydrogens->atoms[i]->atom_name, "H2"));
        assert(strcmp(not_hydrogens->atoms[i]->atom_name, "H3"));
        assert(strcmp(not_hydrogens->atoms[i]->atom_name, "HA"));
        assert(strcmp(not_hydrogens->atoms[i]->atom_name, "HB1"));
        assert(strcmp(not_hydrogens->atoms[i]->atom_name, "HB2"));
        assert(strcmp(not_hydrogens->atoms[i]->atom_name, "HG"));
        assert(strcmp(not_hydrogens->atoms[i]->atom_name, "HD11"));
        assert(strcmp(not_hydrogens->atoms[i]->atom_name, "HD12"));
        assert(strcmp(not_hydrogens->atoms[i]->atom_name, "HD13"));
        assert(strcmp(not_hydrogens->atoms[i]->atom_name, "HD13"));
        assert(strcmp(not_hydrogens->atoms[i]->atom_name, "HD21"));
        assert(strcmp(not_hydrogens->atoms[i]->atom_name, "HD22"));
        assert(strcmp(not_hydrogens->atoms[i]->atom_name, "HD23"));
        assert(strcmp(not_hydrogens->atoms[i]->atom_name, "H"));
    }
    assert(not_hydrogens_from_hydrogens->n_atoms == 0);

    free(manual_hydrogens);
    free(hydrogens);
    free(not_hydrogens);
    free(not_hydrogens_from_hydrogens);


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
    select_t *not_serial = smart_select(all, "! serial 76532 73 6542 9875 23463 2 987654 42653", ndx_groups);
    assert(serial != NULL);
    assert(serial2 != NULL);
    assert(serial3 != NULL);
    assert(serial4 != NULL);
    assert(not_serial != NULL);
    assert(selection_compare_strict(manual_serial, serial));
    assert(selection_compare_strict(manual_serial, serial2));
    assert(selection_compare_strict(manual_serial, serial3));
    assert(selection_compare_strict(manual_serial, serial4));
    assert(not_serial->n_atoms == all->n_atoms - serial4->n_atoms);
    for (size_t i = 0; i < not_serial->n_atoms; ++i) {
        assert(not_serial->atoms[i]->gmx_atom_number != 76532);
        assert(not_serial->atoms[i]->gmx_atom_number != 73);
        assert(not_serial->atoms[i]->gmx_atom_number != 6542);
        assert(not_serial->atoms[i]->gmx_atom_number != 9875);
        assert(not_serial->atoms[i]->gmx_atom_number != 23463);
        assert(not_serial->atoms[i]->gmx_atom_number != 2);
        assert(not_serial->atoms[i]->gmx_atom_number != 987654);
        assert(not_serial->atoms[i]->gmx_atom_number != 42653);
    }

    free(manual_serial);
    free(serial);
    free(serial2);
    free(serial3);
    free(serial4);
    free(not_serial);

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

    select_t *not_wion = smart_select(all, "not W_ION", ndx_groups);
    assert(not_wion != NULL);
    assert(not_wion->n_atoms == all->n_atoms - wion->n_atoms);
    assert(not_wion->atoms[13]->atom_number == 14);
    assert(not_wion->atoms[13]->residue_number == 1);
    assert(!strcmp(not_wion->atoms[13]->residue_name, "LEU"));
    assert(!strcmp(not_wion->atoms[13]->atom_name, "HD12"));

    free(wion);
    free(not_wion);

    // multi-word index selection
    select_t *normal_ion = smart_select(all, "ION", ndx_groups);
    select_t *multi_ion = smart_select(all, "Interesting Selection", ndx_groups);
    assert(multi_ion != NULL);
    assert(selection_compare_strict(normal_ion, multi_ion));
    free(normal_ion);
    free(multi_ion);

    // select all
    select_t *smart_all = smart_select(all, NULL, NULL);
    select_t *smart_all_query1 = smart_select(all, "all", NULL);
    select_t *smart_all_query2 = smart_select(all, "   all   ", NULL);
    assert(smart_all != NULL);
    assert(smart_all_query1 != NULL);
    assert(smart_all_query2 != NULL);
    assert(selection_compare_strict(all, smart_all));
    assert(selection_compare_strict(all, smart_all_query1));
    assert(selection_compare_strict(all, smart_all_query2));
    free(smart_all);
    free(smart_all_query1);
    free(smart_all_query2);

    free(system);
    free(all);
    dict_destroy(ndx_groups);
    printf("OK\n");
}

static void test_smart_select_fails(void)
{
    printf("%-40s", "smart_select (fails) ");
    fflush(stdout);

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);
    dict_t *ndx_groups = read_ndx(NDX_FILE, system);

    assert(ndx_groups != NULL);

    select_t *selection = NULL;

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

static void test_smart_select_advanced(void)
{
    printf("%-40s", "smart_select (advanced) ");
    fflush(stdout);

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);
    dict_t *ndx_groups = read_ndx(NDX_FILE, system);

    // residue names
    select_t *pope = smart_select(all, "resname POPE", NULL);
    select_t *popg = smart_select(all, "resname POPG", NULL);
    select_t *leu = smart_select(all, "resname LEU", NULL); 
    select_t *manual_pope_or_popg = selection_cat_d(pope, popg);
    select_t *manual_pope_or_popg_or_leu = selection_cat(manual_pope_or_popg, leu);
    select_t *auto_pope_or_popg1 = smart_select(all, "resname POPE or resname POPG   ", NULL);
    select_t *auto_pope_or_popg2 = smart_select(all, "   resname POPE  || resname POPG", NULL);
    select_t *auto_pope_and_popg1 = smart_select(all, " resname POPE and resname POPG", NULL);
    select_t *auto_pope_and_popg2 = smart_select(all, "resname POPE &&      resname POPG  ", NULL);
    select_t *auto_pope_or_popg_or_leu = smart_select(all, "resname POPE        or resname POPG || resname LEU", NULL);
    assert(auto_pope_or_popg1 != NULL);
    assert(auto_pope_or_popg2 != NULL);
    assert(auto_pope_and_popg1 != NULL);
    assert(auto_pope_and_popg2 != NULL);
    assert(auto_pope_or_popg_or_leu != NULL);
    assert(selection_compare_strict(manual_pope_or_popg, auto_pope_or_popg1));
    assert(selection_compare_strict(manual_pope_or_popg, auto_pope_or_popg2));
    assert(selection_compare_strict(auto_pope_or_popg_or_leu, manual_pope_or_popg_or_leu));
    assert(auto_pope_and_popg1->n_atoms == 0);
    assert(auto_pope_and_popg2->n_atoms == 0);

    free(leu);
    free(manual_pope_or_popg);
    free(manual_pope_or_popg_or_leu);
    free(auto_pope_or_popg1);
    free(auto_pope_or_popg2);
    free(auto_pope_and_popg1);
    free(auto_pope_and_popg2);
    free(auto_pope_or_popg_or_leu);

    
    // residue numbers
    select_t *res1 = smart_select(all, "resid 8874 7734 4 5", NULL);
    select_t *res2 = smart_select(all, "resid 6 9207 1", NULL);
    select_t *manual_res12 = selection_cat_d(res1, res2);
    select_t *auto_res12 = smart_select(all, "resid 8874 7734 4 5 or resid 6 9207 1", NULL);
    assert(auto_res12 != NULL);
    assert(selection_compare_strict(manual_res12, auto_res12));

    select_t *not_resid = smart_select(all, "! resid 8874 7734 4 5 6 9207 1 || resid 8874", NULL);
    assert(not_resid != NULL);
    assert(not_resid->n_atoms == all->n_atoms - manual_res12->n_atoms + 3);
    assert(not_resid->atoms[not_resid->n_atoms - 1]->residue_number == 8874);

    free(manual_res12);
    free(auto_res12);
    free(not_resid);

    // atom names
    select_t *hydrogens1 = smart_select(all, "name H1 H2 H3 HA HB1 HB2 HG", NULL);
    select_t *hydrogens2 = smart_select(all, "name HD11 HD12 HD13 HD21 HD22 HD23 H", NULL);
    select_t *hd11 = smart_select(all, "name HD11", NULL);
    select_t *manual_hydrogens12 = selection_cat_d(hydrogens1, hydrogens2);
    select_t *auto_hydrogens12 = smart_select(all, "name H1 H2 H3 HA HB1 HB2 HG or name HD11 HD12 HD13 HD21 HD22 HD23 H", NULL);
    select_t *auto_hydrogens_and1 = smart_select(all, "name H1 H2 H3 HA HB1 HD11 HB2 HG and name HD11 HD12 HD13 HD21 HD22 HD23 H", NULL);
    select_t *auto_hydrogens_and2 = smart_select(all, "name H1 H2 H3 HA HB1 HB2 HG HD11 && name HD11 HD12 HD13 HD21 HD22 HD23 H", NULL);
    select_t *auto_hydrogens_and3 = smart_select(all, "name H1 H2 H3 HA HB1 HB2 HG HD11 && name HD11 HD12 HD13 HD21 HD22 HD23 H and name HD13 HD21 H1 H2 HD11", NULL);
    assert(auto_hydrogens12 != NULL);
    assert(auto_hydrogens_and1 != NULL);
    assert(auto_hydrogens_and2 != NULL);
    assert(selection_compare_strict(manual_hydrogens12, auto_hydrogens12));
    assert(selection_compare_strict(hd11, auto_hydrogens_and1));
    assert(selection_compare_strict(hd11, auto_hydrogens_and2));
    assert(selection_compare_strict(hd11, auto_hydrogens_and3));
    
    free(hd11);
    free(manual_hydrogens12);
    free(auto_hydrogens12);
    free(auto_hydrogens_and1);
    free(auto_hydrogens_and2);
    free(auto_hydrogens_and3);

    // atom numbers
    select_t *serial1 = smart_select(all, "serial 73 6542 9875", NULL);
    select_t *serial2 = smart_select(all, "serial 23463 2 42653", NULL);
    select_t *serial3 = smart_select(all, "! serial 73 23463 2", NULL);
    select_t *serial12 = selection_cat_d(serial1, serial2);
    select_t *manual_serial_complex = selection_intersect_d(serial12, serial3);
    select_t *auto_serial_complex = smart_select(all, "serial 73 6542 9875 or serial 23463 2 42653 && not serial 73 23463 2", NULL);
    assert(auto_serial_complex != NULL);
    assert(selection_compare_strict(manual_serial_complex, auto_serial_complex));
    assert(auto_serial_complex->n_atoms == 3);

    free(manual_serial_complex);
    free(auto_serial_complex);

    // ndx groups
    select_t *sidechain = smart_select(all, "SideChain", ndx_groups);
    select_t *protein = smart_select(all, "Protein", ndx_groups);
    select_t *manual_scprotein = selection_cat(sidechain, protein);
    select_t *manual_scprotein_unique = selection_cat_unique_d(sidechain, protein);
    select_t *auto_scprotein = smart_select(all, "SideChain || Protein", ndx_groups);
    assert(auto_scprotein != NULL);
    assert(manual_scprotein->n_atoms != manual_scprotein_unique->n_atoms);
    assert(selection_compare_strict(manual_scprotein_unique, auto_scprotein));


    free(manual_scprotein);
    free(manual_scprotein_unique);
    free(auto_scprotein);

    // residue names with residue numbers
    select_t *pope2 = smart_select(all, "resname POPE", NULL);
    select_t *resid = smart_select(all, "resid 29 33 38 8643 1315", NULL);
    select_t *manual_poperes = selection_intersect(pope2, resid);
    select_t *auto_poperes = smart_select(all, "resname POPE && resid 29 33 38 8643 1315", NULL);
    select_t *auto_poperes2 = smart_select(all, "resid 29 33 38 8643 1315 && resname POPE", NULL);
    assert(auto_poperes != NULL);
    assert(auto_poperes2 != NULL);
    assert(selection_compare_strict(manual_poperes, auto_poperes));
    assert(selection_compare_strict(manual_poperes, auto_poperes2));

    free(manual_poperes);
    free(auto_poperes);
    free(auto_poperes2);

    // residue names with atom names
    select_t *atoms = smart_select(all, "name HN1 HN2 HN3 C2 HS", NULL);
    select_t *manual_popeatoms = selection_intersect(pope2, atoms);
    select_t *auto_popeatoms = smart_select(all, "resname POPE and name HN1 HN2 HN3 C2 HS", NULL);
    assert(auto_popeatoms != NULL);
    assert(selection_compare_strict(manual_popeatoms, auto_popeatoms));

    free(manual_popeatoms);
    free(auto_popeatoms);

    // residue names with atom numbers
    select_t *serial = smart_select(all, "serial 468 469 470 471 472 473 474", NULL);
    select_t *manual_popeserial = selection_intersect(serial, pope2);
    select_t *auto_popeserial = smart_select(all, "serial 468 469 470 471 472 473 474 && resname POPE", NULL);
    assert(auto_popeserial != NULL);
    assert(selection_compare_strict(manual_popeserial, auto_popeserial));
    for (size_t i = 0; i < auto_popeserial->n_atoms; ++i) {
        assert(auto_popeserial->atoms[i]->residue_number == 23);
    }

    free(manual_popeserial);
    free(auto_popeserial);
    free(pope2);

    // residue names with NDX groups
    select_t *serines = smart_select(all, "resname SER", NULL);
    select_t *backbone = smart_select(all, "Backbone", ndx_groups);
    select_t *manual_bbser = selection_intersect(serines, backbone);
    select_t *manual_bbser_or = selection_cat_unique_d(serines, backbone);
    select_t *auto_bbser = smart_select(all, "resname SER and Backbone", ndx_groups);
    select_t *auto_bbser_or = smart_select(all, "resname SER or Backbone", ndx_groups);

    assert(auto_bbser_or != NULL);
    assert(auto_bbser != NULL);
    assert(selection_compare_strict(manual_bbser, auto_bbser));
    assert(selection_compare_strict(manual_bbser_or, auto_bbser_or));

    free(manual_bbser);
    free(auto_bbser);
    free(manual_bbser_or);
    free(auto_bbser_or);

    // residue numbers with atom names
    select_t *manual_residatoms = selection_intersect(resid, atoms);
    select_t *auto_residatoms = smart_select(all, "   resid 29 33 38 8643 1315  && name HN1 HN2 HN3 C2 HS", NULL);

    assert(auto_residatoms != NULL);
    assert(selection_compare_strict(manual_residatoms, auto_residatoms));

    free(manual_residatoms);
    free(auto_residatoms);

    // residue numbers with atom numbers
    select_t *serial_for_resid = smart_select(all, "serial 1814    1817    1819", NULL);
    select_t *manual_residserial = selection_intersect(resid, serial_for_resid);
    select_t *auto_residserial = smart_select(all, "resid 29        33 38 8643 1315 && serial 1814    1817    1819          ", NULL);

    assert(auto_residserial != NULL);
    assert(auto_residserial->n_atoms == 3);
    assert(selection_compare_strict(manual_residserial, auto_residserial));

    free(manual_residserial);
    free(auto_residserial);
    free(serial_for_resid);

    // residue numbers with NDX groups
    select_t *nonprotein = smart_select(all, "non-Protein", ndx_groups);
    select_t *manual_residindex = selection_intersect_d(resid, nonprotein);
    select_t *auto_residindex = smart_select(all, "resid 29 33 38 8643 1315 and non-Protein  ", ndx_groups);

    assert(auto_residindex != NULL);
    assert(selection_compare_strict(manual_residindex, auto_residindex));

    free(manual_residindex);
    free(auto_residindex);

    // atom names with atom numbers
    select_t *atoms2 = smart_select(all, "name P O13 C1", NULL);
    select_t *manual_atomsserial = selection_intersect(atoms2, serial);
    select_t *auto_atomsserial = smart_select(all, "name P O13 C1 and serial 468 469 470 471 472 473 474", NULL);

    assert(auto_atomsserial != NULL);
    assert(auto_atomsserial->n_atoms == 3);
    assert(selection_compare_strict(manual_atomsserial, auto_atomsserial));

    free(atoms2);
    free(manual_atomsserial);
    free(auto_atomsserial);

    // atom names with NDX groups
    select_t *membrane = smart_select(all, "Membrane", ndx_groups);
    select_t *manual_atomsindex = selection_intersect_d(membrane, atoms);
    select_t *auto_atomsindex = smart_select(all, "name HN1 HN2 HN3 C2 HS && Membrane", ndx_groups);

    assert(auto_atomsindex != NULL);
    assert(selection_compare_strict(manual_atomsindex, auto_atomsindex));

    free(manual_atomsindex);
    free(auto_atomsindex);

    // atom numbers with NDX groups
    select_t *pope_ndx = smart_select(all, "POPE", ndx_groups);
    select_t *manual_serialindex = selection_intersect_d(pope_ndx, serial);
    select_t *auto_serialindex = smart_select(all, "POPE && serial 468 469 470 471 472 473 474", ndx_groups);

    assert(auto_serialindex != NULL);
    assert(selection_compare_strict(manual_serialindex, auto_serialindex));

    free(manual_serialindex);
    free(auto_serialindex);


    // complex operation with AND, OR and NOT
    select_t *lipid_residues = smart_select(all, "resid 31 164 165 168", NULL);
    select_t *membrane_ndx = smart_select(all, "Membrane", ndx_groups);
    select_t *not_popg_resname = smart_select(all, "not resname POPG", ndx_groups);
    select_t *not_p_name = smart_select(all, "! name P", NULL);
    select_t *pope_ndx2 = smart_select(all, "POPE", ndx_groups);
    select_t *protein_serial = smart_select(all, "serial 132 191 150 162", NULL);

    select_t *step1 = selection_intersect_d(lipid_residues, membrane_ndx);
    select_t *step2 = selection_intersect_d(step1, not_popg_resname);
    select_t *step3 = selection_intersect_d(step2, not_p_name);
    select_t *step4 = selection_cat_unique_d(step3, pope_ndx2);
    select_t *step5 = selection_cat_unique_d(step4, protein_serial);
    
    select_t *auto_step5 = smart_select(all, "resid 31 164 165 168 and Membrane && not resname POPG and ! name P or POPE || serial 132 191 150 162", ndx_groups);

    assert(auto_step5);
    assert(selection_compare_strict(step5, auto_step5));

    free(step5);
    free(auto_step5);
    

    free(system);
    free(all);
    dict_destroy(ndx_groups);
    printf("OK\n");
}


static void test_smart_select_to(void)
{
    printf("%-40s", "smart_select (to) ");
    fflush(stdout);

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);
    dict_t *ndx_groups = read_ndx(NDX_FILE, system);

    // success selections

    select_t *manual_selection1 = smart_select(all, "serial 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44", NULL);
    select_t *manual_selection2 = smart_select(all, "serial 1 2 3 4 5 6", NULL);
    select_t *manual_selection3 = smart_select(all, "resid 1 2 3", NULL);
    select_t *manual_selection4 = smart_select(all, "serial 1 2 3 4 5 6 7 8 9 10 11 && Backbone", ndx_groups);
    select_t *manual_selection5 = smart_select(all, "serial 1 2 3 4 5 6 12 13 14 15 1654 1655 1656 1657 1658 1659 1660", NULL);
    select_t *manual_selection6 = smart_select(all, "not resid 14 15 16 17 18 19", NULL);

    select_t *selection1 = smart_select(all, "serial 1 to 44", NULL);
    select_t *selection2 = smart_select(all, "serial   1 -    6", NULL);
    select_t *selection3 = smart_select(all, "resid 1 to 3", NULL);
    select_t *selection4 = smart_select(all, "serial 1 to 11 and Backbone", ndx_groups);
    select_t *selection5 = smart_select(all, "serial 1 - 6 12 13 14 15 1654 to 1660", NULL);
    select_t *selection6 = smart_select(all, "not resid 14 - 19", NULL);
    select_t *selection7 = smart_select(all, "serial 8 to 8 to 8", NULL);

    assert(selection1 != NULL);
    assert(selection2 != NULL);
    assert(selection3 != NULL);
    assert(selection4 != NULL);
    assert(selection5 != NULL);
    assert(selection6 != NULL);
    assert(selection7 != NULL);
    assert(selection_compare_strict(manual_selection1, selection1));
    assert(selection_compare_strict(manual_selection2, selection2));
    assert(selection_compare_strict(manual_selection3, selection3));
    assert(selection_compare_strict(manual_selection4, selection4));
    assert(selection_compare_strict(manual_selection5, selection5));
    assert(selection_compare_strict(manual_selection6, selection6));
    assert(selection7->n_atoms == 1);

    free(manual_selection1);
    free(manual_selection2);
    free(manual_selection3);
    free(manual_selection4);
    free(manual_selection5);
    free(manual_selection6);

    free(selection1);
    free(selection2);
    free(selection3);    
    free(selection4);
    free(selection5);
    free(selection6);
    free(selection7);

    // fail selections

    select_t *selection_f1 = smart_select(all, "serial to", NULL);
    select_t *selection_f2 = smart_select(all, "serial 1 to ", NULL);
    select_t *selection_f3 = smart_select(all, "resid 1 - g", NULL);
    select_t *selection_f4 = smart_select(all, "serial - 2", NULL);
    select_t *selection_f5 = smart_select(all, "resid g - h", NULL);
    select_t *selection_f6 = smart_select(all, "serial 1 - to 5", NULL);
    select_t *selection_f7 = smart_select(all, "serial 8 - 6", NULL);
    select_t *selection_f8 = smart_select(all, "serial 4 to 9 - 5", NULL);

    assert(selection_f1 == NULL);
    assert(selection_f2 == NULL);
    assert(selection_f3 == NULL);
    assert(selection_f4 == NULL);
    assert(selection_f5 == NULL);
    assert(selection_f6 == NULL);
    assert(selection_f7 == NULL);
    assert(selection_f8 == NULL);



    dict_destroy(ndx_groups);
    free(all);
    free(system);
    printf("OK\n");
}

static void test_smart_select_advanced_fails(void)
{
    printf("%-40s", "smart_select (advanced, fails) ");
    fflush(stdout);

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);
    dict_t *ndx_groups = read_ndx(NDX_FILE, system);

    select_t *selection1 = smart_select(all, "resid 766 43 234 && 1 to 54", NULL);
    assert(selection1 == NULL);
    select_t *selection2 = smart_select(all, "resid 766 43 24 || 1 to 54", NULL);
    assert(selection2 == NULL);
    select_t *selection3 = smart_select(all, "&& 1 to 44 || resname POPE", NULL);
    assert(selection3 == NULL);
    select_t *selection4 = smart_select(all, "name HD11 HD12 && resname && POPE", NULL);
    assert(selection4 == NULL);
    select_t *selection5 = smart_select(all, "resname POPG && x > 50", NULL);
    assert(selection5 == NULL);
    select_t *selection6 = smart_select(all, "resid 1 to 4 || Nonexistent", ndx_groups);
    assert(selection6 == NULL);
    select_t *selection7 = smart_select(all, "serial 1 - 64 &&", NULL);
    assert(selection7 == NULL);
    select_t *selection8 = smart_select(all, "resname POPE POPG ||", NULL);
    assert(selection8 == NULL);

    free(system);
    free(all);
    dict_destroy(ndx_groups);
    printf("OK\n");
}

static void test_smart_select_parentheses(void)
{
    printf("%-40s", "smart_select (parentheses) ");
    fflush(stdout);

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);
    dict_t *ndx_groups = read_ndx(NDX_FILE, system);

    select_t *selection1 = smart_select(all, "resid 65 to 67", NULL);
    select_t *selection1_par = smart_select(all, "(resid 65 to 67)", NULL);
    assert(selection1_par != NULL);
    assert(selection_compare_strict(selection1, selection1_par));

    free(selection1);
    free(selection1_par);

    select_t *selection2 = smart_select(all, "POPE && serial 500 to 600", ndx_groups);
    select_t *selection2_par = smart_select(all, "(POPE or POPG) and serial 500 to 600", ndx_groups);
    assert(selection2_par != NULL);
    assert(selection_compare_strict(selection2, selection2_par));

    free(selection2);
    free(selection2_par);

    select_t *selection3_1 = smart_select(all, "serial 1 to 13 && name H1 H2 H3", NULL);
    select_t *selection3_2 = smart_select(all, "resname POPE && name P && resid 147 - 149", ndx_groups);
    select_t *selection3 = selection_cat_d(selection3_1, selection3_2);
    select_t *selection3_par = smart_select(all, "(serial 1 to 13 && name H1 H2 H3) || (resname POPE && name P && resid 147 - 149)", ndx_groups);
    select_t *selection3_par_adv = smart_select(all, "(   serial 1 to 13 && name H1 H2 H3 )  || (    resname POPE && name P && resid 147 - 149 )   ", ndx_groups);
    assert(selection3_par != NULL);
    assert(selection3_par_adv != NULL);
    assert(selection_compare_strict(selection3, selection3_par));
    assert(selection_compare_strict(selection3, selection3_par_adv));

    free(selection3);
    free(selection3_par);
    free(selection3_par_adv);

    select_t *selection4 = smart_select(all, "Backbone && serial 1 to 5", ndx_groups);
    select_t *selection4_par = smart_select(all, "(Backbone) and serial 1 to 5", ndx_groups);
    assert(selection4_par != NULL);
    assert(selection_compare_strict(selection4, selection4_par));

    free(selection4);
    free(selection4_par);

    select_t *selection5_1 = smart_select(all, "Protein and resid 1 to 5", ndx_groups);
    select_t *selection5_2 = smart_select(all, "resname POPE && name P", NULL);
    select_t *selection5_12 = selection_cat_unique_d(selection5_1, selection5_2);
    select_t *selection5_3 = smart_select(all, "resid 1 to 33", NULL);
    select_t *selection5_123 = selection_intersect_d(selection5_12, selection5_3);
    select_t *selection5_4 = smart_select(all, "ION or resid 9088 9089", ndx_groups);
    select_t *selection5 = selection_cat_unique_d(selection5_123, selection5_4);
    select_t *selection5_par = smart_select(all, "( ( (Protein and resid 1 to 5) || (resname POPE && name P) ) && resid 1 to 33 ) || (ION or resid 9088 to 9089)", ndx_groups);
    select_t *selection5_par_alt = smart_select(all, "(((Protein and resid 1 to 5) || (resname POPE && name P)) && resid 1 to 33) || (ION or resid 9088 to 9089)", ndx_groups);
    assert(selection5_par != NULL);
    assert(selection5_par_alt != NULL);
    assert(selection_compare_strict(selection5, selection5_par));
    assert(selection_compare_strict(selection5, selection5_par_alt));
    free(selection5);
    free(selection5_par);
    free(selection5_par_alt);

    select_t *selection6 = smart_select(all, "resname LEU && name CA", NULL);
    select_t *selection6_par = smart_select(all, "(resname LEU and name CA ) || (resname LEU and name CA )", NULL);
    assert(selection6_par != NULL);
    assert(selection_compare_strict(selection6, selection6_par));
    free(selection6);
    free(selection6_par);

    select_t *selection7_1 = smart_select(all, "resname SER and not serial 20 to 30", NULL);
    select_t *selection7_2 = smart_select(all, "! resname POPE && name P", NULL);
    select_t *selection7 = selection_cat_unique_d(selection7_1, selection7_2);
    select_t *selection7_par = smart_select(all, "(resname SER and not serial 20 to 30) || (! resname POPE && name P)", NULL);
    assert(selection7_par != NULL);
    assert(selection_compare_strict(selection7, selection7_par));
    free(selection7);
    free(selection7_par);

    select_t *selection8 = smart_select(all, "not resname SOL NA CL && not resid 1 to 15", NULL);
    select_t *selection8_par1 = smart_select(all, "! (resname SOL NA CL || resid 1 to 15)", NULL);
    select_t *selection8_par2 = smart_select(all, "not (resname SOL NA CL or resid 1 to 15)", NULL);
    select_t *selection8_par3 = smart_select(all, "not ( ! (not resname SOL NA CL && not resid 1 to 15))", NULL);
    select_t *selection8_par4 = smart_select(all, "not (! (not resname SOL NA CL && not resid 1 to 15))", NULL);
    assert(selection8_par1 != NULL);
    assert(selection8_par2 != NULL);
    assert(selection8_par3 != NULL);
    assert(selection8_par4 != NULL);
    assert(selection_compare_strict(selection8, selection8_par1));
    assert(selection_compare_strict(selection8, selection8_par2));
    assert(selection_compare_strict(selection8, selection8_par3));
    assert(selection_compare_strict(selection8, selection8_par4));
    free(selection8);
    free(selection8_par1);
    free(selection8_par2);
    free(selection8_par3);
    free(selection8_par4);

    select_t *selection9_1 = smart_select(all, "not resname SOL NA CL && not resid 1 to 15", NULL);
    select_t *selection9 = smart_select(selection9_1, "resname POPE", NULL);
    free(selection9_1);
    select_t *selection9_par = smart_select(all, "resname POPE && ! (resname SOL NA CL or resid 1 to 15)", NULL);
    select_t *selection9_par2 = smart_select(all, "! (resname SOL NA CL or resid 1 to 15) && resname POPE", NULL);
    assert(selection9_par != NULL);
    assert(selection9_par2 != NULL);
    assert(selection_compare_strict(selection9, selection9_par));
    assert(selection_compare_strict(selection9, selection9_par2));
    free(selection9);
    free(selection9_par);
    free(selection9_par2);

    select_t *selection10 = smart_select(all, "resname SOL NA CL or resid 1 to 15 or not resname POPE", NULL);
    select_t *selection10_par = smart_select(all, "! (not (resname SOL NA CL or resid 1 to 15) && resname POPE)", NULL);
    assert(selection10_par != NULL);
    assert(selection_compare(selection10, selection10_par));
    free(selection10);
    free(selection10_par);

    free(system);
    free(all);
    dict_destroy(ndx_groups);
    printf("OK\n");
}

static void test_smart_select_parentheses_fails(void)
{
    printf("%-40s", "smart_select (parentheses, fails) ");
    fflush(stdout);

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);
    dict_t *ndx_groups = read_ndx(NDX_FILE, system);

    // incorrect number of parentheses
    select_t *selection1 = smart_select(all, "((resname SOL)", NULL);
    assert(selection1 == NULL);

    select_t *selection2 = smart_select(all, "(resname SOL))", NULL);
    assert(selection2 == NULL);

    select_t *selection3 = smart_select(all, "((((resname SOL)))", NULL);
    assert(selection3 == NULL);

    select_t *selection4 = smart_select(all, "resname SOL )", NULL);
    assert(selection4 == NULL);

    select_t *selection5 = smart_select(all, " ( ( resname SOL      )  ", NULL);
    assert(selection5 == NULL);

    // incorrect spacing around operators
    select_t *selection6 = smart_select(all, "(resname POPE)&&(name P)", NULL);
    assert(selection6 == NULL);

    // nonsensical queries
    select_t *selection7 = smart_select(all, "resname POPG && ( && 1 to 44 || resname POPE)", NULL);
    assert(selection7 == NULL);

    select_t *selection8 = smart_select(all, "resname POPG && (&& 1 to 44 || resname POPE)", NULL);
    assert(selection8 == NULL);

    select_t *selection9 = smart_select(all, "(POPE && (resxdd 1 to 45 || name P)) || serial 1 to 45", ndx_groups);
    assert(selection9 == NULL);

    // characters before or after parenthesis
    select_t *selection10 = smart_select(all, "resname (POPE && resid 55 to 60) POPG", ndx_groups);
    assert(selection10 == NULL);

    select_t *selection11 = smart_select(all, "resname (POPE && resid 55 to 60)", ndx_groups);
    assert(selection11 == NULL);

    select_t *selection12 = smart_select(all, "(POPE && resid 55 to 60) serial 1 2 5", ndx_groups);
    assert(selection12 == NULL);

    free(system);
    free(all);
    dict_destroy(ndx_groups);
    printf("OK\n");
}

static void test_smart_geometry(void)
{
    printf("%-40s", "smart_geometry ");
    fflush(stdout);

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);
    dict_t *ndx_groups = read_ndx(NDX_FILE, system);

    select_t *protein = smart_select(all, "Protein", ndx_groups);

    // selecting xcylinder
    vec_t centerx = {0., 0., 0.};
    float optionsx[3] = {2.0, -1.1, 1.3};
    select_t *pope = smart_select(all, "resname POPE", NULL);
    select_t *xcylinder_abs_manual = select_geometry(pope, centerx, xcylinder, optionsx, system->box);

    center_of_geometry(protein, centerx, system->box);
    select_t *xcylinder_manual = select_geometry(pope, centerx, xcylinder, optionsx, system->box);

    select_t *xcylinder_abs = smart_geometry(all, "resname POPE", NULL, "xcylinder 2 -1.1-1.3   ", NULL, system->box);
    select_t *xcylinder = smart_geometry(all, "resname POPE", "Protein", "xcylinder 2 -1.1-1.3    ", ndx_groups, system->box);
    assert(xcylinder_abs);
    assert(xcylinder);
    assert(selection_compare_strict(xcylinder_abs_manual, xcylinder_abs));
    assert(selection_compare_strict(xcylinder_manual, xcylinder));

    free(pope);
    free(xcylinder_abs_manual);
    free(xcylinder_manual);
    free(xcylinder_abs);
    free(xcylinder);
    
    // selecting ycylinder
    vec_t centery = {0., 0., 0.};
    float optionsy[3] = {13.0, -0.5, 1.4};
    select_t *popg = smart_select(all, "resname POPG", NULL);
    select_t *ycylinder_abs_manual = select_geometry(popg, centery, ycylinder, optionsy, system->box);

    center_of_geometry(protein, centery, system->box);
    select_t *ycylinder_manual = select_geometry(popg, centery, ycylinder, optionsy, system->box);

    select_t *ycylinder_abs = smart_geometry(all, "resname POPG", NULL, "ycylinder 13.0 -0.5-1.4     ", NULL, system->box);
    select_t *ycylinder = smart_geometry(all, "resname POPG", "Protein", "ycylinder 13.0 -0.5-1.4", ndx_groups, system->box);
    assert(ycylinder_abs);
    assert(ycylinder);
    assert(selection_compare_strict(ycylinder_abs_manual, ycylinder_abs));
    assert(selection_compare_strict(ycylinder_manual, ycylinder));

    free(popg);
    free(ycylinder_abs_manual);
    free(ycylinder_manual);
    free(ycylinder_abs);
    free(ycylinder);

    // selecting zcylinder
    vec_t centerz = {0., 0., 0.};
    float optionsz[3] = {2.3, 0.7, 4.9};
    select_t *water = smart_select(all, "Water", ndx_groups);
    select_t *zcylinder_abs_manual = select_geometry(water, centerz, zcylinder, optionsz, system->box);

    center_of_geometry(protein, centerz, system->box);
    select_t *zcylinder_manual = select_geometry(water, centerz, zcylinder, optionsz, system->box);

    select_t *zcylinder_abs = smart_geometry(all, "Water", NULL, "    zcylinder 2.3 0.7-4.9", ndx_groups, system->box);
    select_t *zcylinder = smart_geometry(all, "resname SOL", "Protein", "zcylinder    2.3    0.7-4.9", ndx_groups, system->box);
    assert(zcylinder_abs);
    assert(zcylinder);
    assert(selection_compare_strict(zcylinder_abs_manual, zcylinder_abs));
    assert(selection_compare_strict(zcylinder_manual, zcylinder));

    free(water);
    free(zcylinder_abs_manual);
    free(zcylinder_manual);
    free(zcylinder_abs);
    free(zcylinder);

    // selecting sphere
    vec_t center_sphere = {0., 0., 0.};
    float options_sphere = {4.23};
    select_t *membrane = smart_select(all, "Membrane", ndx_groups);
    select_t *sphere_abs_manual = select_geometry(membrane, center_sphere, sphere, &options_sphere, system->box);

    center_sphere[0] = 2.4;
    center_sphere[1] = 3.1;
    center_sphere[2] = 7.3;
    select_t *sphere_point_manual = select_geometry(membrane, center_sphere, sphere, &options_sphere, system->box);
    select_t *sphere_point = smart_geometry(all, "Membrane", "point   2.4 3.1       7.3", "sphere 4.23", ndx_groups, system->box);

    center_of_geometry(protein, center_sphere, system->box);
    select_t *sphere_manual = select_geometry(membrane, center_sphere, sphere, &options_sphere, system->box);

    select_t *sphere_abs = smart_geometry(all, "Membrane", NULL, "sphere 4.23", ndx_groups, system->box);
    select_t *sphere = smart_geometry(all, "Membrane", "Protein", "   sphere 4.23", ndx_groups, system->box);
    assert(sphere_abs);
    assert(sphere);
    assert(sphere_point);
    assert(selection_compare_strict(sphere_abs_manual, sphere_abs));
    assert(selection_compare_strict(sphere_manual, sphere));
    assert(selection_compare_strict(sphere_point_manual, sphere_point));

    free(membrane);
    free(sphere_abs_manual);
    free(sphere_manual);
    free(sphere_abs);
    free(sphere);
    free(sphere_point_manual);
    free(sphere_point);

    // selecting box
    vec_t center_box = {0., 0., 0.};
    float options_box[6] = {-4, 4, -1.5, 2, -3.2, 2.2};
    select_t *wion = smart_select(all, "W_ION", ndx_groups);
    select_t *box_abs_manual = select_geometry(wion, center_box, box, options_box, system->box);

    center_of_geometry(protein, center_box, system->box);
    select_t *box_manual = select_geometry(wion, center_box, box, options_box, system->box);

    select_t *box_abs = smart_geometry(all, "W_ION", NULL, "box -4-4 -1.5-2 -3.2-2.2", ndx_groups, system->box);
    select_t *box = smart_geometry(all, "W_ION", "Protein", "box    -4-4   -1.5-2 -3.2-2.2   ", ndx_groups, system->box);
    assert(box_abs);
    assert(box);
    assert(selection_compare_strict(box_abs_manual, box_abs));
    assert(selection_compare_strict(box_manual, box));

    free(wion);
    free(box_abs_manual);
    free(box_manual);
    free(box_abs);
    free(box);

    free(system);
    free(all);
    free(protein);
    dict_destroy(ndx_groups);
    printf("OK\n");
}

static void test_smart_geometry_null(void)
{
    printf("%-40s", "smart_geometry (fail & null) ");
    fflush(stdout);

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);
    dict_t *ndx_groups = read_ndx(NDX_FILE, system);

    // no input selection provided
    select_t *selection1 = smart_geometry(NULL, "resname POPC", "name CA", "sphere 3.2", ndx_groups, system->box);
    assert(!selection1);
    
    // no selection query provided
    vec_t center = {0.};
    float radius = 3.2;
    select_t *ca = smart_select(all, "name CA", ndx_groups);
    center_of_geometry(ca, center, system->box);
    select_t *selection2_manual = select_geometry(all, center, sphere, &radius, system->box);
    select_t *selection2 = smart_geometry(all, NULL, "name CA", "sphere 3.2", ndx_groups, system->box);
    assert(selection2);
    assert(selection_compare_strict(selection2_manual, selection2));
    free(ca);
    free(selection2);
    free(selection2_manual);

    // no geometry query provided
    select_t *popc = smart_select(all, "resname POPC", ndx_groups);
    select_t *selection3 = smart_geometry(all, "resname POPC", "name CA", NULL, ndx_groups, system->box);
    assert(selection3);
    assert(selection_compare_strict(popc, selection3));
    free(popc);
    free(selection3);

    // no selection and geometry query provided
    select_t *selection4 = smart_geometry(all, NULL, "name CA", NULL, ndx_groups, system->box);
    assert(selection4);
    assert(selection_compare_strict(all, selection4));
    free(selection4);

    // no selection, reference query, and geometry query provided
    select_t *selection5 = smart_geometry(all, NULL, NULL, NULL, ndx_groups, system->box);
    select_t *selection5_b = smart_geometry(all, NULL, NULL, NULL, NULL, system->box);
    assert(selection5);
    assert(selection5_b);
    assert(selection_compare_strict(all, selection5));
    assert(selection_compare_strict(all, selection5_b));
    free(selection5);
    free(selection5_b);
    
    // no box
    select_t *selection6 = smart_geometry(all, "resname POPC", "name CA", "sphere 3.2", ndx_groups, NULL);
    assert(!selection6);

    // failed selection query
    select_t *selection_f1 = smart_geometry(all, "resn POPC", "name CA", "sphere 3.2", ndx_groups, system->box);
    assert(!selection_f1);
    select_t *selection_f2 = smart_geometry(all, "Membrane", "name CA", "sphere 3.2", NULL, system->box);
    assert(!selection_f2);
    
    // failed reference query
    select_t *selection_f3 = smart_geometry(all, "resname POPC", "ame CA", "sphere 3.2", ndx_groups, system->box);
    assert(!selection_f3);
    select_t *selection_f4 = smart_geometry(all, "resname POPC", "Protein", "sphere 3.2", NULL, system->box);
    assert(!selection_f4);
    select_t *selection_f5 = smart_geometry(all, "resname POPC", "name XYZ", "sphere 3.2", ndx_groups, system->box);
    assert(!selection_f5);

    // failed cylinder query
    select_t *selection_f6 = smart_geometry(all, "resname POPC", "Protein", "xcylinder 4.2 5-.", ndx_groups, system->box);
    select_t *selection_f7 = smart_geometry(all, "resname POPC", "Protein", "ycylinder 4.2 3.2", ndx_groups, system->box);
    select_t *selection_f8 = smart_geometry(all, "resname POPC", "Protein", "zcylinder 4.2", ndx_groups, system->box);
    select_t *selection_f9 = smart_geometry(all, "resname POPC", "Protein", "xcylinder", ndx_groups, system->box);
    select_t *selection_f10 = smart_geometry(all, "resname POPC", "Protein", "cylinder", ndx_groups, system->box);
    select_t *selection_f11 = smart_geometry(all, "resname POPC", "Protein", "zcylinder 4.2 3.3-4 9.2", ndx_groups, system->box);
    assert(!selection_f6);
    assert(!selection_f7);
    assert(!selection_f8);
    assert(!selection_f9);
    assert(!selection_f10);
    assert(!selection_f11);

    // failed sphere query
    select_t *selection_f12 = smart_geometry(all, "resname POPC", "Protein", "sphere", ndx_groups, system->box);
    //select_t *selection_f13 = smart_geometry(all, "resname POPC", "Protein", "sphere 2-3", ndx_groups, system->box);
    select_t *selection_f14 = smart_geometry(all, "resname POPC", "Protein", "sphere 4.2 1.8", ndx_groups, system->box);
    assert(!selection_f12);
    //assert(!selection_f13);
    assert(!selection_f14);

    // failed box query
    select_t *selection_f15 = smart_geometry(all, "resname POPC", "Protein", "box", ndx_groups, system->box);
    select_t *selection_f16 = smart_geometry(all, "resname POPC", "Protein", "box -3--2", ndx_groups, system->box);
    select_t *selection_f17 = smart_geometry(all, "resname POPC", "Protein", "box -3-2 2-3", ndx_groups, system->box);
    select_t *selection_f18 = smart_geometry(all, "resname POPC", "Protein", "box -3-2 2-3 4-5 4-6", ndx_groups, system->box);
    select_t *selection_f19 = smart_geometry(all, "resname POPC", "Protein", "box 4-5 5 7", ndx_groups, system->box);
    assert(!selection_f15);
    assert(!selection_f16);
    assert(!selection_f17);
    assert(!selection_f18);
    assert(!selection_f19);

    // failed reference query point
    select_t *selection_f20 = smart_geometry(all, "resname POPC", "point", "sphere 3.2", ndx_groups, system->box);
    select_t *selection_f21 = smart_geometry(all, "resname POPC", "point 2.2", "sphere 3.2", ndx_groups, system->box);
    select_t *selection_f22 = smart_geometry(all, "resname POPC", "point 2.2 -4.3", "sphere 3.2", ndx_groups, system->box);
    select_t *selection_f23 = smart_geometry(all, "resname POPC", "point 2.2 4.3 5.1 0.7", "sphere 3.2", ndx_groups, system->box);
    //select_t *selection_f24 = smart_geometry(all, "resname POPC", "point 2.2-4.5 7.0 1.3", "sphere 3.2", ndx_groups, system->box);
    select_t *selection_f25 = smart_geometry(all, "resname POPC", "point 1 2 F", "sphere 3.2", ndx_groups, system->box);
    assert(!selection_f20);
    assert(!selection_f21);
    assert(!selection_f22);
    assert(!selection_f23);
    //assert(!selection_f24);
    assert(!selection_f25);


    free(system);
    free(all);
    dict_destroy(ndx_groups);
    printf("OK\n");
}

void test_selection(void)
{
    test_strsplit_space();
    test_strsplit_spacetab();

    test_strstrip();
    test_strstrip_zero();
    test_strstrip_white();

    test_strremwhite();
    test_strremwhite_white();

    test_load_gro();

    test_match_residue_name();
    test_match_residue_num();
    test_match_atom_name();
    test_match_atom_num();

    test_selection_create();
    test_select_system();
    test_selection_copy();
    test_selection_empty();
    test_selection_add_atom();
    test_selection_add();

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

    test_selection_reverse();
    test_selection_slice();

    test_selection_fixres();
    test_selection_isin();
    test_selection_getnres();
    test_selection_getresnames();
    test_selection_splitbyres();
    test_selection_splitbyres_empty();
    test_selection_splitbyres_broken();

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
    test_smart_select_advanced();
    test_smart_select_to();
    test_smart_select_advanced_fails();
    test_smart_select_parentheses();
    test_smart_select_parentheses_fails();

    test_smart_geometry();
    test_smart_geometry_null();
    
}