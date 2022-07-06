// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#include "tests.h"

/* Simple function wrapping a coordinate into a simulation box. */
static inline void wrap_coordinate(float *x, const float dimension)
{
    while (*x > dimension) *x -= dimension;
    while (*x < 0) *x += dimension;
}

void test_distance1D(void)
{
    printf("%-40s", "distance1D ");

    system_t *system = load_gro(INPUT_GRO_FILE);

    atom_t *a1 = &system->atoms[16094];
    atom_t *a2 = &system->atoms[16095];

    assert(closef(distance1D(a1->position, a2->position, x, system->box),  0.103, 0.000001));
    assert(closef(distance1D(a1->position, a2->position, y, system->box),  0.098 , 0.000001));
    assert(closef(distance1D(a1->position, a2->position, z, system->box), -0.042 , 0.000001));

    a1 = &system->atoms[7712];
    a2 = &system->atoms[7211];

    assert(closef(distance1D(a1->position, a2->position, x, system->box), -1.56125, 0.000001));
    assert(closef(distance1D(a1->position, a2->position, y, system->box),  1.23925 , 0.000001));
    assert(closef(distance1D(a1->position, a2->position, z, system->box),  0.527 , 0.000001));

    a1 = &system->atoms[24569];
    a2 = &system->atoms[42344];

    assert(closef(distance1D(a1->position, a2->position, x, system->box), -3.093, 0.000001));
    assert(closef(distance1D(a1->position, a2->position, y, system->box),  2.85325, 0.000001));
    assert(closef(distance1D(a1->position, a2->position, z, system->box),  3.95112, 0.000001));

    free(system);
    printf("OK\n");
}

void test_distance2D(void)
{
    printf("%-40s", "distance2D ");

    system_t *system = load_gro(INPUT_GRO_FILE);

    atom_t *a1 = &system->atoms[16094];
    atom_t *a2 = &system->atoms[16095];

    /*printf("%f\n", distance2D(a1->position, a2->position, xy, system->box));
    printf("%f\n", distance2D(a1->position, a2->position, yz, system->box));
    printf("%f\n", distance2D(a1->position, a2->position, xz, system->box));*/

    assert(closef(distance2D(a1->position, a2->position, xy, system->box),  0.142172, 0.000001));
    assert(closef(distance2D(a1->position, a2->position, yz, system->box),  0.106621, 0.000001));
    assert(closef(distance2D(a1->position, a2->position, xz, system->box),  0.111234, 0.000001));

    a1 = &system->atoms[7712];
    a2 = &system->atoms[7211];

    /*printf("%f\n", distance2D(a1->position, a2->position, xy, system->box));
    printf("%f\n", distance2D(a1->position, a2->position, yz, system->box));
    printf("%f\n", distance2D(a1->position, a2->position, xz, system->box));*/

    assert(closef(distance2D(a1->position, a2->position, xy, system->box),  1.993299, 0.000001));
    assert(closef(distance2D(a1->position, a2->position, yz, system->box),  1.346651, 0.000001));
    assert(closef(distance2D(a1->position, a2->position, xz, system->box),  1.647795, 0.000001));

    a1 = &system->atoms[24569];
    a2 = &system->atoms[42344];

    /*printf("%f\n", distance2D(a1->position, a2->position, xy, system->box));
    printf("%f\n", distance2D(a1->position, a2->position, yz, system->box));
    printf("%f\n", distance2D(a1->position, a2->position, xz, system->box));*/

    assert(closef(distance2D(a1->position, a2->position, xy, system->box),  4.208050, 0.000001));
    assert(closef(distance2D(a1->position, a2->position, yz, system->box),  4.873641, 0.000001));
    assert(closef(distance2D(a1->position, a2->position, xz, system->box),  5.017768, 0.000001));


    free(system);
    printf("OK\n");
}

void test_distance2D_naive(void)
{
    printf("%-40s", "distance2D_naive ");

    system_t *system = load_gro(INPUT_GRO_FILE);

    atom_t *a1 = &system->atoms[16094];
    atom_t *a2 = &system->atoms[16095];

    /*printf("%f\n", distance2D_naive(a1->position, a2->position, xy));
    printf("%f\n", distance2D_naive(a1->position, a2->position, yz));
    printf("%f\n", distance2D_naive(a1->position, a2->position, xz));*/

    assert(closef(distance2D_naive(a1->position, a2->position, xy),  0.142172, 0.000001));
    assert(closef(distance2D_naive(a1->position, a2->position, yz),  0.106621, 0.000001));
    assert(closef(distance2D_naive(a1->position, a2->position, xz),  0.111234, 0.000001));

    a1 = &system->atoms[7712];
    a2 = &system->atoms[7211];

    /*printf("%f\n", distance2D_naive(a1->position, a2->position, xy));
    printf("%f\n", distance2D_naive(a1->position, a2->position, yz));
    printf("%f\n", distance2D_naive(a1->position, a2->position, xz));*/

    assert(closef(distance2D_naive(a1->position, a2->position, xy),  8.286178, 0.000001));
    assert(closef(distance2D_naive(a1->position, a2->position, yz),  6.041031, 0.000001));
    assert(closef(distance2D_naive(a1->position, a2->position, xz),  5.720327, 0.000001));

    a1 = &system->atoms[24569];
    a2 = &system->atoms[42344];

    /*printf("%f\n", distance2D_naive(a1->position, a2->position, xy));
    printf("%f\n", distance2D_naive(a1->position, a2->position, yz));
    printf("%f\n", distance2D_naive(a1->position, a2->position, xz));*/

    assert(closef(distance2D_naive(a1->position, a2->position, xy),  5.381623, 0.000001));
    assert(closef(distance2D_naive(a1->position, a2->position, yz),  6.714907, 0.000001));
    assert(closef(distance2D_naive(a1->position, a2->position, xz),  5.938131, 0.000001));


    free(system);
    printf("OK\n");
}

void test_distance3D(void)
{
    printf("%-40s", "distance3D ");

    system_t *system = load_gro(INPUT_GRO_FILE);

    atom_t *a1 = &system->atoms[16094];
    atom_t *a2 = &system->atoms[16095];

    //printf("%f\n", distance3D(a1->position, a2->position, system->box));

    assert(closef(distance3D(a1->position, a2->position, system->box), 0.148246, 0.000001));

    a1 = &system->atoms[7712];
    a2 = &system->atoms[7211];

    //printf("%f\n", distance3D(a1->position, a2->position, system->box));

    assert(closef(distance3D(a1->position, a2->position, system->box), 2.061788, 0.000001));

    a1 = &system->atoms[24569];
    a2 = &system->atoms[42344];

    //printf("%f\n", distance3D(a1->position, a2->position, system->box));

    assert(closef(distance3D(a1->position, a2->position, system->box), 5.772264, 0.000001));

    free(system);
    printf("OK\n");
}

void test_distance3D_naive(void)
{
    printf("%-40s", "distance3D_naive ");

    system_t *system = load_gro(INPUT_GRO_FILE);

    atom_t *a1 = &system->atoms[16094];
    atom_t *a2 = &system->atoms[16095];

    //printf("%f\n", distance3D_naive(a1->position, a2->position));

    assert(closef(distance3D_naive(a1->position, a2->position), 0.148246, 0.000001));

    a1 = &system->atoms[7712];
    a2 = &system->atoms[7211];

    //printf("%f\n", distance3D_naive(a1->position, a2->position));

    assert(closef(distance3D_naive(a1->position, a2->position), 8.302918, 0.000001));

    a1 = &system->atoms[24569];
    a2 = &system->atoms[42344];

    //printf("%f\n", distance3D_naive(a1->position, a2->position));

    assert(closef(distance3D_naive(a1->position, a2->position), 7.393012, 0.000001));

    free(system);
    printf("OK\n");
}

void test_selection_translate(void)
{
    printf("%-40s", "selection_translate ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    system_t *original_system = selection_to_system(all, system->box, 0, 0.f);

    select_t *protein = select_atoms(all, "LEU SER", &match_residue_name);
    select_t *membrane = select_atoms(all, "POPE POPG", &match_residue_name);
    select_t *water = select_atoms(all, "SOL", &match_residue_name);

    vec_t translate1  = {74.5, 3.0, -0.5};
    vec_t translateP  = {0.0, -2.2, 8.1};
    vec_t translateM  = {-3.1, 0.45, -4.2};
    vec_t translatePW = {2.1, -0.3, -43.8};

    selection_translate(all, translate1, system->box);
    selection_translate(protein, translateP, system->box);
    selection_translate(membrane, translateM, system->box);
    selection_translate(protein, translatePW, system->box);
    selection_translate(water, translatePW, system->box);

    for (size_t i = 0; i < system->n_atoms; ++i) {
        float predicted_x = original_system->atoms[i].position[0] + translate1[0];
        float predicted_y = original_system->atoms[i].position[1] + translate1[1];
        float predicted_z = original_system->atoms[i].position[2] + translate1[2];

        wrap_coordinate(&predicted_x, system->box[0]);
        wrap_coordinate(&predicted_y, system->box[1]);
        wrap_coordinate(&predicted_z, system->box[2]);

        atom_t *atom = &system->atoms[i];

        
        if (selection_isin(protein, atom)) {
            predicted_x += translateP[0];
            predicted_y += translateP[1];
            predicted_z += translateP[2];

            wrap_coordinate(&predicted_x, system->box[0]);
            wrap_coordinate(&predicted_y, system->box[1]);
            wrap_coordinate(&predicted_z, system->box[2]);

            predicted_x += translatePW[0];
            predicted_y += translatePW[1];
            predicted_z += translatePW[2];

            wrap_coordinate(&predicted_x, system->box[0]);
            wrap_coordinate(&predicted_y, system->box[1]);
            wrap_coordinate(&predicted_z, system->box[2]);
        }

        if (selection_isin(membrane, atom)) {
            predicted_x += translateM[0];
            predicted_y += translateM[1];
            predicted_z += translateM[2];

            wrap_coordinate(&predicted_x, system->box[0]);
            wrap_coordinate(&predicted_y, system->box[1]);
            wrap_coordinate(&predicted_z, system->box[2]);
        }

        if (selection_isin(water, atom)) {

            predicted_x += translatePW[0];
            predicted_y += translatePW[1];
            predicted_z += translatePW[2];

            wrap_coordinate(&predicted_x, system->box[0]);
            wrap_coordinate(&predicted_y, system->box[1]);
            wrap_coordinate(&predicted_z, system->box[2]);
        }

        /*printf("%f vs %f\n", atom->position[0], predicted_x);
        printf("%f vs %f\n", atom->position[1], predicted_y);
        printf("%f vs %f\n", atom->position[2], predicted_z);*/

        assert(closef(atom->position[0], predicted_x, 0.000001f));
        assert(closef(atom->position[1], predicted_y, 0.000001f));
        assert(closef(atom->position[2], predicted_z, 0.000001f));
    } 

    free(protein);
    free(membrane);
    free(water);
    free(all);
    free(system);
    free(original_system);
    printf("OK\n");
    
}

void test_center_of_geometry(void)
{
    printf("%-40s", "center_of_geometry ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *protein = select_atoms(all, "LEU SER", &match_residue_name);
    select_t *membrane = select_atoms(all, "POPE POPG", &match_residue_name);
    select_t *water = select_atoms_d(all, "SOL", &match_residue_name);

    vec_t protein_center = {0.f};
    vec_t membrane_center = {0.f};
    vec_t water_center = {0.f};

    center_of_geometry(protein, protein_center, system->box);
    center_of_geometry(membrane, membrane_center, system->box);
    center_of_geometry(water, water_center, system->box);

    /*printf("%f %f %f\n", protein_center[0], protein_center[1], protein_center[2]);
    printf("%f %f %f\n", membrane_center[0], membrane_center[1], membrane_center[2]);
    printf("%f %f %f\n", water_center[0], water_center[1], water_center[2]);*/

    assert(closef(protein_center[0], 3.443204, 0.000001));
    assert(closef(protein_center[1], 3.718657, 0.000001));
    assert(closef(protein_center[2], 6.074358, 0.000001));

    assert(closef(membrane_center[0], 6.906639, 0.000001));
    assert(closef(membrane_center[1], 1.015586, 0.000001));
    assert(closef(membrane_center[2], 4.253287, 0.000001));

    assert(closef(water_center[0], 3.241589, 0.000001));
    assert(closef(water_center[1], 6.776079, 0.000001));
    assert(closef(water_center[2], 8.794089, 0.000001));

    free(protein);
    free(membrane);
    free(water);
    free(system);
    printf("OK\n");
}

void test_center_of_geometry_translated(void)
{
    printf("%-40s", "center_of_geometry (translated) ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *protein = select_atoms(all, "LEU SER", &match_residue_name);
    select_t *membrane = select_atoms(all, "POPE POPG", &match_residue_name);
    select_t *water = select_atoms(all, "SOL", &match_residue_name);

    vec_t translate1  = {74.5, 3.0, -0.5};
    vec_t translateP  = {0.0, -2.2, 8.1};
    vec_t translateM  = {-3.1, 0.45, -4.2};
    vec_t translatePW = {2.1, -0.3, -43.8};

    selection_translate(all, translate1, system->box);
    selection_translate(protein, translateP, system->box);
    selection_translate(membrane, translateM, system->box);
    selection_translate(protein, translatePW, system->box);
    selection_translate(water, translatePW, system->box);

    vec_t protein_center = {0.f};
    vec_t membrane_center = {0.f};
    vec_t water_center = {0.f};

    center_of_geometry(protein, protein_center, system->box);
    center_of_geometry(membrane, membrane_center, system->box);
    center_of_geometry(water, water_center, system->box);

    /*printf("%f %f %f\n", protein_center[0], protein_center[1], protein_center[2]);
    printf("%f %f %f\n", membrane_center[0], membrane_center[1], membrane_center[2]);
    printf("%f %f %f\n", water_center[0], water_center[1], water_center[2]);*/

    assert(closef(protein_center[0], 0.213467, 0.000001));
    assert(closef(protein_center[1], 4.218657, 0.000001));
    assert(closef(protein_center[2], 5.954837, 0.000001));

    assert(closef(membrane_center[0], 5.734128, 0.000001));
    assert(closef(membrane_center[1], 4.465578, 0.000001));
    assert(closef(membrane_center[2], 8.573410, 0.000001));

    assert(closef(water_center[0], 0.011871, 0.000001));
    assert(closef(water_center[1], 2.218873, 0.000001));
    assert(closef(water_center[2], 0.574566, 0.000001));

    free(all);
    free(protein);
    free(membrane);
    free(water);
    free(system);
    printf("OK\n");
}

void test_center_of_geometry_naive(void)
{
    printf("%-40s", "center_of_geometry_naive ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *protein = select_atoms(all, "LEU SER", &match_residue_name);
    select_t *membrane = select_atoms(all, "POPE POPG", &match_residue_name);
    select_t *water = select_atoms_d(all, "SOL", &match_residue_name);

    vec_t protein_center = {0.f};
    vec_t membrane_center = {0.f};
    vec_t water_center = {0.f};

    center_of_geometry_naive(protein, protein_center);
    center_of_geometry_naive(membrane, membrane_center);
    center_of_geometry_naive(water, water_center);

    /*printf("%f %f %f\n", protein_center[0], protein_center[1], protein_center[2]);
    printf("%f %f %f\n", membrane_center[0], membrane_center[1], membrane_center[2]);
    printf("%f %f %f\n", water_center[0], water_center[1], water_center[2]);*/

    assert(closef(protein_center[0], 3.443604, 0.000001));
    assert(closef(protein_center[1], 3.718559, 0.000001));
    assert(closef(protein_center[2], 6.074302, 0.000001));

    assert(closef(membrane_center[0], 3.649805, 0.000001));
    assert(closef(membrane_center[1], 3.490983, 0.000001));
    assert(closef(membrane_center[2], 4.258356, 0.000001));

    assert(closef(water_center[0], 3.617002, 0.000001));
    assert(closef(water_center[1], 3.628748, 0.000001));
    assert(closef(water_center[2], 4.680107, 0.000001));

    free(protein);
    free(membrane);
    free(water);
    free(system);
    printf("OK\n");
}


void test_analysis_tools(void)
{
    test_distance1D();
    test_distance2D();
    test_distance2D_naive();
    test_distance3D();
    test_distance3D_naive();

    test_selection_translate();

    test_center_of_geometry();
    test_center_of_geometry_translated();
    test_center_of_geometry_naive();
}