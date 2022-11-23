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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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

void test_vector_artificial(void)
{
    printf("%-40s", "calc_vector (artificial) ");
    fflush(stdout);

    // no minimum image convention
    vec_t particle1 = {4., 4., 5.};
    vec_t particle2 = {5., 5., 3.};
    vec_t result = {0.0};

    box_t box = {10.0, 10.0, 10.0, 
                 0.0, 0.0, 0.0,
                 0.0, 0.0, 0.0};
    
    calc_vector(result, particle1, particle2, box);

    assert(closef(result[0], 1.0, 0.000001));
    assert(closef(result[1], 1.0, 0.000001));
    assert(closef(result[2],-2.0, 0.000001));

    // image of the particle is closer in one dimensions
    vec_t particle3 = {3., 0., 7.};
    vec_t particle4 = {1., 2., 1.};

    calc_vector(result, particle3, particle4, box);

    assert(closef(result[0],-2.0, 0.000001));
    assert(closef(result[1], 2.0, 0.000001));
    assert(closef(result[2], 4.0, 0.000001));

    // image of the particle is closer in two dimensions
    vec_t particle5 = {1., 2., 5.};
    vec_t particle6 = {9., 8., 6.};

    calc_vector(result, particle5, particle6, box);

    assert(closef(result[0],-2.0, 0.000001));
    assert(closef(result[1],-4.0, 0.000001));
    assert(closef(result[2], 1.0, 0.000001));

    // image of the particle is closer in all three dimensions
    vec_t particle7 = {8., 9., 2.};
    vec_t particle8 = {1., 3., 9.};

    calc_vector(result, particle7, particle8, box);

    assert(closef(result[0], 3.0, 0.000001));
    assert(closef(result[1], 4.0, 0.000001));
    assert(closef(result[2],-3.0, 0.000001));

    // particles are at the same position
    vec_t particle9  = {0.0, 3.0, 10.0};
    vec_t particle10 = {10.0, 3.0, 0.0};

    calc_vector(result, particle9, particle10, box);

    assert(closef(result[0], 0.0, 0.000001));
    assert(closef(result[1], 0.0, 0.000001));
    assert(closef(result[2], 0.0, 0.000001));

    // real particle and its image are equidistant from the other particle
    vec_t particle11 = {7., 4., 3.};
    vec_t particle12 = {2., 5., 2.};

    calc_vector(result, particle11, particle12, box);

    //printf("%f %f %f\n", result[0], result[1], result[2]);

    assert(closef(result[0], 5.0, 0.000001) || closef(result[0],-5.0, 0.000001));
    assert(closef(result[1], 1.0, 0.000001));
    assert(closef(result[2],-1.0, 0.000001));

    printf("OK\n");
}

void test_vector_system(void)
{
    printf("%-40s", "calc_vector (system) ");
    fflush(stdout);

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    // no effect of the PBC
    select_t *particle1 = smart_select(all, "serial 43243", NULL);
    select_t *particle2 = smart_select(all, "serial 43250", NULL);

    atom_t *atom1 = particle1->atoms[0];
    atom_t *atom2 = particle2->atoms[0];

    //printf("atom1: %f %f %f\n", atom1->position[0], atom1->position[1], atom1->position[2]);
    //printf("atom2: %f %f %f\n", atom2->position[0], atom2->position[1], atom2->position[2]);

    free(particle1);
    free(particle2);

    vec_t result = {0.};

    calc_vector(result, atom1->position, atom2->position, system->box);

    assert(closef(result[0], atom2->position[0] - atom1->position[0], 0.000001));
    assert(closef(result[1], atom2->position[1] - atom1->position[1], 0.000001));
    assert(closef(result[2], atom2->position[2] - atom1->position[2], 0.000001));

    // PBC do have an effect
    select_t *particle3 = smart_select(all, "serial 48263", NULL);

    atom_t *atom3 = particle3->atoms[0];

    //printf("atom3: %f %f %f\n", atom3->position[0], atom3->position[1], atom3->position[2]);
    //printf("box: %f %f %f\n", system->box[0], system->box[1], system->box[2]);

    free(particle3);

    calc_vector(result, atom1->position, atom3->position, system->box);

    //printf("result: %f %f %f\n", result[0], result[1], result[2]);

    assert(closef(result[0], 2.58525, 0.000001));
    assert(closef(result[1],-0.289, 0.000001));
    assert(closef(result[2], 2.169119, 0.000001));

    free(all);
    free(system);

    printf("OK\n");
}

void test_selection_translate(void)
{
    printf("%-40s", "selection_translate ");
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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
    fflush(stdout);

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

void test_rotate_point(void)
{
    printf("%-40s", "rotate_point ");
    fflush(stdout);

    vec_t reference = {2.0, 1.0, 3.0};

    // rotation around the x axis
    vec_t point_rx = {4.0, 3.0, 2.0};
    rotate_point(point_rx, reference, 44, x);

    assert(closef(point_rx[0], 4.0000, 0.0001));
    assert(closef(point_rx[1], 3.1333, 0.0001));
    assert(closef(point_rx[2], 3.6699, 0.0001));

    // rotation around the y axis
    vec_t point_ry = {4.0, 3.0, 2.0};
    rotate_point(point_ry, reference, 63, y);

    assert(closef(point_ry[0], 3.7989, 0.0001));
    assert(closef(point_ry[1], 3.0000, 0.0001));
    assert(closef(point_ry[2], 4.3280, 0.0001));

    // rotation around the z axis
    vec_t point_rz = {4.0, 3.0, 2.0};
    rotate_point(point_rz, reference, -17, z);

    assert(closef(point_rz[0], 4.4973, 0.0001));
    assert(closef(point_rz[1], 2.3278, 0.0001));
    assert(closef(point_rz[2], 2.0000, 0.0001));

    printf("OK\n");
}

void test_selection_rotate_naive(void)
{
    printf("%-40s", "selection_rotate_naive ");
    fflush(stdout);

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    system_t *new_system = selection_to_system(all, system->box, system->step, system->time);
    select_t *new_all = select_system(new_system);

    vec_t origin1 = {2.4, 3.2, 1.8};
    vec_t origin2 = {1.8, 3.1, -1.4};
    selection_rotate_naive(new_all, origin1, 543, x);
    selection_rotate_naive(new_all, origin2, 13, y);
    selection_rotate_naive(new_all, origin1, -32, z);

    for (size_t i = 0; i < new_all->n_atoms; ++i) {
        vec_t opos = {all->atoms[i]->position[0], 
                      all->atoms[i]->position[1],
                      all->atoms[i]->position[2]};

        rotate_point(opos, origin1, 543, x);
        rotate_point(opos, origin2, 13, y);
        rotate_point(opos, origin1, -32, z);
        
        assert(closef(opos[0], new_all->atoms[i]->position[0], 0.000001));
        assert(closef(opos[1], new_all->atoms[i]->position[1], 0.000001));
        assert(closef(opos[2], new_all->atoms[i]->position[2], 0.000001));
    }


    free(system);
    free(all);
    free(new_system);
    free(new_all);
    printf("OK\n");
}

void test_selection_rotate(void)
{
    printf("%-40s", "selection_rotate ");
    fflush(stdout);

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    system_t *new_system = selection_to_system(all, system->box, system->step, system->time);
    select_t *new_all = select_system(new_system);

    vec_t origin1 = {2.4, 3.2, 1.8};
    vec_t origin2 = {1.8, 3.1, -1.4};
    selection_rotate(new_all, origin1, 543, x, system->box);
    selection_rotate(new_all, origin2, 13, y, system->box);
    selection_rotate(new_all, origin1, -32, z, system->box);

    for (size_t i = 0; i < new_all->n_atoms; ++i) {
        vec_t opos = {all->atoms[i]->position[0], 
                      all->atoms[i]->position[1],
                      all->atoms[i]->position[2]};

        rotate_point(opos, origin1, 543, x);
        wrap_coordinate(&opos[0], system->box[0]);
        wrap_coordinate(&opos[1], system->box[1]);
        wrap_coordinate(&opos[2], system->box[2]);
        rotate_point(opos, origin2, 13, y);
        wrap_coordinate(&opos[0], system->box[0]);
        wrap_coordinate(&opos[1], system->box[1]);
        wrap_coordinate(&opos[2], system->box[2]);
        rotate_point(opos, origin1, -32, z);
        wrap_coordinate(&opos[0], system->box[0]);
        wrap_coordinate(&opos[1], system->box[1]);
        wrap_coordinate(&opos[2], system->box[2]);
        
        assert(closef(opos[0], new_all->atoms[i]->position[0], 0.000001));
        assert(closef(opos[1], new_all->atoms[i]->position[1], 0.000001));
        assert(closef(opos[2], new_all->atoms[i]->position[2], 0.000001));
    }


    free(system);
    free(all);
    free(new_system);
    free(new_all);
    printf("OK\n");
}

void test_calc_angle(void)
{
    printf("%-40s", "calc_angle ");
    fflush(stdout);

    vec_t vec1 = {5.0, 7.0, 3.0};
    vec_t vec2 = {3.0, 4.0, 2.0};
    assert(closef(calc_angle(vec1, vec2), 2.8618, 0.0001));

    vec_t vec3 = {5.3, 7.1, 3.8};
    vec_t vec4 = {3.1, 4.9, 2.1};
    assert(closef(calc_angle(vec3, vec4), 5.2744, 0.0001));

    vec_t vec5 = {4.2, -3.5, 6.0};
    vec_t vec6 = {-2.0, 1.12, 4.7};
    assert(closef(calc_angle(vec5, vec6), 68.0303, 0.0001));

    vec_t vec7 = {1.0, 2.0, 5.0};
    vec_t vec8 = {2.0, 4.0, 10.0};
    assert(closef(calc_angle(vec7, vec8), 0.0000, 0.0001));

    printf("OK\n");
}

void test_selection_sort_by_dist(void)
{
    printf("%-40s", "selection_sort_by_dist ");
    fflush(stdout);

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    select_t *phosphates = smart_select(all, "name P", NULL);

    vec_t center = {0.0};
    center_of_geometry_naive(phosphates, center);

    // XYZ
    select_t *phosphates_sorted_xyz = selection_copy(phosphates);
    selection_sort_by_dist(phosphates_sorted_xyz, center, dimensionality_xyz, system->box);

    assert(selection_compare(phosphates, phosphates_sorted_xyz));
    for (size_t i = 0; i < phosphates->n_atoms - 1; ++i) {
        assert( distance3D(phosphates_sorted_xyz->atoms[i]->position, center, system->box)
                <= distance3D(phosphates_sorted_xyz->atoms[i + 1]->position, center, system->box));
    }

    free(phosphates_sorted_xyz);

    // XY
    select_t *phosphates_sorted_xy = selection_copy(phosphates);
    selection_sort_by_dist(phosphates_sorted_xy, center, dimensionality_xy, system->box);
    assert(selection_compare(phosphates, phosphates_sorted_xy));
    for (size_t i = 0; i < phosphates->n_atoms - 1; ++i) {
        assert( distance2D(phosphates_sorted_xy->atoms[i]->position, center, xy, system->box)
                <= distance2D(phosphates_sorted_xy->atoms[i + 1]->position, center, xy, system->box));
    }

    free(phosphates_sorted_xy);

    // XZ
    select_t *phosphates_sorted_xz = selection_copy(phosphates);
    selection_sort_by_dist(phosphates_sorted_xz, center, dimensionality_xz, system->box);
    assert(selection_compare(phosphates, phosphates_sorted_xz));
    for (size_t i = 0; i < phosphates->n_atoms - 1; ++i) {
        assert( distance2D(phosphates_sorted_xz->atoms[i]->position, center, xz, system->box)
                <= distance2D(phosphates_sorted_xz->atoms[i + 1]->position, center, xz, system->box));
    }

    free(phosphates_sorted_xz);

    // YZ
    select_t *phosphates_sorted_yz = selection_copy(phosphates);
    selection_sort_by_dist(phosphates_sorted_yz, center, dimensionality_yz, system->box);
    assert(selection_compare(phosphates, phosphates_sorted_yz));
    for (size_t i = 0; i < phosphates->n_atoms - 1; ++i) {
        assert( distance2D(phosphates_sorted_yz->atoms[i]->position, center, yz, system->box)
                <= distance2D(phosphates_sorted_yz->atoms[i + 1]->position, center, yz, system->box));
    }

    free(phosphates_sorted_yz);

    // X
    select_t *phosphates_sorted_x = selection_copy(phosphates);
    selection_sort_by_dist(phosphates_sorted_x, center, dimensionality_x, system->box);
    assert(selection_compare(phosphates, phosphates_sorted_x));
    for (size_t i = 0; i < phosphates->n_atoms - 1; ++i) {
        assert( fabsf(distance1D(phosphates_sorted_x->atoms[i]->position, center, x, system->box))
                <= fabsf(distance1D(phosphates_sorted_x->atoms[i + 1]->position, center, x, system->box)));
    }

    free(phosphates_sorted_x);

    // Y
    select_t *phosphates_sorted_y = selection_copy(phosphates);
    selection_sort_by_dist(phosphates_sorted_y, center, dimensionality_y, system->box);
    assert(selection_compare(phosphates, phosphates_sorted_y));
    for (size_t i = 0; i < phosphates->n_atoms - 1; ++i) {
        assert( fabsf(distance1D(phosphates_sorted_y->atoms[i]->position, center, y, system->box))
                <= fabsf(distance1D(phosphates_sorted_y->atoms[i + 1]->position, center, y, system->box)));
    }

    free(phosphates_sorted_y);

    // Z
    select_t *phosphates_sorted_z = selection_copy(phosphates);
    selection_sort_by_dist(phosphates_sorted_z, center, dimensionality_z, system->box);
    assert(selection_compare(phosphates, phosphates_sorted_z));
    for (size_t i = 0; i < phosphates->n_atoms - 1; ++i) {
        assert( fabsf(distance1D(phosphates_sorted_z->atoms[i]->position, center, z, system->box))
                <= fabsf(distance1D(phosphates_sorted_z->atoms[i + 1]->position, center, z, system->box)));
    }

    free(phosphates_sorted_z);

    free(phosphates);
    free(all);
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

    test_vector_artificial();
    test_vector_system();

    test_selection_translate();

    test_center_of_geometry();
    test_center_of_geometry_translated();
    test_center_of_geometry_naive();

    test_rotate_point();
    test_selection_rotate_naive();
    test_selection_rotate();

    test_calc_angle();
    
    test_selection_sort_by_dist();
}