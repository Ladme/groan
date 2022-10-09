// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

/* This file contains test functions for xtc and trr input/output. */

#include "tests.h"

/* Simple function wrapping a coordinate into a simulation box. */
static inline void wrap_coordinate(float *x, const float dimension)
{
    while (*x > dimension) *x -= dimension;
    while (*x < 0) *x += dimension;
}

void test_box_xtc2gro(void)
{
    printf("%-40s", "box_xtc2gro ");

    float box[3][3] = {{13.43, 0.0, 0.0}, {0.0, 17.2, 0.0}, {0.0, 0.0, -3.465} };

    box_t gro_box = {0.0};

    box_xtc2gro(box, gro_box);

    assert(closef(gro_box[0], 13.43, 0.00001));
    assert(closef(gro_box[1], 17.2, 0.00001));
    assert(closef(gro_box[2], -3.465, 0.00001));

    printf("OK\n");
}

void test_box_gro2xtc(void)
{
    printf("%-40s", "box_gro2xtc ");

    float box[3][3] = {{0.0}};

    box_t gro_box = {-3, 0.42, 176.256};

    box_gro2xtc(gro_box, box);

    assert(closef(box[0][0], -3, 0.00001));
    assert(closef(box[0][1], 0.0, 0.00001));
    assert(closef(box[0][2], 0.0, 0.00001));

    assert(closef(box[1][0], 0.0, 0.00001));
    assert(closef(box[1][1], 0.42, 0.00001));
    assert(closef(box[1][2], 0.0, 0.00001));

    assert(closef(box[2][0], 0.0, 0.00001));
    assert(closef(box[2][1], 0.0, 0.00001));
    assert(closef(box[2][2], 176.256, 0.00001));

    printf("OK\n");
}

void test_reset_velocities(void)
{
    printf("%-40s", "reset_velocities ");

    system_t *system = load_gro(INPUT_GRO_FILE);

    assert(vec_len(system->atoms[0].velocity) > 0.0);

    reset_velocities(system);

    assert(system->atoms[0].velocity[0] == 0.0);
    assert(system->atoms[0].velocity[1] == 0.0);
    assert(system->atoms[0].velocity[2] == 0.0);

    assert(system->atoms[10004].velocity[0] == 0.0);
    assert(system->atoms[10004].velocity[1] == 0.0);
    assert(system->atoms[10004].velocity[2] == 0.0);

    assert(system->atoms[48283].velocity[0] == 0.0);
    assert(system->atoms[48283].velocity[1] == 0.0);
    assert(system->atoms[48283].velocity[2] == 0.0);

    free(system);
    printf("OK\n");
}

void test_validate_xtc(void)
{
    printf("%-40s", "validate_xtc ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    system_t *system2 = load_gro(SMALL_GRO_FILE);

    assert(validate_xtc(INPUT_XTC_FILE, system->n_atoms));
    assert(!validate_xtc(INPUT_XTC_FILE, system2->n_atoms));

    free(system);
    free(system2);
    printf("OK\n");
}

void test_read_xtc_step_first(char *xtc_file)
{
    printf("%-40s", "read_xtc_step (first frame) ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);
    system_t *system_copy = selection_to_system_d(all, system->box, system->step, system->time);

    XDRFILE *xtc = xdrfile_open(xtc_file, "r");

    assert(read_xtc_step(xtc, system) == 0);

    assert(system->step == 0);
    assert(closef(system->time, 0.0, 0.00001));
    assert(closef(system->precision, 100.0, 0.00001));
    assert(closef(system->box[0], 7.25725, 0.00001));
    assert(closef(system->box[1], 7.25725, 0.00001));
    assert(closef(system->box[2], 9.02012, 0.00001));

    // compare with data from gro file (this is the first step)
    for (size_t i = 0; i < system->n_atoms; ++i) {

        vec_t dist = {0.0};
        calc_vector(dist, system->atoms[i].position, system_copy->atoms[i].position, system->box);

        assert(closef(dist[0], 0.0, 0.01));
        assert(closef(dist[1], 0.0, 0.01));
        assert(closef(dist[2], 0.0, 0.01));
    }

    free(system);
    free(system_copy);
    xdrfile_close(xtc);
    printf("OK\n");
}

void test_read_xtc_step_last(char *xtc_file)
{
    printf("%-40s", "read_xtc_step (last frame) ");

    system_t *system = load_gro(INPUT_GRO_FILE);

    XDRFILE *xtc = xdrfile_open(xtc_file, "r");

    while (read_xtc_step(xtc, system) == 0);

    assert(system->step == 20000);
    assert(closef(system->time, 40.0, 0.00001));
    assert(closef(system->precision, 100.0, 0.00001));
    assert(closef(system->box[0], 7.25505, 0.00001));
    assert(closef(system->box[1], 7.25505, 0.00001));
    assert(closef(system->box[2], 9.03129, 0.00001));

    assert(closef(system->atoms[0].position[0], 5.05, 0.00001));
    assert(closef(system->atoms[0].position[1], 3.82, 0.00001));
    assert(closef(system->atoms[0].position[2], 6.32, 0.00001));

    assert(closef(system->atoms[10004].position[0], 2.64, 0.00001));
    assert(closef(system->atoms[10004].position[1], 2.03, 0.00001));
    assert(closef(system->atoms[10004].position[2], 3.72, 0.00001));

    assert(closef(system->atoms[48283].position[0], 1.78, 0.00001));
    assert(closef(system->atoms[48283].position[1], 5.47, 0.00001));
    assert(closef(system->atoms[48283].position[2], 6.57, 0.00001));

    free(system);
    xdrfile_close(xtc);
    printf("OK\n");
}

void test_write_xtc_step_full(void)
{
    printf("%-40s", "write_xtc_step (full)");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    XDRFILE *xtc = xdrfile_open(INPUT_XTC_FILE, "r");
    XDRFILE *output = xdrfile_open("temporary.xtc", "w");
    while (read_xtc_step(xtc, system) == 0) {
        assert(write_xtc_step(output, all, system->step, system->time, system->box, system->precision) == 0); 
    }

    xdrfile_close(xtc);
    xdrfile_close(output);
    free(all);
    free(system);

    // test the new file
    printf("\n   >>> ");
    test_read_xtc_step_first("temporary.xtc");
    printf("   >>> ");
    test_read_xtc_step_last("temporary.xtc");

    // remove the output file
    remove("temporary.xtc");
}

void test_validate_trr(void)
{
    printf("%-40s", "validate_trr ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    system_t *system2 = load_gro(SMALL_GRO_FILE);

    assert(validate_trr(INPUT_TRR_FILE, system->n_atoms));
    assert(!validate_trr(INPUT_TRR_FILE, system2->n_atoms));

    free(system);
    free(system2);
    printf("OK\n");
}

void test_read_trr_step_first4(char *trr_file)
{
    printf("%-40s", "read_trr_step (first 4 frames) ");

    system_t *system = load_gro(INPUT_GRO_FILE);

    XDRFILE *trr = xdrfile_open(trr_file, "r");

    // general properties

    int steps[4] = {0, 2000, 2500, 4000};
    float time[4] = {0.0, 4.0, 5.0, 8.0};

    float box_x[4] = {7.257250, 7.244047, 7.243149, 7.234244};
    float box_y[4] = {7.257250, 7.244047, 7.243149, 7.234244};
    float box_z[4] = {9.020120, 9.068537, 9.039181, 9.095878};

    float lambda = 0.0;

    // properties of atom 0
    float positions_x_a0[4]  = {5.028, 5.046815, 0.0, 5.049396};
    float positions_y_a0[4]  = {3.864, 3.764930, 0.0, 3.804045};
    float positions_z_a0[4]  = {6.219, 6.355316, 0.0, 6.291129};

    float velocities_x_a0[4] = {-0.2376, 0.0, 0.0, 0.0};
    float velocities_y_a0[4] = {0.2655, 0.0, 0.0, 0.0};
    float velocities_z_a0[4] = {-0.5587, 0.0, 0.0, 0.0};

    float forces_x_a0[4]     = {23.145859, 0.0, 176.134796, 0.0};
    float forces_y_a0[4]     = {336.384064, 0.0, -243.962662, 0.0};
    float forces_z_a0[4]     = {591.736694, 0.0, -233.449295, 0.0};

    // properties of atom 10004
    float positions_x_a10004[4]  = {2.654, 2.731320, 0.0, 2.695636};
    float positions_y_a10004[4]  = {2.122, 2.109899, 0.0, 2.069323};
    float positions_z_a10004[4]  = {3.822, 3.788733, 0.0, 3.812672};

    float velocities_x_a10004[4] = {-0.053, 0.0, 0.0, 0.0};
    float velocities_y_a10004[4] = {0.0798, 0.0, 0.0, 0.0};
    float velocities_z_a10004[4] = {-0.1935, 0.0, 0.0, 0.0};

    float forces_x_a10004[4]     = {-528.195251, 0.0, 289.247955, 0.0};
    float forces_y_a10004[4]     = {-169.680695, 0.0, -994.645020, 0.0};
    float forces_z_a10004[4]     = {-474.444824, 0.0, 319.326111, 0.0};

    // properties of atom 48283
    float positions_x_a48283[4]  = {1.593, 1.573664, 0.0, 1.663715};
    float positions_y_a48283[4]  = {5.569, 5.623546, 0.0, 5.594219};
    float positions_z_a48283[4]  = {6.361, 6.352405, 0.0, 6.413816};

    float velocities_x_a48283[4] = {0.1499, 0.0, 0.0, 0.0};
    float velocities_y_a48283[4] = {-0.0234, 0.0, 0.0, 0.0};
    float velocities_z_a48283[4] = {-0.3537, 0.0, 0.0, 0.0};

    float forces_x_a48283[4]     = {-23.251177, 0.0, 610.640747, 0.0};
    float forces_y_a48283[4]     = {523.358887, 0.0, -784.378357, 0.0};
    float forces_z_a48283[4]     = {152.823563, 0.0, -352.684357, 0.0};

    for (int i = 0; i < 4; ++i) {
        assert(read_trr_step(trr, system) == 0);

        assert(system->step == steps[i]);
        assert(closef(system->time, time[i], 0.00001));
        assert(closef(system->lambda, lambda, 0.00001));

        assert(closef(system->box[0], box_x[i], 0.00001));
        assert(closef(system->box[1], box_y[i], 0.00001));
        assert(closef(system->box[2], box_z[i], 0.00001));

        atom_t *atom = &(system->atoms[0]);
        assert(closef(atom->position[0], positions_x_a0[i], 0.00001));
        assert(closef(atom->position[1], positions_y_a0[i], 0.00001));
        assert(closef(atom->position[2], positions_z_a0[i], 0.00001));

        assert(closef(atom->velocity[0], velocities_x_a0[i], 0.00001));
        assert(closef(atom->velocity[1], velocities_y_a0[i], 0.00001));
        assert(closef(atom->velocity[2], velocities_z_a0[i], 0.00001));

        assert(closef(atom->force[0], forces_x_a0[i], 0.0001));
        assert(closef(atom->force[1], forces_y_a0[i], 0.0001));
        assert(closef(atom->force[2], forces_z_a0[i], 0.0001));

        atom = &(system->atoms[10004]);
        assert(closef(atom->position[0], positions_x_a10004[i], 0.00001));
        assert(closef(atom->position[1], positions_y_a10004[i], 0.00001));
        assert(closef(atom->position[2], positions_z_a10004[i], 0.00001));

        assert(closef(atom->velocity[0], velocities_x_a10004[i], 0.00001));
        assert(closef(atom->velocity[1], velocities_y_a10004[i], 0.00001));
        assert(closef(atom->velocity[2], velocities_z_a10004[i], 0.00001));

        assert(closef(atom->force[0], forces_x_a10004[i], 0.0001));
        assert(closef(atom->force[1], forces_y_a10004[i], 0.0001));
        assert(closef(atom->force[2], forces_z_a10004[i], 0.0001));

        atom = &(system->atoms[48283]);
        assert(closef(atom->position[0], positions_x_a48283[i], 0.00001));
        assert(closef(atom->position[1], positions_y_a48283[i], 0.00001));
        assert(closef(atom->position[2], positions_z_a48283[i], 0.00001));

        assert(closef(atom->velocity[0], velocities_x_a48283[i], 0.00001));
        assert(closef(atom->velocity[1], velocities_y_a48283[i], 0.00001));
        assert(closef(atom->velocity[2], velocities_z_a48283[i], 0.00001));

        assert(closef(atom->force[0], forces_x_a48283[i], 0.0001));
        assert(closef(atom->force[1], forces_y_a48283[i], 0.0001));
        assert(closef(atom->force[2], forces_z_a48283[i], 0.0001));
    }

    free(system);
    xdrfile_close(trr);
    printf("OK\n");
}

void test_read_trr_step_last(char *trr_file)
{
    printf("%-40s", "read_trr_step (last frame) ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    XDRFILE *trr = xdrfile_open(trr_file, "r");

    float box[3] = {7.255046, 7.255046, 9.031292};

    while (read_trr_step(trr, system) == 0);

    assert(system->step == 20000);
    assert(closef(system->time, 40.0, 0.00001));
    assert(closef(system->lambda, 0.0, 0.00001));

    assert(closef(system->box[0], box[0], 0.00001));
    assert(closef(system->box[1], box[1], 0.00001));
    assert(closef(system->box[2], box[2], 0.00001));

    assert(closef(system->atoms[0].position[0], 5.047436, 0.00001));
    assert(closef(system->atoms[0].position[1], 3.817306, 0.00001));
    assert(closef(system->atoms[0].position[2], 6.319234, 0.00001));

    assert(closef(system->atoms[0].velocity[0], -0.133033, 0.00001));
    assert(closef(system->atoms[0].velocity[1], 0.503775, 0.00001));
    assert(closef(system->atoms[0].velocity[2], -0.521044, 0.00001));

    assert(closef(system->atoms[0].force[0], 1351.114624, 0.0001));
    assert(closef(system->atoms[0].force[1], 237.777359, 0.0001));
    assert(closef(system->atoms[0].force[2], -363.947998, 0.0001));


    assert(closef(system->atoms[10004].position[0], 2.642371, 0.00001));
    assert(closef(system->atoms[10004].position[1], 2.033608, 0.00001));
    assert(closef(system->atoms[10004].position[2], 3.723695, 0.00001));

    assert(closef(system->atoms[10004].velocity[0], -0.521758, 0.00001));
    assert(closef(system->atoms[10004].velocity[1], 0.176894, 0.00001));
    assert(closef(system->atoms[10004].velocity[2], -0.014081, 0.00001));

    assert(closef(system->atoms[10004].force[0], 123.242401, 0.0001));
    assert(closef(system->atoms[10004].force[1], 366.518158, 0.0001));
    assert(closef(system->atoms[10004].force[2], 617.743469, 0.0001));


    assert(closef(system->atoms[48283].position[0], 1.778632, 0.00001));
    assert(closef(system->atoms[48283].position[1], 5.469856, 0.00001));
    assert(closef(system->atoms[48283].position[2], 6.566182, 0.00001));

    assert(closef(system->atoms[48283].velocity[0], 0.129547, 0.00001));
    assert(closef(system->atoms[48283].velocity[1], 0.349285, 0.00001));
    assert(closef(system->atoms[48283].velocity[2], 0.072853, 0.00001));

    assert(closef(system->atoms[48283].force[0], 421.829041, 0.0001));
    assert(closef(system->atoms[48283].force[1], 330.728302, 0.0001));
    assert(closef(system->atoms[48283].force[2], -322.925568, 0.0001));

    free(system);
    xdrfile_close(trr);
    printf("OK\n");
}

void test_read_trr_step_compare(char *trr_file)
{
    printf("%-40s", "read_trr_step (full, compare with xtc) ");

    system_t *system_trr = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system_trr);
    system_t *system_xtc = selection_to_system_d(all, system_trr->box, system_trr->step, system_trr->time);

    XDRFILE *trr = xdrfile_open(trr_file, "r");
    XDRFILE *xtc = xdrfile_open(INPUT_XTC_FILE, "r");

    read_xtc_step(xtc, system_xtc);
    while (read_trr_step(trr, system_trr) == 0) {
        if (system_trr->step % 2000 == 0) {
            while (system_xtc->step != system_trr->step) {
                read_xtc_step(xtc, system_xtc);
            }

            for (size_t i = 0; i < system_trr->n_atoms; ++i) {
                assert(closef(system_trr->atoms[i].position[0], system_xtc->atoms[i].position[0], 0.01));
                assert(closef(system_trr->atoms[i].position[1], system_xtc->atoms[i].position[1], 0.01));
                assert(closef(system_trr->atoms[i].position[2], system_xtc->atoms[i].position[2], 0.01));
            }
        }
        
    }

    free(system_trr);
    free(system_xtc);
    xdrfile_close(trr);
    xdrfile_close(xtc);
    printf("OK\n");
}

void test_write_trr_step_full(void)
{
    printf("%-40s", "write_trr_step (full) ");

    system_t *system = load_gro(INPUT_GRO_FILE);
    select_t *all = select_system(system);

    XDRFILE *trr = xdrfile_open(INPUT_TRR_FILE, "r");
    XDRFILE *output = xdrfile_open("temporary.trr", "w");
    while (read_trr_step(trr, system) == 0) {
        assert(write_trr_step(output, all, system->step, system->time, system->box, system->lambda) == 0); 
    }

    xdrfile_close(trr);
    xdrfile_close(output);
    free(all);
    free(system);

    // test the new file
    printf("\n   >>> ");
    test_read_trr_step_first4("temporary.trr");
    printf("   >>> ");
    test_read_trr_step_last("temporary.trr");
    printf("   >>> ");
    test_read_trr_step_compare("temporary.trr");

    // remove the output file
    remove("temporary.trr");
}

void test_xdr(void)
{
    test_box_xtc2gro();
    test_box_gro2xtc();
    test_reset_velocities();

    test_validate_xtc();
    test_read_xtc_step_first(INPUT_XTC_FILE);
    test_read_xtc_step_last(INPUT_XTC_FILE);
    test_write_xtc_step_full();

    test_validate_trr();
    test_read_trr_step_first4(INPUT_TRR_FILE);
    test_read_trr_step_last(INPUT_TRR_FILE);
    test_read_trr_step_compare(INPUT_TRR_FILE);
    test_write_trr_step_full();
}