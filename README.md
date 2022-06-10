# GROmacs ANalaysis (groan)

C library for analysis of Gromacs simulations.

## How to install

Just run `make`.

## How to use in your C programs

Include `groan.h` in your code and link with `-lm -lgroan`.

## Example

`
#include <groan.h>
#include <stdio.h>

int main(void)
{
    // read gro file
    system_t *system = load_gro("md.gro");
    if (system == NULL) return 1;
    
    // select atoms of membrane lipids
    size_t *membrane = NULL;
    size_t n_membrane = select_atoms(system, NULL, 0, &membrane, "POPC", &match_residue_name);
    
    // calculate membrane center
    vec_t membrane_center = {0, 0, 0};
    center_of_geometry(system, membrane, n_membrane, membrane_center, NULL);

    // print membrane center
    printf("Membrane center: %f %f %f\n", membrane_center[0], membrane_center[1], membrane_center[2]);
    
    free(membrane);
    free(system);
    
    return 0;
}
`

Compile with `gcc {YOUR_CODE}.c -L{PATH_TO_GROAN} -I{PATH_TO_GROAN} -lm -lgroan -o {YOUR_PROGRAM}`.