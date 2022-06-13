
/* This file shows how to work with atom selections in groan library. 
 * Compile this code as 'gcc example.c -L.. -I.. -lm -lgroan -o example'.
 * The output of this example code is in 'example_output.gro'.
 */

#ifdef CREATEEXAMPLE
    #define OUTPUT "example_output.gro"
#else
    #define OUTPUT "output.gro" 
#endif

#include <groan.h>

int main(void)
{   
    // open and read gro file
    system_t *system = load_gro("example.gro");
    if (system == NULL) return 1;

    // if we want to work with atom selections, we must always first select all atoms of the system
    atom_selection_t *all_atoms = select_system(system);

    // select all atoms of POPE and POPG lipids
    atom_selection_t *membrane = select_atoms(all_atoms, "POPE POPG", &match_residue_name);
    // ideally, we should check, if any atoms have been selected, but we will not bother with that here
    // if (membrane->n_atoms == 0) ..DEALLOCATE AND RETURN NON-ZERO EXIT CODE..

    // select peptide backbone
    atom_selection_t *backbone = select_atoms(all_atoms, "N CA C", &match_atom_name);

    // calculate peptide backbone center of geometry
    vec_t center_prot = {0};
    if (center_of_geometry(backbone, center_prot) != 0) return 1;

    // define properties of a cylinder for geometric selection of atoms
    // {radius, bottom of the cylinder, top of the cylinder}
    float cylinder_definition[3] = {2.5, -1, 1.5};
    // from membrane, select atoms located inside a cylinder positioned in the peptide center of geometry
    atom_selection_t *membrane_selection = select_geometry(membrane, center_prot, zcylinder, cylinder_definition);

    // we will not be needing "membrane" or "backbone" selections anymore
    free(membrane);
    membrane = NULL;
    free(backbone);
    backbone = NULL;

    // now we select oxygen atoms of water
    atom_selection_t *water_oxygens = select_atoms(all_atoms, "OW", &match_atom_name);

    // we will not be needing "all_atoms" anymore, so we free them
    free(all_atoms);
    all_atoms = NULL;

    // here we select only those oxygens that are located inside a sphere around peptide backbone center
    // sphere definition is just a radius
    float sphere_definition = 3.0;
    // note that we always provide pointer to the geometry definition in select_geometry()         vvvvvvvvvvvvvvvvvv
    atom_selection_t *water_selection_sphere = select_geometry(water_oxygens, center_prot, sphere, &sphere_definition);

    // for some reason, we also want to select different water oxygens: those positioned in an arbitrary box
    float box_definition[9] = {2, 4, 0, 5, 6, 8};
    // we use coordinates x = 0, y = 0, z = 0 as an absolute reference for select_geometry()
    vec_t absolute_center = {0, 0, 0};
    atom_selection_t *water_selection_box = select_geometry(water_oxygens, absolute_center, box, box_definition);
    
    free(water_oxygens);
    water_oxygens = NULL;

    // now we want to join these two selections of water oxygens
    // as we do not want multiple identical water oxygens in our final selection, we have to use selection_cat_unique()
    atom_selection_t *water_selection = selection_cat_unique(water_selection_sphere, water_selection_box);

    free(water_selection_sphere);
    water_selection_sphere = NULL;
    free(water_selection_box);
    water_selection_box = NULL;

    // now we want to join our membrane selection and our water oxygen selection
    // here, we can use simple selection_cat() as there can be no atom overlap 
    // (selection_cat() is faster than selection_cat_unique())
    atom_selection_t *final_selection = selection_cat(water_selection, membrane_selection);

    free(membrane_selection);
    free(water_selection);

    // now, if there is more than 0 atoms in the final selection, we write the selection into a gro file (output.gro)
    if (final_selection->n_atoms > 0) {
        FILE *output = fopen(OUTPUT, "w");
        // note that the numbering of the original gro file is kept in the output gro file
        write_gro(output, final_selection, system->box, velocities, "Arbitrary example selection");
        fclose(output);
    }
    

    // finally, we deallocate system_t structure and the final_selection
    free(system);
    free(final_selection);
    
    return 0;

}