
/* This file shows how to work with atoms and atom selections in groan library. 
 * Compile this code as 'gcc atoms_selections.c -L.. -I.. -lgroan -lm -o atoms_selections'.
 */

#include <groan.h>

int main(void)
{   
    // open and read a gro file
    system_t *system = load_gro("example.gro");
    if (system == NULL) return 1;

    // we can loop through the atoms in the system and print some information about them
    // here, we will loop just through the first five atoms in the system
    for (size_t i = 0; i < 5; ++i) {
        // each atom contains the following information:
        // a) residue_number: index of the residue containing this atom
        printf("> residue_number: %d\n", system->atoms[i].residue_number);
        // b) residue_name: name of the residue containing this atom (max 5 characters)
        printf("> residue_name: %s\n", system->atoms[i].residue_name);
        // c) atom name: name of the atom
        printf("> atom_name: %s\n", system->atoms[i].atom_name);
        // d) atom number: atom number as written in the gro file
        printf("> atom_number: %d\n", system->atoms[i].atom_number);
        // e) gmx_atom_number: atom number as used by gromacs
        printf("> gmx_atom_number: %ld\n", system->atoms[i].gmx_atom_number);
        // f) position: xyz coordinates of the atom stored in vec_t (array of 3 floats)
        printf("> position: %f %f %f\n", system->atoms[i].position[0], system->atoms[i].position[1], system->atoms[i].position[2]);
        // g) velocity: velocity of the atom
        printf("> velocity: %f %f %f\n", system->atoms[i].velocity[0], system->atoms[i].velocity[1], system->atoms[i].velocity[2]);
        // h) force: force acting on the atom (gro files never contain this information, you have to use trr file for this)
        printf("> force: %f %f %f\n\n", system->atoms[i].force[0], system->atoms[i].force[1], system->atoms[i].force[2]);
    }

    // if we want to work with atom selections, we must always first select all atoms of the system
    atom_selection_t *all_atoms = select_system(system);

    // the all_atoms structure then contains:
    // a) the number of atoms in the selection
    printf("Number of atoms: %ld\n", all_atoms->n_atoms);
    // b) an array of pointers to atoms in the selection that you can loop through:
    //for (size_t i = 0; i < all_atoms->n_atoms; ++i) {
    //    printf("ATOM %d, name %s\n", all_atoms->atoms[i]->atom_number, all_atoms->atoms[i]->atom_name);
    //}
    
    // note that the atoms in any atom selection are just pointers to the corresponding atoms stored in the system structure
    // if you change any property of any atom in any selection, this change will propagate through the system_t structure
    // into every other atom selection containing this atom

    // select atoms of POPE and POPG lipids using smart_select function
    // this function uses groan selection language to specify the selection of atoms,
    // see https://github.com/Ladme/groan#groan-selection-language for more information
    atom_selection_t *membrane = smart_select(all_atoms, "resname POPE POPG", NULL);

    // we should always check whether atoms have been successfully selected
    if (membrane == NULL || membrane->n_atoms == 0) {
        free(membrane);
        free(all_atoms);
        free(system);
        return 1;
    }

    // to select peptide backbone, we will first read an ndx file
    dict_t *ndx_groups = read_ndx("index.ndx", system);
    if (ndx_groups == NULL) {
        free(membrane);
        free(all_atoms);
        free(system);
        return 1;
    }

    // the ndx_groups dictionary contains pairs of ndx group names and atom selections
    // we can use this dictionary in the smart_select function
    select_t *backbone = smart_select(all_atoms, "Backbone", ndx_groups);
    // this will select all atoms corresponding to the ndx group 'Backbone'
    // note that we are using `select_t` instead of `atom_selection_t`
    // `select_t` is just typedef for `atom_selection_t`
    if (backbone == NULL || backbone->n_atoms == 0) {
        free(backbone);
        dict_destroy(ndx_groups);
        free(membrane);
        free(all_atoms);
        free(system);
        return 1;
    }

    // you can then: 
    // a) concatenate the atom selections
    select_t *membrane_backbone = selection_cat(membrane, backbone);
    free(membrane_backbone);
    // b) find intersections of the atom selections
    select_t *intersect = selection_intersect(membrane, backbone); // (this will be empty)
    free(intersect);
    // c) remove atoms that are part of selection from some other selection
    size_t removed = selection_remove(membrane, backbone);
    // removed corresponds to the number of atoms removed from `membrane` 
    // here, removed == 0
    // d) use all sorts of other functions, as described in selection.h 

    // finally, we must deallocate all atom selections, the system and the ndx_groups dictionary
    // you can deallocate the atom selections and the system using a simple free
    // dictionaries must be deallocated by calling the dict_destroy function
    free(backbone);
    free(membrane);
    free(all_atoms);
    free(system);
    dict_destroy(ndx_groups);
    
    return 0;

}