// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#ifndef GRO_IO_H
#define GRO_IO_H

#include <ctype.h>
#include <string.h>
#include "gro.h"

/*
 * Used for gro file format printing. 
 * If set to "velocities", velocities of all particles will be printed even if they are zero.
 * If set to "no_velocities", velocities of all particules will NOT be printed.
 */
typedef enum write_mode {
    no_velocities = 0, 
    velocities = 1
} write_mode_t;

/* ! \brief Checks whether the string can be converted into a valid decimal integer.
 * 
 * \param string    string to check
 * 
 * \return Returns 1 if string only contains decimal digits, 
 * whitespaces and at most one +/- sign. Else returns 0.
 */
int isdecimal(const char *string);

/* ! \brief Checks whether the string can be converted into a valid decimal float.
 * 
 * \param string    string to check
 * 
 * \return One if string only contains decimal digits, 
 * whitespaces, at most one '.' and at most one +/- sign. Else returns 0.
 */
int isdecimalf(const char *string);

/* ! \brief Slices a string into a substring of target length starting from target position.
 *
 * Does not check array boundaries. The substring is terminated with '\0'.
 * 
 * \param src       string to be sliced
 * \param dest      string to which the output string will be saved
 * \param start     index from which the slicing will start
 * \param len       length of the final substring
 */ 
void get_fragment(const char *src, char *dest, const size_t start, const size_t len);

/* ! \brief Reads a target part of line and parses it into a gro integer (groint_t).
 * 
 * Uses get_fragment() function to slice the line.
 * 
 * \param line      line to parse
 * \param start     index from which the line will be parsed
 * \param len       length of the parsed segment
 * \param element   pointer to memory into which the parsed gro integer will be saved
 * 
 * \return Zero if parsing was successful. Else returns non-zero.
 */
int parse_int(const char *line, const size_t start, const size_t len, groint_t *element);

/* ! \brief Reads a target part of line, trims it and copies it to the provided pointer.
 * 
 * Uses get_fragment() function to slice the line.
 * 
 * \param line      line to parse
 * \param start     index from which the line will be parsed
 * \param len       length of the parsed segment
 * \param element   pointer to memory into which the parsed string will be saved
 * 
 * \return Zero if parsing was succesful. Else returns non-zero.
 */
int parse_string(const char *line, const size_t start, const size_t len, char *element);

/* ! \brief Reads a target part of line and parses it into a float.
 * 
 * Uses get_fragment() function to slice the line.
 * 
 * \param line      line to parse
 * \param start     index from which the line will be parsed
 * \param len       length of the parsed segment
 * \param element   pointer to memory into which the parsed float will be saved
 * 
 * \return Zero if parsing was successful. Else returns non-zero.
 */
int parse_float(const char *line, const size_t start, const size_t len, float *element);

/* ! \brief Parses a single line of a gro file, reading information about a single atom.
 * 
 * \param line      line to parse
 * \param atom      pointer to an atom which information is read
 * 
 * \return Zero if parsing was successful. Else returns non-zero.
 */
int parse_gro_line(const char *line, atom_t *atom);

/* ! \brief Reads a gro file and loads information about the atoms.
 * 
 * Note on atom numbering:
 * Gro files only support atom numbers of <100,000.
 * If there are more atoms, the atom number wraps, starting from 1 again.
 * Atom_number in the atom_t structure represents the real atom_number from the gro file.
 * In other words, we do not remove the wrapping because there is no way to achieve consistency
 * in case of broken gro files (with some atom numbers missing, for example).
 * 
 * Gromacs, internally, does not even care about atom numbers in the gro file and assigns
 * the atom numbers based on the real number of atoms.
 * To also capture this information in groan, we therefore introduce gmx_atom_number
 * to the atom_t structure which can be >99,999.
 *
 * \param filename  path to the gro file
 * 
 * \return Pointer to a system_t structure, if successful.
 * Else returns NULL and prints an error message to stderr.
 */
system_t *load_gro(const char *filename);

/* ! \brief Prints information about the selected atoms in gro format into stream.
 * 
 * \param stream            output stream for printing (e.g. stdout)
 * \param atoms             selection of atoms to be printed
 * \param boxsize           size of the simulation cell
 * \param write_mode        should the velocities be printed (velocities) or not (no_velocities)
 * \param comment           string that will be printed as the first line of gro file
 * 
 * \return Zero if successful, else non-zero.
 */
int write_gro(
        FILE *stream, 
        const atom_selection_t *atoms,
        const box_t boxsize,
        const write_mode_t write_mode,
        const char *comment);
        
#endif /* GRO_IO_H */