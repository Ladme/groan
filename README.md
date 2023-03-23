# GROmacs ANalysis (groan)

Pure C library for analysis of Gromacs simulations.

## Installing

### Linux

1) Clone or download this repository.

2) Run `make` to compile the library.

### Windows / Mac OS

Sorry, no idea. Good luck.

## Tests

Validate the installation by running `./tests` in the `tests` directory. Do not remove the `examples` directory as the tests use several files located in there.

## Usage

Include `groan.h` in your code and link with `-lgroan -lm`.

## Groan-associated programs

- [center](https://github.com/Ladme/center): center simulation trajectory using Bai & Breen algorithm
- [contact](https://github.com/Ladme/contact): calculate contact matrix for selected atoms
- [cylinder](https://github.com/Ladme/cylinder): calculate density of atoms inside a cylinder
- [gndx](https://github.com/Ladme/gndx): create an ndx group from a selection using groan selection language
- [gselect](https://github.com/Ladme/gselect): select a group of atoms using groan selection language
- [leaflets2ndx](https://github.com/Ladme/leaflets2ndx): create ndx groups for lipids based on their type and membrane leaflet they occupy
- [memdian](https://github.com/Ladme/memdian): collection of several programs for analyzing membrane disruption
	- memthick: calculate average membrane thickness across the membrane
	- wdcalc: calculate average water defect in a cylinder
	- wdmap: calculate average water defect across the membrane
- [order & ordermap](https://github.com/Ladme/order): calculate lipid order parameters from coarse-grained trajectories
- [posdist](https://github.com/Ladme/posdist): calculate positions of atoms and their centers as well as distances between atoms or groups of atoms
- [scramblyzer](https://github.com/Ladme/scramblyzer): toolbox for analyzing lipid scrambling

Note that while groan library should work on pretty much any modern machine, all of the groan-associated programs require at least a basic level of POSIX-compliance and thus will most likely not work on Windows.

## Groan selection language

Groan library and its associated programs use "Groan selection language" to specify selections of atoms. Groan selection language is similar to the selection language used by VMD. You can play around with the groan selection language using [gselect](https://github.com/Ladme/gselect).

### Basic queries
Select atoms based on their:
1) residue name using `resname XYZ`. For example, `resname POPE` will select all atoms corresponding to residues named POPE.
2) residue number using `resid XYZ`. For example, `resid 17` will select all atoms corresponding to residue with number 17.
3) atom name using `name XYZ`. For example, `name P` will select all atoms which name corresponds to P.
4) atom number using `serial XYZ`. For example, `serial 256` will select an atom with atom number 256.
5) ndx group using `NDX_GROUP_NAME`. For example, `Protein` will select all atoms belonging to the ndx group Protein (ndx file must be provided). 

You can specify multiple identifiers in your query. By using `resname POPE POPC`, you will select all atoms corresponding to residues named POPE, as well as all atoms corresponding to residues named POPC.
See examples of similar queries below:
1) `resid 13 15 16 17` will select all atoms corresponding to residues with number 13, 15, 16, or 17.
2) `name P CA HA` will select all atoms with atom name P, CA, or HA.
3) `serial 245 267 269 271` will select atoms with atom number 245, 267, 269, or 271.

You can't use multiple identifiers for ndx groups. In that case, you have to use the `or` operator (see below).

You can also select all atoms by using `all`.

The length of the query is only limited by the size of your computer's memory. Note that longer queries will take longer to parse.

### Ranges
Instead of writing residue or atom numbers explicitly, you can use keyword `to` or `-` to specify a range. For example, instead of writing `resid 14 15 16 17 18 19 20`, you can use `resid 14 to 20` or `resid 14 - 20`. This will select all atoms corresponding to residues with residue numbers 14, 15, 16, 17, 18, 19, and 20. Note that both `to` and `-` must always be separated from the rest of the query by (at least one) whitespace.

You can also specify multiple ranges in a single query or combine ranges with explicitly provided numbers. For example, `serial 1 3 to 6 10 12 to 14 17` will expand to `serial 1 3 4 5 6 10 12 13 14 17`.

### Negations
Using keyword `not` or `!` in front of a query will negate the query. For example, query `not name CA` or `! name CA` will select all atoms which name does NOT correspond to CA. Similarly, `not resname POPE POPG` will select all atoms that correspond to residues with names other than POPE or POPG. Note again that both `not` and `!` must always be separated from the rest of the query by a whitespace.

### 'And' and 'or' operations
You can combine basic queries by using `and` (`&&`) and `or` (`||`) operators. 

Joining two queries by `and` will select only atoms that were selected by BOTH of the queries. For example, `resname POPE and name P` will select all atoms that belong to residues named POPE and that have the name P. Similarly, `resid 17 18 && serial 256 to 271` will select only atoms corresponding to residue 17 or 18 and with atom numbers between 256 and 271 (including 271).

Joining two queries by `or` will select atoms that were selected by AT LEAST ONE of the queries. For example, `resname POPE or name P` will select all atoms that belong to residues named POPE as well as all atoms with the name P. Similarly, `resid 17 18 || serial 256 to 271` will select all atoms corresponding to residue 17 or 18 as well as all atoms with atom numbers between 256 and 271.

In case multiple `and` and/or `or` operators are used in a single query, they are evaluated from left to right. For example, `resname POPE or name CA and not Protein` will select all atoms belonging to residues named POPE or having the atom name CA but all these atoms can not belong to the ndx group called Protein.

Note again that `and`, `&&`, `or`, and `||` must all be separated from the rest of the query by a whitespace. Also note that there can be no more than 50 individual sub-queries connected by operators in a single query.

### Parentheses
You can change the order in which the individual sub-queries and operations are evaluated by using parentheses `(` and `)`. Expressions enclosed in parentheses are evaluated first (think math). For example, `resname POPE or (name CA and not resid 18 to 21)` will select all atoms belonging to residues named POPE along with all atoms that a) have the atom name P and b) do not correspond to residues numbered 18 to 21. Meanwhile `(resname POPE or name CA) and not resid 18 to 21` is equivalent to `resname POPE or name CA and not resid 18 to 21` (i.e. this will select all atoms belonging to residues named POPE or having the atom name CA but all of these atoms can't belong to residues 18 to 21).

You can place parenthetical expressions into other parenthetical expressions. For example `serial 1 to 6 or ( name CA and resname POPE || (resid 1 to 7 or serial 123 to 128) ) and Protein` is a valid query, albeit possibly way too convoluted.

You can also place `not` (`!`) operator in front of a parenthetical expression. For example, `! (serial 1 to 6 && name P)` will select all atoms that do NOT have atom number between 1 and 6 while also having the atom name P.

Note that parentheses are allowed to be but do not have to be separated from the rest of the query by a whitespace.

## Acknowledgments
The groan library uses `xdrfile` library developed by David van der Speol and Erik Lindahl and published under the BSD License. Thank you!

To calculate center of geometry/mass, the groan library and all its associated programs use Bai & Breen algorithm for calculating center of mass in periodic systems (https://doi.org/10.1080/2151237X.2008.10129266).

## Disclaimer
Please note that while groan library is an attempt at creating a fast and reasonably robust library for analyzing Gromacs simulations, it is mostly a tool for improving my rather pathetic C programming skills. In other words, even though I'm doing my best, there WILL be bugs and memory leaks.

Also, the amount of time that went into designing the library to be user-friendly or even understandable to anyone other than me is... zero. There is some documentation and many of the functions do have names that actually make at least some sense but if you want to actually start using groan in your code... well, good luck with that.

At the same time, I'm not saying that the library is completely worthless as I quite successfully use many of the groan-associated programs on a daily basis. And I do believe that THOSE programs could actually be useful to other people.
