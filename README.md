# GROmacs ANalysis (groan)

C library for analysis of Gromacs simulations.

## How to install

Just run `make`.

## Tests

Validate the installation by running `./tests` in the `tests` directory. Do not remove the `examples` directory as the tests use several files located in there.

## How to use in your C programs

Include `groan.h` in your code and link with `-lm -lgroan`.

## Examples

Example code is included in the `examples` directory.