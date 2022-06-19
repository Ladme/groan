// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#ifndef DYN_ARRAY_H
#define DYN_ARRAY_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>


/*! \brief Two-dimensional array that dynamically expands its size when needed.
 * 
 * \par Basic usage
 * 
 * 1) Create an array with `dyn_array2D *array = array_create(n_rows, n_columns)`.
 * 2) Assign value to the array with `array_add(array, row, column, value, '=')`.
 * 3) Select value from the array with `array[row][column]`.
 * 4) Free memory that was allocated to the array with `array_destroy(array)`.
 * 
 * \par Example
 * 
 * Using array_create(2, 5) you create an array with 2 rows, each row containing 5 columns.
 * Using array_add(arr, 5, 7, 4.5, '=') you then assign a float 4.5 to the 6th row and 8th column of the array.
 * The array automatically expands to be able to accomodate this value.
 * The new array has 6 rows, each with 8 columns.
 */
typedef struct dynarray {
    size_t n_rows;    // number of rows in an array
    size_t n_cols;    // number of columns in an array
    float **arr;      // the array itself
} dyn_array2D;


/*! \brief Frees memory allocated for a dynamic 2D array.
 * 
 * \param arr       pointer to a dynamic 2D array
 * 
 */
void array_destroy(dyn_array2D *arr);


/*! \brief Creates a dynamic 2D array of initial size n_rows x n_cols.
 *
 * \par Choosing initial size
 * 
 * Performance-wise it is a good idea to set the initial size of you array 
 * (close) to its final size, so fewer memory reallocations have to be performed.
 * If you don't care about performace, you can use any non-negative numbers as the
 * initial dimensions of the array (i.e. including zero).
 * 
 * \par Initial values
 * 
 * Initial values of all cells in the array are always zero.
 * This also holds true for new cells when the array expands.
 * 
 * \param n_rows        initial number of rows in the array
 * \param n_cols        initial number of columns in the array
 * 
 * \return pointer to the dynamic 2D array
 */
dyn_array2D *array_create(const size_t n_rows, const size_t n_cols);


/*! \brief Adds value to a dynamic 2D array and expands the array if needed.
 *
 * \par Operations
 * 
 * This function does not necessarily assign the value to the array.
 * Instead, it performs a specified operation:
 * '=': value is assigned to the selected value in the array
 * '+': value is added to the selected value in the array
 * '-': value is subtracted from the selected value in the array
 * '*': value is used to multiply the selected value in the array
 * '/': value is used to divide the selected value in the array
 * 
 * \par Initial values
 * 
 * Initial values of all newly added cells in the array are zero.
 * 
 * \param arr       pointer to a dynamic 2D array
 * \param row       index of a row to which the value shall be added
 * \param col       index of a column to which the value shall be added
 * \param value     value to add to the array
 * \param operation operation to perform (=, +, -, *, /) 
 * 
 * \return Zero, if value has been successfully added into the array. Else non-zero.
 */
int array_add(
        dyn_array2D *arr, 
        const size_t row, 
        const size_t col, 
        const float value, 
        const char operation);

#endif /* DYN_ARRAY_H */