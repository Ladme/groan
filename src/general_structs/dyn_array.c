// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#include "dyn_array.h"

void array_destroy(dyn_array2D *arr)
{
    // loop through all the rows and deallocate them
    for (size_t i = 0; i < arr->n_rows; ++i) {
        free(arr->arr[i]);
    }
    // deallocate the array itself
    free(arr->arr);
    // deallocate the wrapping structure
    free(arr);
}

dyn_array2D *array_create(const size_t n_rows, const size_t n_cols)
{
    // allocate memory for the dynamic array structure
    dyn_array2D *arr = calloc(1, sizeof(dyn_array2D));
    if (arr == NULL) {
        return NULL;
    }

    arr->n_rows = n_rows;
    arr->n_cols = n_cols;

    // allocate memory for the array itself
    arr->arr = calloc(n_rows, sizeof(float *));
    if (arr->arr == NULL) {
        free(arr);
        return NULL;
    }
    
    // allocate memory for the individual rows
    for (size_t i = 0; i < n_rows; ++i) {
        arr->arr[i] = calloc(n_cols, sizeof(float));

        if (arr->arr[i] == NULL) {
            array_destroy(arr);
            return NULL;
        }
    }

    return arr;
}


int array_add(
        dyn_array2D *arr, 
        const size_t row, 
        const size_t col, 
        const float value, 
        const char operation)
{

    // save the initial sizes of the array
    size_t old_n_rows = arr->n_rows;
    size_t old_n_cols = arr->n_cols;

    // if the requested row index is larger than the array size
    if (row >= arr->n_rows) {
        // reallocate memory for the entire array increasing the number of rows
        float **new_arr = realloc(arr->arr, sizeof(float *) * (row + 1));
        if (new_arr == NULL) {
            array_destroy(arr);
            return 1;
        }
        arr->arr = new_arr;

        arr->n_rows = row + 1;

        // set the new row pointers to NULL
        // this is necessary if we want to correctly destroy the array in case of an error
        memset(arr->arr + old_n_rows, 0, (arr->n_rows - old_n_rows) * sizeof(float *));

        // increase the number of columns in array, if needed
        if (col >= arr->n_cols) {
            arr->n_cols = col + 1;
        }

        // allocate memory for the NEW rows
        // if the number of columns increased, memory for the new columns
        // will be included in memory allocated to these new rows
        for (size_t i = old_n_rows; i < arr->n_rows; ++i) {
            arr->arr[i] = calloc(arr->n_cols, sizeof(float));

            if (arr->arr[i] == NULL) {
                array_destroy(arr);
                return 1;
            }
        }
    }

    // if the requested columns index is larger than array size
    // (we have to compare with the original size because the current
    // size could have been changed in the above block)
    if (col >= old_n_cols) {
        arr->n_cols = col + 1;
        // reallocate memory for the individual rows, expanding them
        // to accomodate larger number of columns
        // note that we only reallocate rows which were in the original array
        // if there were any new rows added, memory for these has already been allocated in
        // the previous block (new columns have already been included in this allocation)
        for (size_t i = 0; i < old_n_rows; ++i) {
            float *new_row = realloc(arr->arr[i], arr->n_cols * sizeof(float));
            
            if (new_row == NULL) {
                array_destroy(arr);
                return 1;
            }

            arr->arr[i] = new_row;

            // set values to 0
            memset(arr->arr[i] + old_n_cols, 0, (arr->n_cols - old_n_cols) * sizeof(float));
        }
    }

    // perform the requested operation
    switch (operation) {
        case '=': arr->arr[row][col]  = value; break;
        case '+': arr->arr[row][col] += value; break;
        case '-': arr->arr[row][col] -= value; break;
        case '*': arr->arr[row][col] *= value; break;
        case '/': arr->arr[row][col] /= value; break;
    }
    return 0;
}

/* debugging function; prints array to stdout */
/*static void array_print(const dyn_array2D *arr)
{
    printf("N_ROWS: %lu\n", arr->n_rows);
    printf("N_COLS: %lu\n", arr->n_cols);
    for (size_t i = 0; i < arr->n_rows; ++i) {
        for (size_t j = 0; j < arr->n_cols; ++j) {
            printf("% 7.3f ", arr->arr[i][j]);
        }
        printf("\n");
    }
}*/

/*int main(void)
{
    dyn_array2D *array = array_create(0, 0);
    if (array == NULL) {
        fprintf(stderr, "Could not allocate memory.\n");
        return 1;
    }

    array_add(array, 2, 3, 6.5, '=');

    if (array_add(array, 7, 21, 2, '+') != 0) {
        fprintf(stderr, "Error adding item into array.\n");
        return 1;
    }

    array_add(array, 5, 4, 7, '-');
    array_add(array, 7, 21, 8, '=');
    array_add(array, 5, 4, 3, '/');
    array_add(array, 7, 21, 4, '*');
    array_add(array, 7, 27, 3.2, '=');
    array_add(array, 11, 0, 0.4, '+');

    array_print(array);

    array_destroy(array);
    return 0;
}*/