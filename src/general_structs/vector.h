// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#ifndef VECTOR_H
#define VECTOR_H

/*! \brief An array of three floats. */
typedef float vec_t[3];

/*! \brief Sums two vectors by adding 'add' to 'result'.
 *
 * \param result    vector which will be modified
 * \param add       vector which is being added
 * 
 */
void vec_sum(vec_t result, const vec_t add);

/*! \brief Subtracts two vectors by subtracting 'subtract' from 'result'.
 *
 * \param result    vector which will be modified
 * \param subtract  vector which is being subtracted
 * 
 */
void vec_sub(vec_t result, const vec_t subtract);

/*! \brief Multiplies a vector by scalar.
 *
 * \param vec       vector to be multiplied
 * \param scalar    scalar to be used for multiplication
 *
 */
void vec_mul(vec_t result, const float scalar);

/*! \brief Divides a vector by scalar.
 *
 * \param vec       vector to be divided
 * \param scalar    scalar to be used for division
 *
 */
void vec_div(vec_t result, const float scalar);

#endif /* VECTOR_H */