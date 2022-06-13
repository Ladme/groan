// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#include "vector.h"

void vec_sum(vec_t result, const vec_t add)
{
    result[0] += add[0];
    result[1] += add[1];
    result[2] += add[2];
}

void vec_sub(vec_t result, const vec_t subtract)
{
    result[0] -= subtract[0];
    result[1] -= subtract[1];
    result[2] -= subtract[2];
}

void vec_mult(vec_t result, const float scalar)
{
    result[0] *= scalar;
    result[1] *= scalar;
    result[2] *= scalar;
}

void vec_div(vec_t result, const float scalar)
{
    result[0] /= scalar;
    result[1] /= scalar;
    result[2] /= scalar;
}