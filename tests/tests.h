// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#ifndef TESTS_H
#define TESTS_H

#include <assert.h>
#include "../groan.h"

#define INPUT_GRO_FILE "../examples/example.gro"
#define NDX_FILE "../examples/index.ndx"
#define EMPTY_NDX_FILE "../examples/empty.ndx"

/*! @brief Tests whether float is close to some specified value. */
inline int closef(const float a, const float b, const float limit)
{
    if (a > b - limit && a < b + limit) return 1;
    return 0;
}

/*! @brief Collection of unit tests for selection.h. */ 
void test_selection(void);

/*! @brief Collection of unit tests for analysis_tools.h. */ 
void test_analysis_tools(void);

#endif /* TESTS_H */