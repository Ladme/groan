// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#ifndef TESTS_H
#define TESTS_H

#include <assert.h>
#include "../groan.h"

#define INPUT_GRO_FILE "../examples/example.gro"
#define SMALL_GRO_FILE "../examples/micro.gro"
#define NDX_FILE "../examples/index.ndx"
#define EMPTY_NDX_FILE "../examples/empty.ndx"
#define INPUT_XTC_FILE "../examples/example.xtc"
#define INPUT_TRR_FILE "../examples/example.trr"

/*! @brief Collection of unit tests for selection.h. */ 
void test_selection(void);

/*! @brief Collection of unit tests for analysis_tools.h. */ 
void test_analysis_tools(void);

/*! @brief Collection of unit tests for xtc_io.h and trr_io.h. */ 
void test_xdr(void);


#endif /* TESTS_H */