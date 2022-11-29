// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#include "tests.h"

int main(int argc, char **argv)
{
    if (argc > 1) {
        for (int i = 1; i < argc; ++i) {
            if (!strcmp(argv[i], "selection")) {
                test_selection();
            } else if (!strcmp(argv[i], "tools")) {
                test_analysis_tools();
            } else if (!strcmp(argv[i], "xtc") || !strcmp(argv[i], "trr") || !strcmp(argv[i], "xdr")) {
                test_xdr();
            } else if (!strcmp(argv[i], "gro")) {
                test_gro_io();
            }
        }   
    } else {
        test_gro_io();
        test_selection();
        test_analysis_tools();
        test_xdr();
    }
    
    return 0;
}