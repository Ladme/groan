
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
            }
        }   
    } else {
        test_selection();
        test_analysis_tools();
    }
    
    return 0;
}