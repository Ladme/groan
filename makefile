groan: src/xdrfile.o src/xdrfile_xtc.o src/dyn_array.o src/list.o src/dict.o src/gro_io.o src/xtc_io.o src/analysis_tools.o src/vector.o src/selection.o
	ar -rcs libgroan.a src/xdrfile.o src/xdrfile_xtc.o src/dyn_array.o src/list.o src/dict.o src/gro_io.o src/xtc_io.o src/vector.o src/selection.o src/analysis_tools.o
	make tests

src/xdrfile.o: src/xdrfile/xdrfile.c
	gcc -c src/xdrfile/xdrfile.c -o src/xdrfile.o -std=c99 -pedantic -Wall -O3 -march=native

src/xdrfile_xtc.o: src/xdrfile/xdrfile_xtc.c
	gcc -c src/xdrfile/xdrfile_xtc.c -o src/xdrfile_xtc.o -std=c99 -pedantic -Wall -O3 -march=native

src/dyn_array.o: src/general_structs/dyn_array.c
	gcc -c src/general_structs/dyn_array.c -o src/dyn_array.o -std=c99 -pedantic -Wall -Wextra -O3 -march=native

src/list.o: src/general_structs/list.c
	gcc -c src/general_structs/list.c -o src/list.o -std=c99 -pedantic -Wall -Wextra -O3 -march=native

src/dict.o: src/general_structs/dict.c
	gcc -c src/general_structs/dict.c -o src/dict.o -std=c99 -pedantic -Wall -Wextra -O3 -march=native

src/vector.o: src/general_structs/vector.c
	gcc -c src/general_structs/vector.c -o src/vector.o -std=c99 -pedantic -Wall -Wextra -O3 -march=native

src/gro_io.o: src/gro_io.c	
	gcc -c src/gro_io.c -o src/gro_io.o -std=c99 -pedantic -Wall -Wextra -O3 -march=native

src/xtc_io.o: src/xtc_io.c
	gcc -c src/xtc_io.c -o src/xtc_io.o -std=c99 -pedantic -Wall -Wextra -O3 -march=native

src/selection.o: src/selection.c
	gcc -c src/selection.c -o src/selection.o -std=c99 -pedantic -Wall -Wextra -O3 -march=native

src/analysis_tools.o: src/analysis_tools.c
	gcc -c src/analysis_tools.c -o src/analysis_tools.o -std=c99 -pedantic -Wall -Wextra -O3 -march=native

clean:
	rm -f *.a *.o src/*.a src/*.o

example: examples/example.c
	gcc examples/example.c -L. -I. -lgroan -lm -std=c99 -pedantic -Wall -Wextra -DCREATEEXAMPLE -o examples/example

tests: tests/tests.c tests/selection_tests.c tests/analysis_tools_tests.c libgroan.a groan.h
	gcc tests/tests.c tests/selection_tests.c tests/analysis_tools_tests.c -L. -I. -lgroan -lm -g -std=c99 -pedantic -Wall -Wextra -O3 -march=native -o tests/tests

install: libgroan.a groan.h
	cp -a libgroan.a /usr/local/lib/
	cp -a groan.h /usr/local/include/