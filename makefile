groan: xdrfile.o xdrfile_xtc.o gro_io.o xtc_io.o analysis_tools.o
	ar -rcs libgroan.a src/xdrfile.o src/xdrfile_xtc.o src/gro_io.o src/xtc_io.o src/analysis_tools.o

xdrfile.o: src/xdrfile/xdrfile.c
	gcc -c src/xdrfile/xdrfile.c -o src/xdrfile.o -std=c99 -pedantic -Wall -O3 -march=native

xdrfile_xtc.o: src/xdrfile/xdrfile_xtc.c
	gcc -c src/xdrfile/xdrfile_xtc.c -o src/xdrfile_xtc.o -std=c99 -pedantic -Wall -O3 -march=native

gro_io.o: src/gro_io.c	
	gcc -c src/gro_io.c -o src/gro_io.o -std=c99 -pedantic -Wall -Wextra -O3 -march=native

xtc_io.o: src/xtc_io.c
	gcc -c src/xtc_io.c -o src/xtc_io.o -std=c99 -pedantic -Wall -Wextra -O3 -march=native

analysis_tools.o: src/analysis_tools.c
	gcc -c src/analysis_tools.c -o src/analysis_tools.o -std=c99 -pedantic -Wall -Wextra -O3 -march=native
	
clean:
	rm -f *.a *.o src/*.a src/*.o

dev: src/analysis_tools.c src/xtc_io.c src/gro_io.c src/xdrfile/xdrfile_xtc.c src/xdrfile/xdrfile.c
	rm -f *.a *.o src/*.a src/*.o examples/example examples/example_output.gro
	gcc -c src/analysis_tools.c -o src/analysis_tools.o -std=c99 -pedantic -Wall -O0 -g
	gcc -c src/xtc_io.c -o src/xtc_io.o -std=c99 -pedantic -Wall -O0 -g
	gcc -c src/gro_io.c -o src/gro_io.o -std=c99 -pedantic -Wall -Wextra -O0 -g
	gcc -c src/xdrfile/xdrfile_xtc.c -o src/xdrfile_xtc.o -std=c99 -pedantic -Wall -Wextra -O0 -g
	gcc -c src/xdrfile/xdrfile.c -o src/xdrfile.o -std=c99 -pedantic -Wall -Wextra -O0 -g
	ar -rcs libgroan.a src/xdrfile.o src/xdrfile_xtc.o src/gro_io.o src/xtc_io.o src/analysis_tools.o

example: examples/example.c
	gcc examples/example.c -L. -I. -lm -lgroan -std=c99 -pedantic -Wall -Wextra -DCREATEEXAMPLE -o examples/example
