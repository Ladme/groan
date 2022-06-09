gro_analysis: xdrfile/xdrfile.c xdrfile/xdrfile_xtc.c gro_io.c xtc_io.c analysis_tools.c
	gcc -c xdrfile/xdrfile.c -o xdrfile.o -O3 -march=native
	gcc -c xdrfile/xdrfile_xtc.c -o xdrfile_xtc.o -O3 -march=native
	gcc -c gro_io.c -o gro_io.o -O3 -march=native
	gcc -c xtc_io.c -o xtc_io.o -O3 -march=native
	gcc -c analysis_tools.c -o analysis_tools.o -O3 -march=native
	ar -rcs libgroan.a xdrfile.o xdrfile_xtc.o gro_io.o xtc_io.o analysis_tools.o
	rm -f *.o
	cp libgroan.a ~/.local/lib