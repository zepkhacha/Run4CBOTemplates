fitplot: fitplot.o
	g++ -o fitplot fitplot.o $(shell root-config --libs) \
        -L/cvmfs/gm2.opensciencegrid.org/prod/g-2/gm2util/v10_16_00/slf7.x86_64.e15.prof/lib -lgm2util_blinders -lMinuit -L${TBB_LIB} -ltbb

fitplot.o: fitplot.cc Makefile
	g++ -c -Wall -Wextra fitplot.cc $(shell root-config --libs --cflags) -I /cvmfs/gm2.opensciencegrid.org/prod/g-2/gm2util/v10_16_00/include -ffast-math -O2 \
	-I ${TBB_INC}
