slidingwindowfitplot: slidingwindowfitplot.o
	g++ -o slidingwindowfitplot slidingwindowfitplot.o $(shell root-config --libs) -L/cvmfs/gm2.opensciencegrid.org/prod9/g-2/gm2util/v9_52_00/slf7.x86_64.e15.prof/lib -lgm2util_blinders -lMinuit

slidingwindowfitplot.o: slidingwindowfitplot.cc 
	g++ -c -Wall -Wextra slidingwindowfitplot.cc $(shell root-config --libs --cflags) -I /cvmfs/gm2.opensciencegrid.org/prod9/g-2/gm2util/v9_52_00/include -ffast-math -O2

clean:
	rm slidingwindowfitplot.o slidingwindowfitplot
