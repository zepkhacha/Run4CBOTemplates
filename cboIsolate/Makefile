dataDivFcn_noCBO: dataDivFcn_noCBO.o
	g++ -o dataDivFcn_noCBO dataDivFcn_noCBO.o $(shell root-config --libs) -L/cvmfs/gm2.opensciencegrid.org/prod9/g-2/gm2util/v9_52_00/slf7.x86_64.e15.prof/lib -lgm2util_blinders -lMinuit

dataDivFcn_noCBO.o: dataDivFcn_noCBO.cc 
	g++ -c -Wall -Wextra dataDivFcn_noCBO.cc $(shell root-config --libs --cflags) -I /cvmfs/gm2.opensciencegrid.org/prod9/g-2/gm2util/v9_52_00/include -ffast-math -O2
