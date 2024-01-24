createTMethod: createTMethod.o
	g++ -o createTMethod createTMethod.o $(shell root-config --libs) -L/cvmfs/gm2.opensciencegrid.org/prod/g-2/gm2util/v9_52_00/slf7.x86_64.e15.prof/lib -lMinuit -L${TBB_LIB} -ltbb

createTMethod.o: createTMethod.cc 
	g++ -c -Wall -Wextra createTMethod.cc $(shell root-config --libs --cflags) -I /cvmfs/gm2.opensciencegrid.org/prod/g-2/gm2util/v9_52_00/include -ffast-math -O2 -I ${TBB_INC}

percalofitplot: percalofitplot.o
	g++ -o percalofitplot percalofitplot.o $(shell root-config --libs) \
        -L/cvmfs/gm2.opensciencegrid.org/prod9/g-2/gm2util/v9_52_00/slf7.x86_64.e15.prof/lib -lgm2util_blinders -lMinuit -L${TBB_LIB} -ltbb

percalofitplot.o: percalofitplot.cc Makefile
	g++ -c -Wall -Wextra percalofitplot.cc $(shell root-config --libs --cflags) -I /cvmfs/gm2.opensciencegrid.org/prod9/g-2/gm2util/v9_52_00/include -ffast-math -O2 \
	-I ${TBB_INC}

slidingCBO: slidingCBO.o
	g++ -o slidingCBO slidingCBO.o $(shell root-config --libs) -L/cvmfs/gm2.opensciencegrid.org/prod9/g-2/gm2util/v9_52_00/slf7.x86_64.e15.prof/lib -lgm2util_blinders -lMinuit

slidingCBO.o: slidingCBO.cc 
	g++ -c -Wall -Wextra slidingCBO.cc $(shell root-config --libs --cflags) -I /cvmfs/gm2.opensciencegrid.org/prod9/g-2/gm2util/v9_52_00/include -ffast-math -O2

slidingwindowfitplot: slidingwindowfitplot.o
	g++ -o slidingwindowfitplot slidingwindowfitplot.o $(shell root-config --libs) -L/cvmfs/gm2.opensciencegrid.org/prod9/g-2/gm2util/v9_52_00/slf7.x86_64.e15.prof/lib -lgm2util_blinders -lMinuit

slidingwindowfitplot.o: slidingwindowfitplot.cc 
	g++ -c -Wall -Wextra slidingwindowfitplot.cc $(shell root-config --libs --cflags) -I /cvmfs/gm2.opensciencegrid.org/prod9/g-2/gm2util/v9_52_00/include -ffast-math -O2

cboResidual: cboResidual.o
	g++ -o cboResidual cboResidual.o $(shell root-config --libs) -L/cvmfs/gm2.opensciencegrid.org/prod9/g-2/gm2util/v9_52_00/slf7.x86_64.e15.prof/lib -lgm2util_blinders -lMinuit

cboResidual.o: cboResidual.cc 
	g++ -c -Wall -Wextra cboResidual.cc $(shell root-config --libs --cflags) -I /cvmfs/gm2.opensciencegrid.org/prod9/g-2/gm2util/v9_52_00/include -ffast-math -O2

dataDivFcn_noCBO: dataDivFcn_noCBO.o
	g++ -o dataDivFcn_noCBO dataDivFcn_noCBO.o $(shell root-config --libs) -L/cvmfs/gm2.opensciencegrid.org/prod9/g-2/gm2util/v9_52_00/slf7.x86_64.e15.prof/lib -lgm2util_blinders -lMinuit

dataDivFcn_noCBO.o: dataDivFcn_noCBO.cc 
	g++ -c -Wall -Wextra dataDivFcn_noCBO.cc $(shell root-config --libs --cflags) -I /cvmfs/gm2.opensciencegrid.org/prod9/g-2/gm2util/v9_52_00/include -ffast-math -O2

buildTemplate: buildTemplate.o
	g++ -o buildTemplate buildTemplate.o $(shell root-config --libs) -L/cvmfs/gm2.opensciencegrid.org/prod9/g-2/gm2util/v9_52_00/slf7.x86_64.e15.prof/lib -lgm2util_blinders -lMinuit -L${TBB_LIB} -ltbb

buildTemplate.o: buildTemplate.cc Makefile
	g++ -c -Wall -Wextra buildTemplate.cc $(shell root-config --libs --cflags) -I /cvmfs/gm2.opensciencegrid.org/prod9/g-2/gm2util/v9_52_00/include -ffast-math -O2 -I ${TBB_INC}

cboAmplitude: cboAmplitude.o
	g++ -o cboAmplitude cboAmplitude.o $(shell root-config --libs) -L/cvmfs/gm2.opensciencegrid.org/prod9/g-2/gm2util/v9_52_00/slf7.x86_64.e15.prof/lib -lgm2util_blinders -lMinuit

cboAmplitude.o: cboAmplitude.cc 
	g++ -c -Wall -Wextra cboAmplitude.cc $(shell root-config --libs --cflags) -I /cvmfs/gm2.opensciencegrid.org/prod9/g-2/gm2util/v9_52_00/include -ffast-math -O2

plots: plots.o
	g++ -o plots plots.o $(shell root-config --libs) -L/cvmfs/gm2.opensciencegrid.org/prod9/g-2/gm2util/v9_52_00/slf7.x86_64.e15.prof/lib -lgm2util_blinders -lMinuit

plots.o: plots.cc 
	g++ -c -Wall -Wextra plots.cc $(shell root-config --libs --cflags) -I /cvmfs/gm2.opensciencegrid.org/prod9/g-2/gm2util/v9_52_00/include -ffast-math -O2

windowFFT: windowFFT.o
	g++ -o windowFFT windowFFT.o $(shell root-config --libs) -L/cvmfs/gm2.opensciencegrid.org/prod9/g-2/gm2util/v9_52_00/slf7.x86_64.e15.prof/lib -lgm2util_blinders -lMinuit

windowFFT.o: windowFFT.cc 
	g++ -c -Wall -Wextra windowFFT.cc $(shell root-config --libs --cflags) -I /cvmfs/gm2.opensciencegrid.org/prod9/g-2/gm2util/v9_52_00/include -ffast-math -O2

compareCalo: compareCalo.o
	g++ -o compareCalo compareCalo.o $(shell root-config --libs) -L/cvmfs/gm2.opensciencegrid.org/prod9/g-2/gm2util/v9_52_00/slf7.x86_64.e15.prof/lib -lgm2util_blinders -lMinuit

compareCalo.o: compareCalo.cc
	g++ -c -Wall -Wextra compareCalo.cc $(shell root-config --libs --cflags) -I /cvmfs/gm2.opensciencegrid.org/prod9/g-2/gm2util/v9_52_00/include -ffast-math -O2

makeTemplate: makeTemplate.o
	g++ -o makeTemplate makeTemplate.o $(shell root-config --libs) -L/cvmfs/gm2.opensciencegrid.org/prod9/g-2/gm2util/v9_52_00/slf7.x86_64.e15.prof/lib -lgm2util_blinders -lMinuit

makeTemplate.o: makeTemplate.cc
	g++ -c -Wall -Wextra makeTemplate.cc $(shell root-config --libs --cflags) -I /cvmfs/gm2.opensciencegrid.org/prod9/g-2/gm2util/v9_52_00/include -ffast-math -O2

clean:
	rm cboAmplitudeBinned.o cboAmplitudeBinned cboResidual.o cboResidual dataDivFcn_noCBO.o dataDivFcn_noCBO windowFFT.o windowFFT compareCalo.o compareCalo plots.o plots

cboAmplitudeBinned: cboAmplitudeBinned.o
	g++ -o cboAmplitudeBinned cboAmplitudeBinned.o $(shell root-config --libs) -L/cvmfs/gm2.opensciencegrid.org/prod9/g-2/gm2util/v9_52_00/slf7.x86_64.e15.prof/lib -lgm2util_blinders -lMinuit

cboAmplitudeBinned.o: cboAmplitudeBinned.cc 
	g++ -c -Wall -Wextra cboAmplitudeBinned.cc $(shell root-config --libs --cflags) -I /cvmfs/gm2.opensciencegrid.org/prod9/g-2/gm2util/v9_52_00/include -ffast-math -O2

