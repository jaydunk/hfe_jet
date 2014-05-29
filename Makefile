PROGRAM  = ucla_ee
PYTHIAPATH   = /Users/jaydunkelberger/PYTHIA/pythia8162
LHAPDFPATH   = /Users/jaydunkelberger/PYTHIA/lhapdf-5.8.7/lib
CXX      =  g++
CXXFLAGS = -O  -W -Wall
CPPFLAGS = -I$(PYTHIAPATH)/include -I$(ROOTSYS)/include
LDFLAGS  = -L$(PYTHIAPATH)/lib/archive -L$(ROOTSYS)/lib -L$(LHAPDFPATH) -lLHAPDF -lpythia8 -llhapdfdummy -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lpthread -lm -ldl -rdynamic -lThread -lMathCore

$(PROGRAM):	$(PROGRAM).cpp 
			$(CXX) $(CXXFLAGS)  $(PROGRAM).cpp $(CPPFLAGS) $(LDFLAGS) -o $(PROGRAM) 

clean:
			rm -f $(PROGRAM) core *.o

