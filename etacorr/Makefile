CXX = g++
#INCLUDES := $(root-config --incdir)
ROOTCFLAGS := $(root-config --cflags)
ROOTLIBS := $(root-config --libs)

etacorr: etacorr.o
	$(CXX) -L/wsu/home/fr/fr45/fr4526/build/root_v5.34.34/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic etacorr.o -o etacorr 

etacorr.o: etacorr.cxx
	$(CXX) -c -pthread -m64 -I/wsu/home/fr/fr45/fr4526/build/root_v5.34.34/include etacorr.cxx

