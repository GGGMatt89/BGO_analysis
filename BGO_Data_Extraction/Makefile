#---------------------------------------
# Resources section - Symbol definitions
#---------------------------------------
ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --libs)
#THREAD = -lpthread

CXX := g++
#CXX := clang++
CXXFLAGS := -std=c++11 -g #-stdlib=libstdc++
#CXXFLAGS := -std=c++11 -stdlib=libc++


C++ = clang -O
LIB = BGO_Data_extraction.a
AR  = ar -r

#----------------------
# Makefile head section
#----------------------

all:
	make BGO_Data_extraction

BGO_Data_extraction: BGO_Data_extraction.o $(LIB)
	$(CXX) $(CXXFLAGS) -o BGO_Data_extraction BGO_Data_extraction.o $(ROOTFLAGS) $(ROOTLIBS)

clean:
	rm -f *.o *.a

cleanall:
	rm -f *.o *.a BGO_Data_extraction


#----------------------------------
# Dependencies for the source files
#----------------------------------

ANALYSIS_FROMBIN_CPP_DEP = \
	functions_new.h \

#--------------------------------
# Compilation of the source files
#--------------------------------

BGO_Data_extraction.o : BGO_Data_extraction.cxx $(CREATE_DATA_CPP_DEP)
	$(CXX) $(CXXFLAGS) -c BGO_Data_extraction.cxx -o BGO_Data_extraction.o $(ROOTFLAGS)

#----------------------------------
# Library dependencies and creation
#----------------------------------

LIB_DEP = \
	BGO_Data_extraction.o \

$(LIB) : $(LIB_DEP)
	$(AR) $(LIB) $(LIB_DEP) $(POST_AR)
