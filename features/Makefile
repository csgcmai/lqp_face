OPTIMIZATION_FLAG=-O3

CFLAGS =$(INCLUDEDIR) $(OPTIMIZATION_FLAG) -std=c++0x -fmessage-length=0 -Wno-write-strings -DONLINE_MAPPING  #-Wall# 
#ONLINE_MAPPING:  uses unordered_maps for storing a code's mapping to nearest cluster center (computed on fly); if absent by default lookup table will be used for storing the code's cluster information, where the lookup table is populated offline. LUT gives constant access time while unordered_maps give constant access time on average.

# Image Magick Include and LIBS, these can be found by calling 'Magick++-config' with arguments --cppflags and --ldflags  

IMAGEMAGICK_INCLUDE=`Magick++-config --cppflags` #/data/image/hussain/softwares/ImageMagick/include/ImageMagick/ #
IMAGEMAGICK_LIBS=`Magick++-config --ldflags --libs` # -lMagickCore -lMagick++  

LIBDIR =-L/usr/lib/ # -L/data/image/hussain/softwares/ImageMagick/lib   

INCLUDEDIR =-I ../util/  $(IMAGEMAGICK_INCLUDE) #-I ../Kmeans/  
LIBS = -lboost_program_options  -lboost_system -lboost_filesystem  -lboost_filesystem-mt $(IMAGEMAGICK_LIBS) #-lboost_thread-mt

TARGET=mainFeatures
OUTPUTDIR=../build/
BUILDDIR=$(OUTPUTDIR)features$(OPTIMIZATION_FLAG)/


#Some Reference Directories...
UTILDIR=../util/
FEATDIR=./
MPI_KMEANSDIR=../MPI_KMeans/
MPI_KMEANSOBJS=$(BUILDDIR)/mpi_kmeans.o

UTIL = npyramid.cpp ncoord.cpp processImage.cpp fileIO.cpp hogInfo.cpp windowInfo.cpp util.cpp classInfo.cpp
TOOLS= lbpMap.cpp lbpFeaturesRGB.cpp ltpFeaturesRGB.cpp \
lbpFeaturesRGBBlock.cpp ltpFeaturesRGBCell.cpp onlyLTPFeaturesRGBCell.cpp encodedPatchFeatures.cpp\
encodeFeatures.cpp multipleCBQuantizedLQPFeatures.cpp lbpFeatures.cpp ltpFeatures.cpp onlyLTPFeatures.cpp customValidators.cpp
 
TOOLSOBJS=$(TOOLS:.cpp=.o)
INCLUDEOBJS=$(INCLUDE:.cpp=.o)
UTILOBJS=$(UTIL:.cpp=.o)

UTILOBJS := $(addprefix $(BUILDDIR)/,$(UTILOBJS))
TOOLSOBJS := $(addprefix $(BUILDDIR)/,$(TOOLSOBJS))
INCLUDEOBJS := $(addprefix $(BUILDDIR)/,$(INCLUDEOBJS))

checkdir:
	mkdir -p $(BUILDDIR); echo $(BUILDDIR);

mainFeatures:$(TARGET).cpp $(UTILOBJS)  $(TOOLSOBJS)  $(MPI_KMEANSOBJS)
	$(CXX) $(CFLAGS) -o $(OUTPUTDIR)$(TARGET) $(TARGET).cpp $(UTILOBJS) $(LIBS) $(TOOLSOBJS)  $(MPI_KMEANSOBJS)

all: checkdir mainFeatures #$(TARGET)

#TOOLS OBJS

#%.o: %.cpp
#	$(CXX) -c $(CFLAGS) $< -o $@
#    
    
#UTILOBJS
$(BUILDDIR)/util.o: $(UTILDIR)/util.cpp
	$(CXX) -c $(CFLAGS)  $< -o $@ 
$(BUILDDIR)/npyramid.o: $(UTILDIR)/npyramid.cpp
	$(CXX) -c $(CFLAGS)  $< -o $@
$(BUILDDIR)/ncoord.o: $(UTILDIR)/ncoord.cpp
	$(CXX) -c $(CFLAGS)  $< -o $@
$(BUILDDIR)/windowInfo.o: $(UTILDIR)/windowInfo.cpp
	$(CXX) -c $(CFLAGS)  $< -o $@
$(BUILDDIR)/processImage.o: $(UTILDIR)/processImage.cpp
	$(CXX) -c $(CFLAGS)  $< -o $@
$(BUILDDIR)/fileIO.o: $(UTILDIR)/fileIO.cpp
	$(CXX) -c $(CFLAGS)  $< -o $@	
$(BUILDDIR)/hogInfo.o: $(UTILDIR)/hogInfo.cpp
	$(CXX) -c $(CFLAGS)  $< -o $@
$(BUILDDIR)/classInfo.o:$(UTILDIR)/classInfo.cpp $(UTILDIR)/classInfo.h
	$(CXX) -c $(CFLAGS) $< -o $@
	
$(BUILDDIR)/mpi_kmeans.o: $(MPI_KMEANSDIR)/mpi_kmeans.cxx $(MPI_KMEANSDIR)/mpi_kmeans.h
	$(CXX) -c $(CFLAGS)  $< -o $@

#FEATURES
$(BUILDDIR)/lbpMap.o:$(FEATDIR)/lbpMap.cpp $(FEATDIR)/lbpMap.h
	$(CXX) -c $(CFLAGS) $< -o $@
$(BUILDDIR)/lbpFeatures.o:$(FEATDIR)/lbpFeatures.cpp $(FEATDIR)/lbpFeatures.h
	$(CXX) -c $(CFLAGS) $< -o $@
$(BUILDDIR)/ltpFeatures.o:$(FEATDIR)/ltpFeatures.cpp $(FEATDIR)/ltpFeatures.h
	$(CXX) -c $(CFLAGS) $< -o $@
$(BUILDDIR)/onlyLTPFeatures.o:$(FEATDIR)/onlyLTPFeatures.cpp $(FEATDIR)/onlyLTPFeatures.h
	$(CXX) -c $(CFLAGS) $< -o $@	
$(BUILDDIR)/lbpFeaturesRGB.o:$(FEATDIR)/lbpFeaturesRGB.cpp $(FEATDIR)/lbpFeaturesRGB.h
	$(CXX) -c $(CFLAGS) $< -o $@	
$(BUILDDIR)/ltpFeaturesRGB.o:$(FEATDIR)/ltpFeaturesRGB.cpp $(FEATDIR)/ltpFeaturesRGB.h
	$(CXX) -c $(CFLAGS) $< -o $@
$(BUILDDIR)/ltpFeaturesRGBCell.o:$(FEATDIR)/ltpFeaturesRGBCell.cpp $(FEATDIR)/ltpFeaturesRGBCell.h
	$(CXX) -c $(CFLAGS) $< -o $@

$(BUILDDIR)/lbpFeaturesRGBBlock.o:$(FEATDIR)/lbpFeaturesRGBBlock.cpp $(FEATDIR)/lbpFeaturesRGBBlock.h
	$(CXX) -c $(CFLAGS) $< -o $@
$(BUILDDIR)/onlyLTPFeaturesRGBCell.o:$(FEATDIR)/onlyLTPFeaturesRGBCell.cpp $(FEATDIR)/onlyLTPFeaturesRGBCell.h
	$(CXX) -c $(CFLAGS) $< -o $@
$(BUILDDIR)/encodeFeatures.o:$(FEATDIR)/encodeFeatures.cpp $(FEATDIR)/encodeFeatures.h  
	$(CXX) -c $(CFLAGS) $< -o $@	
$(BUILDDIR)/encodedPatchFeatures.o:$(FEATDIR)/encodedPatchFeatures.cpp $(FEATDIR)/encodedPatchFeatures.h  $(BUILDDIR)/encodeFeatures.o
	$(CXX) -c $(CFLAGS) $< -o $@
$(BUILDDIR)/multipleCBQuantizedLQPFeatures.o:$(FEATDIR)/multipleCBQuantizedLQPFeatures.cpp $(FEATDIR)/multipleCBQuantizedLQPFeatures.h $(BUILDDIR)/encodedPatchFeatures.o  $(BUILDDIR)/encodeFeatures.o 
	$(CXX) -c $(CFLAGS) $< -o $@	
$(BUILDDIR)/customValidators.o:$(FEATDIR)/customValidators.cpp
	$(CXX) -c $(CFLAGS) $< -o $@	

clean:
	rm -f $(OBJS) $(TARGET)  $(TOOLSOBJS)  *.o *.gch $(UTILOBJS)


