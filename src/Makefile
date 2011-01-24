DEFINES = 
LIBS  = 
INCLUDES = -I../include
CFLAGS = -O3

SRC := staar.cpp        \
      Atom.cpp          \
      Chain.cpp         \
      Options.cpp       \
      PDB.cpp           \
      Seqres.cpp        \
      Utils.cpp         \
      Coordinates.cpp   \
      AminoAcid.cpp     \
      Geometry.cpp

OBJ := ${SRC:.cpp=.o}

all: staar

%.o : %.cpp
	g++ $(CFLAGS) -c $< -o $@ $(DEFINES) $(INCLUDES)

staar: $(OBJ)
	g++ $(CFLAGS) -o $@ $^ $(LIBS) $(DEFINES) $(INCLUDES)

clean:
	rm staar *.o