DEFINES = -DDEBUG
LIBS  = z gzstream m
LIB_PATH = ../lib
INCLUDES = ../include
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

OBJ      := ${SRC:.cpp=.o}
LIBS     := ${addprefix -l, $(LIBS)}
LIB_PATH := ${addprefix -L, $(LIB_PATH)}
INCLUDES := ${addprefix -I, $(INCLUDES)}

all: staar

%.o : %.cpp
	g++ $(CFLAGS) -c $< -o $@ $(DEFINES) $(INCLUDES)

staar: $(OBJ)
	g++ $(CFLAGS) -o $@ $^ $(LIB_PATH) $(LIBS) $(DEFINES) $(INCLUDES)

clean:
	rm staar *.o