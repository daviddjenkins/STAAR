DEFINES = 
LIBS  = 
CFLAGS = -O3

SRC := staar.cpp

OBJ := ${SRC:.cpp=.o}

all: staar

%.o : %.cpp
	g++ $(CFLAGS) -c $< -o $@ $(DEFINES)

staar: $(OBJ)
	g++ $(CFLAGS) -o $@ $^ $(LIBS) $(DEFINES)

clean:
	rm staar *.o