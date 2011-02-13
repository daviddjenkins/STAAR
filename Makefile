NAME := staar

CC := gcc
CPP := g++
LINK := $(CPP)

MODULES := src include

CFLAGS := -O3
LFLAGS := -O3
LIBS := -lm -lz -lgzstream -Llib

# You shouldn't have to go below here

DIRNAME = `dirname $1`
MAKEDEPS = gcc -MM -MG $2 $3 | sed -e "s@^\(.*\)\.o:@.dep/$1/\1.d obj/$1/\1.o:@"

.PHONY : all

all : $(NAME)

# look for include files in each of the modules
INCLUDEFLAGS := $(patsubst %, -I%, $(MODULES))

CFLAGS += $(INCLUDEFLAGS)
CPPFLAGS += $(INCLUDEFLAGS)

# each module will add to this
SRC := $(wildcard $(patsubst %, %/*.cpp, $(MODULES))) \
        $(wildcard $(patsubst %, %/*.c, $(MODULES)))

# determine the object files
OBJ := $(patsubst %.cpp, obj/%.o, $(filter %.cpp, $(SRC))) \
         $(patsubst %.c, obj/%.o, $(filter %.c, $(SRC))) \
         $(patsubst %.y, obj/%.o, $(filter %.y, $(SRC)))

# link the program
$(NAME) : $(OBJ)
	$(LINK) $(LFLAGS) -o $@ $(OBJ) $(LIBS)

# calculate C include dependencies
.dep/%.d : %.cpp
	@mkdir -p `echo '$@' | sed -e 's|/[^/]*.d$$||'`
	$(call MAKEDEPS,$(call DIRNAME, $<), $(CFLAGS), $<) > $@

obj/%.o : %.cpp
	@mkdir -p `echo '$@' | sed -e 's|/[^/]*.o$$||'`
	$(CPP) $(CPPFLAGS) -c -o $@ $<

# include the C include dependencies
DEP := $(patsubst obj/%.o, .dep/%.d, $(OBJ))

include $(DEP)

clean :
	-@rm $(NAME) $(OBJ) $(DEP)
