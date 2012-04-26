NAME := staar

CC := gcc
CPP := g++
LINK := $(CPP)

MODULES := src include

CFLAGS   := -O3
CPPFLAGS := -O3
LFLAGS   := -O3
LIBS := -lm -lz -lgzstream -lopenbabel -Llib 

# Edit these!
BABEL_LIBDIR := /lustre/AQ/openbabel-2.3.0/build/lib/
BABEL_INCDIR := /lustre/AQ/openbabel-2.3.0/include /data/AQ/openbabel-2.3.0/build/include

ifdef BABEL_LIBDIR
LIBS += -L$(BABEL_LIBDIR)
endif


# You shouldn't have to go below here

DIRNAME = `dirname $1`
MAKEDEPS = $(CC) -MM -MG $2 $3 | sed -e "s@^\(.*\)\.o:@.dep/$1/\1.d obj/$1/\1.o:@"

.PHONY : all

all : $(NAME)

# look for include files in each of the modules
INCLUDEFLAGS := $(patsubst %, -I%, $(MODULES)) $(addprefix -I,$(BABEL_INCDIR))

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
$(NAME) : $(OBJ) gzip
	$(LINK) $(LFLAGS) -o $@ $(OBJ) $(LIBS)

# calculate C include dependencies
.dep/%.d : %.cpp
	@mkdir -p `echo '$@' | sed -e 's|/[^/]*.d$$||'`
	$(call MAKEDEPS,$(call DIRNAME, $<), $(CFLAGS), $<) > $@

.dep/%.d : %.C
	@mkdir -p `echo '$@' | sed -e 's|/[^/]*.d$$||'`
	$(call MAKEDEPS,$(call DIRNAME, $<), $(CFLAGS), $<) > $@

obj/%.o : %.cpp
	@mkdir -p `echo '$@' | sed -e 's|/[^/]*.o$$||'`
	$(CPP) $(CPPFLAGS) -c -o $@ $<

obj/%.o : %.C
	@mkdir -p `echo '$@' | sed -e 's|/[^/]*.o$$||'`
	$(CPP) $(CPPFLAGS) -c -o $@ $<

gzip : gzstream/gzstream.C
	mkdir -p lib
	$(MAKE) -C gzstream install

# include the C include dependencies
DEP := $(patsubst obj/%.o, .dep/%.d, $(OBJ))

ifneq ($(MAKECMDGOALS),clean)
-include $(DEP)
endif

clean :
	-@rm $(NAME) $(OBJ) $(DEP) lib/libgzstream.a
