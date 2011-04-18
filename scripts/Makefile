EXE    := parseGamess
EXEDIR := $(addsuffix .d,$(EXE))

all:
	for dir in $(EXEDIR); do \
		$(MAKE) -C $$dir; \
	done

clean:
	for dir in $(EXEDIR); do \
		$(MAKE) clean -C $$dir; \
	done
	rm $(EXE)
