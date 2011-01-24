DIRECTORIES=src

all:
	for dir in $(DIRECTORIES); do \
                $(MAKE) -C $$dir; \
        done
	cp src/staar .

clean:
	for dir in $(DIRECTORIES); do \
                $(MAKE) clean -C $$dir; \
        done
	rm *~ staar