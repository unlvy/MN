CXX = gcc
EXEC = main
CXXFLAGS += -pedantic -Wall -MMD -MP
LFLAGS += -lm

SRC=$(wildcard *.c)
OBJ=$(SRC:.c=.o)
DEP=$(SRC:.c=.d)

$(EXEC): $(OBJ) 
	$(CXX) $^ -o $@ $(LFLAGS)

.PHONY: clean

clean:
	@rm -f $(EXEC) $(OBJ) $(DEP)

-include $(DEP)