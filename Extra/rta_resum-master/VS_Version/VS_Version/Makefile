
UNAME := $(shell uname)

DIR_MAIN = ./
DIR_SRC = $(DIR_MAIN)
DIR_OBJ = $(DIR_MAIN)

CFLAGS= $(OPTIMIZATION) $(OPTIONS)

COMPILER = gcc
OPTIMIZATION = -O3
CFLAGS := $(CFLAG) -Wno-comment

# Library options
ifeq ($(UNAME), Linux)
LIBS := -lm -L/home/ominusliticus/gsl/lib
endif
ifeq ($(UNAME), Darwin)
LIBS := -lc++ -lm -lgsl
endif

INCLUDES = -I include -I/home/ominusliticus/gsl/include

CPP := $(shell find $(DIR_SRC) -name '*.cpp')
OBJ =$(CPP:$(DIR_SRC)%.cpp=$(DIR_OBJ)%.o)

EXE = rta

# $@ = target
# $< = dependency
# $^ = + everything else

$(EXE): $(OBJ)
	echo "\nLinking: $@ ($(COMPILER))"
	$(COMPILER) -o $@ $^ $(LIBS) $(INCLUDES)
	echo "\nRunning...\n"
	$(DIR_MAIN)$(EXE)
	rm *.o

$(DIR_OBJ)%.o: $(DIR_SRC)%.cpp
	echo "Compiling $< ($(COMPILER))"
	$(COMPILER) $(CFLAGS) $(INCLUDES) -c -o $@ $<

# run:
# 	echo "\nRunning...\n"
# 	$(DIR_MAIN)$(EXE)

clean:
	echo "\nDeleting executable and object files\n"
	rm -rf $(EXE)

.SILENT:

