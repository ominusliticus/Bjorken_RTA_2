DIR_MAIN       = ./
DIR_SRC        = $(DIR_MAIN)build/src
DIR_H          = $(DIR_MAIN)build/include
DIR_BUILD      = $(DIR_MAIN)build
DIR_OBJ        = $(DIR_BUILD)/obj

DEBUG =
OPTIMIZATION = -O2
FLOWTRACE =
OPTIONS = --relocatable-device-code=true
#LINK_OPTIONS = --cudart static --relocatable-device-code=true -link
LINK_OPTIONS = 
CFLAGS = $(DEBUG) $(OPTIMIZATION) $(FLOWTRACE) $(OPTIONS)
COMPILER = nvcc
LIBS = -lm -lgsl -lgslcblas -lconfig
INCLUDES = -I $(DIR_H)

CPP := $(shell find $(DIR_SRC) -name '*.cpp' -not -name '.*')
CU := $(shell find $(DIR_SRC) -name '*.cu' -not -name '.*')
CPP_OBJ  = $(CPP:$(DIR_SRC)%.cpp=$(DIR_OBJ)%.o)
CU_OBJ  = $(CU:$(DIR_SRC)%.cu=$(DIR_OBJ)%.o)
OBJ = $(CPP_OBJ) $(CU_OBJ)

EXE =\
rta-cuda

$(EXE): $(OBJ)
	echo "Linking:   $@ ($(COMPILER))"
	$(COMPILER) $(LINK_OPTIONS) -o $@ $^ $(LIBS) $(INCLUDES)

$(DIR_OBJ)%.o: $(DIR_SRC)%.cpp
	@[ -d $(DIR_OBJ) ] || mkdir -p $(DIR_OBJ)
	@echo "Compiling: $< ($(COMPILER))"
	$(COMPILER) $(CFLAGS) $(INCLUDES) -c -o $@ $<

$(DIR_OBJ)%.o: $(DIR_SRC)%.cu
	@[ -d $(DIR_OBJ) ] || mkdir -p $(DIR_OBJ) 
	@echo "Compiling: $< ($(COMPILER))"
	$(COMPILER) $(CFLAGS) $(INCLUDES) -c -o $@ $<

clean:
	@echo "Object files and executable deleted"
	if [ -d "$(DIR_OBJ)" ]; then rm -rf $(EXE) $(DIR_OBJ)/*; rmdir $(DIR_OBJ); fi

cleandata:
	@echo "Output files deleted"
	rm -f output/*.*
	rm -f output/snapshot/*.*

.SILENT:
