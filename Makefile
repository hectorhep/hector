#    * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#  *                                                         *
# *                   --<--<--  A fast simulator --<--<--     *
# *                 / --<--<--     of particle   --<--<--     *
# *  ----HECTOR----<                                          *
# *                 \ -->-->-- transport through -->-->--     *
# *                   -->-->-- generic beamlines -->-->--     *
# *                                                           *
# * JINST 2:P09005 (2007)                                     *
# *      X Rouby, J de Favereau, K Piotrzkowski (CP3)         *
# *       http://www.fynu.ucl.ac.be/hector.html               *
# *                                                           *
# * Center for Cosmology, Particle Physics and Phenomenology  *
# *              Universite catholique de Louvain             *
# *                 Louvain-la-Neuve, Belgium                 *
#  *                                                         *
#    * * * * * * * * * * * * * * * * * * * * * * * * * * * *


# Makefile - HECTOR
# make : creates the library libHector.dll / libHector.so
# X. Rouby 01.03.2006
# latest update: 15.02.2009

VERSION=$(shell cat VERSION)
# Version number for the tarball. No blank space before and after !

# Warning : by default, the file extension for libraries is .so . However, for Cygwin (MS Windows systems) libraries are called .dll instead of .so . For compatibility, when running on cygwin, just change the LEXT variable value from "so" to "dll" (without the quote symbols).
#
# ----- Usage -----
#  "make what" : prints the possibilities 
#  "make all" or "make" or "make libHector.so": compile libHector library from the sources. (See LIBRARY and LEXT variables)
#  "make clean" : deletes useless files
#  "make depend" or "make Dependencies"  : creates the dependency list for this makefile. (See DEPENDENCIES variable) Called by "make all".
#  "make myfile" : if "myfile.cpp" exists, it will build the executable "myfile". This requires a "main(...)" function in "myfile.cpp"
#  "make myfile.o" : if "myfile.cc" exists, it will build the object file "myfile.o".
#  "make tar" : build the tarball of the sources
#
# In order to use the full power of "make", the list of prerequisites (the dependency files) should be included. 
# With this, "make" sees directly whether a file has to be reprocessed or not. This saves time for the compilation. To avoid any bug, the list of prerequistes is automatically created with "make depend". It is included by the line "-include $(DEPENDECIES)" here below. If the file does not exist when trying to include it, it will be created automatically.

.PHONY: all clean cleanall debug depend what
# avoids any misunderstanding if a file in the current directory exists with one of these names

# ----- Directories -----
# OBJ  contains the object files .o, needed for building of libHector. Is emptied by "make clean"
OBJ= obj/
# SRC  contains the source files .cc .cpp .C with routines and class implementation files.
SRC= src/
# INC  contains the header files .h
INC= include/
# LIB  depository for libHector
LIB= lib/
# ROUTINES contains the routine files .cpp 
ROUTINES= routines/
ROUTINESOBJ= $(OBJ)
ROUTINESSRC= $(ROUTINES)
ROUTINESINC= $(ROUTINES)
#DAT contains some data files
DAT= data/
# vpath  Search path for all prerequisites.
VPATH = .:$(SRC):$(OBJ):$(INC):$(LIB):$(ROUTINES):$(ROUTINESOBJ):$(ROUTINESSRC):$(ROUTINESINC)
#In other words, VPATH is needed by "make" in order to find the dependency files into the file arborescence.

# ----- Library file extensions (dll or so) -----
# for Cygwin, use LEXT= dll
LEXT= .so 

# ----- Filename related variables -----
# LIBRARY  final library name (See LIBFULLNAME and LEXT variables)
LIBRARY= Hector
ROUTINESLIBRARY= $(LIBRARY)_routines
# LIBFULLNAME  Returns the full name of the library (e.g. libHector.so)
LIBFULLNAME = $(addprefix lib,$(addsuffix $(LEXT),$(LIBRARY))) 
ROUTINESLIBFULLNAME = $(addprefix lib,$(addsuffix $(LEXT),$(ROUTINESLIBRARY))) 
# ROOTCFLAGS & ROOTLIBS  flags needed from ROOT for the compilation command
ROOTCFLAGS= -fPIC $(shell root-config --cflags) 
ROOTLIBS= $(shell root-config --libs --glibs) #-lEGPythia8 #-lPythia6 -lEG -lEGPythia6
# HEADERS  List of .h files, without path
HEADERS= $(notdir $(wildcard $(INC)*h))
ROUTINESHEADERS = $(notdir $(wildcard $(ROUTINESINC)*h))
# SOURCES  List of .cc files, without path ; this syntax insures there is one .cc file per .h file
SOURCES= $(addsuffix .cc,$(basename $(notdir $(HEADERS))))
ROUTINESSOURCES= $(addsuffix .cpp,$(basename $(notdir $(ROUTINESHEADERS))))
ROUTINESSOURCESFULLNAME= $(addprefix $(ROUTINESSRC),$(addsuffix .cpp,$(basename $(notdir $(ROUTINESHEADERS)))))
# OBJECTS  List of .o  files, without path ; this syntax insures there is one .o  file per .h file
OBJECTS= $(addsuffix .o, $(basename $(notdir $(HEADERS))))
ROUTINESOBJECTS= $(addsuffix .o, $(basename $(notdir $(ROUTINESHEADERS))))
# DEPENDENCIES  File containing the dependencies for each .o file. It is automaticaly generated.
DEPENDENCIES= Dependencies
ROUTINESDEPENDENCIES= $(DEPENDENCIES)_routines
# TARBALL tar.bz file containing the source code.
TARBALL= $(LIBRARY)$(VERSION).tbz
# WARNINGS  Warning flags for g++
WARNINGS= -Wall -Wno-deprecated -Woverloaded-virtual -Wno-misleading-indentation
# OPTIMIZE= -msse2 -mfpmath=sse -O3 
OPTIMIZE= 
#OPTIMIZE= -pg -g -fprofile-arcs -ftest-coverage



# ----- make all  -----
#all: $(ROUTINESLIBFULLNAME)
all: $(LIBFULLNAME)
# "make all" generates the library. Same as "make" or "make libHector.so"
# Should be left ahead, in order to be the default case (i.e. : "make " = "make all")
# rem : should even be before "include $(DEPENDENCIES)"

# ----- include dependency file -----
ifneq ($(MAKECMDGOALS),clean)
-include $(DEPENDENCIES)
-include $(ROUTINESDEPENDENCIES)
else 
ifneq ($(MAKECMDGOALS),depend)
-include $(DEPENDENCIES)
-include $(ROUTINESDEPENDENCIES)
endif
endif
# Includes the dependencies for each object file. If the file $(DEPENDENCIES) does not exists, "make" creates it. 
# ifneq checks if you are not doing "make clean". If so, you do not need to build $(DEPENDENCIES) as you'll delete it just afterwards.
# ifneq also checks if your are not doing "make depend", as you do not have to build twice the dependence file.

$(OBJ):
	@mkdir $(OBJ)
$(LIB):
	@mkdir $(LIB)
###### Suffix rules (Tells how to build a.o file from a.cc file)
# Building the object files from .cc files
%.o : %.cc | $(OBJ)
	@echo Making $@
	@g++ $< -c $(ROOTCFLAGS) -I$(INC) $(WARNINGS) -g -o $@ $(OPTIMIZE)
	@mv -f $@ $(OBJ)
	@rm -f $(addsuffix .d, $(basename $@))
# This suffix rules requires that the first dependency file for a.o is a.cc
# e.g. : H_Aperture.o: H_Aperture.cc ..and then the other dependencies..
# rem : the macro $< returns the first dependency
# rem : the last line removes the a.d file created with a.o

# Building the object files from .cpp files, which will be included in $(ROUTINESLIBRARY)
%.o : %.cpp | $(OBJ)
	@echo Making $@
#	@g++ $< -c $(ROOTCFLAGS) -I$(INC) -I$(ROUTINESINC) $(WARNINGS) -o $@ $(OPTIMIZE)
	@g++ $< -c $(ROOTCFLAGS) -I$(INC) $(WARNINGS) -o $@ $(OPTIMIZE)
	@mv -f $@ $(OBJ)
	@rm -f $(addsuffix .d, $(basename $@))


# Building the executable files from .cpp files
% : %.cpp $(LIBFULLNAME) $(ROUTINESLIBFULLNAME)
	@echo Compiling $@
#	@g++ $< -o $@ $(ROOTLIBS) $(WARNINGS) $(ROOTCFLAGS) -L$(LIB) -l$(LIBRARY) -l$(ROUTINESLIBRARY) -I$(INC) -I$(ROUTINESINC) $(OPTIMIZE)
	@g++ $< -o $@ $(ROOTLIBS) $(WARNINGS) $(ROOTCFLAGS) -L$(LIB) -l$(LIBRARY) -I$(INC) $(OPTIMIZE) -Wl,-Rlib/
	@cp $(LIB)$(LIBFULLNAME) .
	@rm -f $@.d

# ----- make libHector.so -----
$(LIBFULLNAME): $(OBJECTS) | $(LIB)
	@echo Making $@ 
	@g++ -shared $(ROOTLIBS) $(WARNINGS) -o $@ -Wl,-soname,$@ $(addprefix $(OBJ),$(OBJECTS)) $(OPTIMIZE)
	@strip -s $@
	@mv -f $@ $(LIB) 
	@echo Done : `ls $(LIB)$@`
	@echo

$(ROUTINESLIBFULLNAME): $(ROUTINESOBJECTS) $(LIBFULLNAME) $(ROUTINESSOURCESFULLNAME)
	@echo Making $@
	@g++ -shared $(ROOTLIBS) $(WARNINGS) -o $@ -Wl,-soname,$@ $(addprefix $(ROUTINESOBJ),$(ROUTINESOBJECTS)) $(OPTIMIZE)
	@strip -s $@
	@mv -f $@ $(LIB)
	@echo Done : `ls $(LIB)lib$(ROUTINESLIBRARY)*`
	@echo

# When "undefined reference to" errors occur, check the order of the object files.
# Depends on the object files
# With Mac OS X use instead
#@g++ -dynamiclib $(ROOTLIBS) $(WARNINGS) -o $@ $(addprefix $(OBJ),$(OBJECTS)) $(OPTIMIZE)

# ----- make debug ------
debug : $(OBJECTS)
	@echo Making $(LIBFULLNAME) with debug information
	@g++ -shared $(ROOTLIBS) $(WARNINGS) -o $(LIBFULLNAME) -Wl,-soname,$(LIBFULLNAME) $(addprefix $(OBJ),$(OBJECTS)) $(OPTIMIZE)
	@mv -f $(LIBFULLNAME) $(LIB) 
	@echo Done : `ls $(LIB)lib*`
	@echo
#@g++ -dynamiclib $(ROOTLIBS) $(WARNINGS) -o $@ $(addprefix $(OBJ),$(OBJECTS)) $(OPTIMIZE)

# ----- make depend -----
$(DEPENDENCIES) depend : $(HEADERS) $(SOURCES) 
	@echo Making the dependency file : $(DEPENDENCIES) 
	@g++ -MM $(ROOTCFLAGS) $(SRC)/H_*.cc $(WARNINGS) -I$(INC) -I$(ROOTSYS)/include > $(DEPENDENCIES) 

$(ROUTINESDEPENDENCIES) : $(ROUTINESHEADERS) $(HEADERS)
	@echo Making the dependency file : $(ROUTINESDEPENDENCIES) 
	@echo $(ROUTINESSOURCESFULLNAME)
	@g++ -MM $(ROOTCFLAGS) $(ROUTINESSOURCESFULLNAME) $(WARNINGS) -I$(INC) -I$(ROUTINESINC) -I$(ROOTSYS)/include > $(ROUTINESDEPENDENCIES) 
#Creates the file $(DEPENDENCIES)
#"g++ -MM " produces the same output as "makedepend" but with the $@.cc file added to the list, which is needed here

# ----- make clean -----
clean:
	@echo Deleting object and temporary files
	@-rm -f *.d *~ core $(OBJ)*.o Doxywarn.txt $(DEPENDENCIES) $(ROUTINESDEPENDENCIES) 
# the dash "-" in front avoids an warning if $(DEPENDENCIES) does not exist

cleanall: clean
	@-rm -f ./H_* $(ROUTINES)/H_*.so $(LIB)$(ROUTINESLIBFULLNAME) $(LIB)$(LIBFULLNAME)
	@-rm -f chromaticity_grid_*.txt gmon.out doxywarn.txt RP_acceptance_* m7000.txt $(LIBFULLNAME) output.valgrind 
	@-rm -f *.gcov *.gcda *.gcno
	@echo

distclean: cleanall
	
# ----- make tar -----
$(TARBALL) tar : Makefile README rootlogon.C VERSION $(HEADERS) $(SOURCES) $(notdir $(wildcard $(ROUTINES)H_*cpp)) 
	@echo Building tarball of sources
	@mkdir $(LIBRARY)
	@cp Makefile README rootlogon.C VERSION $(LIBRARY)/
	@mkdir $(LIBRARY)/$(INC)
	@cp $(INC)/*h $(LIBRARY)/$(INC)/
	@mkdir $(LIBRARY)/$(SRC)
	@cp $(SRC)/*cc $(LIBRARY)/$(SRC)/
	@mkdir $(LIBRARY)/$(DAT)
	@cp $(DAT)/*.tfs $(LIBRARY)/$(DAT)/
	@cp $(DAT)/aah_gammae.root $(DAT)/lpair_mumu_2gev.root $(LIBRARY)/$(DAT)/
	@mkdir $(LIBRARY)/$(ROUTINES)
	@cp $(ROUTINES)/H_*cpp $(LIBRARY)/$(ROUTINES)/
	@cp $(ROUTINES)/H_*h $(LIBRARY)/$(ROUTINES)/
	@mkdir $(LIBRARY)/$(OBJ)
	@mkdir $(LIBRARY)/$(LIB)
	@tar cjf $(TARBALL) $(LIBRARY)
	@rm -rf $(LIBRARY)
	@echo Done : `ls $(TARBALL)`

# ----- make what -----
what:
	@echo
	@echo "Hector version = $(VERSION)"
	@echo "Usage : "
	@echo " 'make' or 'make all' or 'make $(LIBFULLNAME)': creates the library file ($(LIBFULLNAME))"
	@echo " 'make debug' : creates $(LIBFULLNAME), including  the debug information"
	@echo " 'make clean' or 'make cleanall' : deletes the useless files"
	@echo " 'make myfile' : if 'myfile.cpp' exists, it will build the executable 'myfile'. This requires a 'main(...)' function in 'myfile.cpp'"
	@echo " 'make tar' : builds the tarball of the sources"
	@echo " 'make what' : prints this help"
	@echo
#	@echo Others, but you do not need to call them
#	@echo " 'make depend' or 'make $(DEPENDENCIES)' : creates the dependency list for this makefile. (done by 'make all')"
#	@echo "	'make myfile.o' : if 'myfile.cc' exists, it will build 'myfile.o'"
