CXX=g++
CXXFLAGS=-g -std=c++11 -Wall -pedantic -O3
LDFLAGS=-g -O3

TARGET=main

all: main

SRCDIR=src
BINDIR=bin
OBJDIR=obj

SRCFILES=main.cpp SmileiCrossSection.cpp SmileiIonizationTables.cpp MyIonizationTables.cpp

OBJFILES:=$(SRCFILES:.cpp=.o)

FULLTARGET=$(BINDIR)/$(TARGET)

# Path to look for source files
VPATH=$(SRCDIR):$(OBJDIR)

main: $(FULLTARGET)


# Rule to build the cpp files
%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) -o $(OBJDIR)/$@ $< $(INCLUDE)

$(FULLTARGET): $(OBJFILES)
	mkdir -p $(BINDIR)
	$(CXX) -o $@ $(addprefix $(OBJDIR)/,$(OBJFILES)) $(LDFLAGS) $(LDLIBS)

clean:
	$(RM) $(OBJDIR)

cleanall: clean
	$(RM) $(BINDIR)

$(OBJFILES): | $(OBJDIR)

$(OBJDIR):
	mkdir -p $(OBJDIR)

.PHONY: clean cleanall main


# All the dependencies
main.o: main.cpp SmileiIonizationTables.h SmileiCrossSection.h MyIonizationTables.h
	$(CXX) -c $(CXXFLAGS) -o $(OBJDIR)/$@ $< $(INCLUDE)

SmileiIonizationTables.o: SmileiIonizationTables.cpp SmileiIonizationTables.h
MyIonizationTables.o: MyIonizationTables.cpp MyIonizationTables.h
SmileiCrossSection.o: SmileiCrossSection.cpp SmileiIonizationTables.h MyIonizationTables.h
