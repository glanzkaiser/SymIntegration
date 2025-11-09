
CXX = g++

SYMINTEGRATION_DIR = ../src
LIB = libsymintegration.so
SOURCES = array.cpp
SOURCES += cloning.cpp constants.cpp derive.cpp divgradcurl.cpp dsolve.cpp dynamicalsystem.cpp equation.cpp functions.cpp
SOURCES += geometryandvectors.cpp geompack.cpp identity.cpp integrate.cpp linearalgebra.cpp 
SOURCES += mathematicalfinance.cpp matnorm.cpp matrix.cpp multinomial.cpp navierstokes2d.cpp number.cpp numericalmethod.cpp
SOURCES += operationsresearch.cpp polynomial.cpp product.cpp quatern.cpp rational.cpp simplify.cpp
SOURCES += solve.cpp specialfunctions.cpp statistics.cpp sum.cpp symbol.cpp symbolic.cpp symerror.cpp 
SOURCES += symintegrationc++.cpp symmatrix.cpp vecnorm.cpp 
SOURCES += vector.cpp verylong.cpp
SOURCES += asa241.cpp

OBJS = $(addsuffix .o, $(basename $(notdir $(SOURCES))))
UNAME_S := $(shell uname -s)

CXXFLAGS = -fPIC
CXXFLAGS += -g -Wall -Wformat

##---------------------------------------------------------------------
## BUILD FLAGS PER PLATFORM
##---------------------------------------------------------------------

ifeq ($(UNAME_S), Linux) #LINUX
	ECHO_MESSAGE = "SymIntegration from GFreya OS"
	
	CFLAGS = $(CXXFLAGS)
endif
##---------------------------------------------------------------------
## BUILD RULES
##---------------------------------------------------------------------

%.o:$(SYMINTEGRATION_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

all: $(LIB)
	@echo Build complete for $(ECHO_MESSAGE)

$(LIB): $(OBJS)
	$(CXX) -shared -o $@ $^ $(CXXFLAGS) 

clean:
	rm -f $(OBJS)
