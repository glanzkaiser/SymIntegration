
CXX = g++

SYMINTEGRATION_DIR = ../src
LIB = libsymintegration.so
SOURCES = array.cpp
SOURCES += cloning.cpp constants.cpp derive.cpp dsolve.cpp equation.cpp functions.cpp
SOURCES += identity.cpp integrate.cpp matnorm.cpp matrix.cpp multinomial.cpp number.cpp
SOURCES += polynomial.cpp product.cpp quatern.cpp rational.cpp solve.cpp sum.cpp symbol.cpp
SOURCES += symbolic.cpp symerror.cpp symintegrationc++.cpp symmatrix.cpp vecnorm.cpp 
SOURCES += vector.cpp verylong.cpp


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
