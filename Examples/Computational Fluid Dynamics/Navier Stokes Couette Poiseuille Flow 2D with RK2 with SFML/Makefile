CFLAGS = -ggdb
DEFINES = -DDEBUGGA
INCLUDES = 
LIBS = -lstdc++ -lsymintegration -lsfml-graphics -lsfml-window -lsfml-system
MAIN = main.o
CC=g++
# Source files
SRC=main.cpp

main::	
	$(CC) $(CFLAGS) $(SRC) $(INCLUDES) -o $@ $(EXE) $(LFLAGS) $(LIBS)

clean: 
	rm -f $(MAIN) main

