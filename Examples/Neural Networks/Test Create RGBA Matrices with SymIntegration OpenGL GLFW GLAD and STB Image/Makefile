CFLAGS = -ggdb
DEFINES = -DDEBUGGA
INCLUDES = 
LIBS = -lstdc++ -lGLEW -lGL -lglfw /root/SourceCodes/CPP/src/glad.c -lsymintegration
MAIN = main.o
CC=g++

.cc.o:
	$(CC) -c $(CFLAGS) $(DEFINES) $(INCLUDES) $<

all::	main

gnuplot_i.o:	
main.o:	main.cpp

main: $(MAIN)
	$(CC) -o $@ $(CFLAGS) $(MAIN) $(LIBS) $(INCLUDES)

clean: 
	rm -f $(MAIN) main
