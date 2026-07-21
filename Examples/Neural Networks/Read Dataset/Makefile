CFLAGS = -ggdb 
DEFINES = -DDEBUGGA
INCLUDES = 
LIBS = -lstdc++ -lsymintegration
MAIN = predict.o
MAIN2 = train.o
CC=g++

.cc.o:
	$(CC) -c $(CFLAGS) $(DEFINES) $(INCLUDES) $<

all::	train predict

train.o:	train.cpp

predict.o:	predict.cpp 

train: $(MAIN2)
	$(CC) -o $@ $(CFLAGS) $(MAIN2) $(LIBS)

predict: $(MAIN)
	$(CC) -o $@ $(CFLAGS) $(MAIN) $(LIBS)

clean: 
	rm -f $(MAIN) predict $(MAIN2) train
