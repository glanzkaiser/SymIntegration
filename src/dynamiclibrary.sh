g++ -fPIC -c *.cpp
g++ -shared -o libsymintegration.so *.o  # create *.so shared library then move it to /usr/lib
