# SymIntegration

SymIntegration is a C++ library that is branching out from SymbolicC++3.35.
The main idea is to improve its symbolic integration codes. Our main focus is to make it able to compute:

1. All kinds of standard integral form (trigonometry, inverse trigonometry, polynomial, transcendental, hyperbolic)
2. The sum, product and divide combination of the standard functions
3. To be able to compute improper integrals with cases (e.g. computing mean and variance for exponential distribution)

## About SymbolicC++
SymbolicC++ was originally written as a collection of header files
for C++. Many of the functions and classes provided are template
functions and classes and in general cannot be compiled in as part
of a library.

This project attempts to extract the parts of SymbolicC++ that can
be compiled as part of a library and so create the include / library
infrastructure. The src and include directories are populated by
scripts from the SymbolicC++ header files.

# Learning the Code

The C++ codes are only located in the `/src` directory , it is designed to contain all the `.cpp` files, then inside the `/include/` we have header files and another folder `/include/symintegral/` is also a folder that contain header files too.

You can open them one by one, read them and making sense of it based on the function that you are looking for.

There are total of 27 `.cpp` files and 26 `.h` / header files. So it won't take a long time for someone to learn it.

# Installation
Assuming you are using Linux ( we are using GFreya OS based on LFS version 11.0 System V) then you should have no problem in following these methods.

## Move all .h / Header Files to /usr/include

Open terminal and from the current working directory / this repository main directory:

```
	cd include
	cp -r * /usr/include
```

When we want to create the shared library it will look for this header files in the default path where the include files are usually located. In Linux OS `usr/include` is the basic / default path.

## Create the Shared Library
We will create the shared library from all the `.cpp` files in `/src` directory.

Open terminal and from the current working directory / this repository main directory:

```
	cd src
	g++ -fPIC -c *.cpp
	g++ -shared -o libsymintegration.so *.o
```

Then move / copy it into `/usr/lib/` since it is the default path for libraries for Linux OS.

The size of libsymbolicc++.so is about 5 MB compared to libsymintegration.so that is only 1.1 MB.

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/1.png" width="83%">

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/2.png" width="83%">

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/3.png" width="83%">

# Test Example

Open terminal and from the current working directory / this repository main directory:

```
	cd Test
	make
	./main
```

If you prefer the old way then you can compile the old way / type `g++ -o main main.cpp -lsymintegration` instead of typing `make`.

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/5.png" width="83%">

# Status

:sunflower: = Done

:writing_hand: = On Progress

| Status | Name | Details |
| -------------     | ------------- | ------------- | 
|:sunflower:   | Compute integral of sine and cosine  					| Done
|:writing_hand:| Compute integral for all trigonometric	       				| Not yet
|:writing_hand:| Compute integral for all hyperbolic	        			| Not yet
|:writing_hand:| Compute integral for product, divide and sum of basic functions	| Not yet
|:writing_hand:| Compute integral for inverse trigonometric        			| Not yet
|:writing_hand:| Compute integral with Integration by parts				| Not yet
|:writing_hand:| Compute definite integral with improper integrals			| Not yet

# Remark

We are using less, minimal amount of code to create the library, if you compare it with the original SymbolicC++ that has `.configure` and `Makefile` that will be spawn after you configure it, we don' t use any of that.

Two simple commands with `g++` can already make a shared symbolic integration computation library, with limitation still.

We have found how to be able to compute `int (sin(2x))` correctly.

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/5.png" width="83%">

By modifying the source code `src/functions.cpp` we can compute the integration for sine and cosine correctly.
By April 10th, 2025: The implementation of sine and cosine have been modified.

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/6.png" width="83%">

For comparison: The basic codes that we are using SymbolicC++, cannot compute the integral of sine and cosine correctly because they only write `return Integral(*this,s)`
<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/7.png" width="83%">

