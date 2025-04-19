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

# Available Features

It is able to compute

1. The integral of `sin(ax+b), cos(ax+b), tan(ax+b), cot(ax+b), sec(ax+b), csc(ax+b)` with `ax + b` is a polynomial of order 1.
2. The integral of `1/(ax+b)`
1. The integral and derivative of `asin(ax+b), acos(ax+b), atan(ax+b), acot(ax+b), asec(ax+b), acsc(ax+b)` with `ax + b` is a polynomial of order 1. Only `asec(ax+b), acsc(ax+b)` integration that have no analytic solution yet since it has cases output.

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
|:sunflower:   | Compute integral for all trigonometric	       				| Done
|:writing_hand:| Compute integral for all hyperbolic	        			| Not yet
|:writing_hand:| Compute integral for product, divide and sum of basic functions	| Not yet
|:writing_hand:| Compute integral for inverse trigonometric        			| Not yet
|:writing_hand:| Compute integral with Integration by parts				| Not yet
|:writing_hand:| Compute definite integral with improper integrals			| Not yet

# Milestone

We are using less, minimal amount of code to create the library, if you compare it with the original SymbolicC++ that has `.configure` and `Makefile` that will be spawn after you configure it, we don' t use any of that.

Two simple commands with `g++` can already make a shared symbolic integration computation library, with limitation still.

We have found how to be able to compute `int (sin(2x))` correctly.

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/5.png" width="83%">

By modifying the source code `src/functions.cpp` we can compute the integration for sine and cosine correctly.
By April 10th, 2025: The implementation of sine and cosine have been modified.

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/6.png" width="83%">

For comparison: The basic codes that we are using SymbolicC++, cannot compute the integral of sine and cosine correctly because they only write `return Integral(*this,s)`
<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/7.png" width="83%">

By April 13th, 2025: The integral of `1/(ax+b)` can be computed nicely and fraction / decimal power computation e.g. `(-1)^(a/b)` can be computed too.
<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/8.png" width="83%">

The code for the decimal power computation is located here (the if statement with 4 conditions):
<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/9.png" width="83%">

By April 17th, 2025: We add the implementation for `tan(x)` and `cot(x)` so it can compute the derivative and integral of all basic trigonometry functions (`sin(ax+b), cos(ax+b), tan(ax+b), cot(ax+b), sec(ax+b), csc(ax+b)`) for polynomial of order 1.

The codes that we modified are located in: 

* src/functions.cpp
* src/symintegrationc++.cpp
* include/symintegral/symintegrationc++.h
* include/symintegral/functions.h (to add the Class of `Tan` and `Cot`)

The test code can be located in the folder `Test/Trigonometry/main.cpp`

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/10.png" width="83%">

Comparing the result with JULIA
<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/11.png" width="83%">

The results produce the same output number at the end
<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/12.png" width="83%">

By April 19th, 2025: We can now compute the integral and derivative of `asin(ax+b), acos(ax+b), atan(ax+b), acot(ax+b), asec(ax+b), acsc(ax+b)`, exception for the integration of `asec(ax+b), acsc(ax+b)` they still have no analytic solution yet.

These are the basic integration formula related to the inverse trigonometry


$$\int \frac{1}{\sqrt{1 - x^{2}}} dx = \sin^{-1} (x) + C, -1 < x < 1$$

$$\int - \frac{1}{\sqrt{1 - x^{2}}} dx = \cos^{-1} (x) + C, -1 < x < 1$$

$$\int \frac{1}{1 + x^{2}} dx = \tan^{-1} (x) + C $$

$$\int \frac{1}{|x| \sqrt{x^{2} - 1}} dx = \sec^{-1} (x) + C, |x| > 1 $$

Remember that

$$\sin^{-1} (x) \neq \frac{1}{\sin (x)} $$

$$\sin^{-1} (x) = asin(x)$$


The codes that we modified are located in: 

* src/functions.cpp
* src/symintegrationc++.cpp
* include/symintegral/symintegrationc++.h
* include/symintegral/functions.h (to add the Class of `Atan, Acot, Asin, Acos, Asec, Acsc, Asinh,` and `Acosh`)

The test code can be located in the folder `Test/Inverse Trigonometry/main.cpp`

Comparing the result with SymPy in JULIA:
<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/13.png" width="83%">

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/14.png" width="83%">

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/15.png" width="83%">

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/16.png" width="83%">

The analytic solution for integral of `asec(ax+b` from SymPy in JULIA:
<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/17.png" width="83%">

