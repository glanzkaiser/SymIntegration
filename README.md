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
2. The integral of `1/(ax+b)`.
3. The integral and derivative of `asin(ax+b), acos(ax+b), atan(ax+b), acot(ax+b), asec(ax+b), acsc(ax+b)` with `ax + b` is a polynomial of order 1. Only `asec(ax+b), acsc(ax+b)` integration that have no analytic solution yet since it has cases output.
4. The integral and derivative of all hyperbolic trigonometry functions`sinh(ax+b), cosh(ax+b), tanh(ax+b), coth(ax+b), sech(ax+b), csch(ax+b)`.
5. The integral with the form of `x^{b} exp(a*x), x^{b} exp(x)` with integration by parts method.

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
|:sunflower:   | Compute integral for all hyperbolic	        			| Done
|:sunflower:   | Compute integral for product, divide and sum of basic functions	| Done
|:sunflower:   | Compute integral for inverse trigonometric        			| Done
|:sunflower:   | Compute integral with Integration by parts				| Done
|:writing_hand:| Compute definite integral with improper integrals			| How to substitute `Inf` and evaluate it directly? e.g. `exp(-Inf*x) = 0`

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

By April 22nd, 2025: We can now compute the integral and derivative of all hyperbolic trigonometry functions`sinh(ax+b), cosh(ax+b), tanh(ax+b), coth(ax+b), sech(ax+b), csch(ax+b)`. We can also compute the numerical result of all hyperbolic trignometry functions, e.g. `sech(2)`.

We are taking note of the implementation for `cot, sec, csc, coth, sech, csch` they are all need to be divided to 2 cases for `integer` and `double` to be able to return numerical computation, e.g. `cot(2)`, because in `src/functions.cpp` for the implementation part `return Number<double>(cot(CastPtr<const Number<double> >(s)->n));` will return `segmentation fault`.

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/18.png" width="83%">

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/19.png" width="83%">

By April 24th, 2025: We can now compute with integration by parts method for `x^{b} exp(a*x), x^{b} exp(x)`.

We modify `src/integration.cpp` and `src/function.cpp` to be able to obtain the correct result for the integral problems above.
<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/20.png" width="83%">
<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/22.png" width="83%">

Previously in SymbolicC++, the integration by parts for the form of `x^{b} exp(a*x), x^{b} exp(x)` have incorrect computation.

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/25.png" width="83%">

Based on SymPy that is called from JULIA for the form of `x^{b} exp(a*x), x^{b} exp(x)` we obtain
<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/21.png" width="51%">
<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/26.png" width="51%">

Thus, after a bit of tinkering the code, we are able to fix the integration by parts for the form of `x^{b} exp(a*x), x^{b} exp(x)`
<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/23.png" width="83%">
<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/24.png" width="83%">

The test code can be located in the folder `Test/Integration by Parts/main.cpp`
