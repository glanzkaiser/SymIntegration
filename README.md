# SymIntegration

SymIntegration is a C++ library that is branching out from SymbolicC++3.35.
The main idea is to improve its symbolic integration codes. Our main focus is to make it able to compute:

1. All kinds of standard integral form (trigonometry, inverse trigonometry, polynomial, transcendental, hyperbolic)
2. The sum, product and divide combination of the standard functions
3. To be able to compute improper integrals with cases (e.g. computing mean and variance for exponential distribution)

But as time goes by we add more than integral from calculus, we add differential equations solver, statistics, and more in the future.

# Manual / Documentation

<a href="https://github.com/glanzkaiser/SymIntegration/blob/main/SymIntegration.pdf">SymIntegration Manual</a>


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
2. The integral of $`\frac{1}{(ax+b)}`$, $`\frac{1}{ax^{2}+bx + c}`$, $`\frac{a_{1}x^2 + b_{1} x + c_{1}}{ax^{2}+bx + c}`$, $`\frac{a_{1}x + b_{1}}{ax^{2}+bx + c}`$
3. The integral and derivative of `asin(ax+b), acos(ax+b), atan(ax+b), acot(ax+b), asec(ax+b), acsc(ax+b)` with `ax + b` is a polynomial of order 1. Only `asec(ax+b), acsc(ax+b)` integration that have no analytic solution yet since it has cases output.
4. The integral and derivative of all hyperbolic trigonometry functions`sinh(ax+b), cosh(ax+b), tanh(ax+b), coth(ax+b), sech(ax+b), csch(ax+b)`.
5. The integral with the form of $`x^{b} e^{ax}, x^{b} e^{x}`$ with integration by parts method.
6. The integral of $`\sin(ax) \cos(bx), \cos(ax) \cos (bx), \sin(ax) \sin (bx)`$.
7. The integral of $`\int \sin^{n} (x) \ dx , \int \cos^{n} (x) \ dx, \int \tan^{n} (x) \ dx, \int \sec^{n} (x) \ dx, \int \csc^{n} (x) \ dx `$, but the speed to show the whole symbolic integral is slow compared to SymPy, needs to be fixed.
8. The solution ($`y(t)`$) of first order linear ordinary differential equation of type $`ay' + ty = b`$, separable equation of type $`\frac{dy}{dx} = \frac{f(x)}{g(y)}`$, and separable equation with homogeneous ratio equations of type $`\frac{dy}{dx} = \frac{f(x,y)}{g(x,y)}`$ along with the initial value problem solution.
9. The computation of divergence, gradient, and curl.
10. The probability mass function, mean, variance, moment generating function of discrete distributions.
11. The probability density function, cumulative distribution function, mean, variance, moment generating function of continuous distributions.
12. The regression line equation along with Pearson's correlation ($`r`$) with input of vector of 2 columns
13. Bond pricing and amortization table for mortgage loan or other credit related loan.
14. Compute direction fields and plot it with gnuplot.
15. Linear programming solution with Simplex Method.
16. Numerical differentiation and Richardson's Extrapolation.

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

## Create the Shared Library (Manual Way)
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

## Create the Shared Library (with Makefile)

We have provide a Makefile that can be used to make the library faster, because if you only making change at, e.g. `src/integrate.cpp`, then instead of creating the object files for all the `.cpp` files it will only create a new object file for the newly modified `src/integrate.cpp` file.

The downside is, the size of the library is now 3.6 MB, maybe it is because we add `-Wall -Wformat` in the CXX flags.

To create the SymIntegration shared library, open terminal and from the current working directory / this repository main directory:

```
	cd src
	make
```
Then move the newly build `libsymintegration.so` to `/usr/lib`, assuming your libraries default location is there.
 
<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/41.png" width="51%">

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/42.png" width="51%">

# Test Example

Open terminal and from the current working directory / this repository main directory:

```
	cd Examples/Integral
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
|:sunflower:   | Compute definite integral with improper integrals			| Done (write `INFINITY` to substitute `Inf`)
|:sunflower:   | Compute the solution of first order linear ODE of several types	| Done
|:sunflower:   | Compute the probability mass function, mean, variance, moment generating function of geometric distribution, binomial distribution, negative binomial distribution, and Poisson distribution	| Done
|:sunflower:   | Compute the probability density function, cumulative distribution function, mean, variance, moment generating function of uniform distribution, normal distribution, gamma, exponential, beta, cauchy, laplace, logistic, chi-squared, students't, and F distribution	| Done
|:sunflower:   | Compute the root / solution of a function of one variable with bisection method and Newton-Raphson method	| Done
|:sunflower:   | Compute the minimization solution with gradient descent and downhill simplex	| Done
|:sunflower:   | Compute the regression line equation along with Pearson's correlation ($`r`$) with input of vector of 2 columns	| Done
|:sunflower:   | Bond pricing computation and amortization table			| Done
|:sunflower:   | Linear Programming: Simplex Method					| Done
|:sunflower:   | Generate Random Number with Mersenne Twister				| Done
|:sunflower:   | Numerical differentiation and Richardson's Extrapolation		| Done

# Milestone

We are using less, minimal amount of code to create the library, if you compare it with the original SymbolicC++ that has `.configure` and `Makefile` that will be spawn after you configure it, we don' t use any of that.

Two simple commands with `g++` can already make a shared symbolic integration computation library, with limitation still.

We have found how to be able to compute $`\int sin(2x) dx`$ correctly.

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/5.png" width="60%">

By modifying the source code `src/functions.cpp` we can compute the integration for sine and cosine correctly.

By April 10th, 2025: The implementation of sine and cosine have been modified.

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/6.png" width="60%">

For comparison: The basic codes that we are using SymbolicC++, cannot compute the integral of sine and cosine correctly because they only write `return Integral(*this,s)`

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/7.png" width="60%">

By April 13th, 2025: The integral of $`\frac{1}{ax+b}`$ can be computed nicely and fraction / decimal power computation e.g. `(-1)^(a/b)` can be computed too.

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/8.png" width="60%">

The code for the decimal power computation is located here (the if statement with 4 conditions):

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/9.png" width="60%">

By April 17th, 2025: We add the implementation for `tan(x)` and `cot(x)` so it can compute the derivative and integral of all basic trigonometry functions (`sin(ax+b), cos(ax+b), tan(ax+b), cot(ax+b), sec(ax+b), csc(ax+b)`) for polynomial of order 1.

The codes that we modified are located in: 

* src/functions.cpp
* src/symintegrationc++.cpp
* include/symintegral/symintegrationc++.h
* include/symintegral/functions.h (to add the Class of `Tan` and `Cot`)

The test code can be located in the folder `Examples/Trigonometry/main.cpp`

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/10.png" width="60%">

Comparing the result with JULIA

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/11.png" width="60%">

The results produce the same output number at the end

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/12.png" width="60%">

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

The test code can be located in the folder `Examples/Inverse Trigonometry/main.cpp`

Comparing the result with SymPy in JULIA:

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/13.png" width="60%">

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/14.png" width="60%">

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/15.png" width="60%">

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/16.png" width="60%">

The analytic solution for integral of `asec(ax+b)` from SymPy in JULIA:
<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/17.png" width="60%">

By April 22nd, 2025: We can now compute the integral and derivative of all hyperbolic trigonometry functions`sinh(ax+b), cosh(ax+b), tanh(ax+b), coth(ax+b), sech(ax+b), csch(ax+b)`. We can also compute the numerical result of all hyperbolic trignometry functions, e.g. `sech(2)`.

We are taking note of the implementation for `cot, sec, csc, coth, sech, csch` they are all need to be divided to 2 cases for `integer` and `double` to be able to return numerical computation, e.g. `cot(2)`, because in `src/functions.cpp` for the implementation part `return Number<double>(cot(CastPtr<const Number<double> >(s)->n));` will return `segmentation fault`.

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/18.png" width="60%">

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/19.png" width="60%">

By April 24th, 2025: We can now compute with integration by parts method for $`x^{b} e^{ax}, x^{b} e^{x}`$.

We modify `src/integration.cpp` and `src/functions.cpp` to be able to obtain the correct result for the integral problems above.

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/20.png" width="60%">

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/22.png" width="60%">

Previously in SymbolicC++, the integration by parts for the form of $`x^{b} e^{ax}, x^{b} e^{x}`$ have incorrect computation.

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/25.png" width="60%">

Based on SymPy that is called from JULIA for the form of $`x^{b} e^{ax}, x^{b} e^{x}`$ we obtain

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/21.png" width="51%">
<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/26.png" width="51%">

Thus, after a bit of tinkering the code, we are able to fix the integration by parts for the form of $`x^{b} e^{ax}, x^{b} e^{x}`$

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/23.png" width="60%">
<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/24.png" width="60%">

The test code can be located in the folder `Examples/Integration by Parts/main.cpp`

By April 25th, 2025: We can now compute the integral of the general form $`\frac{1}{ax^{2} + bx + c}`$.

By modifying the `src/functions.cpp`

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/27.png" width="60%">

et voilà

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/28.png" width="60%">

The test code can be located in the folder `Examples/Integration of Polynomial Order 2 in Denominator/main.cpp`

By April 28th, 2025: We can now compute the integral of the general form $`\cos(ax) \sin (bx), \sin(ax) \sin(bx)`$ and $`\cos(ax) \cos (bx)`$.

By modifying the `src/integrate.cpp`, we can handle the $`\cos(mx) \sin (nx)`$ nicely, but there is a problem for $`\cos(ax) \cos (bx), \sin(ax) \sin(bx)`$ when $`a = b`$ or $`a = -b`$, so we deal with it in `src/functions.cpp` 

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/29.png" width="60%">

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/30.png" width="60%">

The test code can be located in the folder `Examples/Trigonometry Integration Level 2/main.cpp`

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/31.png" width="60%">

By April 29th, 2025: We can now compute the integral of the form $`\cos(x) \cos(x), \sin(x) \sin(x), \tan(x) \tan(x), \cot(x) \cot(x), \sec(x) \sec(x), \csc(x) \csc(x)`$.

By modifying the `src/functions.cpp`

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/33.png" width="60%">

The test code can be located in the folder `Examples/Trigonometry Integration Level 1.5/main.cpp`

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/32.png" width="60%">

By May 1st, 2025: We add a test example to compute Jacobi Polynomials that is using combinatorial and factorial.

The test code can be located in the folder `Examples/Compute Jacobi Polynomials/main.cpp`

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/34.png" width="60%">

We also add a test example to compute integral of the form $`\int \sin^{n} (x) \ dx`$, we can try with any number of `n`.

(this test example is a very manual way to compute  $`\int \sin^{n} (x) \ dx`$, next we will add the implementation in the `src/functions.cpp` so we can use the integrate formula directly)

The test code can be located in the folder `Examples/Compute Hypergeometric Function 2F1/main.cpp`

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/35.png" width="60%">

The formula is based on

```math
\int \sin^{n} (x) \ dx = - \cos(x) \sin^{n+1} (x) (\sin^{2} (x))^{-\frac{n}{2} - \frac{1}{2}} {}_{2}F_{1} \left( \frac{1}{2}, \frac{1-n}{2} ; \frac{3}{2} ; \cos^{2} (x) \right)
```

We obtain the formula from Wolfram Alpha (wolframalpha.com). It is nice knowing that the general indefinite integral for  $`\int \sin^{n} (x) \ dx`$ and other trigonometry like $`\cos, \tan, \cot, \sec, \csc`$ are using Hypergeometric function.

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/sin1.png" width="60%">

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/sin2.png" width="60%">

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/hypergeometric1.png" width="60%">

By May 4th, 2025: We are abandoning the idea to use Hypergeometric function to compute the integral of basic trigonometric function to the power of $`n`$.

We use the Reduction formula instead, it is a recursive method, and we are adding the formula in `src/functions.cpp`

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/36.png" width="60%">

We have done for $`\int \sin^{n} (x) \ dx , \int \cos^{n} (x) \ dx `$ for both the case where $`n`$ is even and odd.

For  $`\int \sec^{n} (x) \ dx , \int \csc^{n} (x) \ dx `$ for the case where $`n`$ is odd.

We will finish the whole in a few days.

By May 9th, 2025: We have finished to add all the formulas to compute $`\int \sin^{n} (x) \ dx , \int \cos^{n} (x) \ dx, \int \tan^{n} (x) \ dx, \int \cot^{n} (x) \ dx, \int \sec^{n} (x) \ dx, \int \csc^{n} (x) \ dx`$

By May 24th, 2025: We have finished to make the C++ codes to compute $`\int \sin^{n} (x) \cos^{m} (x)\ dx `$ for all possible cases of $`n`$ and $`m`$.

The test code can be located in the folder `Examples/Test SymIntegration Trigonometry Integration Level 4 Sin^n Cos^m Even Manual n>m/main.cpp`

The test code can be located in the folder `Examples/Test SymIntegration Trigonometry Integration Level 4 Sin^n Cos^m Even Manual n<m/main.cpp`

The test code can be located in the folder `Examples/Test SymIntegration Trigonometry Integration Level 4 Sin^n Cos^m for n = m/main.cpp`

The test code can be located in the folder `Examples/Test SymIntegration Trigonometry Integration Level 4 Sin^n Cos^m Odd Manual/main.cpp`

We haven't put it in `src/functions` yet.

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/37.png" width="60%">

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/38.png" width="60%">

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/39.png" width="60%">

By June 23rd, 2025: SymIntegration is able to compute the solution ($`y(t)`$) of first order linear ordinary differential equation of type $`ay' + ty = b`$ with method of integration and separable equation.

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/44.png" width="60%">

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/43.png" width="60%">

By August 18th, 2025: SymIntegration is able to compute the div, grad, curl operators.

The test code can be located in the folder `Examples/Test SymIntegration DivGradCurl/main.cpp`

By August 23rd, 2025: SymIntegration is able to compute the pmf, mean, variance, and mgf of some discrete probability distributions: geometric, binomial, negative binomial, Poisson.

The test code can be located in the folder `Examples/Test SymIntegration Discrete Distributions/main.cpp`

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/45.png" width="60%">

By August 28th, 2025: SymIntegration is able to compute the pdf, cdf, mean, variance, and mgf of some continuous probability distributions: uniform, normal, gamma, exponential, beta.

The test code can be located in the folder `Examples/Test SymIntegration Continuous Distributions/main.cpp`

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/46.png" width="60%">

By September 2nd, 2025: We add some continuous probability distributions (cauchy, laplace, logistic, chi-squared, students't, and F distribution)

The test code can be located in the folder `Examples/Test SymIntegration Continuous Distributions Part 2/main.cpp`

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/47.png" width="60%">

By September 20th, 2025: We add some functions to compute bond price, amortization table, regression line, and Pearson's correlation coefficient.

The test code can be located in the folder `Examples/Test SymIntegration Amortization for Mortgage/main.cpp`

The test code can be located in the folder `Examples/Test SymIntegration Pearson's Product Moment Correlation and Regression Line with Armadillo/main.cpp`

We are also able to compute direction fields from input of symbolic function $`f(t,y)`$, then plot it with gnuplot afterwards.

The test code can be located in the folder `Examples/Test SymIntegration Compute Direction Fields/main.cpp`

<img src="https://github.com/glanzkaiser/SymIntegration/blob/main/images/48.png" width="60%">
