/*
    SymbolicC++ : An object oriented computer algebra system written in C++

    Copyright (C) 2008 Yorick Hardy and Willi-Hans Steeb

    This library is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/


// dual.h

#include <iostream>
#include <cmath>
#include <complex>
using namespace std;
using Complex = std::complex<double>;

#ifndef STRUCT_DUAL
#define STRUCT_DUAL

// Defining struct  in .h file (header file) will make it work when we call it from .cpp file from anywhere that use SymIntegration library
// Structure to store dual number structure for Automatic Differentiation with respect to 'r'
struct Dual {
	complex<double> val; // Function value: a_n(r)
	complex<double> der; // Derivative value: a_n'(r)

	Dual(complex<double> v = 0.0, complex<double> d = 0.0) : val(v), der(d) {}

	Dual operator+(const Dual& o) const { return Dual(val + o.val, der + o.der); }
	Dual operator-(const Dual& o) const { return Dual(val - o.val, der - o.der); }
	Dual operator*(const Dual& o) const { return Dual(val * o.val, val * o.der + der * o.val); }
	Dual operator/(const Dual& o) const { return Dual(val / o.val, (der * o.val - val * o.der) / (o.val * o.val)); }

	Dual operator+(Complex c) const { return Dual(val + c, der); }
	Dual operator-(Complex c) const { return Dual(val - c, der); }
	Dual operator*(Complex c) const { return Dual(val * c, der * c); }
 	Dual operator/(Complex c) const { return Dual(val / c, der / c); }
}; 
	inline Dual operator+(Complex c, const Dual& d) { return Dual(c + d.val, d.der); } // must have either zero or one argument
	inline Dual operator-(Complex c, const Dual& d) { return Dual(c - d.val, -d.der); }
	inline Dual operator*(Complex c, const Dual& d) { return Dual(c * d.val, c * d.der); }
	inline Dual operator/(Complex c, const Dual& d) { return Dual(c / d.val, -c * d.der / (d.val * d.val)); }

#endif

#ifndef STRUCT_DUAL1
#define STRUCT_DUAL1
// Dual1 number structure for automatic differentiation over complex numbers
struct Dual1 {
    complex<double> val;
    complex<double> der;

    Dual1(complex<double> v = 0.0, complex<double> d = 0.0) : val(v), der(d) {}

    friend Dual1 operator+(const Dual1& a, const Dual1& b) { return Dual1(a.val + b.val, a.der + b.der); }
    friend Dual1 operator-(const Dual1& a, const Dual1& b) { return Dual1(a.val - b.val, a.der - b.der); }
    friend Dual1 operator*(const Dual1& a, const Dual1& b) { return Dual1(a.val * b.val, a.val * b.der + a.der * b.val); }
    friend Dual1 operator/(const Dual1& a, const Dual1& b) { return Dual1(a.val / b.val, (a.der * b.val - a.val * b.der) / (b.val * b.val));  }

    // Overloads for mixing with complex<double> / double
    friend Dual1 operator+(const Dual1& a, complex<double> b) { return Dual1(a.val + b, a.der); }
    friend Dual1 operator+(complex<double> a, const Dual1& b) { return Dual1(a + b.val, b.der); }
    friend Dual1 operator-(const Dual1& a, complex<double> b) { return Dual1(a.val - b, a.der); }
    friend Dual1 operator-(complex<double> a, const Dual1& b) { return Dual1(a - b.val, -b.der); }
    friend Dual1 operator*(const Dual1& a, complex<double> b) { return Dual1(a.val * b, a.der * b); }
    friend Dual1 operator*(complex<double> a, const Dual1& b) { return Dual1(a * b.val, a * b.der); }
    friend Dual1 operator/(const Dual1& a, complex<double> b) { return Dual1(a.val / b, a.der / b); }
    friend Dual1 operator/(complex<double> a, const Dual1& b) { return Dual1(a / b.val, -a * b.der / (b.val * b.val)); }
};

#endif

#ifndef STRUCT_DUAL2
#define STRUCT_DUAL2
// Second-order Dual number struct to handle automatic differentiation up to O(eps^2)
// Expressed as: a + b*eps + c*eps^2
struct Dual2 {
	complex<double> a, b, c;

	Dual2(complex<double> a = 0.0, complex<double> b = 0.0, complex<double> c = 0.0) : a(a), b(b), c(c) {}

	Dual2 operator+(const Dual2& o) const { return Dual2(a + o.a, b + o.b, c + o.c); }
	Dual2 operator-(const Dual2& o) const { return Dual2(a - o.a, b - o.b, c - o.c); }
    
	Dual2 operator*(const Dual2& o) const { return Dual2(a * o.a, a * o.b + b * o.a, a * o.c + b * o.b + c * o.a);}

	Dual2 operator/(const Dual2& o) const 
	{
		if (std::abs(o.a) > 1e-12) {
			Complex a_out = a / o.a;
			Complex b_out = (b - a_out * o.b) / o.a;
			Complex c_out = (c - a_out * o.c - b_out * o.b) / o.a;
		return Dual2(a_out, b_out, c_out); 
		} 
		else 
		{
			// L'Hopital's rule for 0/0 forms when denominator's value is zero
			// when roots differ by an integer(for Frobenius method) calculating a_{N} (r2) leads to a division by zero
			// due to the indicial function hitting a root at r1
			Complex a_out = b / o.b;
			Complex b_out = (c - a_out * o.c) / o.b;
			return Dual2(a_out, b_out, 0.0); 
		}
	}

	Dual2 operator+(Complex o) const { return Dual2(a + o, b, c); }
	Dual2 operator-(Complex o) const { return Dual2(a - o, b, c); }
	Dual2 operator*(Complex o) const { return Dual2(a * o, b * o, c * o); }
	Dual2 operator/(Complex o) const { return Dual2(a / o, b / o, c / o); }
};

inline Dual2 operator+(Complex g, const Dual2& d) { return d + g; }
inline Dual2 operator-(Complex g, const Dual2& d) { return Dual2(g - d.a, -d.b, -d.c); }
inline Dual2 operator*(Complex g, const Dual2& d) { return d * g; }


#endif

// --- MULTIVARIABLE COMPLEX DUAL NUMBER CLASS FOR AUTOMATIC DIFFERENTIATION ---
#ifndef DUAL_CLASS
#define DUAL_CLASS
class DualClass {
public:
	complex<double> val;           // Function value
	vector<complex<double>> der;   // Gradient / Partial derivatives vector

	DualClass() : val(0.0) {}
	DualClass(complex<double> v, size_t size = 0) : val(v), der(size, 0.0) {}
	DualClass(complex<double> v, const vector<complex<double>>& d) : val(v), der(d) {}
	DualClass(complex<double> v, size_t idx, size_t size) : val(v), der(size, 0.0) 
	{	
		if (idx < size) der[idx] = 1.0;
	}

	// Addition
	DualClass operator+(const DualClass& other) const 
	{
		size_t n = max(der.size(), other.der.size());
 		DualClass res(val + other.val, n);
		for (size_t i = 0; i < n; ++i) 
		{
			complex<double> d1 = (i < der.size()) ? der[i] : 0.0;
			complex<double> d2 = (i < other.der.size()) ? other.der[i] : 0.0;
			res.der[i] = d1 + d2;
		}
		return res;
	}
	DualClass operator+(complex<double> c) const 
	{
		DualClass res = *this; res.val += c; return res;
	}
	friend DualClass operator+(complex<double> c, const DualClass& d) { return d + c; }

	// Subtraction
	DualClass operator-(const DualClass& other) const 
	{
		size_t n = max(der.size(), other.der.size());
		DualClass res(val - other.val, n);
		for (size_t i = 0; i < n; ++i) 
		{
			complex<double> d1 = (i < der.size()) ? der[i] : 0.0;
			complex<double> d2 = (i < other.der.size()) ? other.der[i] : 0.0;
			res.der[i] = d1 - d2;
		}
		return res;
	}
	DualClass operator-(complex<double> c) const 
	{
		DualClass res = *this; res.val -= c; return res;
	}
	friend DualClass operator-(complex<double> c, const DualClass& d) 
	{
		size_t n = d.der.size();
		DualClass res(c - d.val, n);
		for (size_t i = 0; i < n; ++i) res.der[i] = -d.der[i];
		{
			return res;
		}
	}

	// Multiplication
	DualClass operator*(const DualClass& other) const 
	{
		size_t n = max(der.size(), other.der.size());
		DualClass res(val * other.val, n);
		for (size_t i = 0; i < n; ++i) 
		{
			complex<double> d1 = (i < der.size()) ? der[i] : 0.0;
			complex<double> d2 = (i < other.der.size()) ? other.der[i] : 0.0;
			res.der[i] = val * d2 + d1 * other.val;
		}
		return res;
	}
	DualClass operator*(complex<double> c) const 
	{
		DualClass res(val * c, der.size());
		for (size_t i = 0; i < der.size(); ++i) 
		{
			res.der[i] = der[i] * c;
		}
		return res;
	}
	friend DualClass operator*(complex<double> c, const DualClass& d) { return d * c; }

	// Division
	DualClass operator/(const DualClass& other) const 
	{
		size_t n = max(der.size(), other.der.size());
		DualClass res(val / other.val, n);
		for (size_t i = 0; i < n; ++i) 
		{
			complex<double> d1 = (i < der.size()) ? der[i] : 0.0;
			complex<double> d2 = (i < other.der.size()) ? other.der[i] : 0.0;
			res.der[i] = (d1 * other.val - val * d2) / (other.val * other.val);
		}
		return res;
	}

	// Unary Negation
	DualClass operator-() const 
	{
	        DualClass res(-val, der.size());
		for (size_t i = 0; i < der.size(); ++i) 
		{
			res.der[i] = -der[i];
		}
		return res;
	}
};
// IF I COMMENTED FROM 245 TO 275 THEN COMPUTING NUMERICAL SOLUTIONS FOR NONLINEAR SYSTEMS FOR QUASINEWON STEEPEST DESCENT NEWTONS METHOD WONT WORK
// WHAT SHOULD I DO?
// add inline
// --- ELEMENTARY TRANSCENDENTAL FUNCTIONS FOR DUAL NUMBERS ---
inline DualClass cos(const DualClass& d) 
{
	DualClass res(cos(d.val), d.der.size());
	complex<double> s = -sin(d.val);
	for (size_t i = 0; i < d.der.size(); ++i) 
	{
		res.der[i] = d.der[i] * s;
	}
	return res;
}

inline DualClass sin(const DualClass& d) {
	DualClass res(sin(d.val), d.der.size());
	complex<double> c = cos(d.val);
	for (size_t i = 0; i < d.der.size(); ++i) 
	{
		res.der[i] = d.der[i] * c;
	}
	return res;
}

inline DualClass exp(const DualClass& d) 
{
	DualClass res(exp(d.val), d.der.size());
	for (size_t i = 0; i < d.der.size(); ++i) 
	{
		res.der[i] = d.der[i] * res.val;
	}
	return res;
}

inline DualClass pow(const DualClass& a, int p) 
{
	if (p == 0) 
	{
		return DualClass(1.0, vector<complex<double>>(a.der.size(), 0.0));
	}
	complex<double> val = pow(a.val, p);
	complex<double> der_mult = double(p) * pow(a.val, p - 1);
	vector<complex<double>> d(a.der.size());
	for (size_t i = 0; i < a.der.size(); ++i) 
	{
		d[i] = a.der[i] * der_mult;
	}
	return DualClass(val, d);
}

#endif

#ifndef DUAL_CLASS_TEMPLATE
#define DUAL_CLASS_TEMPLATE
// Supports complex numbers and multi-variable automatic differentiation
template<typename T>
class DualClassTemplate {
public:
	T val;
	vector<T> der;

	DualClassTemplate() : val(0) {}
	DualClassTemplate(T v) : val(v) {}
	DualClassTemplate(T v, const std::vector<T>& d) : val(v), der(d) {}

	// Overloaded operators
	DualClassTemplate operator+(const DualClassTemplate& other) const 
	{
		size_t n = std::max(der.size(), other.der.size());
		vector<T> d(n, T(0));
		for (size_t i = 0; i < n; ++i) 
		{
			if (i < der.size()) 
			{
				d[i] += der[i];
			}
			if (i < other.der.size()) 
			{
				d[i] += other.der[i];
			}
		}
		return DualClassTemplate(val + other.val, d);
	}

	DualClassTemplate operator+(const T& scalar) const { return DualClassTemplate(val + scalar, der); }
	friend DualClassTemplate operator+(const T& scalar, const DualClassTemplate& d) { return d + scalar; }

	DualClassTemplate operator-(const DualClassTemplate& other) const 
	{
		size_t n = std::max(der.size(), other.der.size());
		vector<T> d(n, T(0));
		for (size_t i = 0; i < n; ++i) 
		{
			if (i < der.size()) 
			{
				d[i] += der[i];
			}
			if (i < other.der.size()) 
			{
				d[i] -= other.der[i];
			}
		}
		return DualClassTemplate(val - other.val, d);
	}

	DualClassTemplate operator-(const T& scalar) const { return DualClassTemplate(val - scalar, der); }
	friend DualClassTemplate operator-(const T& scalar, const DualClassTemplate& d) 
	{
		vector<T> neg_der(d.der.size());
		for (size_t i = 0; i < d.der.size(); ++i) 
		{
			neg_der[i] = -d.der[i];
		}
		return DualClassTemplate(scalar - d.val, neg_der);
	}

	DualClassTemplate operator*(const DualClassTemplate& other) const 
	{
		size_t n = std::max(der.size(), other.der.size());
		vector<T> d(n, T(0));
		for (size_t i = 0; i < n; ++i) 
		{
			if (i < der.size()) 
			{
				d[i] += der[i] * other.val;
			}
			if (i < other.der.size()) 
			{
				d[i] += val * other.der[i];
			}
		}
		return DualClassTemplate(val * other.val, d);
	}

	DualClassTemplate operator*(const T& scalar) const 
	{
		vector<T> d(der.size());
		for (size_t i = 0; i < der.size(); ++i)
		{ 
			d[i] = der[i] * scalar;
		}
		return DualClassTemplate(val * scalar, d);
	
}
	friend DualClassTemplate operator*(const T& scalar, const DualClassTemplate& d) { return d * scalar; }

	DualClassTemplate operator/(const T& scalar) const 
	{
		vector<T> d(der.size());
		for (size_t i = 0; i < der.size(); ++i) 
		{
			d[i] = der[i] / scalar;
		}
		return DualClassTemplate(val / scalar, d);
	}

	// Elementary functions for DualClassTemplate
	friend DualClassTemplate cos(const DualClassTemplate& d) 
	{
		vector<T> der(d.der.size());
		T s = -std::sin(d.val);
		for (size_t i = 0; i < d.der.size(); ++i) 
		{
			der[i] = d.der[i] * s;
		}
		return DualClassTemplate(std::cos(d.val), der);
	}

	friend DualClassTemplate sin(const DualClassTemplate& d) 
	{
		vector<T> der(d.der.size());
		T c = std::cos(d.val);
		for (size_t i = 0; i < d.der.size(); ++i) 
		{
			der[i] = d.der[i] * c;
		}
		return DualClassTemplate(std::sin(d.val), der);
	}

	friend DualClassTemplate exp(const DualClassTemplate& d) 
	{
		T ev = std::exp(d.val);
		vector<T> der(d.der.size());
		for (size_t i = 0; i < d.der.size(); ++i) 
		{
			der[i] = d.der[i] * ev;
		}
		return DualClassTemplate(ev, der);
	}
};

#endif

/*


#ifndef DUAL_H
#define DUAL_H
class Dual {
private:
    double val;     // Value of the function
    double der;     // Derivative value (infinitesimal part)

public:
    // Constructors
    Dual(double v = 0.0, double d = 0.0) : val(v), der(d) {}

    // Getters and Setters
    double value() const { return val; }
    double derivative() const { return der; }
    void setValue(double v) { val = v; }
    void setDerivative(double d) { der = d; }

    // Addition
    Dual operator+(const Dual& other) const {
        return Dual(val + other.val, der + other.der);
    }
    Dual operator+(double scalar) const {
        return Dual(val + scalar, der);
    }
    friend Dual operator+(double scalar, const Dual& d) {
        return Dual(scalar + d.val, d.der);
    }

    // Subtraction
    Dual operator-(const Dual& other) const {
        return Dual(val - other.val, der - other.der);
    }
    Dual operator-(double scalar) const {
        return Dual(val - scalar, der);
    }
    friend Dual operator-(double scalar, const Dual& d) {
        return Dual(scalar - d.val, -d.der);
    }
    Dual operator-() const {
        return Dual(-val, -der);
    }

    // Multiplication (Product Rule: d(u*v) = u'*v + u*v')
    Dual operator*(const Dual& other) const {
        return Dual(val * other.val, der * other.val + val * other.der);
    }
    Dual operator*(double scalar) const {
        return Dual(val * scalar, der * scalar);
    }
    friend Dual operator*(double scalar, const Dual& d) {
        return Dual(scalar * d.val, scalar * d.der);
    }

    // Division (Quotient Rule: d(u/v) = (u'*v - u*v') / v^2)
    Dual operator/(const Dual& other) const {
        return Dual(val / other.val, (der * other.val - val * other.der) / (other.val * other.val));
    }
    Dual operator/(double scalar) const {
        return Dual(val / scalar, der / scalar);
    }
    friend Dual operator/(double scalar, const Dual& d) {
        return Dual(scalar / d.val, (-scalar * d.der) / (d.val * d.val));
    }

    // Compound Assignment Operators
    Dual& operator+=(const Dual& other) {
        val += other.val; der += other.der; return *this;
    }
    Dual& operator-=(const Dual& other) {
        val -= other.val; der -= other.der; return *this;
    }
    Dual& operator*=(const Dual& other) {
        der = der * other.val + val * other.der; val *= other.val; return *this;
    }
    Dual& operator/=(const Dual& other) {
        der = (der * other.val - val * other.der) / (other.val * other.val); val /= other.val; return *this;
    }

    // Elementary Math Functions (Chain Rule implementation)
    friend Dual exp(const Dual& d) {
        double e = std::exp(d.val);
        return Dual(e, e * d.der);
    }

    friend Dual log(const Dual& d) {
        return Dual(std::log(d.val), d.der / d.val);
    }

    friend Dual sin(const Dual& d) {
        return Dual(std::sin(d.val), std::cos(d.val) * d.der);
    }

    friend Dual cos(const Dual& d) {
        return Dual(std::cos(d.val), -std::sin(d.val) * d.der);
    }

    friend Dual pow(const Dual& d, double p) {
        return Dual(std::pow(d.val, p), p * std::pow(d.val, p - 1) * d.der);
    }

    // Stream Output
    friend std::ostream& operator<<(std::ostream& os, const Dual& d) {
        os << d.val << " + " << d.der << "ε";
        return os;
    }
};

#endif
*/