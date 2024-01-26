#ifndef TYPES_H
#define TYPES_H

#include <vector>
#include <iostream>
#include <stack>
#include <functional>
#include <string>
#include <algorithm>

#include <boost/icl/interval.hpp>
#include <boost/icl/interval_set.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <fftw3.h>

typedef int     Exponent;
typedef uint64_t  Prime;
typedef uint64_t  PrimePower;
typedef int  PrimeIndex;
typedef std::vector<int> Tuple;
typedef std::vector<Tuple> Tuples;

// Boost library shortcut typenames.
using interval = boost::icl::interval<int>;
using interval_set = boost::icl::interval_set<int>;
using float_dec_100 = boost::multiprecision::cpp_dec_float_50;
using complex_128 = boost::multiprecision::complex128;

using namespace boost::multiprecision;
//typedef boost::multiprecision::uint128_t uint128_t;

namespace Types
{
	struct Complex
	{
		float_dec_100 re;
		float_dec_100 im;

		Complex(float_dec_100 real = float_dec_100(0), float_dec_100 imag = float_dec_100(0))
		{
			re = real;
			im = imag;
		}

		Complex& operator*=(Complex z)
		{
			//std::cout << std::string(*this) << " " << std::string(z) << std::endl;
			float_dec_100 temp = re * z.re - im * z.im;
			float_dec_100 temp2 = re * z.im + im * z.re;

			re = temp;
			im = temp2;

			//std::cout << std::string(*this) << std::endl;

			return *this;
		}

		Complex operator*(Complex z)
		{
			return Complex(re * z.re - im * z.im, re * z.im + im * z.re);
		}

		Complex operator+(Complex z)
		{
			return Complex(re + z.re, im + z.im);
		}

		Complex operator-(Complex z)
		{
			return Complex(re - z.re, im - z.im);
		}

		operator std::string()
		{
			return re.str(6) + " + " + im.str(6) + "i";
		}
	};

	struct Rectangle
	{
		uint64_t N, M, m0, d0, a, b, a0, a0_inv, q, s;
		double delta;

		Rectangle(uint64_t N, uint64_t M, uint64_t m0, uint64_t d0, uint64_t a, uint64_t b);

		std::tuple<double, uint64_t> beta_r0(uint64_t m);
	};

	// Placeholder for logarithms of integers.
	struct Log
	{
		uint64_t n;

		// Constructor
		Log(uint64_t n);

		// Retrieve approximate value.
		float_dec_100 numerical() const;
	};

	// Structure to store a pair of a prime and an exponent.
	// Also provides implicit conversion to string type.
	struct PrimeFactor
	{
		Prime prime;
		Exponent exponent;

		PrimeFactor(Prime p, Exponent exp);

		operator std::string();
	};

	// Structure to store the primefactorization of an integer.
	// Primefactors can be added individually.
	struct Factorization
	{
	private:
		std::vector<PrimeFactor> primeFactors_;

		// Product of the prime factors which have been added
		// with corresponding parameter 'update_n' equal to true.
		uint64_t n_;

	public:
		// Constructor
		Factorization() : n_(1) {};

		// Returns n_.
		const uint64_t n();

		// Make it possible to read the prime factors
		// without being able to change them out-of-scope.
		const std::vector<PrimeFactor> primeFactors();

		void AddFactor(Prime p, Exponent k, bool update_n = false);

		operator std::string();
	};
}
#endif // TYPES_H