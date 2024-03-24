// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#ifndef TYPES_H
#define TYPES_H

#include <vector>

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>

#define PRECISION 33
#define DOUBLE_ERROR_THRESHOLD 1.0E-8

typedef char Exponent;
typedef int64_t Prime;
typedef std::vector<int> Tuple;
typedef std::vector<Tuple> Tuples;

typedef boost::multiprecision::int128_t int128_t;
typedef boost::multiprecision::cpp_bin_float_quad float128_t;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<PRECISION>> float_dec_T;

namespace Types
{
	struct Complex
	{
		float_dec_T re;
		float_dec_T im;

		Complex(float_dec_T real = 0, float_dec_T imag = 0);

		Complex& operator*=(Complex z);
		Complex operator*(Complex z);
		Complex operator+(Complex z);
		Complex operator-(Complex z);
	};

	struct Fraction
	{
		int128_t num;
		int128_t denom;

		Fraction(int128_t num, int128_t denom = 1) : num(num), denom(denom)
		{ }

		// Return a numerical representation of the fraction.
		float128_t numerical() const;

		// Returns the floor as integer.
		int128_t Floor() const;

		// Rounds the fraction to the nearest integer.
		int128_t Round() const;

		// Returns fractional part as a fraction.
		Fraction FractionalPart() const;

		// Returns the sign of the fraction: 1 if positive, -1 if negative, 0 if 0.
		int Sign() const;

		// Checks if the fraction is an integer.
		bool IsIntegral() const;

		// Checks if the fraction is negative.
		bool IsNegative() const;

		// Checks if the fraction is zero.
		bool IsZero() const;

		// Inverts the fraction.
		void Invert();

		// Negates the fraction.
		void Negate();

		Fraction operator*(Fraction f);
		Fraction operator+(Fraction f);
		Fraction operator-(Fraction f);
		Fraction operator/(Fraction f);
	};

	struct Interval
	{
		int64_t start;
		int64_t end;

		// Returns an interval [start, end].
		// If start > end, should be read as an empty interval.
		Interval(int64_t start = 1, int64_t end = 0);

		// Returns the integer interval enclosed 
		// by the roots of the quadratic ax^2 + bx + c.
		// If a < 0, potential roots are included;
		// if a > 0, we instead disregard them,
		// according to the implementation by Helfgott & Thompson.
		// Requires that a != 0.
		Interval(int64_t a, int64_t b, int64_t c);

		// Shifts 'this' by 'a' units (x -> x + a).
		void Shift(int64_t a);

		// Intersects 'this' with I.
		void Intersect(Interval& I);
	};

	// Placeholder for logarithms of positive integers.
	struct Log
	{
		int64_t n;

		Log(int64_t n) : n(n)
		{ }

		// Returns approximation of 'this' as float_dec_T.
		float_dec_T numerical() const
		{
			if (n > 1)
			{
				return boost::multiprecision::log(float_dec_T(n));
			}
			return 0;
		}

		operator float_dec_T()
		{
			return numerical();
		}
	};

	// Structure to store a pair of a prime and an exponent.
	// Also provides implicit conversion to string type.
	struct PrimeFactor
	{
		Prime prime;
		Exponent exponent;

		PrimeFactor(Prime p, Exponent exp) : prime(p), exponent(exp)
		{ }

		// Returns integer presentation of the prime factor.
		operator int64_t()
		{
			return (int64_t)std::pow(prime, exponent);
		}
	};

	// Structure to store the primefactorization of an integer.
	// Primefactors can be added individually.
	struct Factorization
	{
		private:
		// List of prime factors.
		std::vector<PrimeFactor> primeFactors_;

		// Product of the prime factors.
		int64_t n_;

		public:
		// Constructor
		Factorization() : n_(1)
		{ }

		// Returns n_.
		const int64_t n() const
		{
			return n_;
		}

		// Make it possible to read the prime factors
		// without being able to change them out-of-scope.
		const std::vector<PrimeFactor> primeFactors() const
		{
			return primeFactors_;
		}

		void AddFactor(Prime p, Exponent k, bool update_n = false)
		{
			primeFactors_.push_back(PrimeFactor(p, k));
			if (update_n)
			{
				n_ *= (int64_t)std::pow(p, k);
			}
		}
	};
}
#endif // TYPES_H