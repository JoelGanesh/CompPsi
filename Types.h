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

typedef int     Exponent;
typedef uint64_t  Prime;
typedef uint64_t  PrimePower;
typedef int  PrimeIndex;
typedef std::vector<int> Tuple;
typedef std::vector<Tuple> Tuples;

// Boost library shortcut typenames.
using interval = boost::icl::interval<int>;
using interval_set = boost::icl::interval_set<int>;
using float_dec_100 = boost::multiprecision::cpp_dec_float_100;
using complex_128 = boost::multiprecision::complex128;

//typedef boost::multiprecision::uint128_t uint128_t;

namespace Types
{
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