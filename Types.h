#ifndef TYPES_H
#define TYPES_H

#include <vector>
#include <iostream>
#include <stack>
#include <functional>
#include <string>
#include <algorithm>
#include <numeric>

#include <boost/icl/interval.hpp>
#include <boost/icl/interval_set.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <fftw3.h>

typedef int     Exponent;
typedef int64_t  Prime;
typedef int64_t  PrimePower;
typedef int  PrimeIndex;
typedef std::vector<int> Tuple;
typedef std::vector<Tuple> Tuples;

// Boost library shortcut typenames.
//using interval = boost::icl::interval<int>;
//using interval_set = boost::icl::interval_set<int>;
using float_dec_100 = boost::multiprecision::cpp_dec_float_100;
using complex_128 = boost::multiprecision::complex128;

using namespace boost::multiprecision;
//typedef boost::multiprecision::int128_t int128_t;

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
		int64_t N, M, m0, d0, a, b, a0, a0_inv, q, s;
		double delta;

		Rectangle(int64_t N, int64_t M, int64_t m0, int64_t d0, int64_t a, int64_t b);

		std::tuple<double, int64_t> beta_r0(int64_t m);
	};

	struct Interval
	{
		int64_t start;
		int64_t end;

		// Returns an interval [start, end].
		// If start > end, should be read as an empty interval.
		Interval(int64_t start = 1, int64_t end = 0) : start(start), end(end)
		{
		};

		// Returns the integer interval enclosed 
		// by the roots of the quadratic ax^2 + bx + c.
		// If a < 0, potential roots are included;
		// if a > 0, we instead disregard them,
		// according to the implementation by Helfgott & Thompson.
		// Requires that a != 0.
		Interval(int64_t a, int64_t b, int64_t c) : start(1), end(0)
		{
			int64_t D = b * b - 4 * a * c;
			if (D >= 0)
			{
				int64_t Q = std::sqrt(D);
				if (a < 0)
				{
					start = std::ceil((double)(-b + Q) / (2 * a));
					end = std::floor((double)(-b - Q) / (2 * a));
				}
				else if (a > 0)
				{
					// We have to distinguish the possibilities of D being a square or not;
					// for certain values of a and b this results in a different integer interval.
					if (Q * Q != D)
					{
						start = std::ceil((double)(-b - Q) / (2 * a));
						end = std::floor((double)(-b + Q) / (2 * a));
					}
					else
					{
						start = std::floor((double)(-b - Q) / (2 * a)) + 1;
						end = std::ceil((double)(-b + Q) / (2 * a)) - 1;
					}
				}
			}
		};

		void Shift(int64_t a)
		{
			if (start != LLONG_MIN && start != LLONG_MAX)
			{
				start += a;
			}
			if (end != LLONG_MIN && end != LLONG_MAX)
			{
				end += a;
			}
		}

		// Intersects this interval with I.
		void Intersect(Interval I)
		{
			start = std::max(start, I.start);
			end = std::min(end, I.end);
		};
	};

	// Placeholder for logarithms of positive integers.
	struct Log
	{
		int64_t n;

		// Constructor
		Log(int64_t n);

		// Retrieve approximate value.
		float_dec_100 numerical() const;

		operator float_dec_100()
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

		PrimeFactor(Prime p, Exponent exp);

		operator int64_t();
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
		int64_t n_;

		public:
		// Constructor
		Factorization() : n_(1) {};

		// Returns n_.
		const int64_t n();

		// Make it possible to read the prime factors
		// without being able to change them out-of-scope.
		const std::vector<PrimeFactor> primeFactors();

		void AddFactor(Prime p, Exponent k, bool update_n = false);

		operator std::string();
	};

	struct Fraction
	{
		public:
		int64_t num;
		int64_t denom;

		Fraction(int64_t num, int64_t denom) : num(num), denom(denom)
		{
			//Normalize();
		};

		Fraction(int64_t num) : num(num), denom(1)
		{};

		double numerical() const
		{
			if (denom != 0)
			{
				return (double)num / denom;
			}
			return 1;
		};

		int64_t Floor() const
		{
			if (num < 0 && num % denom != 0)
			{
				return (num / denom) - 1;
			}
			return num / denom;
		}

		bool IsIntegral() const
		{
			return num % denom == 0;
		}

		bool IsNegative() const
		{
			if (denom > 0)
			{
				return num < 0;
			}
			else return num > 0;
		}

		bool IsZero() const
		{
			return num == 0;
		}

		int Sign() const
		{
			if (IsNegative())
			{
				return -1;
			}
			else if (IsZero())
			{
				return 0;
			}
			return 1;
		}

		Fraction operator*(Fraction f)
		{
			return Fraction(num * f.num, denom * f.denom);
		};

		Fraction operator+(Fraction f)
		{
			return Fraction(num * f.denom + f.num * denom,
							denom * f.denom);
		};

		Fraction operator-(Fraction f)
		{
			return Fraction(num * f.denom - f.num * denom, 
							denom * f.denom);
		};

		Fraction operator/(Fraction f)
		{
			return Fraction(num * f.denom, denom * f.num);
		};

		void Invert()
		{
			int64_t num_temp = num;
			num = denom;
			denom = num_temp;

		};

		void Negate()
		{
			num *= -1;
		}

		// Fractional part; i.e., num/denom - floor(num/denom).
		Fraction FractionalPart()
		{
			int64_t temp_num = num % denom;
			if (temp_num < 0)
			{
				temp_num += denom;
			}
			return Fraction(temp_num, denom);
		};

		private:
		// Makes sure that denominator stays positive,
		// also reduces fraction if possible.
		void Normalize()
		{
			if (denom < 0)
			{
				num *= -1;
				denom *= -1;
			}

			int64_t d = std::_Gcd(num, denom);
			num /= d;
			denom /= d;
		}
	};
	/*
	// Structure to store the prime divisors of an integer.
	struct SqFreeFactorization
	{
	protected:
		std::vector<Prime> primes_;

		// Product of the prime factors.
		int64_t n_;

	public:
		// Constructor
		SqFreeFactorization() : n_(1), primes_() {};

		const std::vector<Prime> Primes();
		virtual int64_t n();

		virtual void AddFactor(Prime p);

		virtual operator std::string();
	};


	struct PrimeFactorization : SqFreeFactorization
	{
	private:
		std::vector<Exponent> multiplicities_;

	public:
		// Constructor
		PrimeFactorization() : multiplicities_() {};

		const std::vector<Exponent> Multiplicities();

		void AddFactor(Prime p) override;
		void AddFactor(Prime p, Exponent j, PrimePower q = 0);

		operator std::string() override;
	};*/
}
#endif // TYPES_H