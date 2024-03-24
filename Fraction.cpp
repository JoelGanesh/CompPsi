#include "Types.h"

namespace Types
{
	float128_t Fraction::numerical() const
	{
		if (denom != 0)
		{
			return float128_t(num) / float128_t(denom);
		}
		return 1;
	};

	int128_t Fraction::Floor() const
	{
		if (num < 0 && num % denom != 0)
		{
			return (num / denom) - 1;
		}
		return num / denom;
	}

	bool Fraction::IsIntegral() const
	{
		return num % denom == 0;
	}

	bool Fraction::IsNegative() const
	{
		if (denom > 0)
		{
			return num < 0;
		}
		else return num > 0;
	}

	bool Fraction::IsZero() const
	{
		return num == 0;
	}

	int Fraction::Sign() const
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

	int128_t Fraction::Round() const
	{
		Fraction q(2 * num + denom, 2 * denom);
		return q.Floor();
	}

	Fraction Fraction::operator*(Fraction f)
	{
		return Fraction(num * f.num, denom * f.denom);
	};

	Fraction Fraction::operator+(Fraction f)
	{
		return Fraction(num * f.denom + f.num * denom,
						denom * f.denom);
	};

	Fraction Fraction::operator-(Fraction f)
	{
		return Fraction(num * f.denom - f.num * denom,
						denom * f.denom);
	};

	Fraction Fraction::operator/(Fraction f)
	{
		return Fraction(num * f.denom, denom * f.num);
	};

	void Fraction::Invert()
	{
		int128_t num_temp = num;
		num = denom;
		denom = num_temp;

	};

	void Fraction::Negate()
	{
		num *= -1;
	}

	// Fractional part; i.e., num/denom - floor(num/denom).
	Fraction Fraction::FractionalPart() const
	{
		int128_t temp_num = num % denom;
		if (temp_num < 0)
		{
			temp_num += denom;
		}
		return Fraction(temp_num, denom);
	};
}