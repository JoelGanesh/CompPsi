#include "Types.h"

namespace Types
{
	Complex::Complex(float_dec_T real, float_dec_T imag) : re(real), im(imag)
	{ };

	Complex& Complex::operator*=(Complex z)
	{
		float_dec_T temp = re * z.re - im * z.im;
		float_dec_T temp2 = re * z.im + im * z.re;

		re = temp;
		im = temp2;
		return *this;
	}

	Complex Complex::operator*(Complex z)
	{
		return Complex(re * z.re - im * z.im, re * z.im + im * z.re);
	}

	Complex Complex::operator+(Complex z)
	{
		return Complex(re + z.re, im + z.im);
	}

	Complex Complex::operator-(Complex z)
	{
		return Complex(re - z.re, im - z.im);
	}
}