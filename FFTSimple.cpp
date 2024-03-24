// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#include "Fourier.h"

namespace Fourier
{
	FFTSimple::FFTSimple(int size, int n)
	{
		int k = size * n;
		int N = 1;
		while (N < k)
		{
			N *= 2;
		}
		this->size = N;
	}

	void FFTSimple::InitializeContainer(std::vector<Complex>& Fv)
	{
		Fv = std::vector<Complex>(size, Complex(1, 0));
	}

	void FFTSimple::Multiply(std::vector<Complex>& Fv, std::vector<Complex>& Fv_i)
	{
		for (int j = 0; j < size; j++)
		{
			Fv[j] *= Fv_i[j];
		}
	}

	std::vector<Complex> FFTSimple::FFT(std::vector<float_dec_T>& v)
	{
		int n = v.size();
		std::vector<Complex> Fv(n);
		if (n == 1)
		{
			Fv[0].re = v[0];
			return Fv;
		}

		int m = n >> 1;
		std::vector<float_dec_T> v_even(m), v_odd(m);
		for (int k = 0; k < m; k++)
		{
			v_even[k] = v[2 * k];
			v_odd[k] = v[2 * k + 1];
		}

		std::vector<Complex> Fv_even = FFT(v_even), Fv_odd = FFT(v_odd);

		for (int k = 0; k < m; k++)
		{
			float_dec_T theta = (2 * boost::math::constants::pi<float_dec_T>() * k) / n;
			Complex w(boost::multiprecision::cos(theta), boost::multiprecision::sin(theta));
			Fv[k] = Fv_even[k] + w * Fv_odd[k];
			Fv[k + m] = Fv_even[k] - w * Fv_odd[k];
		}

		return Fv;
	}

	std::vector<float_dec_T> FFTSimple::IFFT(std::vector<Complex>& Fv)
	{
		std::vector<Complex> v = IFFT_(Fv);
		std::vector<float_dec_T> w(v.size());
		for (int i = 0; i < w.size(); i++)
		{
			w[i] = v[i].re;
		}
		return w;
	}

	std::vector<Complex> FFTSimple::IFFT_(std::vector<Complex>& Fv)
	{
		int n = Fv.size();
		std::vector<Complex> v(n);
		if (n == 1)
		{
			return Fv;
		}

		int m = n >> 1;
		std::vector<Complex> Fv_even(m), Fv_odd(m);
		for (int k = 0; k < m; k += 1)
		{
			Fv_even[k] = Fv[2 * k];
			Fv_odd[k] = Fv[2 * k + 1];
		}

		std::vector<Complex> v_even = IFFT_(Fv_even), v_odd = IFFT_(Fv_odd);

		for (int k = 0; k < m; k++)
		{
			float_dec_T theta = (2 * boost::math::constants::pi<float_dec_T>() * k) / n;
			Complex w(boost::multiprecision::cos(theta), -boost::multiprecision::sin(theta));
			v[k] = v_even[k] + w * v_odd[k];
			v[k + m] = v_even[k] - w * v_odd[k];
		}

		return v;
	}
}