/*
#include "Fourier.h"
#include "Utility.h"

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <fftw3.h>

namespace Fourier
{
	FFTW::FFTW(int size, int n)
	{
		// To compute the non-modular convolution on n vectors,
		// the vectors need to be padded with zeros until
		// each vector has size equal to n times the original size.
		this->size = size * n;

		// Initialize arrays for input and output by allocating memory.
		in = (double*)fftw_malloc(sizeof(double) * this->size);
		out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * this->size);

		// Set up plans for both the forward and backward Fourier transform
		// to be executed multiple times.
		forward_plan = fftw_plan_dft_r2c_1d(this->size, in, out, FFTW_ESTIMATE);
		backward_plan = fftw_plan_dft_c2r_1d(this->size, out, in, FFTW_ESTIMATE);
	}

	const fftw_complex* FFTW::FFT(std::vector<double> v)
	{
		in = &v[0]; // We replace 'in' by the data in 'v'.
		fftw_execute(forward_plan);
		return out;
	}

	const double* FFTW::IFFT(std::vector<fftw_complex> Fv)
	{	
		out = &Fv[0]; // We replace 'out' by the data in 'Fv'.
		fftw_execute(backward_plan);
		return in;
	}

	FFTW::~FFTW()
	{
		// Destroy the plans, free the allocated space for in and out, and clean up.
		fftw_destroy_plan(forward_plan);
		fftw_destroy_plan(backward_plan);
		fftw_free(in);
		fftw_free(out);
		fftw_cleanup();
	}

	// Note to self: I'm probably better of using the unsupported version of Eigen's implementation.
	template <typename T_real, typename T_complex>
	class FFTSimple : FFTLibrary<T_real, T_complex>
	{
		std::vector<T_complex> FFT(std::vector<T_real> L, int b = 1) // Assumption: Size of L is a power of 2.
		{
			int n = L.size();
			if (n == 1)
			{
				return L;
			}

			T_real theta = (boost::math::pi<T_real>() * 2 * b) / n;
			T_complex w = boost::math::exp(std::complex<T_real>(0, theta)); // Primitive n-th root of unity.

			std::vector<T_complex> L_even;
			std::vector<T_complex> L_odd;
			for (int k = 0; k < n; k++)
			{
				if (k % 2 == 0)
				{
					L_even.push_back(L[k]);
				}
				else
				{
					L_odd.push_back(L[k]);
				}
			}

			std::vector<T_complex> FFT_even = FFT(L_even, b);
			std::vector<T_complex> FFT_odd = FFT(L_odd, b);

			std::vector<T_complex> FFT(n, 0);
			T_complex v(1);
			for (int k = 0; k < n / 2; k++)
			{
				FFT[k] = FFT_even[k] + v * FFT_odd[k];
				FFT[k + n / 2] = FFT_even[k] - v * FFT_odd[k];
				v *= w;
			}

			return FFT;
		}
	};
}*/