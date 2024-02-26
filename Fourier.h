#ifndef FOURIER_H
#define FOURIER_H

#include <vector>
#include <complex.h>
#include <fftw3.h>

#include "Utility.h"
#include "Types.h"
using namespace Types;

#define RE 0
#define IM 1

namespace Fourier
{
	// Abstract class defining FFT and IFFT operations on vectors.
	// (Real -(FFT)-> Complex -(IFFT)-> Real; RR = Real, CC = Complex)
	template <typename RR, typename CC, typename CC_In, typename CC_Out>
	class FFTLibrary
	{
	public:
		// Computes the modular convolution of a list of vectors.
		std::vector<RR> Convolve(std::vector<std::vector<RR>> vs)
		{
			// Initialize the Fourier transform of the convolution of 'vs'.
			CC_In Fv;
			InitializeContainer(Fv);

			// Compute the DFT of the elements of vs, and
			// store the results in Fv by computing term-wise products.
			for (int i = 0; i < vs.size(); i++)
			{
				// Pad zeros to the end of vs[i] so that the length is size.
				vs[i].insert(vs[i].end(), size - vs[i].size(), RR(0));

				// Apply FFT and process the resulting vector in Fv.
				CC_Out Fv_i = FFT(vs[i]);
				Multiply(Fv, Fv_i);
			}

			// Apply the inverse DFT, remove allocated memory for Fv. The resulting vector should still be normalized.
			std::vector<RR> v = IFFT(Fv);
			DestructContainer(Fv);
			//Utility::IO::Print(v);
			for (size_t i = 0; i < v.size(); i++)
			{
				v[i] /= size;
			}
			return v;
		}

	protected:
		// Universal size of vectors to be processed.
		int size = 0;

	private:
		virtual void InitializeContainer(CC_In &Fv) = 0;

		virtual void DestructContainer(CC_In& Fv) {};
		virtual void Multiply(CC_In &Fv, CC_Out Fv_i) = 0;

		// Pure virtual function that computes the Fourier transform of a vector.
		// (Implementation to be provided by subclass.)
		virtual CC_Out FFT(std::vector<RR> in) = 0;

		// Pure virtual function that computes the inverse Fourier transform of a vector.
		// (Implementation to be provided by subclass.)
		virtual std::vector<RR> IFFT(CC_In in) = 0;
	};

	// Subclass of FFTLibrary using the FFTW library.
	class FFTW : public FFTLibrary<double, fftw_complex, fftw_complex*, fftw_complex*>
	{
	public:
		// Constructor
		FFTW(int size, int n)
		{
			// To compute the non-modular convolution on n vectors,
			// the vectors need to be padded with zeros until
			// each vector has size equal to n times the original size.
			this->size = size * n;

			// Initialize arrays for input and output by allocating memory.
			in = (double*)fftw_malloc(sizeof(double) * this->size);
			out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * this->size);

			inv_in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * this->size);
			inv_out = (double*)fftw_malloc(sizeof(double) * this->size);

			// Set up plans for both the forward and backward Fourier transform
			// to be executed multiple times.
			forward_plan = fftw_plan_dft_r2c_1d(this->size, in, out, FFTW_ESTIMATE);
			backward_plan = fftw_plan_dft_c2r_1d(this->size, inv_in, inv_out, FFTW_ESTIMATE);
		}

		// Destructor
		~FFTW()
		{
			// Destroy the plans, free the allocated space for in and out, and clean up.
			fftw_destroy_plan(forward_plan);
			fftw_destroy_plan(backward_plan);
			fftw_free(in); fftw_free(out);
			fftw_free(inv_in); fftw_free(inv_out);
			fftw_cleanup();
		}

	private:
		// Serves for real input to compute forward/backward DFTs.
		double* in;
		fftw_complex* inv_in;

		// Serves for complex output to store forward/backward DFTs.
		fftw_complex* out;
		double* inv_out;

		// A plan created to compute FFTs with data from 'in' (real) to 'out' (complex).
		fftw_plan forward_plan;
		fftw_plan backward_plan;

		void InitializeContainer(fftw_complex* &Fv) override
		{
			Fv = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * this->size);
			for (int i = 0; i < size; i++)
			{
				Fv[i][RE] = 1.0;
				Fv[i][IM] = 0.0;
			}
		}

		void DestructContainer(fftw_complex*& Fv) override
		{
			fftw_free(Fv);
		}

		void Multiply(fftw_complex* &Fv, fftw_complex* Fv_i) override
		{
			for (int j = 0; j < size; j++)
			{
				double temp = Fv[j][RE] * Fv_i[j][RE] - Fv[j][IM] * Fv_i[j][IM];

				Fv[j][IM] = Fv[j][RE] * Fv_i[j][IM] + Fv[j][IM] * Fv_i[j][RE];
				Fv[j][RE] = temp;
			}
		}

		// Implementation of FFT using the FFTW library.
		fftw_complex* FFT(std::vector<double> v) override
		{
			for (int j = 0; j < size; j++)
			{
				in[j] = v[j];
			}
			fftw_execute(forward_plan);
			return out;
		}

		// Implementation of IFFT using the FFTW library.
		std::vector<double> IFFT(fftw_complex* Fv) override
		{
			for (int j = 0; j < size; j++)
			{
				inv_in[j][RE] = Fv[j][RE];
				inv_in[j][IM] = Fv[j][IM];
			}
			fftw_execute(backward_plan);

			std::vector<double> result(inv_out, inv_out + size);
			return result;
		}

	};

	class FFTSimple : public FFTLibrary<float_dec_100, Complex, std::vector<Complex>, std::vector<Complex>>
	{
	public:
		FFTSimple(int size, int n)
		{
			int k = size * n;
			int N = 1;
			while (N < k)
			{
				N *= 2;
			}
			this->size = N;
		}

	private:
		void InitializeContainer(std::vector<Complex>& Fv) override
		{
			Fv = std::vector<Complex>(size, Complex(1, 0));
		}

		void Multiply(std::vector<Complex>& Fv, std::vector<Complex> Fv_i) override
		{
			//std::function<std::string(Complex)> g = [](Complex z) { return std::string(z); };
			//std::cout << "Fv: ";
			//Utility::IO::Print(Utility::Generic::Map(Fv, g));
			//std::cout << std::endl << "Fv_i: ";
			//Utility::IO::Print(Utility::Generic::Map(Fv_i, g));

			for (int j = 0; j < size; j++)
			{
				Fv[j] *= Fv_i[j]; // Multiplication is defined in struct Complex (in Types.h)
			}
			//std::cout << std::endl << "Fv: ";
			//Utility::IO::Print(Utility::Generic::Map(Fv, g));
			//std::cout << std::endl;
		}

		std::vector<Complex> FFT(std::vector<float_dec_100> v) override
		{
			int n = v.size();
			std::vector<Complex> Fv(n);
			if (n == 1)
			{
				Fv[0].re = v[0];
				return Fv;
			}

			int m = n / 2;
			std::vector<float_dec_100> v_even(m), v_odd(m);
			for (int k = 0; k < m; k++)
			{
				v_even[k] = v[2 * k];
				v_odd[k] = v[2 * k + 1];
			}

			std::vector<Complex> Fv_even = FFT(v_even), Fv_odd = FFT(v_odd);

			for (int k = 0; k < m; k++)
			{
				float_dec_100 theta = (2 * boost::math::constants::pi<float_dec_100>() * k) / n;
				Complex w(boost::multiprecision::cos(theta), boost::multiprecision::sin(theta));
				Fv[k] = Fv_even[k] + w * Fv_odd[k];
				Fv[k + m] = Fv_even[k] - w * Fv_odd[k];
			}

			return Fv;
		}

		std::vector<float_dec_100> IFFT(std::vector<Complex> Fv) override
		{
			std::vector<Complex> v = IFFT_(Fv);
			std::vector<float_dec_100> w(v.size());
			for (int i = 0; i < w.size(); i++)
			{
				w[i] = v[i].re;
			}
			return w;
		}

		std::vector<Complex> IFFT_(std::vector<Complex> Fv)
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
				float_dec_100 theta = (2 * boost::math::constants::pi<float_dec_100>() * k) / n;
				Complex w(boost::multiprecision::cos(theta), -boost::multiprecision::sin(theta));
				v[k] = v_even[k] + w * v_odd[k];
				v[k + m] = v_even[k] - w * v_odd[k];
			}

			return v;
		}
	};
}
#endif // FOURIER_H