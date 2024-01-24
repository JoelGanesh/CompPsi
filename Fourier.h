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
	template <typename T_real, typename T_complex>
	class FFTLibrary
	{
	public:
		// Computes the modular convolution of a list of vectors.
		std::vector<T_real> Convolve(std::vector<std::vector<T_real>> vs)
		{
			// Initialize the Fourier transform of the convolution of 'vs'.
			std::unique_ptr<T_complex[]> Fv = std::make_unique<T_complex[]>(size);
			std::function<int(std::complex<float_dec_100>)> f = [](std::complex<float_dec_100> c) { std::cout << c.real() << "+" << c.imag() << "i" << " "; return 0; };
			for (int i = 0; i < size; i++)
			{
				Fv[i] = { float_dec_100(1), float_dec_100(0) }; // Fv[i][IM] is already set to 0.
				f(Fv[i]);
			}
			std::cout << std::endl;

			// Compute the DFT of the elements of vs, and
			// store the results in Fv by computing term-wise products.
			for (int i = 0; i < vs.size(); i++)
			{
				// Pad zeros to the end of vs[i] so that the length is size.
				T_real zero(0);
				vs[i].insert(vs[i].end(), size - vs[i].size(), zero);

				// Apply FFT 
				T_complex* Fv_i = FFT(vs[i]);

				for (int i = 0; i < size; i++)
				{
					f(Fv_i[i]);
				}
				//Utility::IO::Print(Fv_i);
				for (int j = 0; j < size; ++j)
				{
					Fv[j] *= Fv_i[j];
					// We multiply Fv[j] by Fv_i[j].
					/*float_dec_100 Fv_j_Re = Fv[j].real();
					float_dec_100 Fv_j_Im = Fv[j].imag();
					float_dec_100 Fv_i_j_Re = Fv_i[j].real();
					float_dec_100 Fv_i_j_Im = Fv_i[j].imag();

					float_dec_100 new_Fv_j_Re = Fv_j_Re * Fv_i_j_Re - Fv_j_Im * Fv_i_j_Im;
					float_dec_100 new_Fv_j_Im = Fv_j_Re * Fv_i_j_Im + Fv_j_Im * Fv_i_j_Re;
					Fv[j] = std::complex<float_dec_100>(new_Fv_j_Re, new_Fv_j_Im);*/
					//Fv[j][RE] = Fv_j_Re * Fv_i[j][RE] - Fv_j_Im * Fv_i[j][IM];
					//Fv[j][IM] = Fv_j_Re * Fv_i[j][IM] + Fv_j_Im * Fv_i[j][RE];
				}
			}

			// Apply the inverse DFT, normalize the result and return it.
			std::vector<T_real> v = IFFT(Fv.get());
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
		// Pure virtual function that computes the Fourier transform of a vector.
		// (Implementation to be provided by subclass.)
		virtual T_complex* FFT(std::vector<T_real> in) = 0;

		// Pure virtual function that computes the inverse Fourier transform of a vector.
		// (Implementation to be provided by subclass.)
		virtual std::vector<T_real> IFFT(T_complex* in) = 0;
	};

	// Subclass of FFTLibrary using the FFTW library.
	class FFTW : public FFTLibrary<double, fftw_complex>
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

	class FFTSimple : public FFTLibrary<float_dec_100, std::complex<float_dec_100>>
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
		std::complex<float_dec_100>* FFT(std::vector<float_dec_100> v) override
		{
			std::vector<std::complex<float_dec_100>> Fv = FFT_(v);
			return &Fv[0];
		}

		std::vector<float_dec_100> IFFT(std::complex<float_dec_100>* Fv) override
		{
			std::vector<std::complex<float_dec_100>> Fv_(Fv, Fv + size);
			return IFFT_(Fv_);
		}

		std::vector<std::complex<float_dec_100>> FFT_(std::vector<float_dec_100> v)
		{
			int n = v.size();
			std::vector<std::complex<float_dec_100>> Fv(n, { float_dec_100(0), float_dec_100(0) });
			if (n == 1)
			{
				Fv[0] = { v[0], float_dec_100(0) };
				return Fv;
			}

			int m = n >> 1;
			std::vector<float_dec_100> v_even(m), v_odd(m);
			for (int k = 0; k < m; k += 1)
			{
				v_even[k] = v[2 * k];
				v_odd[k] = v[2 * k + 1];
			}

			std::vector<std::complex<float_dec_100>> Fv_even = FFT_(v_even), Fv_odd = FFT_(v_odd);

			for (int k = 0; k < m; k++)
			{
				float_dec_100 theta = 2 * boost::math::constants::pi<float_dec_100>() * k / n;
				std::complex<float_dec_100> w = { boost::multiprecision::cos(theta), boost::multiprecision::sin(theta) };
				Fv[k] = Fv_even[k] + w * Fv_odd[k];
				Fv[k + m] = Fv_even[k] - w * Fv_odd[k];
			}
		}

		std::vector<float_dec_100> IFFT_(std::vector<std::complex<float_dec_100>> Fv)
		{
			int n = Fv.size();
			std::vector<float_dec_100> v(n);
			if (n == 1)
			{
				v[0] = Fv[0].real();
				return v;
			}

			int m = n >> 1;
			std::vector<std::complex<float_dec_100>> Fv_even(m), Fv_odd(m);
			for (int k = 0; k < m; k += 1)
			{
				Fv_even[k] = Fv[2 * k];
				Fv_odd[k] = Fv[2 * k + 1];
			}

			std::vector<float_dec_100> v_even = IFFT_(Fv_even), v_odd = IFFT_(Fv_odd);

			for (int k = 0; k < m; k++)
			{
				float_dec_100 theta = 2 * boost::math::constants::pi<float_dec_100>() * k / n;
				std::complex<float_dec_100> w = { boost::multiprecision::cos(-theta), boost::multiprecision::sin(-theta) };
				std::complex<float_dec_100> temp = v_even[k] + w * Fv_odd[k];
				v[k] = temp.real();
				temp = v_even[k] - w * Fv_odd[k];
				v[k + m] = temp.real();
			}

			return v;
		}
	};
}
#endif // FOURIER_H