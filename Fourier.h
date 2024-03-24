// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#ifndef FOURIER_H
#define FOURIER_H

#include "Utility.h"
#include "Types.h"

using namespace Types;

namespace Fourier
{
	// Abstract class defining FFT and IFFT operations on vectors.
	// (Real -(FFT)-> Complex -(IFFT)-> Real; RR = Real, CC = Complex)
	template <typename RR, typename CC, typename CC_Container>
	class FFTLibrary
	{
		public:
		// Computes the modular convolution of a list of vectors.
		std::vector<RR> Convolve(std::vector<std::vector<RR>> vs)
		{
			// Initialize the Fourier transform of the convolution of 'vs'.
			CC_Container Fv;
			InitializeContainer(Fv);

			// Compute the DFT of the elements of vs, and
			// store the results in Fv by computing term-wise products.
			for (int i = 0; i < vs.size(); i++)
			{
				// Pad zeros to the end of vs[i] so that the length is size.
				vs[i].insert(vs[i].end(), size - vs[i].size(), RR(0));

				// Apply FFT and process the resulting vector in Fv.
				CC_Container Fv_i = FFT(vs[i]);
				Multiply(Fv, Fv_i);
			}

			// Apply the inverse DFT, remove allocated memory for Fv. The resulting vector should still be normalized.
			std::vector<RR> v = IFFT(Fv);
			DestructContainer(Fv);
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
		virtual void InitializeContainer(CC_Container& Fv) = 0;

		virtual void DestructContainer(CC_Container& Fv) {};
		virtual void Multiply(CC_Container& Fv, CC_Container& Fv_i) = 0;

		// Pure virtual function that computes the Fourier transform of a vector.
		virtual CC_Container FFT(std::vector<RR>& in) = 0;

		// Pure virtual function that computes the inverse Fourier transform of a vector.
		virtual std::vector<RR> IFFT(CC_Container& in) = 0;
	};

	class FFTSimple : public FFTLibrary<float_dec_T, Complex, std::vector<Complex>>
	{
		public:
		FFTSimple(int size, int n);

		private:
		void InitializeContainer(std::vector<Complex>& Fv) override;
		void Multiply(std::vector<Complex>& Fv, std::vector<Complex>& Fv_i) override;
		std::vector<Complex> FFT(std::vector<float_dec_T>& v) override;
		std::vector<float_dec_T> IFFT(std::vector<Complex>& Fv) override;

		std::vector<Complex> IFFT_(std::vector<Complex>& Fv);
	};
}
#endif // FOURIER_H