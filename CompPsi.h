// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#ifndef COMPPSI_H
#define COMPPSI_H

#include "Types.h"
#include "SegmentationArray.h"

#include <vector>

namespace CompPsi
{
	// A brute force (or naive) approach to compute Psi(N).
	class PsiBF
	{
		public:
		// Computes Psi(N) naively, sufficient for "small" values of N.
		static float_dec_T Psi(int64_t N);
	};

	// "Abstract" subclass of CompPsi where an implementation of Psi(N)
	// is based on an underlying implementation of Psi0, to be defined in a subclass.
	class PsiPrep
	{
		public:
		// Computes and returns an approximation of Psi(N).
		float_dec_T Psi(int64_t N);

		private:
		// Pure virtual method implemented in subclasses.
		// Computes and returns an approximation of
		// Psi0(N) := sum_{mdk <= N : m, d <= sqrt(N)} Lambda(m)mu(d).
		virtual float_dec_T Psi0(int64_t N) = 0;

		// Transforms the computation of Psi(N) to computing
		// sum_{mdk <= N : m, d <= sqrt(N)} Lambda(m)mu(d).
		static float_dec_T PsiPreparation(int64_t N);
	};

	// An elementary approach to compute Psi(N) efficiently.
	class PsiElem : public PsiPrep
	{
		private:
		// Computes Psi0(N) using elementary methods.
		// Definition of Psi0 can be found in base class.
		float_dec_T Psi0(int64_t N) override;

		// Computes sum_{mdk <= N, M0 < max(m,d) <= M} Lambda(m)mu(d).
		static float_dec_T DependentVar(int64_t N, int64_t M, int64_t M0);

		// Computes sum_{mdk <= N: m,d <= M0} Lambda(m)mu(d) = sum_{m,d <= M0} Lambda(m)mu(d)floor(N/md).
		static float_dec_T IndependentVar(int64_t N, int64_t M0);
	};

	// An FFT approach to compute Psi(N) efficiently.
	class PsiFFT : public PsiPrep
	{
		private:
		// Computes Psi0(N) using the Fast Fourier Transform.
		// Definition of Psi0 can be found in base class.
		float_dec_T Psi0(int64_t N) override;

		// Applies FFT on segmentation arrays to get a rough approximation of Psi0(N).
		static float_dec_T SegmentationFFT(int64_t N, Elementary::SegmentationArray<float_dec_T>& array);

		// Computes the error in the rough approximation of Psi0(N) using the segmentation arrays.
		static float_dec_T SegmentationError(int64_t N, int32_t S, Elementary::SegmentationArray<float_dec_T>& array);
	};
}
#endif // COMPPSI_H