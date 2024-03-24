// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#include "CompPsi.h"
#include "Elementary.h"

#include <chrono>

Elementary::SegmentedSieve Elementary::sieve = Elementary::SegmentedSieve();

enum COMP_MODE
{
	INVALID,
	BF,
	ELEM,
	FFT
};

int main()
{
	CompPsi::PsiBF psiBF;
	CompPsi::PsiElem psiElem;
	CompPsi::PsiFFT psiFFT;

	COMP_MODE mode = INVALID;
	char c;
	int64_t n;
	int m;
	while (mode == INVALID)
	{
		std::cout << "Usage: <mode> <initial value> <multiplier>\nmode: 'B' <-> Brute Force | 'E' <-> Elementary | 'F' <-> FFT" << std::endl;
		std::cin >> c >> n >> m;
		switch (c)
		{
			case 'B':
			case 'b':
				mode = BF;
				break;
			case 'E':
			case 'e':
				mode = ELEM;
				break;
			case 'F':
			case 'f':
				mode = FFT;
				break;
			default:
				break;
		}
	}

	std::chrono::steady_clock::time_point start, end;
	for (;;)
	{
		start = std::chrono::high_resolution_clock::now();
		std::string s;
		switch (mode)
		{
			case BF:
				s = psiBF.Psi(n).str(PRECISION);
				break;
			case ELEM:
				s = psiElem.Psi(n).str(PRECISION);
				break;
			case FFT:
				s = psiFFT.Psi(n).str(PRECISION);
				break;
		}

		std::cout << "Psi(" << n << ") = " << s << std::endl;
		end = std::chrono::high_resolution_clock::now();
		std::cout << "Computation time = " << std::setprecision(8) << (double)std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000000 << "s" << std::endl;

		n *= m;
	}
}