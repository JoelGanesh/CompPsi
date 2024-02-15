#include "CompPsi.h"
#include "Fourier.h"
#include "Types.h"

using namespace Types;

Elementary::SegmentedSieve Elementary::sieve = Elementary::SegmentedSieve();

int main()
{
	CompPsi::PsiBF psiBF;
	CompPsi::PsiElem psiElem;
	CompPsi::PsiFFT psiFFT;

	while (true)
	{
		uint64_t n, k;
		std::cin >> n;
		//std::cin >> k;

		//std::vector<double> v1(n, 1);
		//std::vector<float_dec_100> v2(n, 1);
		//std::vector<std::vector<double>> w1(k, v1);
		//std::vector<std::vector<float_dec_100>> w2(k, v2);

		//Fourier::FFTSimple fftSimple(n, k);
		//Fourier::FFTW fftw(n, k);

		//std::vector<double> x = fftw.Convolve(w1);
		//std::vector<float_dec_100> y = fftSimple.Convolve(w2);

		//Utility::IO::Print(x);
		//std::function<std::string(float_dec_100)> f = [](float_dec_100 q) { return q.str(6); };
		//Utility::IO::Print(Utility::Generic::Map(y, f));
		
		std::cout << std::setprecision(100) << psiFFT.Psi(n) << std::endl;
		std::cout << std::setprecision(100) << psiBF.Psi(n) << std::endl;
		//std::cout << CompPsi::PsiElem::Psi(n) << std::endl;
	}
}