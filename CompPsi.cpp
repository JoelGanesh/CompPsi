#include "CompPsi.h"
#include "Fourier.h"
#include "Types.h"

using namespace Types;

#include <chrono>

Elementary::SegmentedSieve Elementary::sieve = Elementary::SegmentedSieve();

int main()
{
    CompPsi::PsiBF psiBF;
    CompPsi::PsiElem psiElem;
    CompPsi::PsiFFT psiFFT;

    uint64_t n = (uint64_t)(1 << 26);
    while (true)
    {
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
        
        auto start = std::chrono::high_resolution_clock::now();
        //std::cout << std::setprecision(100) << "Psi(" << n << "): " << psiFFT.Psi(n) << std::endl;
        std::cout << std::setprecision(100) << "Psi(" << n << "): " << psiBF.Psi(n) << std::endl;
        auto end = std::chrono::high_resolution_clock::now();

        std::cout << "Computation took: " << std::setprecision(8) << (double)std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()/1000000 << "s" << std::endl;

        n *= 2;
        //std::cout << CompPsi::PsiElem::Psi(n) << std::endl;
    }
}