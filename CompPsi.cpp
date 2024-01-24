#include "CompPsi.h"

Elementary::SegmentedSieve Elementary::sieve = Elementary::SegmentedSieve();

int main()
{
    CompPsi::PsiBF psiBF;
    CompPsi::PsiElem psiElem;
    CompPsi::PsiFFT psiFFT;

    while (true)
    {
        uint64_t n;
        std::cin >> n;
        std::cout << std::setprecision(100) << psiFFT.Psi(n) << std::endl;
        std::cout << std::setprecision(100) << psiBF.Psi(n) << std::endl;
        //std::cout << CompPsi::PsiElem::Psi(n) << std::endl;
    }
}