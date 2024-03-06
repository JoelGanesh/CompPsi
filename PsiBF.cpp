#include "CompPsi.h"
#include <boost/multiprecision/cpp_dec_float.hpp>

namespace CompPsi
{
	float_dec_100 PsiBF::Psi(int64_t N)
	{
		float_dec_100 sum = 0;
		std::vector<Prime> primes = Elementary::sieve.Primes(N);
		for (Prime p : primes)
		{
			//if (p > N)
			//{
			//	break;
			//}
			//std::cout << std::setprecision(100) << k << " * " << log << std::endl;
			sum += Elementary::log(N, p) * boost::multiprecision::log(float_dec_100(p));
		}
		//std::cout << "Psi(M) = " << sum << std::endl;
		return sum;
	}
}