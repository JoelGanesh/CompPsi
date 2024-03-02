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
			if (p > N)
			{
				break;
			}
			int k = (int)(std::log(N) / std::log(p)); // Represents floor of log_p(N).
			float_dec_100 log = boost::multiprecision::log(float_dec_100(p));
			//std::cout << std::setprecision(100) << k << " * " << log << std::endl;
			sum += k * log;
		}
		//std::cout << "Psi(M) = " << sum << std::endl;
		return sum;
	}
}