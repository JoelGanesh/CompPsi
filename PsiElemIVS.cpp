#include "CompPsi.h"
#include "Elementary.h"

#include <numeric>

#include <boost/multiprecision/cpp_dec_float.hpp>

namespace CompPsi
{
	float_dec_100 PsiElem::IndependentVar(uint64_t N, uint64_t M0) // Returns sum_{mdk <= N, m,d <= M0} Lambda(m)mu(d) = sum_{m,d <= M0} Lambda(m)mu(d)floor(N/md).
	{
		float_dec_100 sum = 0;

		for (uint64_t A = M0 + 1; A >= 1; A /= 2)
		{
			double cbrt = std::cbrt((double)A / (6 * N));
			uint64_t a = (uint64_t)std::floor(A * cbrt);
			for (uint64_t B = M0 + 1; B >= 1; B /= 2)
			{
				uint64_t b = (uint64_t)std::floor((double)B * cbrt);

				// Compute the sum over the region [A, 2A) x [B, 2B) and add to total.
				sum += IndependentVar(N, M0, A, B, a, b);
			}
		}
		return sum;
	}


	// Computes sum_{mdk <= N, m <= d <= M0} Lambda(m)mu(d) = sum_{m <= d <= M0} Lambda(m)mu(d)floor(N/md).
	float_dec_100 PsiElem::IndependentVar_md(uint64_t N, uint64_t M0)
	{
		
	}

	// Computes sum_{mdk <= N, d <= m <= M0} Lambda(m)mu(d) = sum_{d <= m <= M0} Lambda(m)mu(d)floor(N/md).
	float_dec_100 PsiElem::IndependentVar_dm(uint64_t N, uint64_t M0)
	{
		float_dec_100 sum(0);

	}

	// Computes sum_{(m,d) in [A, 2A) x [B, 2B)} Lambda(m)mu(d)floor(N/md) in batches of rectangles (size 2a * 2b)
	float_dec_100 PsiElem::IndependentVar(uint64_t N, uint64_t M0, uint64_t A, uint64_t B, uint64_t a, uint64_t b)
	{
		float_dec_100 sum(0);

		// Make sure that 2A and 2B do not overshoot M0 + 1.
		uint64_t A2 = std::min(2 * A, M0 + 1);
		uint64_t B2 = std::min(2 * B, M0 + 1);
		for (uint64_t m_ = A; m_ < A2; m_ += 2 * a)
		{
			for (uint64_t d_ = B; d_ < B2; d_ += 2 * b)
			{
				// We consider the sum over the rectangle R = [m_, m_ + 2a) x [d_, d_ + 2b),
				// restricted to [A, 2A) x [B, 2B).
				Rectangle R(N, M0, m_, std::min(m_ + 2 * a, A2), d_, std::min(d_ + 2 * b, B2));
				sum += IndependentVar(R);
			}
		}
		return sum;
	}

	float_dec_100 PsiElem::IndependentVar(Rectangle R) // Computes sum_{(m,d) in [m0 - a, m0 + a) x [d0 - b, d0 + b)} Lambda(m)mu(d)floor(N/md)
	{
		return IndependentVar_S0mS1(R) + IndependentVar_S1mS2(R) + IndependentVar_S2(R);
	}

	float_dec_100 PsiElem::IndependentVar_S0mS1(Rectangle R)
	{
		int max_index = 1 + (2 * R.b) / R.q;
		std::vector<std::vector<int>> rho_table(R.q, std::vector<int>());

		std::vector<int> mu = Elementary::sieve.MuSegmented(R.d0 - R.b, R.b);
		for (uint64_t d = R.d0 - R.b; d < std::min(R.d0 + R.b, R.M + 1); d++)
		{
			int residue = R.a0 * (d - R.d0) % R.q;
			rho_table[residue].push_back(mu[d]);
		}
		for (uint64_t r = 0; r < R.q; r++)
		{
			std::vector<int> v = rho_table[r];
			std::partial_sum(v.begin(), v.end(), rho_table[r].begin());
		}

		for (uint64_t m = R.m0 - R.a; m < R.m0 + R.a; m++)
		{
			double beta; uint64_t r0;
			std::tuple<double, uint64_t>(beta, r0) = R.beta_r0(m);
			for (uint64_t d = R.d0 - R.b; d < R.d0 + R.b; d++)
			{

			}
		}
		return 0;
	}

	float_dec_100 PsiElem::IndependentVar_S1mS2(Rectangle R)
	{
		return 0;
	}

	float_dec_100 PsiElem::IndependentVar_S2(Rectangle R) // r0 and beta are dependent on the value of m.
	{
		return 0;
	}

	// Algorithm by Helfgott & Thompson (2023).
	// Returns sum_{(d,m) in [d0 - a, d0 + a) x [m0 - b, m0 + b)} f(d)g(m)floor(N/dm).
	// It is assumed that a and b are such that the error in the linear approximation of N/xy is less than 1/2b.
	template <typename T1, typename T2>
	float_dec_100 PsiElem::SumByLinAppr(std::vector<T1> f, std::vector<T2> g, uint64_t N, uint64_t d0, uint64_t m0, uint64_t a, uint64_t b)
	{
		double alpha0 = (double)N / (d0 * m0), alpha1 = (double)-N / (d0 * d0 * m0), alpha2 = (double)-N / (d0 * m0 * m0);
		
		float_dec_100 S = LinearSum(f, g, a, b, alpha0, alpha1, alpha2);
		uint64_t a0, uint64_t a0_inv, uint64_t q, int sgn_delta;
		(a0, a0_inv, q, sgn_delta) = Elementary::ApprByRedFrac(alpha2, 2 * b);
		double delta = alpha2 - a0 / q;

		float_dec_100 Z = RaySum(g, q, b, s);
		std::vector<float_dec_100> G, rho, sigma;
		(G, rho, sigma) = SumTable(g, q, b, a0);

		for (uint64_t d = 0; d < 2 * a)
		{
			if (f[d] != 0)
			{
				double R0 = alpha0 + alpha1 * d;
				double R0_frac = R0 - std::floor(R0);
				int64_t r0 = std::floor(R0_frac * q + 0.5);
				int64_t d_ = d0 + d;
				double beta = R0_frac - r0 / q;
				int sgn_beta = Sgn(beta);

				int64_t Q(1);
				if (delta != 0)
				{
					Q = beta / delta;
				}
				float_dec_100 T = sigma[r0] + Special0A(G, q, a0, a0_inv, r0, b, Q, sgn_beta, sgn_delta);

				if (q > 1)
				{
					// Account for case a0(d-d0) + r0 = 1 mod q.
					T += Special1(G, N, q, a0, a0_inv, R0, r0, m0, d_, b);

					// Account for case a0(d-d0) + r0 = 0 mod q.
					T += Special0B(G, N, q, a0, a0_inv, R0, r0, m0, d_, b, Q, sgn_beta, sgn_delta);
				}
				else // Case q = 1.
				{
					T += Special00(G, N, q, a0, a0_inv, R0, r0, m0, d_, b, Q, sgn_delta);
				}
				if (0 < r0 < q)
				{
					T += Z;
				}

				S += f[d] * T;
			}
		}
		return S;
	}

	float_dec_100 Special1(std::vector<float_dec_100> G, uint64_t N, uint64_t q, int64_t a, int64_t a_inv,
						   double R0, int64_t r0, int64_t n0, int64_t m, int64_t b)
	{
		int64_t gamma1 = m * (-std::floor(R0) * q - (r0 + 1) + a * n0);
		int64_t r = a_inv * (-1 - r0);
		Interval I(-a * m, gamma1, N * q); I.Shift(-n0);

		return SumInter(G, r, Interval(LLONG_MIN, LLONG_MAX), b, q) - SumInter(G, r, I, b, q);
	}

	float_dec_100 Special0B(std::vector<float_dec_100> G, uint64_t N, uint64_t q, int64_t a, int64_t a_inv,
						   double R0, int64_t r0, int64_t n0, int64_t m, int64_t b, double Q, int s_beta, int s_delta)
	{
		int64_t gamma1 = m * (-std::floor(R0) * q - r0 + a * n0);
		Interval I(-a * m, gamma1, N * q); I.Shift(-n0);
		Interval J;
		if (s_delta > 0)
		{
			J = Interval(LLONG_MIN, -std::floor(Q) - 1);
		}
		else if (s_delta < 0)
		{
			J = Interval(-std::floor(Q) + 1, LLONG_MAX);
		}
		else if (s_beta < 0)
		{
			J = Interval(LLONG_MIN, LLONG_MAX);
		}
		I.Intersect(J);
		return SumInter(G, -r0 * a_inv, J, b, q) - SumInter(G, -r0 * a_inv, I, b, q);
	}

	float_dec_100 Special00(std::vector<float_dec_100> G, uint64_t N, uint64_t q, int64_t a, int64_t a_inv,
							double R0, int64_t r0, int64_t n0, int64_t m, int64_t b, double Q, int s_delta)
	{
		Interval J;
		if (s_delta > 0)
		{
			J = Interval(LLONG_MIN, -std::floor(Q) - 1);
		}
		else if (s_delta < 0)
		{
			J = Interval(-std::floor(Q) + 1, LLONG_MAX);
		}

		std::vector<Interval> I(2);
		for (int j = 0; j <= 1; j++)
		{
			if (a != 0)
			{
				int64_t gamma1 = m * (-std::floor(R0) - (r0 + j) + a * n0);
				I[j] = Interval(-a * m, gamma1, N); I[j].Shift(-n0);
			}
			else
			{
				I[j] = Interval(LLONG_MIN, N / (m * (std::floor(R0) + r0 + j) - n0));
			}
		}

		I[0].Intersect(J);

		// We now consider the complement of J.
		// As this consists of two intervals, we have to split the computation.
		Interval JC1(LLONG_MIN, J.start - 1), JC2(J.end + 1, LLONG_MAX);

		JC1.Intersect(I[1]);
		JC2.Intersect(I[1]);
		return SumInter(G, 0, Interval(LLONG_MIN, LLONG_MAX), b, q) 
			- SumInter(G, 0, I[0], b, q) - SumInter(G, 0, JC1, b, q) - SumInter(G, 0, JC2, b, q);
	}

	// Returns sum_{(d,m) in [d0 - a, d0 + a) x [m0 - b, m0 + b)} f(d)g(m) (floor(alpha0 + alpha1 d) + floor(alpha2 m))
	// by separation of variables.
	template <typename T1, typename T2>
	float_dec_100 PsiElem::LinearSum(std::vector<T1> f, std::vector<T2> g, uint64_t a, uint64_t b, double alpha0, double alpha1, double alpha2)
	{
		float_dec_100 S_f0(0), S_f1(0), S_g0(0), S_g1(0);
		for (uint64_t d = d0 - a; d < d0 + a; d++)
		{
			S_f0 += float_dec_100(f[d]);
			S_f1 += float_dec_100(f[d]) * std::floor(alpha0 + alpha1 * d);
		}
		for (uint64_t m = m0 - b; m < m0 + b; m++)
		{
			S_g0 += float_dec_100(g[m]);
			S_g1 += float_dec_100(g[m]) * std::floor(alpha2 * m);
		}
		return S_f0 * S_g1 + S_f1 * S_g0;
	}

	template <typename T>
	float_dec_100 RaySum(std::vector<T> g, uint64_t q, uint64_t b, double delta)
	{
		float_dec_100 S(0);
		if (delta < 0)
		{
			for (uint64_t m = q; m < b; m += q)
			{
				S += g[b + m];
			}
		}
		else if (delta > 0)
		{
			for (uint64_t m = q; m <= b; m += q)
			{
				S += g[b - m];
			}
		}
	}

	template <typename T>
	std::tuple<std::vector<float_dec_100>, std::vector<float_dec_100>, std::vector<float_dec_100>>
		SumTable(std::vector<T> g, uint64_t q, uint64_t b, uint64_t a0)
	{
		std::vector<float_dec_100> G(2 * b), rho(q), sigma(q + 1);
		for (uint64_t n = 0; n < q; n++)
		{
			G[n] = g[n];
		}
		for (uint64_t n = q; n < b; n++)
		{
			G[n] = G[n - q] + g[n];
		}
		uint64_t r = Mod(a0 * (b - q), q);
		uint64_t a = Mod(a0, q);

		for (uint64_t n = 0; n < q; n++)
		{
			rho[r] = G[b + n];
			r += a;
			if (r >= q)
			{
				r -= q;
			}
		}

		sigma[0] = 0;
		sigma[1] = 0;
		for (uint64_t r = 1; r < q; r++)
		{
			sigma[r + 1] = sigma[r] + rho[q - r];
		}
		return std::make_tuple(G, rho, sigma);
	}

	// Returns the sum of g(n) over n in I with n = r mod q.
	// G represents the 'q' partial sums of g, restricted to the congruence classes modulo q.
	// Algorithm by Helfgott & Thompson, 2023.
	float_dec_100 SumInter(std::vector<float_dec_100> G, int64_t r, Interval I, int64_t b, int64_t q)
	{
		int64_t I0 = I.start, I1 = I.end;
		if (I0 <= I1)
		{
			int64_t r0 = FlCong(I0 - 1, r, q),
					r1 = FlCong(std::min(I1, b - 1), r, q);
			if (r0 > r1 || r1 < -b)
			{
				return 0;
			}
			else if (r0 < -b)
			{
				return G[r1 + b];
			}
			return G[r1 + b] - G[r0 + b];
		}
		return 0;
	}

	// Returns largest integer <= n congruent to a mod q.
	// Algorithm by Helfgott & Thompson, 2023.
	int64_t FlCong(uint64_t n, uint64_t a, uint64_t q)
	{
		return n - Mod(n - a, q);
	}

	uint64_t Mod(int64_t a, uint64_t q)
	{
		// If a < 0, a % q returns the representative of a in {-q+1, -q+2, ..., 0}.
		if (a >= 0)
		{
			return a % q;
		}
		if (a < 0)
		{
			return a % q + q;
		}
	}

	// Returns the sign of a real number.
	int Sgn(double delta)
	{
		if (delta < 0)
		{
			return -1;
		}
		if (delta > 0)
		{
			return 1;
		}
		return 0;
	}
}