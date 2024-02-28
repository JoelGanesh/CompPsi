#include "CompPsi.h"
#include "Elementary.h"

#include <numeric>

#include <boost/multiprecision/cpp_dec_float.hpp>

#define LINEAR_APPR 0
#define BRUTE_FORCE 1

namespace CompPsi
{
	// Returns the sign of a real number.
	static int Sgn(double delta)
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

	static uint64_t Mod(int64_t a, uint64_t q)
	{
		// If a < 0, a % q returns the representative of a in {-q+1, -q+2, ..., 0}.
		if (a >= 0)
		{
			return a % q;
		}
		else return a % q + q;
	}

	// Returns largest integer <= n congruent to a mod q.
	// Algorithm by Helfgott & Thompson, 2023.
	static int64_t FlCong(uint64_t n, uint64_t a, uint64_t q)
	{
		return n - Mod(n - a, q);
	}

	// Returns the sum of g(n) over n in I with n = r mod q.
	// G represents the 'q' partial sums of g, restricted to the congruence classes modulo q.
	// Algorithm by Helfgott & Thompson, 2023.
	static float_dec_100 SumInter(std::vector<float_dec_100> G, int64_t r, Interval I, int64_t b, int64_t q)
	{
		int64_t I0 = I.start, I1 = I.end;
		if (I0 <= I1)
		{
			int64_t r0 = FlCong(I0 - 1, r, q);
			int64_t r1 = FlCong(std::min(I1, b - 1), r, q);
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

	// Returns three specific tables of values of g:
	// * G(n) = sum_{m in I_m : m <= n, m = n mod q} g(m);
	// * rho(r) = sum_{m in I_m : a0(m-m0) = r mod q} g(m);
	// * sigma(r) = sum_{m in I_m : a0(m-m0) mod q > q - r} g(m).
	// It is assumed that T can be cast to float_dec_100.
	// Algorithm by Helfgott & Thompson, 2023.
	template <typename T>
	static std::tuple<std::vector<float_dec_100>, std::vector<float_dec_100>, std::vector<float_dec_100>>
		SumTable(std::vector<T> g, uint64_t q, uint64_t b, uint64_t m_index, uint64_t a0)
	{
		std::vector<float_dec_100> G(2 * b), rho(q), sigma(q + 1, 0);
		for (uint64_t m = 0; m < q; m++)
		{
			G[m] = float_dec_100(g[m_index + m]);
		}
		for (uint64_t m = q; m < b; m++)
		{
			G[m] = G[m - q] + float_dec_100(g[m_index + m]);
		}
		uint64_t r = Mod(a0 * (b - q), q);
		uint64_t a = Mod(a0, q);

		for (uint64_t m = b - q; m < b; m++)
		{
			rho[r] = G[m];
			r += a;
			if (r >= q)
			{
				r -= q;
			}
		}

		for (uint64_t r = 1; r < q; r++)
		{
			sigma[r + 1] = sigma[r] + rho[q - r];
		}
		return std::make_tuple(G, rho, sigma);
	}

	// Returns the sum of values of g(n) for n of the form b + kq, for integers k with sk < 0.
	// It is assumed that T can be cast to float_dec_100.
	// Algorithm by Helfgott & Thompson, 2023.
	template <typename T>
	static float_dec_100 RaySum(std::vector<T> g, uint64_t q, uint64_t b, uint64_t m_index, int s)
	{
		float_dec_100 S(0);
		if (s < 0)
		{
			for (uint64_t m = b + q; m < 2 * b; m += q)
			{
				S += float_dec_100(g[m_index + m]);
			}
		}
		else if (s > 0)
		{
			for (uint64_t m = q; m <= b; m += q)
			{
				S += float_dec_100(g[m_index + b - m]);
			}
		}
		return S;
	}

	// Returns sum_{(d,m) in [d0 - a, d0 + a) x [m0 - b, m0 + b)} f(d)g(m) (floor(alpha0 + alpha1 d) + floor(alpha2 m))
	// by separation of variables.
	// It is assumed that T1 and T2 can be cast to float_dec_100.
	// Algorithm by Helfgott & Thompson, 2023.
	template <typename T1, typename T2>
	static float_dec_100 LinearSum(std::vector<T1> f, std::vector<T2> g,
								   uint64_t a, uint64_t b, uint64_t d_index, uint64_t m_index,
								   double alpha0, double alpha1, double alpha2)
	{
		float_dec_100 S_f0(0), S_f1(0), S_g0(0), S_g1(0);
		for (uint64_t d = 0; d < 2 * a; d++)
		{
			S_f0 += float_dec_100(f[d_index + d]);
			S_f1 += float_dec_100(f[d_index + d]) * std::floor(alpha0 + alpha1 * (d - a));
		}
		for (uint64_t m = 0; m < 2 * b; m++)
		{
			S_g0 += float_dec_100(g[m_index + m]);
			S_g1 += float_dec_100(g[m_index + m]) * std::floor(alpha2 * (m - b));
		}
		return S_f0 * S_g1 + S_f1 * S_g0;
	}

	// Computation of L - L1 for the special case where a0(m-m0) + r0 = -1 mod q, assuming q > 1.
	static float_dec_100 Special1(std::vector<float_dec_100> G, uint64_t N, uint64_t q, int64_t a, int64_t a_inv,
								  double R0, int64_t r0, int64_t m0, int64_t d, int64_t b)
	{
		// Note that a0m + r0 = -1 mod q iff m = r mod q, with r = (-1 - r0) * a_inv.
		int64_t r = a_inv * (-1 - r0);

		int64_t gamma1 = d * (-std::floor(R0) * q - (r0 + 1) + a * m0);
		Interval J(-a * d, gamma1, N * q); J.Shift(-m0);

		// J is actually the complement of the set we are trying to sum over.
		// So we first sum over all values, then subtract the ones belonging to J.
		return SumInter(G, r, Interval(LLONG_MIN, LLONG_MAX), b, q) - SumInter(G, r, J, b, q);
	}

	// Computing L - L1 for a0(m-m0) + r0 = 0 mod q, with q > 1.
	// Algorithm by Helfgott & Thompson, 2023. Some comments are added for clarity.
	static float_dec_100 Special0B(std::vector<float_dec_100> G, uint64_t N, uint64_t q, int64_t a, int64_t a_inv,
								   double R0, int64_t r0, int64_t m0, int64_t d, int64_t b, double Q, int s_beta, int s_delta)
	{
		Interval I;
		if (s_delta > 0)
		{
			I = Interval(LLONG_MIN, std::ceil(-Q) - 1);
		}
		else if (s_delta < 0)
		{
			I = Interval(std::ceil(-Q) + 1, LLONG_MAX);
		}
		else if (s_beta < 0)
		{
			I = Interval(LLONG_MIN, LLONG_MAX);
		}

		int64_t gamma1 = d * (-std::floor(R0) * q - r0 + a * m0);
		Interval J(-a * d, gamma1, N * q); J.Shift(-m0);

		// J is actually the complement of the set J' we want; we try to sum over I \cap J'.
		// Note that I \cap J' = I - (J \cap I).
		J.Intersect(I);
		return SumInter(G, -r0 * a_inv, I, b, q) - SumInter(G, -r0 * a_inv, J, b, q);
	}

	// Computation of L - L1 for the special case where q = 1.
	// Algorithm by Helfgott & Thompson, 2023. Some comments are added for clarity.
	static float_dec_100 Special00(std::vector<float_dec_100> G, uint64_t N, uint64_t q, int64_t a, int64_t a_inv,
								   double R0, int64_t r0, int64_t m0, int64_t d, int64_t b, double Q, int s_delta)
	{
		Interval I, I_c;
		if (s_delta > 0)
		{
			I = Interval(LLONG_MIN, std::ceil(-Q) - 1);
			I_c = Interval(std::ceil(-Q), LLONG_MAX);
		}
		else if (s_delta < 0)
		{
			I = Interval(std::floor(-Q) + 1, LLONG_MAX);
			I_c = Interval(LLONG_MIN, std::floor(-Q));
		}

		std::vector<Interval> J(2);
		for (int j = 0; j <= 1; j++)
		{
			if (a != 0)
			{
				int64_t gamma1 = d * (-std::floor(R0) - (r0 + j) + a * m0);
				J[j] = Interval(-a * d, gamma1, N); J[j].Shift(-m0);
			}
			else
			{
				J[j] = Interval(LLONG_MIN, N / (d * ((int64_t)std::floor(R0) + r0 + j) - m0));
			}
		}

		// The intervals J[0] and J[1] are actually the complements of the sets J'[0] and J'[1] we want;
		// we try to sum over the set (I \cap J'[0]) \sqcup (I_c \cap J'[1]).
		// Note that (I \cap J'[0]) \sqcup (I_c \cap J'[1]) = R - (I \cap J[0]) - (I_c \cap J[1]).
		J[0].Intersect(I);
		J[1].Intersect(I_c);
		return SumInter(G, 0, Interval(LLONG_MIN, LLONG_MAX), b, q)
			- SumInter(G, 0, J[0], b, q) - SumInter(G, 0, J[1], b, q);
	}

	// Computation of L1 - L2 for the special case where a0(m-m0) + r0 = 0 mod q.
	// Algorithm by Helfgott & Thompson, 2023. Some comments are added for clarity.
	static float_dec_100 Special0A(std::vector<float_dec_100> G, uint64_t q, int64_t a, int64_t a_inv,
								   int64_t r0, int64_t b, double Q, int s_beta, int s_delta)
	{
		// Note that am + r0 = 0 mod q iff m = r mod q with r = -r0 * a_inv. 
		int r = -r0 * a_inv;

		Interval I;
		if (0 < r0 && r0 < q)
		{
			// The interval I should be taken s.t. m in I iff beta + delta m >= 0.
			if (s_delta > 0)
			{
				I = Interval(std::ceil(-Q), LLONG_MAX);
			}
			else if (s_delta < 0)
			{
				I = Interval(LLONG_MIN, std::floor(-Q));
			}
			else if (s_beta >= 0) // s_delta = 0.
			{
				I = Interval(LLONG_MIN, LLONG_MAX);
			}
		}
		else
		{
			// The interval I should be taken such that m in I if and only if
			// (beta < 0 AND delta m < 0) OR (beta delta m < 0 AND beta + delta m >= 0).
			if (s_beta < 0)
			{
				if (s_delta < 0)
				{
					return SumInter(G, r, Interval(LLONG_MIN, std::floor(-Q)), b, q)
						+ SumInter(G, r, Interval(1, LLONG_MAX), b, q);
				}
				else if (s_delta > 0)
				{
					return SumInter(G, r, Interval(LLONG_MIN, -1), b, q)
						+ SumInter(G, r, Interval(std::ceil(-Q), LLONG_MAX), b, q);
				}
			}
			else if (s_beta > 0)
			{
				if (s_delta > 0)
				{
					I = Interval(std::ceil(-Q), -1);
				}
				else if (s_delta < 0)
				{
					I = Interval(1, std::floor(-Q));
				}
			}
		}
		return SumInter(G, r, I, b, q);
	}

	// Returns sum_{(d,m) in [d0 - a, d0 + a) x [m0 - b, m0 + b)} f(d)g(m)floor(N/dm).
	// Algorithm by Helfgott & Thompson, 2023.
	// It is assumed that a and b are such that the error in
	// the linear approximation of N/xy is less than 1/2b.
	template <typename T1, typename T2>
	static float_dec_100 SumByLinAppr(std::vector<T1> f, std::vector<T2> g, uint64_t N, uint64_t d0, uint64_t m0,
									  uint64_t a, uint64_t b, uint64_t d_index, uint64_t m_index)
	{
		double alpha0 = (double)N / (d0 * m0), alpha1 = -(double)N / (d0 * d0 * m0), alpha2 = -(double)N / (d0 * m0 * m0);

		float_dec_100 S = LinearSum(f, g, a, b, d_index, m_index, alpha0, alpha1, alpha2);
		std::tuple<uint64_t, uint64_t, uint64_t, int> tuple = Elementary::DiophAppr::ApprByRedFrac(alpha2, 2 * b);
		uint64_t a0 = std::get<0>(tuple);
		uint64_t a0_inv = std::get<1>(tuple);
		uint64_t q = std::get<2>(tuple); 
		int sgn_delta = std::get<3>(tuple);
		double delta = alpha2 - a0 / q;

		float_dec_100 Z = RaySum(g, q, b, m_index, sgn_delta);

		std::tuple<std::vector<float_dec_100>, std::vector<float_dec_100>, std::vector<float_dec_100>> 
			tuple2 = SumTable(g, q, b, m_index, a0);
		std::vector<float_dec_100> G = std::get<0>(tuple2);
		std::vector<float_dec_100> rho = std::get<1>(tuple2);
		std::vector<float_dec_100> sigma = std::get<2>(tuple2);

		for (uint64_t d = 0; d < 2 * a; d++)
		{
			if (float_dec_100(f[d_index + d]) != 0)
			{
				double R0 = alpha0 + alpha1 * (d - a);
				double R0_frac = R0 - std::floor(R0);
				int64_t r0 = std::floor(R0_frac * q + 0.5);
				int64_t d_ = d0 + d - a;
				double beta = R0_frac - r0 / q;
				int sgn_beta = Sgn(beta);

				double Q(1);
				if (delta != 0)
				{
					Q = beta / delta;
				}
				float_dec_100 T = sigma[r0] + Special0A(G, q, a0, a0_inv, r0, b, Q, sgn_beta, sgn_delta);

				if (q > 1)
				{
					// Account for case a0(d-d0) + r0 = -1 mod q.
					T += Special1(G, N, q, a0, a0_inv, R0, r0, m0, d_, b);

					// Account for case a0(d-d0) + r0 = 0 mod q.
					T += Special0B(G, N, q, a0, a0_inv, R0, r0, m0, d_, b, Q, sgn_beta, sgn_delta);
				}
				else // Case q = 1.
				{
					T += Special00(G, N, q, a0, a0_inv, R0, r0, m0, d_, b, Q, sgn_delta);
				}

				if (0 < r0 && r0 < q)
				{
					T += Z;
				}

				S += float_dec_100(f[d_index + d]) * T;
			}
		}
		return S;
	}

	// Computes the sum of f(d)g(m)floor(N/dm) over the rectangle [d0, d1) x [m0, m1) by brute force.
	// It is assumed that T1 and T2 can be cast to the type float_dec_100.
	template <typename T1, typename T2>
	static float_dec_100 BruteDoubleSum(uint64_t d0, uint64_t d1, uint64_t m0, uint64_t m1,
										std::vector<T1> f, std::vector<T2> g, uint64_t N)
	{
		float_dec_100 S(0);
		for (uint64_t d = d0; d < d1; d++)
		{
			for (uint64_t m = m0; m < m1; m++)
			{
				S += float_dec_100(f[d - d0]) * float_dec_100(g[m - m0]) * (N / (d * m));
			}
		}
		return S;
	}

	// Computes the sum of f(d)g(m)floor(N/dm) over the rectangle [d0, d1) x [m0, m1) by
	// linear approximation of the term floor(N/dm) on small rectangles of fixed size a * b.
	template <typename T1, typename T2>
	static float_dec_100 DoubleSum(uint64_t d0, uint64_t d1, uint64_t m0, uint64_t m1, uint64_t a, uint64_t b,
								   std::vector<T1> f, std::vector<T2> g, uint64_t N)
	{
		float_dec_100 S(0);
		for (uint64_t mm = m0; mm < m1; mm += 2 * a)
		{
			uint64_t mp = std::min(mm + 2 * a, m1);
			uint64_t m_mid = (mm + mp) / 2, m_len = (mp - mm) / 2;
			uint64_t m_index = mm - m0;
			for (uint64_t dm = d0; dm < d1; dm += 2 * b)
			{
				uint64_t dp = std::min(dm + 2 * a, d1);
				uint64_t d_mid = (dm + dp) / 2, d_len = (dp - dm) / 2;
				uint64_t d_index = dm - d0;
				S += SumByLinAppr(f, g, N, d_mid, m_mid, d_len, m_len, d_index, m_index);
			}
		}
		return S;
	}

	float_dec_100 DDSum(uint64_t A0, uint64_t A1, uint64_t B0, uint64_t B1,
						uint64_t N, uint64_t D, int mode, uint64_t a, uint64_t b)
	{
		float_dec_100 S(0);
		for (uint64_t d0 = A0; d0 < A1; d0 += D)
		{
			uint64_t d1 = std::min(d0 + D, A1);
			std::vector<Log> Lambda = Elementary::sieve.LambdaSegmented(d0, D);
			for (uint64_t m0 = B0; m0 < B1; m0 += D)
			{
				uint64_t m1 = std::min(m0 + D, B1);
				std::vector<int> mu = Elementary::sieve.MuSegmented(m0, D);
				switch (mode)
				{
					case LINEAR_APPR:
						S += DoubleSum(d0, d1, m0, m1, a, b, Lambda, mu, N);

						if (A1 != B1) // Non-diagonal case; consider region mirrored in diagonal.
						{
							Lambda = Elementary::sieve.LambdaSegmented(m0, D);
							mu = Elementary::sieve.MuSegmented(d0, D);
							S += DoubleSum(d0, d1, m0, m1, a, b, mu, Lambda, N);
						}
						break;
					case BRUTE_FORCE:
						S += BruteDoubleSum(d0, d1, m0, m1, Lambda, mu, N);

						if (A0 != B0 || A1 != B1) // Non-diagonal case.
						{
							Lambda = Elementary::sieve.LambdaSegmented(m0, D);
							mu = Elementary::sieve.MuSegmented(d0, D);
							S += BruteDoubleSum(d0, d1, m0, m1, mu, Lambda, N);
						}
						break;
				}
			}
		}
		return S;
	}

	float_dec_100 PsiElem::IndependentVar(uint64_t N, uint64_t M0) // Returns sum_{mdk <= N, m,d <= M0} Lambda(m)mu(d) = sum_{m,d <= M0} Lambda(m)mu(d)floor(N/md).
	{
		float_dec_100 S(0);
		uint64_t A1 = M0 + 1, B1 = M0 + 1;
		int C = 10, D = 8; // Hand-tuned by Helfgott & Thompson, 2023.

		while (A1 >= 2 * std::pow(6 * C * C * C * N, 0.25) &&
			   A1 >= std::sqrt(M0) + 1 &&
			   A1 >= 2 * D)
		{
			B1 = A1;
			uint64_t A = A1 - 2 * (A1 / (2 * D));
			double cbrt = std::cbrt((double)A / (6 * N));
			uint64_t a = A * cbrt;
			while (B1 >= 2 * std::cbrt(6 * C * C * C * N / A) &&
				   B1 >= std::sqrt(M0) + 1 &&
				   B1 >= 2 * D)
			{
				uint64_t B = B1 - 2 * (B1 / (2 * D));
				uint64_t b = B * cbrt;

				uint64_t D = (1 + M0 / std::max(2 * a, 2 * b)) * std::max(2 * a, 2 * b);

				S += DDSum(A, A1, B, B1, N, D, LINEAR_APPR, a, b);

				B1 = B;
			}

			// Remaining part is done by brute force.
			S += DDSum(A, A1, 1, B1, N, std::sqrt(M0) + 1, BRUTE_FORCE, 0, 0);
			A1 = A;
		}

		// The remaining rectangle is done again by brute force.
		S += DDSum(1, A1, 1, B1, N, std::sqrt(M0) + 1, BRUTE_FORCE, 0, 0);
		return S;
	}
}