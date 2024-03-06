#include "CompPsi.h"
#include "Elementary.h"

#include <numeric>
#include <functional>
#include <boost/multiprecision/cpp_dec_float.hpp>

#define INFTY	LLONG_MAX
#define MINFTY  LLONG_MIN

namespace CompPsi
{
	// Used to distinguish computation on dyadic boxes done
	// using linear approximation or by brute force.
	enum COMP_MODE
	{
		LINEAR_APPR,
		BRUTE_FORCE
	};

	// Used to distinguish diagonal dyadic boxes from non-diagonal ones.
	enum DIAG_MODE
	{
		DIAGONAL,
		NON_DIAG
	};

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

	static int64_t Mod(int64_t a, int64_t q)
	{
		// If a < 0, a % q returns the representative of a in {-q+1, -q+2, ..., 0}.
		if (a >= 0 || a % q == 0)
		{
			return a % q;
		}
		else return a % q + q;
	}

	// Returns largest integer <= n congruent to a mod q.
	// Algorithm by Helfgott & Thompson, 2023.
	static int64_t FlCong(int64_t n, int64_t a, int64_t q)
	{
		if (n == MINFTY || n == INFTY)
		{
			return n;
		}
		return n - Mod(n - a, q);
	}

	// Returns the sum of g(n) over n in I with n = r mod q.
	// G represents the 'q' partial sums of g, restricted to the congruence classes modulo q.
	// Algorithm by Helfgott & Thompson, 2023.
	template <typename T>
	static T SumInter(std::vector<T> G, int64_t r, Interval I, int64_t b, int64_t q)
	{
		int64_t I0 = I.start, I1 = I.end;
		if (I0 <= I1)
		{
			if (I0 != MINFTY && I0 != INFTY)
			{
				I0--;
			}
			int64_t r0 = FlCong(I0, r, q);
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
	static std::tuple<std::vector<T>, std::vector<T>, std::vector<T>>
		SumTable(std::function<T(int64_t)> g, int64_t q, int64_t b, int64_t a0)
	{
		std::vector<T> G(2 * b), rho(q), sigma(q + 1, 0);
		for (int64_t m = -b; m < -b + q; m++)
		{
			G[m + b] = g(m);
		}
		for (int64_t m = -b + q; m < b; m++)
		{
			G[m + b] = G[m + b - q] + g(m);
		}
		int64_t r = Mod(a0 * (b - q), q);
		int64_t a = Mod(a0, q);

		for (int64_t m = b - q; m < b; m++)
		{
			rho[r] = G[m + b];
			r += a;
			if (r >= q)
			{
				r -= q;
			}
		}

		for (int64_t r = 1; r < q; r++)
		{
			sigma[r + 1] = sigma[r] + rho[q - r];
		}
		return std::make_tuple(G, rho, sigma);
	}

	// Returns the sum of values of g(n) for n of the form b + kq, for integers k with sk < 0.
	// It is assumed that T can be cast to float_dec_100.
	// Algorithm by Helfgott & Thompson, 2023.
	template <typename T>
	static T RaySum(std::function<T(int64_t)> g, int64_t q, int64_t b, int s)
	{
		T S(0);
		if (s < 0)
		{
			for (int64_t m = q; m < b; m += q)
			{
				S += g(m);
			}
		}
		else if (s > 0)
		{
			for (int64_t m = q; m <= b; m += q)
			{
				S += g(-m);
			}
		}
		return S;
	}

	// Returns sum_{(d,m) in [d0 - a, d0 + a) x [m0 - b, m0 + b)} f(d)g(m) (floor(alpha0 + alpha1 d) + floor(alpha2 m))
	// by separation of variables.
	// It is assumed that T1 and T2 can be cast to float_dec_100.
	// Algorithm by Helfgott & Thompson, 2023.
	template <typename Tf, typename Tg>
	static float_dec_100 LinearSum(std::function<Tf(int64_t)> f, std::function<Tg(int64_t)> g,
								   int64_t a, int64_t b,
								   Fraction alpha0, Fraction alpha1, Fraction alpha2)
	{
		Tf S_f0(0), S_f1(0);
		Tg S_g0(0), S_g1(0);
		for (int64_t d = -a; d < a; d++)
		{
			S_f0 += f(d);
			S_f1 += f(d) * (alpha0 + alpha1 * Fraction(d)).Floor();
		}
		for (int64_t m = -b; m < b; m++)
		{
			S_g0 += g(m);
			S_g1 += g(m) * (alpha2 * Fraction(m)).Floor();
		}
		return S_f0 * S_g1 + S_f1 * S_g0;
	}

	// Computation of L - L1 for a0(m-m0) + r0 = -1 mod q, with q > 1.
	template <typename T>
	static T Special1(std::vector<T> G, int64_t N, int64_t q, int64_t a, int64_t a_inv,
								  int64_t R0, int64_t r0, int64_t m0, int64_t d, int64_t b)
	{
		// Note that a0m + r0 = -1 mod q iff m = r mod q, with r = (-1 - r0) * a_inv.
		int64_t r = a_inv * (-1 - r0);

		int64_t gamma1 = d * (-R0 * q - (r0 + 1) + a * m0);
		Interval J(-a * d, gamma1, N * q); J.Shift(-m0);

		// J is actually the complement of the set we are trying to sum over.
		// So we first sum over all values, then subtract the ones belonging to J.
		return SumInter<T>(G, r, Interval(MINFTY, INFTY), b, q) - SumInter(G, r, J, b, q);
	}

	// Computing L - L1 for a0(m-m0) + r0 = 0 mod q, with q > 1.
	// Algorithm by Helfgott & Thompson, 2023. Some comments are added for clarity.
	template <typename T>
	static T Special0B(std::vector<T> G, int64_t N, int64_t q, int64_t a, int64_t a_inv,
					   int64_t R0, int64_t r0, int64_t m0, int64_t d, int64_t b, double Q, int s_beta, int s_delta)
	{
		Interval I;
		if (s_delta > 0)
		{
			I = Interval(MINFTY, std::ceil(-Q) - 1);
		}
		else if (s_delta < 0)
		{
			I = Interval(std::floor(-Q) + 1, INFTY);
		}
		else if (s_beta < 0) // delta = 0
		{
			I = Interval(MINFTY, INFTY);
		}

		int64_t gamma1 = d * (-R0 * q - r0 + a * m0);
		Interval J(-a * d, gamma1, N * q); J.Shift(-m0);

		// J is actually the complement of the set J' we want; we try to sum over I \cap J'.
		// Note that I \cap J' = I - (I \cap J).
		J.Intersect(I);
		return SumInter(G, -r0 * a_inv, I, b, q) - SumInter(G, -r0 * a_inv, J, b, q);
	}

	// Computation of L - L1 for the special case where q = 1.
	// Algorithm by Helfgott & Thompson, 2023. Some comments are added for clarity.
	template <typename T>
	static T Special00(std::vector<T> G, int64_t N, int64_t q, int64_t a, int64_t a_inv,
					   int64_t R0, int64_t r0, int64_t m0, int64_t d, int64_t b, double Q, int s_delta)
	{
		Interval I, I_c(MINFTY, INFTY);
		if (s_delta > 0)
		{
			I = Interval(MINFTY, std::ceil(-Q) - 1);
			I_c = Interval(std::ceil(-Q), INFTY);
		}
		else if (s_delta < 0)
		{
			I = Interval(std::floor(-Q) + 1, INFTY);
			I_c = Interval(MINFTY, std::floor(-Q));
		}

		std::vector<Interval> J(2);
		for (int j = 0; j <= 1; j++)
		{
			if (a != 0)
			{
				int64_t gamma1 = d * (-R0 - (r0 + j) + a * m0);
				J[j] = Interval(-a * d, gamma1, N); J[j].Shift(-m0);
			}
			else
			{
				J[j] = Interval(N / (d * (R0 + r0 + j)) - m0 + 1, INFTY);
				//J[j] = Interval(MINFTY, N / (d * ((int64_t)std::floor(R0) + r0 + j)) - m0);
			}
		}

		// The intervals J[0] and J[1] are actually the complements of the sets J'[0] and J'[1] we want;
		// we try to sum over the set (I \cap J'[0]) \sqcup (I_c \cap J'[1]).
		// Note that (I \cap J'[0]) \sqcup (I_c \cap J'[1]) = R - (I \cap J[0]) - (I_c \cap J[1]).
		J[0].Intersect(I); //J[0].Shift(b);
		J[1].Intersect(I_c); //J[1].Shift(b);
		return SumInter(G, 0, Interval(MINFTY, INFTY), b, q)
			- SumInter(G, 0, J[0], b, q) - SumInter(G, 0, J[1], b, q);
	}

	// Computation of L1 - L2 for the case where a0(m-m0) + r0 = 0 mod q.
	// Algorithm by Helfgott & Thompson, 2023. Some comments are added for clarity.
	template <typename T>
	static T Special0A(std::vector<T> G, int64_t q, int64_t a, int64_t a_inv,
					   int64_t r0, int64_t b, double Q, int s_beta, int s_delta)
	{
		// Note that am + r0 = 0 mod q iff m = r mod q with r = -r0 * a_inv. 
		int r = -r0 * a_inv;

		Interval I;
		if (0 < r0 && r0 < q)
		{
			// The interval I should be taken s.t. m in I iff beta + delta m >= 0.
			// Recall that Q = beta / delta.
			if (s_delta > 0)
			{
				I = Interval(std::ceil(-Q), INFTY);
			}
			else if (s_delta < 0)
			{
				I = Interval(MINFTY, std::floor(-Q));
			}
			else if (s_beta >= 0) // delta = 0.
			{
				I = Interval(MINFTY, INFTY);
			}
		}
		else
		{
			// The interval I should be taken such that m in I if and only if
			// (beta < 0 AND delta m < 0) OR (beta delta m < 0 AND beta + delta m >= 0).
			// Recall that Q = beta / delta.
			if (s_beta < 0)
			{
				if (s_delta < 0)
				{
					return SumInter(G, r, Interval(MINFTY, std::floor(-Q)), b, q)
						+ SumInter(G, r, Interval(1, INFTY), b, q);
				}
				else if (s_delta > 0)
				{
					return SumInter(G, r, Interval(MINFTY, -1), b, q)
						+ SumInter(G, r, Interval(std::ceil(-Q), INFTY), b, q);
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
	template <typename Tf, typename Tg>
	static float_dec_100 SumByLinAppr(std::function<Tf(int64_t)> f, std::function<Tg(int64_t)> g, int64_t N, 
									  int64_t d0, int64_t m0, int64_t a, int64_t b)
	{
		Fraction alpha0 = Fraction(N, d0 * m0);
		Fraction alpha1 = Fraction(-N, d0 * d0 * m0);
		Fraction alpha2 = Fraction(-N, d0 * m0 * m0);

		/*
		std::function<int64_t(int64_t, int64_t)> L = [N, d0, m0](int64_t d, int64_t m)
			{
				return N / ((d + d0) * (m + m0));
			};
		std::function<int64_t(int64_t, int64_t)> L1 = [N, d0, m0](int64_t d, int64_t m)
			{
				return std::floor((double)(N * (d0 * m0 - d * m0 - m * d0)) / (d0 * d0 * m0 * m0));
			};
		std::function<int64_t(int64_t, int64_t)> L2 = [N, d0, m0](int64_t d, int64_t m)
			{
				return (N * (d0 * m0 - d * m0)) / (d0 * d0 * m0 * m0) + std::floor((double)(-(N * m)) / (d0 * m0 * m0));
			};*/

		float_dec_100 S = LinearSum(f, g, a, b, alpha0, alpha1, alpha2);

		std::tuple<int64_t, int64_t, int64_t, int> tuple = Elementary::DiophAppr::ApprByRedFrac(alpha2, 2 * b);
		int64_t a0 = std::get<0>(tuple);
		int64_t a0_inv = std::get<1>(tuple);
		int64_t q = std::get<2>(tuple); 
		int sgn_delta = std::get<3>(tuple);
		Fraction delta = alpha2 - Fraction(a0, q);

		Tg raySum = RaySum<Tg>(g, q, b, sgn_delta);

		std::tuple<std::vector<Tg>, std::vector<Tg>, std::vector<Tg>>
			tuple2 = SumTable(g, q, b, a0);
		std::vector<Tg> G = std::get<0>(tuple2);
		std::vector<Tg> rho = std::get<1>(tuple2);
		std::vector<Tg> sigma = std::get<2>(tuple2);

		for (int64_t d = -a; d < a; d++)
		{
			if (f(d) != 0) // I.e., f(d) is non-zero.
			{
				//int64_t R0_num = N * (d0 - d);
				//int64_t R0_denom = d0 * d0 * m0;
				//int64_t R0 = R0_num / R0_denom;
				Fraction R0 = alpha0 + alpha1 * Fraction(d);
				//double R0_frac = R0 - R0.Floor();
				Fraction R0_frac = R0.FractionalPart();

				//int64_t R0_frac_num = R0_num % R0_denom;
				int64_t r0 = (R0_frac * q + Fraction(1, 2)).Floor(); // (2 * q * R0_frac_num + R0_denom) / (2 * R0_denom);
				int64_t d_ = d0 + d;
				Fraction beta = R0_frac - Fraction(r0, q);//Fraction(R0_frac_num, R0_denom) - Fraction(r0, q);
				int sgn_beta = beta.Sign();

				double Q(1);
				if (!delta.IsZero())
				{
					Q = (beta / delta).numerical();
				}

				// Computation of difference g(m) * (L(d, m) - L1(d, m)) for m : a0(m - m0) + r0 = 0, -1 mod q.
				// (For other values of m, the difference is 0.)
				Tg T1(0);
				if (q > 1)
				{
					// Account for case a0(m - m0) + r0 = -1 mod q.
					T1 += Special1<Tg>(G, N, q, a0, a0_inv, R0.Floor(), r0, m0, d_, b);

					// Account for case a0(m - m0) + r0 = 0 mod q.
					T1 += Special0B<Tg>(G, N, q, a0, a0_inv, R0.Floor(), r0, m0, d_, b, Q, sgn_beta, sgn_delta);
				}
				else // Case q = 1.
				{
					T1 += Special00<Tg>(G, N, q, a0, a0_inv, R0.Floor(), r0, m0, d_, b, Q, sgn_delta);
				}

				// Computation of difference g(m) * (L1(d, m) - L2(d, m)).
				Tg T2(0);
				// Contribution of m s.t. a0(m - m0) + r0 != 0 mod q.
				T2 += sigma[r0];
				if (0 < r0 && r0 < q)
				{
					T2 += raySum;
				}

				// Contribution of m s.t. a0(m - m0) + r0 = 0 mod q.
				T2 += Special0A<Tg>(G, q, a0, a0_inv, r0, b, Q, sgn_beta, sgn_delta);

				/*
				Tg sum_LmL1(0), sum_L1mL2(0);
				for (int m = -b; m < b; m++)
				{
					sum_LmL1 += g(m) * (L(d, m) - L1(d, m));
					sum_L1mL2 += g(m) * (L1(d, m) - L2(d, m));
				}
				if (T1 - sum_LmL1 > 0.01 || T1 - sum_LmL1 < -0.01 || T2 - sum_L1mL2 > 0.01 || T2 - sum_L1mL2 < -0.01)
				{
					d--;
					if (T1 - sum_LmL1 > 0.01 || T1 - sum_LmL1 < -0.01)
					{
						//std::cout << "d: " << d0 + d << ", m: " << m0 - b << " -- " << m0 + b - 1 << std::endl;
						std::cout << "T1: " << T1 << ", sum_LmL1: " << sum_LmL1 << std::endl;
					}
					if (T2 - sum_L1mL2 > 0.01 || T2 - sum_L1mL2 < -0.01)
					{
						//std::cout << "d: " << d0 + d << ", m: " << m0 - b << " -- " << m0 + b - 1 << std::endl;
						std::cout << "T2: " << T2 << ", sum_L1mL2: " << sum_L1mL2 << std::endl;
					}
					continue;
				}*/
				S += (T1 + T2) * f(d);
			}
		}
		//std::cout << d0 - a << " " << d0 + a << " " << m0 - b << " " << m0 + b << std::endl;//": " << S << std::endl;
		return S;
	}

	// Computes the sum of f(d)g(m)floor(N/dm) over the rectangle [d0, d1) x [m0, m1) by brute force.
	// It is assumed that T1 and T2 can be cast to the type float_dec_100.
	static float_dec_100 BruteDoubleSum(int64_t d0, int64_t d1, int64_t m0, int64_t m1,
										std::vector<Log> f, std::vector<int> g, int64_t N)
	{
		float_dec_100 S(0);
		for (int64_t d = d0; d < d1; d++)
		{
			int64_t T(0);
			for (int64_t m = m0; m < m1; m++)
			{
				int64_t k = N / (d * m);
				T += k * g[m - m0];
			}
			S += T * f[d - d0].numerical();
		}
		return S;
	}

	// Computes the sum of f(d)g(m)floor(N/dm) over the rectangle [d0, d1) x [m0, m1) by
	// linear approximation of the term floor(N/dm) on small rectangles of fixed size a * b.
	template <typename Tf, typename Tg>
	static float_dec_100 DoubleSum(int64_t d0, int64_t d1, int64_t m0, int64_t m1, int64_t a, int64_t b,
								   std::function<Tf(int64_t)> f, std::function<Tg(int64_t)> g, int64_t N)
	{
		float_dec_100 S(0);
		int64_t mm = m0, dm = d0;
		for (int64_t dm = d0; dm < d1; dm += 2 * a)
		{
			int64_t dp = std::min(dm + 2 * a, d1);
			int64_t d_mid = (dm + dp) / 2, d_len = (dp - dm) / 2;
			std::function<Tf(int64_t)> f_ = [f, d0, dm, d_len](int64_t d)
				{
					return f(d + dm - d0 + d_len);
				};
			for (int64_t mm = m0; mm < m1; mm += 2 * b)
			{
				int64_t mp = std::min(mm + 2 * b, m1);
				int64_t m_mid = (mm + mp) / 2, m_len = (mp - mm) / 2;
				std::function<Tg(int64_t)> g_ = [g, m0, mm, m_len](int64_t m)
					{
						return g(m + mm - m0 + m_len);
					};
				S += SumByLinAppr(f_, g_, N, d_mid, m_mid, d_len, m_len);
			}
		}
		return S;
	}

	float_dec_100 DDSum(int64_t A0, int64_t A1, int64_t B0, int64_t B1, int64_t N, int64_t D, 
						COMP_MODE comp_mode, int64_t a, int64_t b, DIAG_MODE diag_mode)
	{
		float_dec_100 S(0);

		std::vector<int> mu_d, mu_m;
		std::function<int64_t(int64_t)> m_d, m_m;
		std::vector<Log> Lambda_d, Lambda_m;
		std::function<float_dec_100(int64_t)> L_d, L_m;
		for (int64_t d0 = A0; d0 < A1; d0 += D)
		{
			int64_t d1 = std::min(d0 + D, A1);
			Lambda_d = Elementary::sieve.LambdaSegmented(d0, d1 - d0);
			L_d = [Lambda_d](int64_t d)
				{
					return Lambda_d[d].numerical();
				};
			if (diag_mode == NON_DIAG)
			{
				mu_d = Elementary::sieve.MuSegmented(d0, d1 - d0);
				m_d = [mu_d](int64_t d)
					{
						return mu_d[d];
					};
			}

			for (int64_t m0 = B0; m0 < B1; m0 += D)
			{
				int64_t m1 = std::min(m0 + D, B1);
				mu_m = Elementary::sieve.MuSegmented(m0, m1 - m0);
				m_m = [mu_m](int64_t m)
					{
						return mu_m[m];
					};
				if (diag_mode == NON_DIAG)
				{
					Lambda_m = Elementary::sieve.LambdaSegmented(m0, m1 - m0);
					L_m = [Lambda_m](int64_t m)
						{
							return Lambda_m[m].numerical();
						};
				}

				switch (comp_mode)
				{
					case LINEAR_APPR:
						//std::cout << "Lin-appr: " << d0 << " " << d1 << " " << m0 << " " << m1 << std::endl;
						S += DoubleSum(d0, d1, m0, m1, a, b, L_d, m_m, N);
						if (diag_mode == NON_DIAG)
						{
							//std::cout << "Lin-appr: " << m0 << " " << m1 << " " << d0 << " " << d1 << std::endl;
							S += DoubleSum<int64_t, float_dec_100>(d0, d1, m0, m1, a, b, m_d, L_m, N);
						}
						break;
					case BRUTE_FORCE:
						//std::cout << "Brute-force: " << d0 << " " << d1 << " " << m0 << " " << m1 << std::endl;
						S += BruteDoubleSum(d0, d1, m0, m1, Lambda_d, mu_m, N);
						if (diag_mode == NON_DIAG)
						{
							//std::cout << "Brute-force: " << m0 << " " << m1 << " " << d0 << " " << d1 << std::endl;
							S += BruteDoubleSum(m0, m1, d0, d1, Lambda_m, mu_d, N);
						}
						break;
				}
			}
		}
		return S;
	}

	// Returns sum_{mdk <= N, m,d <= M0} Lambda(m)mu(d) = sum_{m,d <= M0} Lambda(m)mu(d)floor(N/md).
	// Based on algorithm by Helfgott & Thompson, 2023.
	float_dec_100 PsiElem::IndependentVar(int64_t N, int64_t M0)
	{
		float_dec_100 S(0);
		int64_t A1 = M0 + 1, B1 = M0 + 1;
		int64_t C = 10, D = 8; // Hand-tuned by Helfgott & Thompson, 2023.

		while (A1 >= 2 * std::pow(6 * C * C * C * N, 0.25) &&
			   A1 >= std::sqrt(M0) + 1 &&
			   A1 >= 2 * D)
		{
			int64_t A = A1 - 2 * (A1 / (2 * D));
			double cbrt = std::cbrt((double)A / (6 * N));
			int64_t a = A * cbrt;
			while (B1 >= 2 * C * std::cbrt(6 * N / A) &&
				   B1 >= std::sqrt(M0) + 1 &&
				   B1 >= 2 * D)
			{
				int64_t B = B1 - 2 * (B1 / (2 * D));
				int64_t b = B * cbrt;

				int64_t D = (1 + M0 / std::max(2 * a, 2 * b)) * std::max(2 * a, 2 * b);

				DIAG_MODE diag_mode = (A == B && A1 == B1) ? DIAGONAL : NON_DIAG;
				//std::cout << "REGION: " << A << " " << A1 << ", " << B << " " << B1 << std::endl;
				S += DDSum(A, A1, B, B1, N, D, LINEAR_APPR, a, b, diag_mode);

				B1 = B;
			}


			// Remaining part is done by brute force.
			//std::cout << "REGION: " << A << " " << A1 << ", " << "1" << " " << B1 << std::endl;
			S += DDSum(A, A1, 1, B1, N, std::sqrt(M0) + 1, BRUTE_FORCE, 0, 0, NON_DIAG);
			A1 = A;
			B1 = A;
		}

		// The remaining rectangle is done again by brute force.
		//std::cout << "REGION: " << "1" << " " << A1 << ", " << "1" << " " << B1 << std::endl;
		S += DDSum(1, A1, 1, B1, N, std::sqrt(M0) + 1, BRUTE_FORCE, 0, 0, DIAGONAL);

		//std::cout << "IV: " << S << std::endl;
		return S;
	}
}