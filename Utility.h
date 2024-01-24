#ifndef UTILITY_H
#define UTILITY_H

#include <vector>
#include <iostream>
#include <stack>
#include <functional>
#include <string>

#include "Types.h"
using namespace Types;

namespace Utility
{
	// Class for generic methods
	class Generic
	{
	public:
		// Apply function f: X -> Y to each of entries of v (which are of type X)
		template <typename X, typename Y>
		static std::vector<Y> Map(std::vector<X> v, std::function<Y(X)> f)
		{
			std::vector<Y> res;
			for (X x : v)
			{
				res.push_back(f(x));
			}
			return res;
		}

		// List all indices with the specified value.
		// Also allows an offset parameter, used to replace an index i by i + offset.
		// T should support comparison (equality) of elements.
		template <typename T>
		static std::vector<uint64_t> IndexAll(std::vector<T> v, T value, uint64_t offset = 0)
		{
			std::vector<uint64_t> res;
			for (uint64_t i = 0; i < v.size(); i++)
			{
				if (v[i] == value)
				{
					res.push_back(i + offset);
				}
			}
			return res;
		}

		// Compute the sum of values of some function f over integers in some interval I.
		// T should support addition.
		//template <typename T>
		//static T Sum(std::function<T(uint64_t)> f, Interval I);

		// Merges objects of a certain type together with a separator.
		// It is assumed that an implict conversion from type S to type T has been implemented.
		//template <typename S, typename T>
		//static T Merge(std::vector<S> v, T sep, std::function<T(T, T)> merge);

		// Finds the last index i such that v[i] <= value, assuming it exists,
		// using a Binary Search approach. If it does not exist, returns 0.
		// Assumes that v[i] <= v[j] for i <= j, with '<' defined on T.
		template <typename T>
		static uint64_t FindIndex(std::vector<T> v, T value)
		{
			uint64_t min_index = 0;
			uint64_t max_index = v.size() - 1;
			while (min_index < max_index)
			{
				uint64_t mean = (min_index + max_index) / 2;
				if (v[mean] < value)
				{
					min_index = mean + 1;
				}
				else
				{
					max_index = mean;
				}
			}
			return max_index;
		}

		// Copies the data referenced by a pointer to a vector.
		// It is assumed that 'size' equals the size of 'data'.
		template <typename T>
		static std::vector<T> Copy(T* data, uint64_t size)
		{
			std::vector<T> v(data, data + size);
			return v;
		}

		// Fills the data a pointer points to with the specified data.
		// It is assumed that the pointer has been allocated enough space.
		template <typename T>
		static void Assign(std::vector<T> data, T* ptr)
		{
			for (uint64_t i = 0; i < data.size(); i++)
			{
				ptr[i] = data[i];
			}
			return;
		}

		// Pads 'n' zeros to the end of the vector 'data'.
		// Assumes that T(0) is defined by T.
		template <typename T>
		static void PadZeros(std::vector<T>& data, int n)
		{
			data.insert(data.end, n, T(0));
		}
	};

	// Class specifically for finding Diophantine approximations.
	class DiophAppr
	{
	public:
		// Returns a tuple (p, p', q, s) of integers so that abs(alpha - p/q) <= 1/(qQ)
		// with gcd(p, q) = 1, q <= Q, and pp' = 1 mod q, while s = sign(alpha - p/q).
		// Algorithm taken from a published paper from 2023 by H. A. Helfgott and L. Thompson.
		static std::tuple<int, int, int, int> ApprByRedFrac(double alpha, int Q);
	};

	// Class helping with indexing certain objects.
	class Indexation
	{
	public:
		// Generates all non-negative integer tuples (k_1, ..., k_d) such that k_i <= a_i,
		// while k_1 + ... + k_d = n. Here 'd' represents the size of 'a'.
		// Note: due to how 'a' is processed, it might be beneficial to order elements in 'a' from small to large. TODO: Figure out if this is the case.
		static std::vector<std::vector<int>> RestrictedTuples(int n, std::vector<int> a);

		// Converts integers to tuples with non-negative integer entries restricted by the array 'a'.
		// Restricted to any integer interval of length P = prod_i a_i,
		// it is a bijection to {(n_1, ..., n_k) : 0 <= n_i < a_i for all i}.
		static std::vector<int> IntToTuple(int n, std::vector<int> a);

	private:
		// Recursive function computing tuples (k_j, ..., k_d) as in 'RestrictedTuples' above, where j = index + 1.
		static std::vector<std::vector<int>> RestrictedTuples(int n, std::vector<int> a, uint64_t index);
	};

	// Class for functions related to input and output.
	class IO
	{
	public:
		// Shortcut for printing out arrays nicely.
		// T should support conversion to string.
		template <typename T>
		static void Print(std::vector<T> v, std::string sep = " ")
		{
			std::cout << v[0];
			for (uint64_t n = 1; n < v.size(); n++)
			{
				std::cout << sep << v[n];
			}
			std::cout << std::endl;
		}
	};


	/*class Interval
	{
	public:
		// Create interval_set I so that x >= 1 lies in I iff ax^2 + bx + c >= 0.
		static interval_set FromQuadIneq(int a, int b, int c);

		// Create interval_set I such that x >= 1 in I iff ax + b >= 0.
		static interval_set FromLinIneq(int a, int b);
	};*/
};
#endif