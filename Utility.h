// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#ifndef UTILITY_H
#define UTILITY_H

#include "Types.h"

namespace Utility
{
	class Generic
	{
		public:
		// Apply function f: X -> Y to each entry of v (which are of type X).
		template <typename X, typename Y>
		static std::vector<Y> Map(std::vector<X> v, std::function<Y(X)> f)
		{
			std::vector<Y> res;
			for (const X x : v)
			{
				res.push_back(f(x));
			}
			return res;
		}

		// List all indices with the specified value.
		// Also allows an offset parameter, used to replace an index i by i + offset.
		// T should support comparison (equality) of elements.
		template <typename T>
		static std::vector<int64_t> IndexAll(std::vector<T> v, T value, int64_t offset = 0)
		{
			std::vector<int64_t> res;
			for (int64_t i = 0; i < v.size(); i++)
			{
				if (v[i] == value)
				{
					res.push_back(i + offset);
				}
			}
			return res;
		}

		// Finds the last index i such that v[i] <= value, assuming it exists,
		// using Binary Search. If it does not exist, returns 0.
		// Assumes that v[i] <= v[j] for i <= j, with '<' defined on T.
		template <typename T>
		static int64_t FindIndex(std::vector<T> v, T value)
		{
			int64_t min_index = 0;
			int64_t max_index = v.size() - 1;
			while (min_index < max_index)
			{
				int64_t mean = (min_index + max_index) / 2;
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
		static std::vector<T> Copy(T* data, int64_t size)
		{
			std::vector<T> v(data, data + size);
			return v;
		}

		// Fills the data a pointer points to with the specified data.
		// It is assumed that the pointer has been allocated enough space.
		template <typename T>
		static void Assign(std::vector<T> data, T* ptr)
		{
			for (int64_t i = 0; i < data.size(); i++)
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
	
	class Indexation
	{
		public:
		// Generates all non-negative integer tuples (k_1, ..., k_d) such that k_i <= a_i,
		// while k_1 + ... + k_d = n. Here 'd' represents the size of 'a'.
		static std::vector<std::vector<int>> RestrictedTuples(int n, std::vector<int> a);

		// Converts integers to tuples with non-negative integer entries restricted by the array 'a'.
		// Restricted to any integer interval of length P = prod_i a_i,
		// it is a bijection to {(n_1, ..., n_k) : 0 <= n_i < a_i for all i}.
		static std::vector<int> IntToTuple(int n, std::vector<int> a);

		private:
		// Recursive function computing tuples (k_j, ..., k_d) as in 'RestrictedTuples' above, where j = index + 1.
		static std::vector<std::vector<int>> RestrictedTuples(int n, std::vector<int> a, int64_t index);
	};
}
#endif