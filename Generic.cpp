/*#include "Utility.h"

namespace Utility
{
	template <typename X, typename Y>
	std::vector<Y> Generic::Map(std::vector<X> v, std::function<Y(X)> f)

	template <typename T>
	T Generic::Sum(std::function<T(int64_t)> f, Interval I)
	{
		T sum = 0;
		for (int64_t n = I.start; n <= I.end; n++)
		{
			sum += f(n);
		}
		return sum;
	}

	template <typename T>
	std::vector<int64_t> Generic::IndexAll(std::vector<T> v, T value, int64_t offset)
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

	template <typename S, typename T>
	T Generic::Merge(std::vector<S> v, T sep, std::function<T(T, T)> merge)
	{
		if (v.empty())
		{
			throw std::invalid_argument("The vector 'v' should not be empty.");
		}

		T res(v[0]);
		for (int64_t i = 1; i < v.size(); i++)
		{
			res = merge(res, sep);
			res = merge(res, T(v[i]));
		}

		return res;
	}

	template <typename T>
	int64_t Generic::FindIndex(std::vector<T> v, T value)
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

	template <typename T>
	std::vector<T> Generic::Copy(T* data, int64_t size)
	{
		std::vector<T> v(size);
		for (int i = 0; i < size; i++)
		{
			v[i] = *(data + i);
		}
		return v;
	}

	template <typename T>
	void Generic::Assign(std::vector<T> data, T* ptr)
	{
		for (int64_t i = 0; i < data.size(); i++)
		{
			ptr[i] = data[i];
		}
		return;
	}

	template <typename T>
	void Generic::PadZeros(std::vector<T>& data, int n)
	{
		data.insert(data.end, n, T(0));
	}
}*/
