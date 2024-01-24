#include "Utility.h"

namespace Utility
{
	std::vector<std::vector<int>> Indexation::RestrictedTuples(int n, std::vector<int> a)
	{
		return RestrictedTuples(n, a, a.size() - 1);
	}

	std::vector<std::vector<int>> Indexation::RestrictedTuples(int n, std::vector<int> a, uint64_t index)
	{
		std::vector<std::vector<int>> tuples;
		if (n >= 0 && index >= 0)
		{
			// Base case: one element.
			// In this case there can only be one tuple, namely { n }.
			if (index == 0)
			{
				if (n <= a[index])
				{
					std::vector<int> tuple{ n };
					tuples.push_back(tuple);
				}
				return tuples;
			}
			else // General case.
			{
				for (int k = 0; k <= std::min(a[index], n); k++)
				{
					std::vector<std::vector<int>> prevTuples = RestrictedTuples(n - k, a, index - 1);
					for (std::vector<int> tuple : prevTuples)
					{
						tuple.push_back(k);
						tuples.push_back(tuple);
					}
				}
			}
		}
		return tuples;
	}

	std::vector<int> Indexation::IntToTuple(int n, std::vector<int> a)
	{
		std::vector<int> v;
		for (int k = 0; k < a.size(); k++)
		{
			v.push_back(n % a[k]);
			n /= a[k];
		}
		return v;
	}
}