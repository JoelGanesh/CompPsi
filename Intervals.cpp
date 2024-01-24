/*#include "Types.h"

#include <vector>
#include <map>
#include <algorithm>

namespace Types
{
	// Stores multiple integer intervals.
	// Allows operations on intervals, such as union/intersection with an interval.
	class Intervals
	{
	private:
		std::map<uint64_t, uint64_t> ranges;

		void RemoveKeys(std::vector<int> keys)
		{
			// We need to take care of the fact that removing elements in the array changes 
			// indices of other elements. We do this by ordering the elements from high to low.
			std::sort(keys.begin(), keys.end(), std::greater<>());
			for (uint64_t key : keys)
			{
				ranges.erase(key);
			}
			return;
		}
	public:
		std::vector<Interval> getIntervals()
		{
			std::vector<Interval> intervals;
			for (const auto& range : ranges)
			{
				intervals.push_back(IntervalLR(range.first, range.second));
			}
			return intervals;
		}

		// Constructor from one interval
		Intervals(Interval I)
		{
			ranges[I.start] = I.end;
		}

		void Union(Interval I)
		{
			auto afterRange = ranges.upper_bound(I.start);
			std::vector<uint64_t> removeIndices;
			for (uint64_t n = 0; n < intervals.size(); n++)
			{
				Interval J = intervals[n];

				if (J.Contains(I))
				{
					// I is already part of 'intervals'.
					return;
				}
				if (!J.Empty())
				{
					if (I.Contains(J))
					{
						removeIndices.push_back(n);
						continue;
					}
					if (I.Contains(J.start))
					{
						// I extends J from the left (and not from the right)
						// We merge I and J by removing J and changing I.end accordingly.
						I.end = J.end;
						removeIndices.push_back(n);
					}
					if (I.Contains(J.end - 1))
					{
						// I extends J from the right (and not from the left)
						// We merge I and J by removing J and changing I.start accordingly.
						I.start = J.start;
						removeIndices.push_back(n);
					}
				}
			}

			RemoveIntervals(removeIndices);
			intervals.push_back(I);
			return;
		}

		void Intersection(Interval I)
		{
			// 'newIntervals' replaces 'intervals' at the end.
			std::vector<Interval> newIntervals;
			for (
			
			
			
			n = 0; n < intervals.size(); n++)
			{
				Interval J = intervals[n];
				if (J.Contains(I))
				{
					// I is contained in J, so we can just shrink intervals to I.
					newIntervals = { I };
					break;
				}
				if (I.Contains(J) && !J.Empty())
				{
					newIntervals.push_back(J);
					continue;
				}
				if (I.Contains(J.start))
				{
					// I extends J from the left (but not from the right).
					Interval newInterval = IntervalLR(J.start, I.end);
					newIntervals.push_back(newInterval);
				}
				if (I.Contains(J.end - 1))
				{
					// I extends J from the right (and not from the left).
					Interval newInterval = IntervalLR(I.start, J.end);
					newIntervals.push_back(newInterval);
				}
			}

			intervals = newIntervals;
			return;
		}

		void SetMinus(Interval I)
		{

		}

		void Complement()
		{

		}
	};

	// Structure to store (non-negative) integer intervals of the form [a, b)
	// Provides multiple possibilities to create instances via different subclasses.
	class Interval
	{
	public:
		uint64_t start;
		uint64_t end;

		Interval() {}

		bool Contains(Interval I)
		{
			if (I.Empty())
			{
				return true;
			}
			return (start <= I.start && I.end <= end);
		}

		bool Contains(uint64_t a)
		{
			return (start <= a && a < end);
		}
		
		bool Empty()
		{
			return (start == end);
		}

		bool operator<(Interval I)
		{
			return end < I.start;
		}
	};

	class IntervalMid : public Interval
	{
	public:
		// Constructor: interval [m - r, m + r).
		IntervalMid(uint64_t m, uint64_t r)
		{
			start = m - r;
			end = m + r;
		};
	};

	class IntervalLR : public Interval
	{
	public:
		IntervalLR(uint64_t a, uint64_t b)
		{
			start = std::min(a, b);
			end = std::max(a, b);
		}

	};

	class IntervalIneq : public Interval
	{
	public:
		IntervalIneq(intmax_t a, intmax_t b)
		{

		}

		IntervalIneq(intmax_t a, intmax_t b, intmax_t c)
		{

		}

	};
}*/