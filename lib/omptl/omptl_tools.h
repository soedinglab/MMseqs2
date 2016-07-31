// Copyright (C) 2006-2011 Fokko Beekhof
// Email contact: Fokko.Beekhof@unige.ch

// The OMPTL library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef OMPTL_TOOLS_H
#define OMPTL_TOOLS_H 1

#include <utility>
#include <vector>
#include <cassert>
#include <algorithm>
#include <climits>
#include <iterator>

#include <tr1/cmath>

namespace omptl
{

namespace detail
{


// Log of the number of operations that is expected to run faster in a single
// thread.
const unsigned C = 12;

template <typename T>
T log2N_(T n)
{
	assert(n > 0);
	const std::size_t N = CHAR_BIT*sizeof(T);

	T result = 0;
	for (std::size_t i = 1; i < N; ++i)
	{
		const std::size_t M = N-i;
		if ( n >= (std::size_t(1) << M) )
		{
			n >>= M;
			result |= M;
		}
	}

	return result;
}

template<typename Iterator>
bool _linear_serial_is_faster(Iterator first, Iterator last,
			     const unsigned P)
{
	assert(P > 0u);
	assert(std::distance(first, last) >= 0);
	const std::size_t N = std::distance(first, last);

	return (N < 4u*P) || (log2N_(N) < C);
}

template<typename Iterator>
bool _logn_serial_is_faster(Iterator first, Iterator last,
			    const unsigned P)
{
	assert(P > 0u);
	assert(std::distance(first, last) >= 0);
	const std::size_t N = std::distance(first, last);

	return (N < 4u*P) || (log2N_(N) < (1 << C));
}

template<typename Iterator>
bool _nlogn_serial_is_faster(Iterator first, Iterator last,
			    const unsigned P)
{
	assert(P > 0u);
	assert(std::distance(first, last) >= 0);
	const std::size_t N = std::distance(first, last);

	return (N < 4u*P) || (N*log2N_(N) < (1 << C));
}

template<typename Iterator1, typename Iterator2>
void _copy_partitions(const std::vector< std::pair<Iterator1, Iterator1> >
			&source_partitions, Iterator2 first,
		std::vector<Iterator2> &dest_partitions, const unsigned P)
{
	assert(source_partitions.size() == P);
	assert(dest_partitions.size() == P);
	for (unsigned i = 0; i < P; ++i)
	{
		dest_partitions[i] = first;

		// The last "advance" is very important, it may create space
		// if it is an InsertIterator or something like that.
		std::advance(first, std::distance(
						source_partitions[i].first,
						source_partitions[i].second) );
	}
}

// Divide a given range into P partitions
template<typename Iterator>
void _partition_range(const Iterator first, const Iterator last,
		std::vector< std::pair<Iterator, Iterator> > &partitions,
		const unsigned P)
{
	assert(partitions.size() == P);

	typedef std::pair<Iterator, Iterator> Partition;

	const std::size_t N = std::distance(first, last);

	// All but last partition have same range
	Iterator currentLast = first;
	for (unsigned i = 0; i < P - 1; ++i)
	{
		const Iterator prev = currentLast;
		currentLast = first;
		std::advance(currentLast, (i+1)*N/P);
		partitions[i] = Partition(prev, currentLast);
	}
	assert(std::distance(currentLast, last) >= 0);

	// Last range may be shorter
	partitions[P - 1] = Partition(currentLast, last);
}

// Given a range, re-arrange the items such that all elements smaller than
// the pivot precede all other values. Returns an Iterator to the first
// element not smaller than the pivot.
template<typename Iterator, class StrictWeakOrdering>
Iterator _stable_pivot_range(Iterator first, Iterator last,
	const typename std::iterator_traits<Iterator>::value_type pivot,
	StrictWeakOrdering comp = std::less<
		typename std::iterator_traits<Iterator>::value_type>())
{
	Iterator pivotIt = last;
	while (first < last)
	{
		if (comp(*first, pivot))
			++first;
		else
		{
			Iterator high = first;
			while ( (++high < last) && !comp(*high, pivot) )
				/* nop */;
			if (high < last)
				std::iter_swap(first, last);
			first = pivotIt = ++high;
		}
	}

	return pivotIt;
}

template<typename Iterator>
void _partition_range_stable_by_pivots(Iterator first, Iterator last,
	const std::vector<typename
			std::iterator_traits<Iterator>::value_type> &pivots,
	std::vector< std::pair<Iterator, Iterator> > &partitions,
	std::less<typename std::iterator_traits<Iterator>::value_type> comp,
	const unsigned P)
{
	assert(partitions.size() == P);
	assert(pivots.size() == P);
	typedef std::pair<Iterator, Iterator> Partition;

	Iterator start = first;
	for (unsigned i = 0; i < P - 1; ++i)
	{
		Iterator low = start;

		while (low < last)
		{
			// Find a value not lower than the pivot.
			while( (*low < pivots[i]) && (low < last) )
				std::advance(low, 1);

			// Entire range scanned ?
			if (low == last) break;

			// Find a value lower than the pivot, starting from
			// low, working our way up.
			Iterator high = low;
			std::advance(high, 1);
			while( !(*high < pivots[i]) && (high < last) )
				std::advance(high, 1);

			// Entire range scanned ?
			if (high == last) break;

			// Swap values
			assert( !(*low<pivots[i]) && (*high<pivots[i]) );
			std::iter_swap(low, high);
		}

		partitions[i] = Partition(start, low);
		start = low;
	}
	partitions[P - 1] = Partition(start, last);
}

template<typename RandomAccessIterator, class StrictWeakOrdering>
void _find_pivots(RandomAccessIterator first, RandomAccessIterator last,
	std::vector<typename
	std::iterator_traits<RandomAccessIterator>::value_type> &pivots,
	StrictWeakOrdering comp, const unsigned P)
{
	const std::size_t N = std::distance(first, last);

	assert(N > P);

	pivots.clear();
	pivots.reserve(P - 1);

	typedef typename
	    std::iterator_traits<RandomAccessIterator>::value_type value_type;

	/*
	 * The sample ratio of 3 is used to sample more data. This way, the pivots can be
	 * chosen more wisely, which is our only guarantee we can generate partitions
	 * of equal size.
	 */
	const std::size_t NSAMPLES = std::min( 3u*std::size_t(P), N);
	std::vector<value_type> samples;
	samples.reserve(NSAMPLES);

	for (std::size_t i = 0; i < NSAMPLES; ++i)
	{
		const std::size_t index = i * (N-1) / (NSAMPLES - 1);
		assert(index < N);
		samples.push_back(*(first + index));
// std::cout << "index: " << index << " sample: " << samples[i] << std::endl;
	}
	assert(samples.size() == NSAMPLES);

	// Sort samples to create relative ordering in data
	std::sort(samples.begin(), samples.end(), comp );

	// Take pivots from sampled data
	for (std::size_t i = 0; i < P-1; ++i)
	{
		pivots.push_back(samples[std::min(1+3*i, N-1)]);
/*std::cout << "pivot: " << i << " idx: " << (i * samples.size() / P)
	<< " " << pivots[i-1] << std::endl;*/
	}
	assert(pivots.size() == P - 1);
}

}  // namespace detail

}  // namespace omptl

#endif /* OMPTL_TOOLS_H */
