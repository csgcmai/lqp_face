/*! \file
 * \brief Sorting algo.
 *
 * sort_elements(x.begin(), x.end()) <=> [sortedValues,sortIndex]=sort(x) of octave and matlab.
 * see Sorter-example.C for simple examples.
 * 
 * A container/carray of any type/class can be sorted efficiently, 
 * only pointers are copied/assigned.
 *
 * @author Eric Nowak (http://lear.inrialpes.fr/people/nowak)
 */

#ifndef __PASS_RAZSORTER
#define __PASS_RAZSORTER

#include <vector>
#include <algorithm>

/// Used in output of sort_elements.
///	- contains ID and { pointer | iterator } on value
/// - T <=> iterators on object or pointers on objects
template<class T> struct SorterElement {
	SorterElement(int ind, T e) :
		_ind(ind), _e(e) {
	}
	;
	int _ind;
	T _e;
};

/// Compare SorterElement <=> compare their elements.
/// - T <=> iterators on object or pointers on objects
/// - BinPredCmp_Element => binary predicate for 2 *T
template<class T, class BinPredCmp_Element> struct BinPredCmp_SorterElement {
	bool operator()(const SorterElement<T>& s1, const SorterElement<T>& s2) {
		BinPredCmp_Element p;
		return p(*s1._e, *s2._e);
	}
};

/// Sort - output contains info about value and rank in sequence begin, end.
/// Sort_elements(x.begin(), x.end()) <=> [sortedValues,sortIndex]=sort(x) of octave and matlab..
/// - T <=> iterators on object or pointers on objects
/// - CMP => binary predicate for 2 T*
template<class T, class CMP>
std::vector<SorterElement<T> > sort_elements(T begin, T end, CMP cmp) {
	std::vector<SorterElement<T> > res;
	BinPredCmp_SorterElement<T, CMP> cmpSE;
	for (int id = 0; begin != end; ++begin, ++id)
		res.push_back(SorterElement<T> (id, begin));
	sort(res.begin(), res.end(), cmpSE);
	return res;
}
;

#endif // __PASS_RAZSORTER
