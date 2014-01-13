#ifndef __MPI_KMEANS_MEX_H__
#define __MPI_KMEANS_MEX_H__
#include<float.h>
#include<vector>
#include "util.hpp"
using namespace std;
#ifndef KMEANS_VERBOSE 
#define KMEANS_VERBOSE 0
#endif

/* Double precision is default*/
#ifndef INPUT_TYPE
#define INPUT_TYPE 1
#endif

#if INPUT_TYPE==0
#define PREC double
#define PREC_MAX DBL_MAX
#elif INPUT_TYPE==1
#define PREC float
#define PREC_MAX FLT_MAX
#endif

#ifndef BOUND_PREC
#define BOUND_PREC float
#endif

#ifndef BOUND_PREC_MAX
#define BOUND_PREC_MAX FLT_MAX
#endif

#define BOUND_EPS 1e-6

typedef unsigned long uInt;

// Modifications for the codes kmeans where there is extra entry to have the counts of that code...
extern "C" {
PREC kmeans(PREC *CXp, const PREC *X, uInt *c, uInt dim, uInt npts, uInt nclus,
		uInt maxiter, uInt nr_restarts, uInt *pcounts);
}
//PREC compute_distance(const PREC *vec1, const PREC *vec2, const uInt dim);
//uInt assign_point_to_cluster_ordinary(const PREC *px, const PREC *CX, uInt dim,uInt nclus);
//uInt assign_point_to_cluster_ordinary(const int *px, const PREC *CX, uInt dim,uInt nclus);
void randperm(uInt *order, uInt npoints);

template<class T>
PREC compute_distance(const T *vec1, const PREC *vec2, const uInt dim) {
	PREC d = 0.0;
	for (uInt k = 0; k < dim; k++) {
		PREC df = (vec1[k] - vec2[k]);
		d += df * df;
	}
	assert(d >= 0.0);
	d = sqrt(d);
	return d;
}
template<class T>
uInt assign_point_to_cluster_ordinary(const T *px, const PREC *CX, uInt dim,
		uInt nclus) {
	uInt assignment = nclus;
	PREC mind = PREC_MAX;
	const PREC *pcx = CX;
	for (uInt j = 0; j < nclus; j++, pcx += dim) {
		PREC d = compute_distance(px, pcx, dim);
		if (d < mind) {
			mind = d;
			assignment = j;
		}
	}
	assert(assignment < nclus);
	return (assignment);
}
template<class T>
// computes the mean of distances and map top
// clusters measured by percentage of mean...
void assign_point_to_cluster_soft(const T *px, const PREC *CX, uInt dim,
		uInt nclus, float *map, double mratio = 1) {
	fill(map, map + nclus, 0);
	uInt assignment = nclus;
	PREC mind = PREC_MAX;
	const PREC *pcx = CX;
	double mean = 0, mini = 1000000, midx = -1;
	for (uInt j = 0; j < nclus; j++, pcx += dim) {
		map[j] = compute_distance(px, pcx, dim); // store cluster distance;
		mean += map[j]; // mean of distance...
		if (map[j] < mini) {
			mini = map[j];
			midx = j;
		}
	}
	mean /= nclus;
	if (mratio < 0) { // sort and select top with difference to mean distance...
		mratio *= -1;
		typedef float* T4;
		vector<SorterElement<T4> > res4 = sort_elements(map, map + nclus,
				CmpAbsDescend<REAL> ());
		uInt count = 0;
		for (vector<SorterElement<T4> >::reverse_iterator i = res4.rbegin()
				+ mratio; i != res4.rend(); ++i, ++count) {
			map[i->_ind] = 0;
			//			spweights[count] = *i->_e;
		}
		count = 0;
		double diff, nsum = 0;
		for (vector<SorterElement<T4> >::reverse_iterator i = res4.rbegin(); i
				!= res4.rend() && count < mratio; ++i, ++count) {
			diff = mean - map[i->_ind];
			map[i->_ind] = diff > 0 ? diff : 0;
			nsum += map[i->_ind];
			//			spweights[count] = *i->_e;
		}
		count=0;
		// normalize to have them as probabilties
		for (vector<SorterElement<T4> >::reverse_iterator i = res4.rbegin(); i
				!= res4.rend() && count < mratio; ++i, ++count) {
			map[i->_ind] /= nsum;
			//			nsum += map[i->_ind];
			//			spweights[count] = *i->_e;
		}

		//		assert(count == nselements);
	} else {
		mean /= mratio;
		mean = max(mini + 1e-4, mean);
		//	mean /= mratio; // to select top features ...
		float diff;
		for (uInt j = 0; j < nclus; j++, pcx += dim) {
			diff = (mean - map[j]);
			map[j] = diff > 0 ? diff : 0;
		}
	}
}
#endif

