#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <memory.h>
#include <math.h>
#include <assert.h>
#include "mpi_kmeans.h"

#if KMEANS_VERBOSE>1
uInt saved_two=0,saved_three_one=0,saved_three_two=0,saved_three_three=0,saved_three_b=0;
#endif

void kmeans_error(char *msg) {
	printf("%s", msg);
	exit(-1);
}

int comp_randperm(const void * a, const void * b) {
	return ((int) (*(double*) a - *(double*) b));
}

void randperm(uInt *order, uInt npoints) {
	double *r = (double*) malloc(2 * npoints * sizeof(double));
	for (uInt i = 0; i < 2 * npoints; i++, i++) {
		r[i] = rand();
		r[i + 1] = i / 2;
	}
	qsort(r, npoints, 2 * sizeof(double), comp_randperm);

	for (uInt i = 1; i < 2 * npoints; i++, i++)
		order[i / 2] = (uInt) r[i];

	free(r);
}
/* Remove the template function and use it instead
PREC compute_distance(const PREC *vec1, const PREC *vec2, const uInt dim) {
	PREC d = 0.0;
	for (uInt k = 0; k < dim; k++) {
		PREC df = (vec1[k] - vec2[k]);
		d += df * df;
	}
	assert(d>=0.0);
	d = sqrt(d);
	return d;
}
PREC compute_distance(const int *vec1, const PREC *vec2, const uInt dim) {
	PREC d = 0.0;
	for (uInt k = 0; k < dim; k++) {
		PREC df = (vec1[k] - vec2[k]);
		d += df * df;
	}
	assert(d>=0.0);
	d = sqrt(d);

	return d;
}
*/
PREC compute_sserror(const PREC *CX, const PREC *X, const uInt *c, uInt dim,
		uInt npts, uInt *ptcount) {
	PREC sse = 0.0;
	const PREC *px = X;
	const uInt *ptrcount = ptcount;
	for (uInt i = 0; i < npts; i++, px += dim, ++ptrcount) {
		const PREC *pcx = CX + c[i] * dim;
		PREC d = compute_distance(px, pcx, dim);
		sse += d * d * *ptrcount;
	}
	assert(sse>=0.0);
	return (sse);
}

void remove_point_from_cluster(uInt cluster_ind, PREC *CX, const PREC *px,
		uInt *nr_points, uInt dim, uInt ptcount) {
	PREC *pcx = CX + cluster_ind * dim;

	/* empty cluster after or before removal */
	if (nr_points[cluster_ind] < 2 || nr_points[cluster_ind] <= ptcount) {
		for (uInt k = 0; k < dim; k++)
			pcx[k] = 0.0;
		nr_points[cluster_ind] = 0;
	} else {
		/* pgehler: remove PREC here */
		PREC nr_old, nr_new;
		nr_old = (PREC) nr_points[cluster_ind];
		nr_points[cluster_ind] -= ptcount;
		nr_new = (PREC) nr_points[cluster_ind];

		for (uInt k = 0; k < dim; k++)
			pcx[k] = (nr_old * pcx[k] - ptcount * px[k]) / nr_new;
	}
}

void add_point_to_cluster(uInt cluster_ind, PREC *CX, const PREC *px,
		uInt *nr_points, uInt dim, uInt ptcount) {

	PREC *pcx = CX + cluster_ind * dim;

	/* first point in cluster */
	if (nr_points[cluster_ind] == 0) {
		nr_points[cluster_ind] += ptcount;
		for (uInt k = 0; k < dim; k++)
			pcx[k] = px[k];
	} else {
		/* remove PREC here */
		PREC nr_old = (PREC) (nr_points[cluster_ind]);
		nr_points[cluster_ind] += ptcount;
		PREC nr_new = (PREC) (nr_points[cluster_ind]);
		for (uInt k = 0; k < dim; k++)
			pcx[k] = (nr_old * pcx[k] + ptcount * px[k]) / nr_new;
	}
}

bool remove_identical_clusters(PREC *CX, BOUND_PREC *cluster_distance,
		const PREC *X, uInt *cluster_count, uInt *c, uInt dim, uInt nclus,
		uInt npts, uInt *ptcount) {
	bool stat = false;
	for (uInt i = 0; i < (nclus - 1); i++) {
		for (uInt j = i + 1; j < nclus; j++) {
			if (cluster_distance[i * nclus + j] <= BOUND_EPS) {
#if KMEANS_VERBOSE>1
				printf("found identical cluster : %d\n",j);
#endif
				stat = true;
				/* assign the points from j to i */
				const PREC *px = X;
				for (uInt n = 0; n < npts; n++, px += dim) {
					if (c[n] != j)
						continue;
					remove_point_from_cluster(c[n], CX, px, cluster_count, dim,
							ptcount[n]);
					c[n] = i;
					add_point_to_cluster(c[i], CX, px, cluster_count, dim,
							ptcount[n]);
				}
			}
		}
	}
	return (stat);
}

void compute_cluster_distances(BOUND_PREC *dist, BOUND_PREC *s, const PREC *CX,
		uInt dim, uInt nclus, const bool *cluster_changed) {
	for (uInt j = 0; j < nclus; j++)
		s[j] = BOUND_PREC_MAX;

	const PREC *pcx = CX;
	for (uInt i = 0; i < nclus - 1; i++, pcx += dim) {
		const PREC *pcxp = CX + (i + 1) * dim;
		uInt cnt = i * nclus + i + 1;
		for (uInt j = i + 1; j < nclus; j++, cnt++, pcxp += dim) {
			if (cluster_changed[i] || cluster_changed[j]) {
				dist[cnt] = (BOUND_PREC) (0.5
						* compute_distance(pcx, pcxp, dim));
				dist[j * nclus + i] = dist[cnt];

				if (dist[cnt] < s[i])
					s[i] = dist[cnt];

				if (dist[cnt] < s[j])
					s[j] = dist[cnt];
			}
		}
	}
}

uInt init_point_to_cluster(uInt point_ind, const PREC *px, const PREC *CX,
		uInt dim, uInt nclus, PREC *mindist, BOUND_PREC *low_b,
		const BOUND_PREC *cl_dist) {
	bool use_low_b = true;

	if (low_b == NULL)
		use_low_b = false;
	uInt bias = point_ind * nclus;

	const PREC *pcx = CX;
	PREC mind = compute_distance(px, pcx, dim);
	if (use_low_b)
		low_b[bias] = (BOUND_PREC) mind;
	uInt assignment = 0;
	pcx += dim;
	for (uInt j = 1; j < nclus; j++, pcx += dim) {
		if (mind + BOUND_EPS <= cl_dist[assignment * nclus + j])
			continue;

		PREC d = compute_distance(px, pcx, dim);
		if (use_low_b)
			low_b[j + bias] = (BOUND_PREC) d;

		if (d < mind) {
			mind = d;
			assignment = j;
		}
	}
	mindist[point_ind] = mind;
	return (assignment);
}
/*
uInt assign_point_to_cluster_ordinary(const int *px, const PREC *CX, uInt dim,
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
uInt assign_point_to_cluster_ordinary(const PREC *px, const PREC *CX, uInt dim,
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
*/
uInt assign_point_to_cluster(uInt point_ind, const PREC *px, const PREC *CX,
		uInt dim, uInt nclus, uInt old_assignment, PREC *mindist,
		BOUND_PREC *s, BOUND_PREC *cl_dist, BOUND_PREC *low_b) {
	bool up_to_date = false, use_low_b = true;
	;

	uInt bias = point_ind * nclus;
	if (low_b == NULL)
		use_low_b = false;

	PREC mind = mindist[point_ind];

	if (mind + BOUND_EPS <= s[old_assignment]) {
#ifdef KMEANS_VEBOSE
		saved_two++;
#endif
		return (old_assignment);
	}

	uInt assignment = old_assignment;
	uInt counter = assignment * nclus;
	const PREC *pcx = CX;
	for (uInt j = 0; j < nclus; j++, pcx += dim) {
		if (j == old_assignment) {
#if KMEANS_VERBOSE>1
			saved_three_one++;
#endif
			continue;
		}

		if (use_low_b && (mind + BOUND_EPS <= low_b[j + bias])) {
#if KMEANS_VERBOSE>1
			saved_three_two++;
#endif
			continue;
		}

		if (mind + BOUND_EPS <= cl_dist[counter + j]) {
#if KMEANS_VERBOSE>1
			saved_three_three++;
#endif
			continue;
		}

		PREC d = 0.0;
		if (!up_to_date) {
			d = compute_distance(px, CX + assignment * dim, dim);
			mind = d;
			if (use_low_b)
				low_b[assignment + bias] = (BOUND_PREC) d;
			up_to_date = true;
		}

		if (!use_low_b)
			d = compute_distance(px, pcx, dim);
		else if ((mind > BOUND_EPS + low_b[j + bias]) || (mind > BOUND_EPS
				+ cl_dist[counter + j])) {
			d = compute_distance(px, pcx, dim);
			low_b[j + bias] = (BOUND_PREC) d;
		} else {
#if KMEANS_VERBOSE>1
			saved_three_b++;
#endif
			continue;
		}

		if (d < mind) {
			mind = d;
			assignment = j;
			counter = assignment * nclus;
			up_to_date = true;
		}
	}
	mindist[point_ind] = mind;

	return (assignment);
}

PREC kmeans_run(PREC *CX, const PREC *X, uInt *c, uInt dim, uInt npts,
		uInt nclus, uInt maxiter, uInt *ptcount) {
	PREC *tCX = (PREC *) calloc(nclus * dim, sizeof(PREC));
	if (tCX == NULL)
		kmeans_error((char*) "Failed to allocate mem for Cluster points");

	/* number of points per cluster */
	uInt *CN = (uInt *) calloc(nclus, sizeof(uInt));
	if (CX == NULL)
		kmeans_error((char*) "Failed to allocate mem for assignment");

	/* old assignement of points to cluster */
	uInt *old_c = (uInt *) malloc(npts * sizeof(uInt));
	if (old_c == NULL)
		kmeans_error((char*) "Failed to allocate mem for temp assignment");

	/* assign to value which is out of range */
	for (uInt i = 0; i < npts; i++)
		old_c[i] = nclus;

#if KMEANS_VERBOSE>0
	printf("compile without setting the KMEANS_VERBOSE flag for no output\n");
#endif

	BOUND_PREC *low_b = (BOUND_PREC *) calloc(npts * nclus, sizeof(BOUND_PREC));
	bool use_low_b = false;
	if (low_b == NULL) {
#if KMEANS_VERBOSE>0
		printf("not enough memory for lower bound, will compute without\n");
#endif
		use_low_b = false;
	} else {
		use_low_b = true;
		assert(low_b);
	}

	BOUND_PREC *cl_dist = (BOUND_PREC *) calloc(nclus * nclus,
			sizeof(BOUND_PREC));
	if (cl_dist == NULL)
		kmeans_error(
				(char*) "Failed to allocate mem for cluster-cluster distance");

	BOUND_PREC *s = (BOUND_PREC *) malloc(nclus * sizeof(BOUND_PREC));
	if (s == NULL)
		kmeans_error((char*) "Failed to allocate mem for assignment");

	BOUND_PREC *offset = (BOUND_PREC *) malloc(nclus * sizeof(BOUND_PREC)); /* change in distance of a cluster mean after a iteration */
	if (offset == NULL)
		kmeans_error(
				(char*) "Failed to allocate mem for bound points-nearest cluster");

	PREC *mindist = (PREC *) malloc(npts * sizeof(PREC));
	if (mindist == NULL)
		kmeans_error((char*) "Failed to allocate mem for bound points-clusters");

	for (uInt i = 0; i < npts; i++)
		mindist[i] = PREC_MAX;

	bool *cluster_changed = (bool *) malloc(nclus * sizeof(bool)); /* did the cluster changed? */
	if (cluster_changed == NULL)
		kmeans_error(
				(char*) "Failed to allocate mem for variable cluster_changed");
	for (uInt j = 0; j < nclus; j++)
		cluster_changed[j] = true;

	uInt iteration = 0;
	uInt nchanged = 1;
	while (iteration < maxiter || maxiter == 0) {

		/* compute cluster-cluster distances */
		compute_cluster_distances(cl_dist, s, CX, dim, nclus, cluster_changed);

		/* assign all points from identical clusters to the first occurence of that cluster */
		remove_identical_clusters(CX, cl_dist, X, CN, c, dim, nclus, npts,
				ptcount);

		/* find nearest cluster center */
		if (iteration == 0) {

			const PREC *px = X;
			for (uInt i = 0; i < npts; i++, px += dim) {
				c[i] = init_point_to_cluster(i, px, CX, dim, nclus, mindist,
						low_b, cl_dist);
				add_point_to_cluster(c[i], tCX, px, CN, dim, ptcount[i]);
			}
			nchanged = npts;
		} else {
			for (uInt j = 0; j < nclus; j++)
				cluster_changed[j] = false;

			nchanged = 0;
			const PREC *px = X;
			for (uInt i = 0; i < npts; i++, px += dim) {
				c[i] = assign_point_to_cluster(i, px, CX, dim, nclus, old_c[i],
						mindist, s, cl_dist, low_b);

#ifdef KMEANS_DEBUG
				{
					/* If the assignments are not the same, there is still the BOUND_EPS difference 
					 which can be the reason of this*/
					uInt tmp = assign_point_to_cluster_ordinary(px,CX,dim,nclus);
					if (tmp != c[i])
					{
						PREC d1 = compute_distance(px,CX+(tmp*dim),dim);
						PREC d2 = compute_distance(px,CX+(c[i]*dim),dim);
						assert( (d1>d2)?((d1-d2)<BOUND_EPS):((d2-d1)<BOUND_EPS) );
					}
				}
#endif

				if (old_c[i] == c[i])
					continue;

				nchanged++;

				cluster_changed[c[i]] = true;
				cluster_changed[old_c[i]] = true;

				remove_point_from_cluster(old_c[i], tCX, px, CN, dim,
						ptcount[i]);
				add_point_to_cluster(c[i], tCX, px, CN, dim, ptcount[i]);
			}

		}

		/* fill up empty clusters */
		for (uInt j = 0; j < nclus; j++) {
			if (CN[j] > 0)
				continue;
			uInt *rperm = (uInt*) malloc(npts * sizeof(uInt));
			if (cluster_changed == NULL)
				kmeans_error((char*) "Failed to allocate mem for permutation");

			randperm(rperm, npts);
			uInt i = 0;
			while (rperm[i] < npts && CN[c[rperm[i]]] < 2)
				i++;
			if (i == npts)
				continue;
			i = rperm[i];
#if KMEANS_VERBOSE>0
			printf("empty cluster [%d], filling it with point [%d]\n",j,i);
#endif
			cluster_changed[c[rperm[i]]] = true;
			cluster_changed[j] = true;
			const PREC *px = X + i * dim;
			remove_point_from_cluster(c[i], tCX, px, CN, dim, ptcount[i]);
			c[i] = j;
			add_point_to_cluster(j, tCX, px, CN, dim, ptcount[i]);
			/* void the bounds */
			s[j] = (BOUND_PREC) 0.0;
			mindist[i] = 0.0;
			if (use_low_b)
				for (uInt k = 0; k < npts; k++)
					low_b[k * nclus + j] = (BOUND_PREC) 0.0;

			nchanged++;
			free(rperm);
		}

		/* no assignment changed: done */
		if (nchanged == 0)
			break;

		/* compute the offset */

		PREC *pcx = CX;
		PREC *tpcx = tCX;
		for (uInt j = 0; j < nclus; j++, pcx += dim, tpcx += dim) {
			offset[j] = (BOUND_PREC) 0.0;
			if (cluster_changed[j]) {
				offset[j] = (BOUND_PREC) compute_distance(pcx, tpcx, dim);
				memcpy(pcx, tpcx, dim * sizeof(PREC));
			}
		}

		/* update the lower bound */
		if (use_low_b) {
			for (uInt i = 0, cnt = 0; i < npts; i++)
				for (uInt j = 0; j < nclus; j++, cnt++) {
					low_b[cnt] -= offset[j];
					if (low_b[cnt] < (BOUND_PREC) 0.0)
						low_b[cnt] = (BOUND_PREC) 0.0;
				}
		}

		for (uInt i = 0; i < npts; i++)
			mindist[i] += (PREC) offset[c[i]];

		memcpy(old_c, c, npts * sizeof(uInt));

#if KMEANS_VERBOSE>0
		PREC sse = compute_sserror(CX,X,c,dim,npts,ptcount);
		printf("iteration %4d, #(changed points): %4d, sse: %4.2f\n",(int)iteration,(int)nchanged,sse);
#endif

#if KMEANS_VERBOSE>1
		printf("saved at 2) %d\n",saved_two);
		printf("saved at 3i) %d\n",saved_three_one);
		printf("saved at 3ii) %d\n",saved_three_two);
		printf("saved at 3iii) %d\n",saved_three_three);
		printf("saved at 3b) %d\n",saved_three_b);
		saved_two=0;
		saved_three_one=0;
		saved_three_two=0;
		saved_three_three=0;
		saved_three_b=0;
#endif

		iteration++;

	}

#ifdef KMEANS_DEBUG
	for ( uInt j=0;j<nclus;j++)
	assert(CN[j]!=0); /* Empty cluster after all */
#endif

	/* find nearest cluster center if iteration reached maxiter */
	if (nchanged > 0) {
		const PREC *px = X;
		for (uInt i = 0; i < npts; i++, px += dim)
			c[i] = assign_point_to_cluster_ordinary(px, CX, dim, nclus);
	}
	PREC sse = compute_sserror(CX, X, c, dim, npts, ptcount);

#if KMEANS_VERBOSE>0
	printf("iteration %4d, #(changed points): %4d, sse: %4.2f\n",(int)iteration,(int)nchanged,sse);
#endif

	if (low_b)
		free(low_b);
	free(cluster_changed);
	free(mindist);
	free(s);
	free(offset);
	free(cl_dist);
	free(tCX);
	free(CN);
	free(old_c);

	return (sse);
}

PREC kmeans(PREC *CX, const PREC *X, uInt *assignment, uInt dim, uInt npts,
		uInt nclus, uInt maxiter, uInt restarts, uInt *ptcount) {
	//pt count refers to number of count each point...
	if (npts < nclus) {
		CX = (PREC*) calloc(nclus * dim, sizeof(PREC));
		memcpy(CX, X, dim * nclus * sizeof(PREC));
		PREC sse = 0.0;
		return (sse);
	} else if (npts == nclus) {
		memcpy(CX, X, dim * nclus * sizeof(PREC));
		PREC sse = 0.0;
		return (sse);
	} else if (nclus == 0) {
		printf("Error: Number of clusters is 0\n");
		exit(-1);
	}

	/*
	 * No starting point is given, generate a new one
	 */
	if (CX == NULL) {
		uInt *order = (uInt*) malloc(npts * sizeof(uInt));

		CX = (PREC*) calloc(nclus * dim, sizeof(PREC));
		/* generate new starting point */
		randperm(order, npts);
		for (uInt i = 0; i < nclus; i++)
			for (uInt k = 0; k < dim; k++)
				CX[(i * dim) + k] = X[order[i] * dim + k];
		free(order);

	}
	assert(CX != NULL);
	PREC sse =
			kmeans_run(CX, X, assignment, dim, npts, nclus, maxiter, ptcount);

	uInt res = restarts;
	if (res > 0) {
		PREC minsse = sse;
		uInt *order = (uInt*) malloc(npts * sizeof(uInt));
		PREC *bestCX = (PREC*) malloc(dim * nclus * sizeof(PREC));
		uInt *bestassignment = (uInt*) malloc(npts * sizeof(uInt));

		memcpy(bestCX, CX, dim * nclus * sizeof(PREC));
		memcpy(bestassignment, assignment, npts * sizeof(uInt));

		while (res > 0) {

			/* generate new starting point */
			randperm(order, npts);
			for (uInt i = 0; i < nclus; i++)
				for (uInt k = 0; k < dim; k++)
					CX[(i * dim) + k] = X[order[i] * dim + k];

			sse = kmeans_run(CX, X, assignment, dim, npts, nclus, maxiter,
					ptcount);
			if (sse < minsse) {
#if KMEANS_VERBOSE>1
				printf("found a better clustering with sse = %g\n",sse);
#endif
				minsse = sse;
				memcpy(bestCX, CX, dim * nclus * sizeof(PREC));
				memcpy(bestassignment, assignment, npts * sizeof(uInt));
			}
			printf("\n For Round # = %d , SSE = %f, Min SSE = %f", restarts
					- res, sse, minsse);
			res--;

		}
		memcpy(CX, bestCX, dim * nclus * sizeof(PREC));
		memcpy(assignment, bestassignment, npts * sizeof(uInt));
		sse = minsse;
		free(bestassignment);
		free(bestCX);
		free(order);
	}
	assert(CX != NULL);

	return (sse);

}
