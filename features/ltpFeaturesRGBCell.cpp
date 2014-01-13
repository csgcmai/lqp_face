/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/

/*
 * features.cpp
 *
 *  Created on: Apr 1, 2009
 *      Author: shussain
 */
#define EPS 0.0001
#include "ltpFeaturesRGBCell.h"
void LTPFeaturesRGBCell::InitalizeMaps(Image &image) {
	sspace = NoPyramid;
	UINT xbmax, ybmax;
	lbpfeatures.Initialize(1, 0, nbins * 3 * dimmult);
	GetXYBlocks(image, xbmax, ybmax);
	lbpfeatures.SetFeature(xbmax, ybmax, 1.0);
	ComputeLBPFeatures(image, 0, cell, lbpstride);
}
void LTPFeaturesRGBCell::InitalizeMaps(Pyramid &pyobj, PyramidType sspace_) {
	UINT xbmax, ybmax, pycount = 0;
	sspace = sspace_;
	// Read the image....
	if (sspace == SingleRes) {
		pycount = pyobj.GetNLevels();
		lbpfeatures.Initialize(pycount, 0, nbins * 3 * dimmult);
		for (UINT i = 0; i < pycount; i++) {
			Image& imgref = pyobj.GetImageRef(i);
			// Compute the Gradient Image
			GetXYBlocks(imgref, xbmax, ybmax);
			lbpfeatures.SetFeature(xbmax, ybmax, pyobj.GetScale(i));
			ComputeLBPFeatures(imgref, i, cell, lbpstride);
		}
	} else /*if (sspace == DoubleRes)*/{

		pycount = pyobj.GetNLevels();
		UINT nlpyramid = pyobj.GetInitIndex();
		lbpfeatures.Initialize(pycount + nlpyramid, nlpyramid,
				nbins * 3 * dimmult);
		// Number of level in 2x...
		for (UINT i = 0; i < pycount; i++) {
			Image &imgref = pyobj.GetImageRef(i);

			if (i < nlpyramid) {
				Get2XYBlocks(imgref, xbmax, ybmax);
				lbpfeatures.SetFeature(i, xbmax, ybmax, pyobj.GetScale(i) / 2);
				ComputeLBPFeatures(imgref, i, cell / 2, lbpstride / 2);
			}
			// Compute the features...
			//		imgref=pyobj.GetImageRef(i);
			GetXYBlocks(imgref, xbmax, ybmax);
			lbpfeatures.SetFeature(i + nlpyramid, xbmax, ybmax,
					pyobj.GetScale(i));
			ComputeLBPFeatures(imgref, i + nlpyramid, cell, lbpstride);
		}
	}
}
void LTPFeaturesRGBCell::InitalizeMaps(Image&image, PyramidType sspace_) {
	UINT xbmax, ybmax, pycount = 0;
	sspace = sspace_;
	// Read the image....
	if (sspace == SingleRes) {
		pycount = pyobj.GeneratePyramid(image);
		lbpfeatures.Initialize(pycount, 0, nbins * 3 * dimmult);
		for (UINT i = 0; i < pycount; i++) {
			Image& imgref = pyobj.GetImageRef(i);
			// Compute the Gradient Image
			GetXYBlocks(imgref, xbmax, ybmax);
			lbpfeatures.SetFeature(xbmax, ybmax, pyobj.GetScale(i));
			ComputeLBPFeatures(imgref, i, cell, lbpstride);
		}
	} else if (sspace == DoubleRes) {

		pycount = pyobj.GeneratePyramid(image);
		UINT nlpyramid = pyobj.GetInitIndex();
		lbpfeatures.Initialize(pycount + nlpyramid, nlpyramid,
				nbins * 3 * dimmult);
		// Number of level in 2x...
		for (UINT i = 0; i < pycount; i++) {
			Image &imgref = pyobj.GetImageRef(i);

			if (i < nlpyramid) {
				Get2XYBlocks(imgref, xbmax, ybmax);
				lbpfeatures.SetFeature(i, xbmax, ybmax, pyobj.GetScale(i) / 2);
				ComputeLBPFeatures(imgref, i, cell / 2, lbpstride / 2);
			}
			// Compute the features...
			//		imgref=pyobj.GetImageRef(i);
			GetXYBlocks(imgref, xbmax, ybmax);
			lbpfeatures.SetFeature(i + nlpyramid, xbmax, ybmax,
					pyobj.GetScale(i));
			ComputeLBPFeatures(imgref, i + nlpyramid, cell, lbpstride);
		}
	} else {
		lbpfeatures.Initialize(1, 0, nbins * 3 * dimmult);
		GetXYBlocks(image, xbmax, ybmax);
		lbpfeatures.SetFeature(xbmax, ybmax, 1.0);
		ComputeLBPFeatures(image, 0, cell, lbpstride);
	}
}
void LTPFeaturesRGBCell::ComputeLBPFeatures(Image &image, UINT index,
		UINT sbin, UINT lbpstride) {
	LBPMap tlbpmap[3], tposmap[3], tnegmap[3];
	ComputeLTPMap(image, tlbpmap, tposmap, tnegmap);

	if (hmethod == HM_Bilinear) {
		vector<REAL> &tvec = lbpfeatures.GetFeatureRef(index);
		REAL *feat = &tvec[0];
		ComputeHistBilinear(image, sbin, tlbpmap, tposmap, tnegmap, feat);
	} else
		ComputeHistDiscrete(image, index, sbin, tlbpmap, tposmap, tnegmap);

}
void LTPFeaturesRGBCell::ComputeHistDiscrete(Image& imgref, UINT index,
		UINT winstride, LBPMap *tlbpmap, LBPMap *tposmap, LBPMap *tnegmap) {
	// for the speed modularization is sacrificed...

	UINT x, y, xbmax = lbpfeatures.GetXBlock(index) * winstride, ybmax =
			lbpfeatures.GetYBlock(index) * winstride, offset = 0, endx, endy;

	vector<REAL> &features = lbpfeatures.GetFeatureRef(index);
	vector<vector<REAL> > lbpfeat(3, vector<REAL> (nbins)), pltpfeat(3,
			vector<REAL> (nbins)), nltpfeat(3, vector<REAL> (nbins));
	REAL lbpsum[3], possum[3], negsum[3]; // normalization constants...
	//	fill(features.begin(), features.end(), 0);

	for (UINT rind = 0; rind < ybmax; rind += winstride) {
		y = rind + winstride/*foffset*/- radius;
		for (UINT colind = 0; colind < xbmax; colind += winstride) {
			x = colind + winstride /*foffset*/- radius;

			endx = MIN(x+winstride,tlbpmap[0].nxpoints);
			endy = MIN(y+winstride,tlbpmap[0].nypoints);

			for (UINT yiter = y; yiter < endy; ++yiter)
				for (UINT xiter = x; xiter < endx; ++xiter) {
					lbpfeat[0][tlbpmap[0].GetValue(xiter, yiter)]++;// red
					lbpfeat[1][tlbpmap[1].GetValue(xiter, yiter)]++;// green
					lbpfeat[2][tlbpmap[2].GetValue(xiter, yiter)]++;// blue

					pltpfeat[0][tposmap[0].GetValue(xiter, yiter)]++;
					pltpfeat[1][tposmap[1].GetValue(xiter, yiter)]++;
					pltpfeat[2][tposmap[2].GetValue(xiter, yiter)]++;

					nltpfeat[0][tnegmap[0].GetValue(xiter, yiter)]++;
					nltpfeat[1][tnegmap[1].GetValue(xiter, yiter)]++;
					nltpfeat[2][tnegmap[2].GetValue(xiter, yiter)]++;

				}
			if (!add) { // don't add but concatenate the RGB channels....
				// normalize each channel separately and then Concatenate....
				lbpsum[0] = lbpsum[1] = lbpsum[2] = EPS;
				possum[0] = possum[1] = possum[2] = EPS;
				negsum[0] = negsum[1] = negsum[2] = EPS;

				// Currently only doing L1-Sqrt Normalization further to be added...
				for (UINT biter = 0; biter < nbins; ++biter) {
					lbpsum[0] += lbpfeat[0][biter];
					lbpsum[1] += lbpfeat[1][biter];
					lbpsum[2] += lbpfeat[2][biter];

					possum[0] += pltpfeat[0][biter];
					possum[1] += pltpfeat[1][biter];
					possum[2] += pltpfeat[2][biter];

					negsum[0] += nltpfeat[0][biter];
					negsum[1] += nltpfeat[1][biter];
					negsum[2] += nltpfeat[2][biter];
				}
				REAL *lbpptr = &features[offset], *posptr = &features[offset
						+ 3 * nbins], *negptr = &features[offset + 6 * nbins];
				for (UINT biter = 0; biter < nbins; ++biter) {
					lbpptr[biter] = sqrt(lbpfeat[0][biter] / lbpsum[0]); //r for lbp features....
					lbpptr[biter + nbins] = sqrt(lbpfeat[1][biter] / lbpsum[1]);//g
					lbpptr[biter + nbins * 2] = sqrt(
							lbpfeat[2][biter] / lbpsum[2]);//b

					posptr[biter] = sqrt(pltpfeat[0][biter] / possum[0]);
					posptr[biter + nbins]
							= sqrt(pltpfeat[1][biter] / possum[1]);
					posptr[biter + nbins * 2] = sqrt(
							pltpfeat[2][biter] / possum[2]);

					negptr[biter] = sqrt(nltpfeat[0][biter] / negsum[0]);
					negptr[biter + nbins]
							= sqrt(nltpfeat[1][biter] / negsum[1]);
					negptr[biter + nbins * 2] = sqrt(
							nltpfeat[2][biter] / negsum[2]);

				}
			} else {

				if (nseparate) {
					// normalize each channel separately and then add....
					lbpsum[0] = lbpsum[1] = lbpsum[2] = EPS;
					possum[0] = possum[1] = possum[2] = EPS;
					negsum[0] = negsum[1] = negsum[2] = EPS;

					// Currently only doing L1-Sqrt Normalization further to be added...
					for (UINT biter = 0; biter < nbins; ++biter) {
						lbpsum[0] += lbpfeat[0][biter];
						lbpsum[1] += lbpfeat[1][biter];
						lbpsum[2] += lbpfeat[2][biter];

						possum[0] += pltpfeat[0][biter];
						possum[1] += pltpfeat[1][biter];
						possum[2] += pltpfeat[2][biter];

						negsum[0] += nltpfeat[0][biter];
						negsum[1] += nltpfeat[1][biter];
						negsum[2] += nltpfeat[2][biter];
					}

					for (UINT biter = 0; biter < nbins; ++biter) {
						lbpfeat[0][biter] = sqrt(lbpfeat[0][biter] / lbpsum[0]);
						lbpfeat[1][biter] = sqrt(lbpfeat[1][biter] / lbpsum[1]);
						lbpfeat[2][biter] = sqrt(lbpfeat[2][biter] / lbpsum[2]);

						pltpfeat[0][biter] = sqrt(
								pltpfeat[0][biter] / possum[0]);
						pltpfeat[1][biter] = sqrt(
								pltpfeat[1][biter] / possum[1]);
						pltpfeat[2][biter] = sqrt(
								pltpfeat[2][biter] / possum[2]);

						nltpfeat[0][biter] = sqrt(
								nltpfeat[0][biter] / negsum[0]);
						nltpfeat[1][biter] = sqrt(
								nltpfeat[1][biter] / negsum[1]);
						nltpfeat[2][biter] = sqrt(
								nltpfeat[2][biter] / negsum[2]);

					}

					for (UINT biter = 0; biter < nbins; ++biter) {
						features[offset + biter] = lbpfeat[0][biter]
								+ lbpfeat[1][biter] + lbpfeat[2][biter];
						features[nbins + offset + biter] = pltpfeat[0][biter]
								+ pltpfeat[1][biter] + pltpfeat[2][biter];
						features[2 * nbins + offset + biter]
								= nltpfeat[0][biter] + nltpfeat[1][biter]
										+ nltpfeat[2][biter];
					}
				} else {

					REAL *lbpptr = &features[offset], *posptr =
							&features[offset + nbins], *negptr =
							&features[offset + 2 * nbins];

					lbpsum[0] = EPS;
					possum[0] = EPS;
					negsum[0] = EPS;

					for (UINT biter = 0; biter < nbins; ++biter) { // add the three channels
						lbpptr[biter] = lbpfeat[0][biter] + lbpfeat[1][biter]
								+ lbpfeat[2][biter];
						posptr[biter] = pltpfeat[0][biter] + pltpfeat[1][biter]
								+ pltpfeat[2][biter];
						negptr[biter] = nltpfeat[0][biter] + nltpfeat[1][biter]
								+ nltpfeat[2][biter];

						lbpsum[0] += lbpptr[biter];
						possum[0] += posptr[biter];
						negsum[0] += negptr[biter];
					}

					for (UINT biter = 0; biter < nbins; ++biter) {
						lbpptr[biter] = sqrt(lbpptr[biter] / lbpsum[0]);
						posptr[biter] = sqrt(posptr[biter] / possum[0]);
						negptr[biter] = sqrt(negptr[biter] / negsum[0]);
					}

				}
			}

			for (UINT i = 0; i < 3; ++i) {// Zero out the temporary  histograms...
				fill(lbpfeat[i].begin(), lbpfeat[i].end(), 0);
				fill(pltpfeat[i].begin(), pltpfeat[i].end(), 0);
				fill(nltpfeat[i].begin(), nltpfeat[i].end(), 0);
			}
			offset += 3 * nbins;
		}

	}

}
void LTPFeaturesRGBCell::ComputeHistBilinear(Image &image, UINT sbin,
		LBPMap *tlbpmap, LBPMap *tposmap, LBPMap *tnegmap, REAL *feat) {
	int dims[2];
	dims[0] = image.rows();
	dims[1] = image.columns();

	// memory for caching histograms & their norms
	int blocks[2], histdim, bdim;
	blocks[0] = (int) round((double) dims[0] / (double) sbin);
	blocks[1] = (int) round((double) dims[1] / (double) sbin);
	histdim = blocks[0] * blocks[1] * nbins;
	bdim = blocks[0] * blocks[1];
	double *rhist = new double[(histdim)], *ghist = new double[(histdim)],
			*bhist = new double[(histdim)], *hist = new double[(histdim)];

	double *prhist = new double[(histdim)], *pghist = new double[(histdim)],
			*pbhist = new double[(histdim)], *phist = new double[(histdim)];
	double *nrhist = new double[(histdim)], *nghist = new double[(histdim)],
			*nbhist = new double[(histdim)], *nhist = new double[(histdim)];
	double *lbpnorm = new double[(bdim)], *pnorm = new double[(bdim)], *nnorm =
			new double[(blocks[0] * blocks[1])];
	fill(lbpnorm, lbpnorm + bdim, 0);
	fill(pnorm, pnorm + bdim, 0);
	fill(nnorm, nnorm + bdim, 0);
	fill(rhist, rhist + histdim, 0);
	fill(ghist, ghist + histdim, 0);
	fill(bhist, bhist + histdim, 0);
	fill(hist, hist + histdim, 0); // for pooling of RGB channels....

	fill(prhist, prhist + histdim, 0);
	fill(pghist, pghist + histdim, 0);
	fill(pbhist, pbhist + histdim, 0);
	fill(phist, phist + histdim, 0); // for pooling of RGB channels....

	fill(nrhist, nrhist + histdim, 0);
	fill(nghist, nghist + histdim, 0);
	fill(nbhist, nbhist + histdim, 0);
	fill(nhist, nhist + histdim, 0); // for pooling of RGB channels....

	// memory for LTP features
	int out[2];
	out[0] = max(blocks[0] - 2, 0);
	out[1] = max(blocks[1] - 2, 0);

	int visible[2];
	visible[0] = blocks[0] * sbin;
	visible[1] = blocks[1] * sbin;
#ifdef INVERSE
	for (int x = radius; x < visible[1] - radius; x++) {
		for (int y = radius; y < visible[0] - radius; y++) {

			// first color channel
			/*because the value of map is offset by radius, value of x at x-radius */
			/*-1 because of < ; MIN(x, dims[1] - radius-1) for which lbp map is computed;
			 * final -radius the offset in lbpmap*/
			tx = MIN(x, dims[1] - radius-1) - radius;

			ty = MIN(y, dims[0] - radius-1) - radius;
#else
	UINT nxpoints = tposmap[0].nxpoints, nypoints = tposmap[0].nypoints;
	for (int tx = 0; tx < nxpoints; ++tx) {
		for (int ty = 0; ty < nypoints; ++ty) {

			// first color channel
			/*because the value of map is offset by radius, value of x at x-radius */
			/*-1 because of < ; MIN(x, dims[1] - radius-1) for which lbp map is computed;
			 * final -radius the offset in lbpmap*/
			int x = tx + radius;
			int y = ty + radius;
#endif
			int ro = tlbpmap[0].GetValue(tx, ty), // red, GREEN & BLUE
					go = tlbpmap[1].GetValue(tx, ty), //
					bo = tlbpmap[2].GetValue(tx, ty); //
			int pro = tposmap[0].GetValue(tx, ty), // red, GREEN & BLUE
					pgo = tposmap[1].GetValue(tx, ty), //
					pbo = tposmap[2].GetValue(tx, ty); //
			int nro = tnegmap[0].GetValue(tx, ty), // red, GREEN & BLUE
					ngo = tnegmap[1].GetValue(tx, ty), //
					nbo = tnegmap[2].GetValue(tx, ty); //
			// add to 4 histograms around pixel using linear interpolation
			double xp = ((double) x + 0.5) / (double) sbin - 0.5;
			double yp = ((double) y + 0.5) / (double) sbin - 0.5;
			int ixp = (int) floor(xp);
			int iyp = (int) floor(yp);
			double vx0 = xp - ixp;
			double vy0 = yp - iyp;
			double vx1 = 1.0 - vx0;
			double vy1 = 1.0 - vy0;
			double v = 1; // pooling of gradient magnitude can be an option


			if (ixp >= 0 && iyp >= 0) {
				*(rhist + ixp * blocks[0] + iyp + ro * bdim) += vx1 * vy1 * v;
				*(ghist + ixp * blocks[0] + iyp + go * bdim) += vx1 * vy1 * v;
				*(bhist + ixp * blocks[0] + iyp + bo * bdim) += vx1 * vy1 * v;

				*(prhist + ixp * blocks[0] + iyp + pro * bdim) += vx1 * vy1 * v;
				*(pghist + ixp * blocks[0] + iyp + pgo * bdim) += vx1 * vy1 * v;
				*(pbhist + ixp * blocks[0] + iyp + pbo * bdim) += vx1 * vy1 * v;

				*(nrhist + ixp * blocks[0] + iyp + nro * bdim) += vx1 * vy1 * v;
				*(nghist + ixp * blocks[0] + iyp + ngo * bdim) += vx1 * vy1 * v;
				*(nbhist + ixp * blocks[0] + iyp + nbo * bdim) += vx1 * vy1 * v;

			}

			if (ixp + 1 < blocks[1] && iyp >= 0) {
				*(rhist + (ixp + 1) * blocks[0] + iyp + ro * bdim) += vx0 * vy1
						* v;
				*(ghist + (ixp + 1) * blocks[0] + iyp + go * bdim) += vx0 * vy1
						* v;
				*(bhist + (ixp + 1) * blocks[0] + iyp + bo * bdim) += vx0 * vy1
						* v;

				*(prhist + (ixp + 1) * blocks[0] + iyp + pro * bdim) += vx0
						* vy1 * v;
				*(pghist + (ixp + 1) * blocks[0] + iyp + pgo * bdim) += vx0
						* vy1 * v;
				*(pbhist + (ixp + 1) * blocks[0] + iyp + pbo * bdim) += vx0
						* vy1 * v;

				*(nrhist + (ixp + 1) * blocks[0] + iyp + nro * bdim) += vx0
						* vy1 * v;
				*(nghist + (ixp + 1) * blocks[0] + iyp + ngo * bdim) += vx0
						* vy1 * v;
				*(nbhist + (ixp + 1) * blocks[0] + iyp + nbo * bdim) += vx0
						* vy1 * v;
			}

			if (ixp >= 0 && iyp + 1 < blocks[0]) {
				*(rhist + ixp * blocks[0] + (iyp + 1) + ro * bdim) += vx1 * vy0
						* v;
				*(ghist + ixp * blocks[0] + (iyp + 1) + go * bdim) += vx1 * vy0
						* v;
				*(bhist + ixp * blocks[0] + (iyp + 1) + bo * bdim) += vx1 * vy0
						* v;

				*(prhist + ixp * blocks[0] + (iyp + 1) + pro * bdim) += vx1
						* vy0 * v;
				*(pghist + ixp * blocks[0] + (iyp + 1) + pgo * bdim) += vx1
						* vy0 * v;
				*(pbhist + ixp * blocks[0] + (iyp + 1) + pbo * bdim) += vx1
						* vy0 * v;

				*(nrhist + ixp * blocks[0] + (iyp + 1) + nro * bdim) += vx1
						* vy0 * v;
				*(nghist + ixp * blocks[0] + (iyp + 1) + ngo * bdim) += vx1
						* vy0 * v;
				*(nbhist + ixp * blocks[0] + (iyp + 1) + nbo * bdim) += vx1
						* vy0 * v;

			}

			if (ixp + 1 < blocks[1] && iyp + 1 < blocks[0]) {
				*(rhist + (ixp + 1) * blocks[0] + (iyp + 1) + ro * bdim) += vx0
						* vy0 * v;
				*(ghist + (ixp + 1) * blocks[0] + (iyp + 1) + go * bdim) += vx0
						* vy0 * v;
				*(bhist + (ixp + 1) * blocks[0] + (iyp + 1) + bo * bdim) += vx0
						* vy0 * v;

				*(prhist + (ixp + 1) * blocks[0] + (iyp + 1) + pro * bdim)
						+= vx0 * vy0 * v;
				*(pghist + (ixp + 1) * blocks[0] + (iyp + 1) + pgo * bdim)
						+= vx0 * vy0 * v;
				*(pbhist + (ixp + 1) * blocks[0] + (iyp + 1) + pbo * bdim)
						+= vx0 * vy0 * v;

				*(nrhist + (ixp + 1) * blocks[0] + (iyp + 1) + nro * bdim)
						+= vx0 * vy0 * v;
				*(nghist + (ixp + 1) * blocks[0] + (iyp + 1) + ngo * bdim)
						+= vx0 * vy0 * v;
				*(nbhist + (ixp + 1) * blocks[0] + (iyp + 1) + nbo * bdim)
						+= vx0 * vy0 * v;

			}
		}
	}
	// compute energy in each block by summing over orientations
	for (int o = 0; o < nbins; o++) {
		double *srcr = rhist + o * bdim, *srcg = ghist + o * bdim, *srcb =
				bhist + o * bdim, *src = hist + o * bdim;

		double *psrcr = prhist + o * bdim, *psrcg = pghist + o * bdim, *psrcb =
				pbhist + o * bdim, *psrc = phist + o * bdim;

		double *nsrcr = nrhist + o * bdim, *nsrcg = nghist + o * bdim, *nsrcb =
				nbhist + o * bdim, *nsrc = nhist + o * bdim;

		double *dst = lbpnorm, *pdst = pnorm, *ndst = nnorm;
		double *end = lbpnorm + blocks[1] * blocks[0];
		while (dst < end) {
			*src = *srcr + *srcg + *srcb; // pooling of rgb channels...
			*psrc = *psrcr + *psrcg + *psrcb; // pooling of rgb channels...
			*nsrc = *nsrcr + *nsrcg + *nsrcb; // pooling of rgb channels...
			if (norm == L2) {
				*(dst++) += *src + *src;
				*(pdst++) += *psrc + *psrc;
				*(ndst++) += *nsrc + *nsrc;
			} else {
				*(dst++) += *src;
				*(pdst++) += *psrc;
				*(ndst++) += *nsrc;
			}
			srcr++;
			srcg++;
			srcb++;
			src++;

			psrcr++;
			psrcg++;
			psrcb++;
			psrc++;

			nsrcr++;
			nsrcg++;
			nsrcb++;
			nsrc++;
		}
	}
	// compute features
	REAL sc1, sc2, sc3, sc4;
	for (int x = 0; x < out[1]; x++) {
		for (int y = 0; y < out[0]; y++) {
			REAL *dst = feat + nbins * 3 * (x + y * out[1]); // write in column major order
			REAL *pdst = dst + nbins, *ndst = pdst + nbins;
			double *p, n1, n2, n3, n4;
			double *src = hist + (x + 1) * blocks[0] + (y + 1);
			double *psrc = phist + (x + 1) * blocks[0] + (y + 1);
			double *nsrc = nhist + (x + 1) * blocks[0] + (y + 1);
			// single cell normalization...
			if (norm == L2) {
				REAL normval = 1 / sqrt(
						*(lbpnorm + (x + 1) * blocks[0] + y + 1) + EPS);
				REAL pnormval = 1 / sqrt(
						*(pnorm + (x + 1) * blocks[0] + y + 1) + EPS);
				REAL nnormval = 1 / sqrt(
						*(nnorm + (x + 1) * blocks[0] + y + 1) + EPS);
				for (int o = 0; o < nbins; o++) {
					*dst = *src * normval;
					*pdst = *psrc * pnormval;
					*ndst = *nsrc * nnormval;
					src += bdim;
					psrc += bdim;
					nsrc += bdim;
					++dst;
					++pdst;
					++ndst;
				}
			} else if (norm == LOG2) {

				for (int o = 0; o < nbins; o++) {
					*dst = log2(*src + 1);
					*pdst = log2(*psrc + 1);
					*ndst = log2(*nsrc + 1);
					src += bdim;
					psrc += bdim;
					nsrc += bdim;
					++dst;
					++pdst;
					++ndst;
				}
			} else {
				REAL normval = *(lbpnorm + (x + 1) * blocks[0] + y + 1) + EPS;
				REAL pnormval = *(pnorm + (x + 1) * blocks[0] + y + 1) + EPS;
				REAL nnormval = *(nnorm + (x + 1) * blocks[0] + y + 1) + EPS;
				for (int o = 0; o < nbins; o++) {
					*dst = sqrt(*src / normval);
					*pdst = sqrt(*psrc / pnormval);
					*ndst = sqrt(*nsrc / nnormval);
					src += bdim;
					psrc += bdim;
					nsrc += bdim;
					++dst;
					++pdst;
					++ndst;
				}
			}
		}
	}
	delete[] hist;
	delete[] rhist;
	delete[] ghist;
	delete[] bhist;
	delete[] lbpnorm;

	delete[] phist;
	delete[] prhist;
	delete[] pghist;
	delete[] pbhist;
	delete[] pnorm;

	delete[] nhist;
	delete[] nrhist;
	delete[] nghist;
	delete[] nbhist;
	delete[] nnorm;
}
void LTPFeaturesRGBCell::ComputeTextureFeatures(Image &imgref, UINT xcellsize,
		UINT ycellsize, vector<REAL>&features) {
	ComputeTextureFeatures(imgref, xcellsize, ycellsize, 0, 0,
			imgref.columns() - 1, imgref.rows() - 1, features);
}
void LTPFeaturesRGBCell::ComputeTextureFeatures(Image &imgref, int sx, int sy,
		int ex, int ey, UINT xcellsize, UINT ycellsize, vector<REAL>&features) {
	LBPMap tlbpmap[3], tposmap[3], tnegmap[3];
	ComputeLTPMap(imgref, tlbpmap, tposmap, tnegmap);
	UINT width = ex - sx + 1, height = ey - sy + 1;
	int xoffset = -radius/*+MAX(ceil((REAL) (width % xcellsize) / 2),radius)*/;
	int yoffset = -radius/*+MAX(ceil((REAL) (height % ycellsize) / 2),radius)*/;
	UINT xbmax = width / xcellsize, ybmax = height / ycellsize;
	UINT offset = 0;
	features.resize(xbmax * ybmax * nbins * dimmult * 3, 0);

	vector<vector<REAL> > lbpfeat(3, vector<REAL> (nbins)), pltpfeat(3,
			vector<REAL> (nbins)), nltpfeat(3, vector<REAL> (nbins));
	REAL lbpsum[3], possum[3], negsum[3]; // normalization constants...

	for (UINT rind = 0; rind < ybmax; rind += ycellsize) {
		for (UINT colind = 0; colind < xbmax; colind += xcellsize) {

			for (int yiter = rind; yiter < rind + ycellsize; ++yiter) {
				int y = MIN(MAX(yiter+yoffset+sy,0),tposmap[0].nypoints-1);
				for (int xiter = colind; xiter < colind + xcellsize; ++xiter) {
					int x = MIN(MAX(xiter+xoffset+sx,0),tposmap[0].nxpoints-1);
					assert(tposmap[0].GetValue(x, y)!=-1);

					lbpfeat[0][tlbpmap[0].GetValue(x, y)]++;// red
					lbpfeat[1][tlbpmap[1].GetValue(x, y)]++;// green
					lbpfeat[2][tlbpmap[2].GetValue(x, y)]++;// blue

					pltpfeat[0][tposmap[0].GetValue(x, y)]++;
					pltpfeat[1][tposmap[1].GetValue(x, y)]++;
					pltpfeat[2][tposmap[2].GetValue(x, y)]++;

					nltpfeat[0][tnegmap[0].GetValue(x, y)]++;
					nltpfeat[1][tnegmap[1].GetValue(x, y)]++;
					nltpfeat[2][tnegmap[2].GetValue(x, y)]++;

				}
			}
			{// normalization

				REAL *lbpptr = &features[offset], *posptr = &features[offset
						+ nbins], *negptr = &features[offset + 2 * nbins];

				lbpsum[0] = EPS;
				possum[0] = EPS;
				negsum[0] = EPS;

				if (norm == L2) {
					for (UINT biter = 0; biter < nbins; ++biter) { // add the three channels
						lbpptr[biter] = lbpfeat[0][biter] + lbpfeat[1][biter]
								+ lbpfeat[2][biter];
						posptr[biter] = pltpfeat[0][biter] + pltpfeat[1][biter]
								+ pltpfeat[2][biter];
						negptr[biter] = nltpfeat[0][biter] + nltpfeat[1][biter]
								+ nltpfeat[2][biter];

						lbpsum[0] += lbpptr[biter] * lbpptr[biter];
						possum[0] += posptr[biter] * posptr[biter];
						negsum[0] += negptr[biter] * negptr[biter];
					}
					lbpsum[0] = sqrt(lbpsum[0]);
					possum[0] = sqrt(possum[0]);
					negsum[0] = sqrt(negsum[0]);
					for (UINT biter = 0; biter < nbins; ++biter) {
						lbpptr[biter] = (lbpptr[biter] / lbpsum[0]);
						posptr[biter] = (posptr[biter] / possum[0]);
						negptr[biter] = (negptr[biter] / negsum[0]);
					}
				} else if (norm == NT_Average) {

					for (UINT biter = 0; biter < nbins; ++biter) { // add the three channels
						lbpptr[biter] = (lbpfeat[0][biter] + lbpfeat[1][biter]
								+ lbpfeat[2][biter]) / 3.0;
						posptr[biter] = (pltpfeat[0][biter]
								+ pltpfeat[1][biter] + pltpfeat[2][biter])
								/ 3.0;
						negptr[biter] = (nltpfeat[0][biter]
								+ nltpfeat[1][biter] + nltpfeat[2][biter])
								/ 3.0;
					}

				} else {
					for (UINT biter = 0; biter < nbins; ++biter) { // add the three channels
						lbpptr[biter] = lbpfeat[0][biter] + lbpfeat[1][biter]
								+ lbpfeat[2][biter];
						posptr[biter] = pltpfeat[0][biter] + pltpfeat[1][biter]
								+ pltpfeat[2][biter];
						negptr[biter] = nltpfeat[0][biter] + nltpfeat[1][biter]
								+ nltpfeat[2][biter];

						lbpsum[0] += lbpptr[biter];
						possum[0] += posptr[biter];
						negsum[0] += negptr[biter];
					}
					if (norm == L1) {
						for (UINT biter = 0; biter < nbins; ++biter) {
							lbpptr[biter] = (lbpptr[biter] / lbpsum[0]);
							posptr[biter] = (posptr[biter] / possum[0]);
							negptr[biter] = (negptr[biter] / negsum[0]);
						}
					} else {
						for (UINT biter = 0; biter < nbins; ++biter) {
							lbpptr[biter] = sqrt(lbpptr[biter] / lbpsum[0]);
							posptr[biter] = sqrt(posptr[biter] / possum[0]);
							negptr[biter] = sqrt(negptr[biter] / negsum[0]);
						}
					}
				}
			}
			for (UINT i = 0; i < 3; ++i) {// Zero out the temporary  histograms...
				fill(lbpfeat[i].begin(), lbpfeat[i].end(), 0);
				fill(pltpfeat[i].begin(), pltpfeat[i].end(), 0);
				fill(nltpfeat[i].begin(), nltpfeat[i].end(), 0);
			}
			offset += 3 * nbins;
		}
	}
}
//#else

void LTPFeaturesRGBCell::GetFeatures(UINT index, int x, int y, int width_,
		int height_, vector<REAL>&feat) {// There is offset of one cell so, input x=0 means x=8 by taking
	//	into account offset and the removed border
	int xcell = x / cell, ycell = y / cell;
	width_ /= cell;
	height_ /= cell;
	lbpfeatures.GetSkippedHOGCells(index, xcell, ycell, width_, height_, skip,
			&feat[0]);
}
void LTPFeaturesRGBCell::GetFeatures(UINT index, int x, int y, int width_,
		int height_, REAL*feat) {// There is offset of one cell so, input x=0 means x=8 by taking
	//	into account offset and the removed border
	int xcell = x / cell, ycell = y / cell;
	width_ /= cell;
	height_ /= cell;
	lbpfeatures.GetSkippedHOGCells(index, xcell, ycell, width_, height_, skip,
			feat);
}
/*
 * Contains only folded HOG's....*/
void LTPFeaturesRGBCell::GetFoldedFeatures(UINT index, int x, int y,
		int twidth, int theight, vector<REAL>&feat) {// There is offset of one cell so, input x=0 means x=8 by taking
	//	into account offset and the removed border
	int xcell = x / cell, ycell = y / cell;
	twidth /= cell;
	theight /= cell;
	lbpfeatures.GetSkippedFoldedHOGCells(index, xcell, ycell, twidth, theight,
			skip, &feat[0]);
}

/*
 * Contains only folded HOG's....*/
void LTPFeaturesRGBCell::GetFoldedFeatures(UINT index, int x, int y,
		int twidth, int theight, REAL*feat) {// There is offset of one cell so, input x=0 means x=8 by taking
	//	into account offset and the removed border
	int xcell = x / cell, ycell = y / cell;
	twidth /= cell;
	theight /= cell;
	lbpfeatures.GetSkippedFoldedHOGCells(index, xcell, ycell, twidth, theight,
			skip, feat);
}

void LTPFeaturesRGBCell::GetFlippedFeatures(UINT index, int x, int y,
		int twidth, int theight, vector<REAL>&feat) {// There is offset of one cell so, input x=0 means x=8 by taking
	//	into account offset and the removed border
	int xcell = x / cell, ycell = y / cell;
	twidth /= cell;
	theight /= cell;
	lbpfeatures.GetSkippedFlippedFeatures(index, xcell, ycell, twidth, theight,
			skip, &feat[0]);

}

void LTPFeaturesRGBCell::GetFlippedFeatures(int twidth, int theight,
		REAL *ifeat, REAL *ofeat) {// There is offset of one cell so, input x=0 means x=8 by taking
	//	into account offset and the removed border
	twidth /= cell;
	theight /= cell;
	HInfo::FlipFeatures(twidth / skip, theight / skip, nbins * 3, ifeat, ofeat);
}

void LTPFeaturesRGBCell::UnFoldWeights(REAL *input, UINT w, UINT h,
		vector<REAL>& weights) {

	UINT tcwidth = w / cell, tcheight = h / cell;
	HInfo hinfo((UINT) ceil((double) tcwidth / (2.0 * skip)),
			(tcheight / skip), 1, nbins * 3, input);
	hinfo.GetUnFoldedHOGCells(tcwidth / skip, &weights[0]);// Flip the features
}

void LTPFeaturesRGBCell::PadFeatureMap(UINT index, UINT padx, UINT pady) { // does padding in the pixel space
	// first pad the hogFeatures map...
	padx /= cell;
	pady /= cell;
	lbpfeatures.PadFeatureMap(index, padx, pady);
}

void LTPFeaturesRGBCell::DotProduct(UINT index, UINT fwidth, UINT fheight,
		vector<REAL>&filter, Store&response) {
	fwidth /= cell;
	fheight /= cell;
	lbpfeatures.SkippedDotProduct(index, fwidth, fheight, skip, &filter[0],
			response);
}
void LTPFeaturesRGBCell::DotProduct(UINT index, UINT fwidth, UINT fheight,
		REAL *filter, Store&response) {
	fwidth /= cell;
	fheight /= cell;
	lbpfeatures.SkippedDotProduct(index, fwidth, fheight, skip, filter,
			response);
}
