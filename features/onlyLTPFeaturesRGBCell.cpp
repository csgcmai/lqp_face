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
#include "onlyLTPFeaturesRGBCell.h"

void OnlyLTPFeaturesRGBCell::InitalizeMaps(Image &image) {
	sspace = NoPyramid;
	UINT xbmax, ybmax;
	lbpfeatures.Initialize(1, 0, celldim * dimmult);
	GetXYBlocks(image, xbmax, ybmax);
	lbpfeatures.SetFeature(xbmax, ybmax, 1.0);
	ComputeLBPFeatures(image, 0, cell, lbpstride);
}
void OnlyLTPFeaturesRGBCell::InitalizeMaps(Pyramid &pyobj, PyramidType sspace_) {
	UINT xbmax, ybmax, pycount = 0;
	sspace = sspace_;
	// Read the image....
	if (sspace == SingleRes) {
		pycount = pyobj.GetNLevels();
		lbpfeatures.Initialize(pycount, 0, celldim * dimmult);
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
				celldim * dimmult);
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
void OnlyLTPFeaturesRGBCell::InitalizeMaps(Image&image, PyramidType sspace_) {
	UINT xbmax, ybmax, pycount = 0;
	sspace = sspace_;
	// Read the image....
	if (sspace == SingleRes) {
		pycount = pyobj.GeneratePyramid(image);
		lbpfeatures.Initialize(pycount, 0, celldim * dimmult);
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
				celldim * dimmult);
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
		lbpfeatures.Initialize(1, 0, celldim * dimmult);
		GetXYBlocks(image, xbmax, ybmax);
		lbpfeatures.SetFeature(xbmax, ybmax, 1.0);
		ComputeLBPFeatures(image, 0, cell, lbpstride);
	}
}
void OnlyLTPFeaturesRGBCell::ComputeLTPMap(Image &image, LBPMap tpmap[],
		LBPMap tnmap[]) {
	UINT xdim = image.columns(), ydim = image.rows();
	int nypoints = ydim - 2 * radius, reg = 1, fx, fy, cx, cy, nxpoints = xdim
			- 2 * radius, cenx, ceny, count, tval, value[3], posvalue[3],
			negvalue[3];
	REAL tmpx, tmpy, fracx, fracy, pixval[3], cenval[3];
	REAL x, y, w1, w2, w3, w4;

	for (int i = 0; i < 3; ++i)// for 3 channels RGB
	{
		tpmap[i].Init(nxpoints, nypoints);
		tnmap[i].Init(nxpoints, nypoints);
	}

	vector<REAL> rimpix(xdim * ydim, 0), gimpix(xdim * ydim, 0), bimpix(
			xdim * ydim, 0);
	pim.Process(image, rimpix, gimpix, bimpix);
	for (int i = 0; i < nypoints; ++i) {
		for (int j = 0; j < nxpoints; ++j) {
			cenx = j + radius;
			ceny = i + radius;

			value[0] = value[1] = value[2] = 0;
			posvalue[0] = posvalue[1] = posvalue[2] = 0;
			negvalue[0] = negvalue[1] = negvalue[2] = 0;

			count = 0;

			cenval[0] = rimpix[ceny * xdim + cenx];
			cenval[1] = gimpix[ceny * xdim + cenx];
			cenval[2] = bimpix[ceny * xdim + cenx];
			for (vector<Point<REAL> >::iterator iter = spoints.begin(); iter
					!= spoints.end(); ++iter, ++count) {
				iter->GetCoord(x, y);
				tmpy = ceny + y;
				tmpx = cenx + x;

				if (ptype == PT_Circular && gridtype == Circular) { // do bilinear interpolation

					fx = floor(tmpx);
					fy = floor(tmpy);
					cx = ceil(tmpx);
					cy = ceil(tmpy);
					fracx = tmpx - fx;
					fracy = tmpy - fy;

					if (ABS(fracx) < 1e-6 && ABS(fracy) < 1e-6) {
						// no interpolation needed
						pixval[0] = rimpix[fy * xdim + fx];
						pixval[1] = gimpix[fy * xdim + fx];
						pixval[2] = bimpix[fy * xdim + fx];
					} else {

						w1 = (1 - fracx) * (1 - fracy);
						w2 = fracx * (1 - fracy);
						w3 = (1 - fracx) * fracy;
						w4 = fracx * fracy;

						pixval[0] = w1 * rimpix[fy * xdim + fx] + w2
								* rimpix[fy * xdim + cx] + w3 * rimpix[cy
								* xdim + fx] + w4 * rimpix[cy * xdim + cx];

						pixval[1] = w1 * gimpix[fy * xdim + fx] + w2
								* gimpix[fy * xdim + cx] + w3 * gimpix[cy
								* xdim + fx] + w4 * gimpix[cy * xdim + cx];

						pixval[2] = w1 * bimpix[fy * xdim + fx] + w2
								* bimpix[fy * xdim + cx] + w3 * bimpix[cy
								* xdim + fx] + w4 * bimpix[cy * xdim + cx];
					}
					// posvalue and negvalue contains the thresholded values for r,g & b
					// tcount is used to run over the values...
				} else {
					pixval[0] = rimpix[tmpy * xdim + tmpx];
					pixval[1] = gimpix[tmpy * xdim + tmpx];
					pixval[2] = bimpix[tmpy * xdim + tmpx];
				}

				tval = reg << count;
				for (int iter = 0; iter < 3; ++iter) {
					posvalue[iter]
							+= pixval[iter] >= cenval[iter] + tolerance ? tval
									: 0;
					negvalue[iter]
							+= pixval[iter] <= cenval[iter] - tolerance ? tval
									: 0;
				}
			}
#ifdef TMPDEBUG
			//			cout << " x =" << j << " y =" << i << " pr =" << posvalue[0]
			//					<< " pg =" << posvalue[1] << " pb =" << posvalue[2]
			//					<< " nr =" << negvalue[0] << " ng =" << negvalue[1]
			//					<< " nb =" << negvalue[2];
#endif

			for (int iter = 0; iter < 3; ++iter) {
				tpmap[iter](j, i) = map[posvalue[iter]];
				tnmap[iter](j, i) = map[negvalue[iter]];
#ifdef COUNT_CODES // count uniform codes...
				if (map[posvalue[iter]] == nbins - 1)
				++pnunicodes;
				else
				++punicodes;

				if (map[negvalue[iter]] == nbins - 1)
				++nnunicodes;
				else
				++nunicodes;
#endif
#ifdef TMPDEBUG
				//				cout << " pmap [" << iter << "]=" << tpmap[iter](j, i);
				//				cout << " nmap [" << iter << "]=" << tnmap[iter](j, i) << endl
				//						<< flush;
#endif
			}
		}
	}
}

void OnlyLTPFeaturesRGBCell::ComputeLBPFeatures(Image &image, UINT index,
		UINT sbin, UINT lbpstride) {
	LBPMap tposmap[3], tnegmap[3];
	ComputeLTPMap(image, tposmap, tnegmap);

	if (hmethod == HM_Bilinear) {
		vector<REAL> &tvec = lbpfeatures.GetFeatureRef(index);
		REAL *feat = &tvec[0];
		ComputeHistBilinear(image, sbin, tposmap, tnegmap, feat);
	} else
		ComputeHistDiscrete(image, index, sbin, tposmap, tnegmap);
}

void OnlyLTPFeaturesRGBCell::ComputeHistDiscrete(Image& imgref, UINT index,
		UINT winstride, LBPMap *tposmap, LBPMap *tnegmap) {
	// for the speed modularization is sacrificed...

	UINT x, y, xbmax = lbpfeatures.GetXBlock(index) * winstride, ybmax =
			lbpfeatures.GetYBlock(index) * winstride, offset = 0, endx, endy;

	vector<REAL> &features = lbpfeatures.GetFeatureRef(index);
	vector<vector<REAL> > pltpfeat(3, vector<REAL> (nbins)), nltpfeat(3,
			vector<REAL> (nbins));
	REAL possum[3], negsum[3]; // normalization constants...
	//	fill(features.begin(), features.end(), 0);

	for (UINT rind = 0; rind < ybmax; rind += winstride) {
		y = rind + winstride/*foffset*/- radius;
		for (UINT colind = 0; colind < xbmax; colind += winstride) {
			x = colind + winstride /*foffset*/- radius;

			endx = MIN(x+winstride,tposmap[0].nxpoints);
			endy = MIN(y+winstride,tposmap[0].nypoints);

			for (UINT yiter = y; yiter < endy; ++yiter)
				for (UINT xiter = x; xiter < endx; ++xiter) {

					pltpfeat[0][tposmap[0].GetValue(xiter, yiter)]++;
					pltpfeat[1][tposmap[1].GetValue(xiter, yiter)]++;
					pltpfeat[2][tposmap[2].GetValue(xiter, yiter)]++;

					nltpfeat[0][tnegmap[0].GetValue(xiter, yiter)]++;
					nltpfeat[1][tnegmap[1].GetValue(xiter, yiter)]++;
					nltpfeat[2][tnegmap[2].GetValue(xiter, yiter)]++;

				}
			if (!add) { // don't add but concatenate the RGB channels....
				// normalize each channel separately and then Concatenate....
				possum[0] = possum[1] = possum[2] = EPS;
				negsum[0] = negsum[1] = negsum[2] = EPS;

				// Currently only doing L1-Sqrt Normalization further to be added...
				for (UINT biter = 0; biter < nbins; ++biter) {

					possum[0] += pltpfeat[0][biter];
					possum[1] += pltpfeat[1][biter];
					possum[2] += pltpfeat[2][biter];

					negsum[0] += nltpfeat[0][biter];
					negsum[1] += nltpfeat[1][biter];
					negsum[2] += nltpfeat[2][biter];
				}
				REAL *posptr = &features[offset], *negptr = &features[offset
						+ 3 * nbins];
				for (UINT biter = 0; biter < nbins; ++biter) {

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
					possum[0] = possum[1] = possum[2] = EPS;
					negsum[0] = negsum[1] = negsum[2] = EPS;

					// Currently only doing L1-Sqrt Normalization further to be added...
					for (UINT biter = 0; biter < nbins; ++biter) {

						possum[0] += pltpfeat[0][biter];
						possum[1] += pltpfeat[1][biter];
						possum[2] += pltpfeat[2][biter];

						negsum[0] += nltpfeat[0][biter];
						negsum[1] += nltpfeat[1][biter];
						negsum[2] += nltpfeat[2][biter];
					}

					for (UINT biter = 0; biter < nbins; ++biter) {

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
						features[offset + biter] = pltpfeat[0][biter]
								+ pltpfeat[1][biter] + pltpfeat[2][biter];
						features[nbins + offset + biter] = nltpfeat[0][biter]
								+ nltpfeat[1][biter] + nltpfeat[2][biter];
					}
				} else {

					REAL *posptr = &features[offset], *negptr =
							&features[offset + nbins];

					possum[0] = EPS;
					negsum[0] = EPS;

					for (UINT biter = 0; biter < nbins; ++biter) { // add the three channels
						posptr[biter] = pltpfeat[0][biter] + pltpfeat[1][biter]
								+ pltpfeat[2][biter];
						negptr[biter] = nltpfeat[0][biter] + nltpfeat[1][biter]
								+ nltpfeat[2][biter];

						possum[0] += posptr[biter];
						negsum[0] += negptr[biter];
					}

					for (UINT biter = 0; biter < nbins; ++biter) {
						posptr[biter] = sqrt(posptr[biter] / possum[0]);
						negptr[biter] = sqrt(negptr[biter] / negsum[0]);
					}

				}
			}

			for (UINT i = 0; i < 3; ++i) {// Zero out the temporary  histograms...
				fill(pltpfeat[i].begin(), pltpfeat[i].end(), 0);
				fill(nltpfeat[i].begin(), nltpfeat[i].end(), 0);
			}
			offset += 2 * nbins;
		}
	}

}
void OnlyLTPFeaturesRGBCell::ComputeHistBilinear(Image &image, UINT sbin,
		LBPMap tposmap[3], LBPMap tnegmap[3], REAL *feat) {
	int dims[2];
	dims[0] = image.rows();
	dims[1] = image.columns();

	// memory for caching orientation histograms & their norms
	int blocks[2], histdim, bdim;
	blocks[0] = (int) round((double) dims[0] / (double) sbin);
	blocks[1] = (int) round((double) dims[1] / (double) sbin);
	histdim = blocks[0] * blocks[1] * clbpdim;
	bdim = blocks[0] * blocks[1];
	double *prhist = new double[histdim], *pghist = new double[(histdim)],
			*pbhist = new double[histdim], *phist = new double[(histdim)];
	double *nrhist = new double[histdim], *nghist = new double[(histdim)],
			*nbhist = new double[histdim], *nhist = new double[(histdim)];
	double *pnorm = new double[(bdim)], *nnorm = new double[(bdim)];
	fill(pnorm, pnorm + bdim, 0);
	fill(nnorm, nnorm + bdim, 0);

	fill(prhist, prhist + histdim, 0);
	fill(pghist, pghist + histdim, 0);
	fill(pbhist, pbhist + histdim, 0);
	fill(phist, phist + histdim, 0); // for pooling of RGB channels....

	fill(nrhist, nrhist + histdim, 0);
	fill(nghist, nghist + histdim, 0);
	fill(nbhist, nbhist + histdim, 0);
	fill(nhist, nhist + histdim, 0); // for pooling of RGB channels....

	// memory for HOG features
	int out[2];
	out[0] = max(blocks[0] - 2, 0);
	out[1] = max(blocks[1] - 2, 0);

	int visible[2];
	visible[0] = blocks[0] * sbin;
	visible[1] = blocks[1] * sbin;

	UINT width = image.columns(), height = image.rows(), mx, my;
	const PixelPacket *pix = image.getConstPixels(0, 0, width, height);
#ifdef INVERSE
	for (int x = radius; x < visible[1] - radius; x++) {
		for (int y = radius; y < visible[0] - radius; y++) {

			// first color channel
			/*because the value of map is offset by radius, value of x at x-radius */
			/*-1 because of < ; MIN(x, dims[1] - radius-1) for which lbp map is computed;
			 * final -radius the offset in lbpmap*/
			mx = MIN(x, dims[1] - radius-1);
			my = MIN(y, dims[0] - radius-1);
			tx = mx - radius;

			ty = my - radius;

#else
	UINT nxpoints = tposmap[0].nxpoints, nypoints = tposmap[0].nypoints;
	for (int tx = 0; tx < nxpoints; ++tx) {
		for (int ty = 0; ty < nypoints; ++ty) {

			// first color channel
			/*because the value of map is offset by radius, value of x at x-radius */
			/*-1 because of < ; MIN(x, dims[1] - radius-1) for which lbp map is computed;
			 * final -radius the offset in lbpmap*/
			int x = tx + radius, mx = x;
			int y = ty + radius, my = y;
#endif
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
			if (usegradmag) {
				// first color channel
				double dy = (pix[mx + (my + 1) * width].red - pix[mx + (my - 1)
						* width].red), dx = (pix[(mx + 1) + my * width].red
						- pix[(mx - 1) + my * width].red),

				dy2 = (pix[mx + (my + 1) * width].green - pix[mx + (my - 1)
						* width].green), dx2 =
						(pix[(mx + 1) + my * width].green - pix[(mx - 1) + my
								* width].green),

				dy3 = (pix[mx + (my + 1) * width].blue - pix[mx + (my - 1)
						* width].blue), dx3 = (pix[(mx + 1) + my * width].blue
						- pix[(mx - 1) + my * width].blue), v2 = dx2 * dx2
						+ dy2 * dy2, v3 = dx3 * dx3 + dy3 * dy3;
				v = dx * dx + dy * dy;
				// pick channel with strongest gradient
				if (v2 > v) {
					v = v2;
					dx = dx2;
					dy = dy2;
				}
				if (v3 > v) {
					v = v3;
					dx = dx3;
					dy = dy3;
				}
				v = sqrt(v);
			}
#ifdef TMPDEBUG
			ixp = x / sbin;
			iyp = y / sbin;
			*(prhist + ixp * blocks[0] + iyp + pro * bdim) += v;
			*(pghist + ixp * blocks[0] + iyp + pgo * bdim) += v;
			*(pbhist + ixp * blocks[0] + iyp + pbo * bdim) += v;

			*(nrhist + ixp * blocks[0] + iyp + nro * bdim) += v;
			*(nghist + ixp * blocks[0] + iyp + ngo * bdim) += v;
			*(nbhist + ixp * blocks[0] + iyp + nbo * bdim) += v;
#else
			if (ixp >= 0 && iyp >= 0) {
				*(prhist + ixp * blocks[0] + iyp + pro * bdim) += vx1 * vy1 * v;
				*(pghist + ixp * blocks[0] + iyp + pgo * bdim) += vx1 * vy1 * v;
				*(pbhist + ixp * blocks[0] + iyp + pbo * bdim) += vx1 * vy1 * v;

				*(nrhist + ixp * blocks[0] + iyp + nro * bdim) += vx1 * vy1 * v;
				*(nghist + ixp * blocks[0] + iyp + ngo * bdim) += vx1 * vy1 * v;
				*(nbhist + ixp * blocks[0] + iyp + nbo * bdim) += vx1 * vy1 * v;

			}

			if (ixp + 1 < blocks[1] && iyp >= 0) {

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
#endif
		}
	}
	// compute energy in each block by summing over orientations
	for (int o = 0; o < clbpdim; o++) {

		double *psrcr = prhist + o * bdim, *psrcg = pghist + o * bdim, *psrcb =
				pbhist + o * bdim, *psrc = phist + o * bdim;

		double *nsrcr = nrhist + o * bdim, *nsrcg = nghist + o * bdim, *nsrcb =
				nbhist + o * bdim, *nsrc = nhist + o * bdim;

		double *pdst = pnorm, *ndst = nnorm;
		double *end = pnorm + blocks[1] * blocks[0];
		while (pdst < end) {
			*psrc = *psrcr + *psrcg + *psrcb; // pooling of rgb channels...
			*nsrc = *nsrcr + *nsrcg + *nsrcb; // pooling of rgb channels...
			if (norm == L2) {
				*(pdst++) += *psrc * *psrc;
				*(ndst++) += *nsrc * *nsrc;
			} else {
				*(pdst++) += *psrc;
				*(ndst++) += *nsrc;
			}

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
	int citer = useblocknorm == true ? clbpdim / 4 : clbpdim;
	for (int x = 0; x < out[1]; x++) {
		for (int y = 0; y < out[0]; y++) {
			REAL *pdst = feat + celldim * (x + y * out[1]); // write in column major order
			//			REAL *dst = feat + celldim * dimmult * (x + y * out[1]); // write in column major order
			REAL *ndst = pdst + clbpdim; //	clbpdim = nbins * (useblocknorm == true ? 4 : 1);
			double *p, pn1, pn2, pn3, pn4, nn1, nn2, nn3, nn4;
			double *psrc = phist + (x + 1) * blocks[0] + (y + 1);
			double *nsrc = nhist + (x + 1) * blocks[0] + (y + 1);

			if (useblocknorm) {/*normlaize by lbp energy of four neighbouring blocks*/

				// positive normalization
				p = pnorm + (x + 1) * blocks[0] + y + 1;
				pn1 = 1.0 / (*p + *(p + 1) + *(p + blocks[0]) + *(p + blocks[0]
						+ 1) + EPS);
				p = pnorm + (x + 1) * blocks[0] + y;
				pn2 = 1.0 / (*p + *(p + 1) + *(p + blocks[0]) + *(p + blocks[0]
						+ 1) + EPS);
				p = pnorm + x * blocks[0] + y + 1;
				pn3 = 1.0 / (*p + *(p + 1) + *(p + blocks[0]) + *(p + blocks[0]
						+ 1) + EPS);
				p = pnorm + x * blocks[0] + y;
				pn4 = 1.0 / (*p + *(p + 1) + *(p + blocks[0]) + *(p + blocks[0]
						+ 1) + EPS);

				// negative normalization
				p = nnorm + (x + 1) * blocks[0] + y + 1;
				nn1 = 1.0 / (*p + *(p + 1) + *(p + blocks[0]) + *(p + blocks[0]
						+ 1) + EPS);
				p = nnorm + (x + 1) * blocks[0] + y;
				nn2 = 1.0 / (*p + *(p + 1) + *(p + blocks[0]) + *(p + blocks[0]
						+ 1) + EPS);
				p = nnorm + x * blocks[0] + y + 1;
				nn3 = 1.0 / (*p + *(p + 1) + *(p + blocks[0]) + *(p + blocks[0]
						+ 1) + EPS);
				p = nnorm + x * blocks[0] + y;
				nn4 = 1.0 / (*p + *(p + 1) + *(p + blocks[0]) + *(p + blocks[0]
						+ 1) + EPS);

				for (int o = 0; o < citer; o++) {// 3rd dimension processing
					double psum = *psrc, nsum = *nsrc;
					double pf1 = sqrt(psum * pn1);
					double pf2 = sqrt(psum * pn2);
					double pf3 = sqrt(psum * pn3);
					double pf4 = sqrt(psum * pn4);
					*(pdst + o) = pf4;
					*(pdst + citer + o) = pf2;
					*(pdst + citer * 2 + o) = pf3;
					*(pdst + citer * 3 + o) = pf1;

					double nf1 = sqrt(nsum * nn1);
					double nf2 = sqrt(nsum * nn2);
					double nf3 = sqrt(nsum * nn3);
					double nf4 = sqrt(nsum * nn4);
					*(ndst + o) = nf4;
					*(ndst + citer + o) = nf2;
					*(ndst + citer * 2 + o) = nf3;
					*(ndst + citer * 3 + o) = nf1;

					psrc += bdim;
					nsrc += bdim;
					//					++pdst;
					//					++ndst;
				}

				//				for (int o = 0; o < citer; o++) {
				//					double sum = *src;
				//					double h1 = sqrt(sum * n1);
				//					double h2 = sqrt(sum * n2);
				//					double h3 = sqrt(sum * n3);
				//					double h4 = sqrt(sum * n4);
				//					*(dst + o) = h4;
				//					*(dst + celldim + o) = h2;
				//					*(dst + celldim * 2 + o) = h3;
				//					*(dst + celldim * 3 + o) = h1;
				//					src += bdim;
				//				}
			} else {
				// single cell normalization...
				if (norm == L2) {
					REAL pnormval = 1 / sqrt(
							*(pnorm + (x + 1) * blocks[0] + y + 1) + EPS);
					REAL nnormval = 1 / sqrt(
							*(nnorm + (x + 1) * blocks[0] + y + 1) + EPS);
					for (int o = 0; o < clbpdim; o++) {
						*pdst = *psrc * pnormval;
						*ndst = *nsrc * nnormval;

						psrc += bdim;
						nsrc += bdim;

						++pdst;
						++ndst;
					}
				} else if (norm == LOG2) {

					for (int o = 0; o < nbins; o++) {
						*pdst = log2(*psrc + 1);
						*ndst = log2(*nsrc + 1);
						psrc += bdim;
						nsrc += bdim;
						++pdst;
						++ndst;
					}
				} else {
					REAL pnormval = *(pnorm + (x + 1) * blocks[0] + y + 1)
							+ EPS;
					REAL nnormval = *(nnorm + (x + 1) * blocks[0] + y + 1)
							+ EPS;
					for (int o = 0; o < clbpdim; o++) {
						*pdst = sqrt(*psrc / pnormval);
						*ndst = sqrt(*nsrc / nnormval);
						psrc += bdim;
						nsrc += bdim;
						++pdst;
						++ndst;
					}
				}
			}
		}
	}

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
void OnlyLTPFeaturesRGBCell::GetFlippedFeatures(int twidth, int theight,
		REAL *ifeat, REAL *ofeat) {// There is offset of one cell so, input x=0 means x=8 by taking //	into account offset and the removed border
	twidth /= cell;
	theight /= cell;
	HInfo::FlipFeatures(twidth / skip, theight / skip, celldim, ifeat, ofeat);
}

void OnlyLTPFeaturesRGBCell::UnFoldWeights(REAL *input, UINT w, UINT h,
		vector<REAL>& weights) {

	UINT tcwidth = w / cell, tcheight = h / cell;
	HInfo hinfo((UINT) ceil((double) tcwidth / (2.0 * skip)),
			(tcheight / skip), 1, celldim, input);
	hinfo.GetUnFoldedHOGCells(tcwidth / skip, &weights[0]);// Flip the features
}

void OnlyLTPFeaturesRGBCell::ComputeTextureFeatures(Image &imgref,
		UINT xcellsize, UINT ycellsize, vector<REAL>&features) {
	ComputeTextureFeatures(imgref, xcellsize, ycellsize, 0, 0,
			imgref.columns() - 1, imgref.rows() - 1, features);
}

void OnlyLTPFeaturesRGBCell::ComputeTextureFeatures(Image &imgref, int sx,
		int sy, int ex, int ey, UINT xcellsize, UINT ycellsize,
		vector<REAL>&features) {
	LBPMap tposmap[3], tnegmap[3];
	ComputeLTPMap(imgref, tposmap, tnegmap);

	UINT width = ex - sx + 1, height = ey - sy + 1;
	int xoffset = -radius/*+MAX(ceil((REAL) (width % xcellsize) / 2),radius)*/;
	int yoffset = -radius/*+MAX(ceil((REAL) (height % ycellsize) / 2),radius)*/;

	UINT xbmax = width / xcellsize, ybmax = height / ycellsize;
	UINT x, y, offset = 0, endx, endy;

	features.resize(xbmax * ybmax * celldim * dimmult, 0);

	vector<vector<REAL> > pltpfeat(3, vector<REAL> (nbins)), nltpfeat(3,
			vector<REAL> (nbins));
	REAL possum[3], negsum[3]; // normalization constants...

	for (UINT rind = 0; rind < ybmax; rind += ycellsize) {
		for (UINT colind = 0; colind < xbmax; colind += xcellsize) {

			for (int yiter = rind; yiter < rind + ycellsize; ++yiter) {
				int y = MIN(MAX(yiter+yoffset+sy,0),tposmap[0].nypoints-1);
				for (int xiter = colind; xiter < colind + xcellsize; ++xiter) {
					int x = MIN(MAX(xiter+xoffset+sx,0),tposmap[0].nxpoints-1);
					assert(tposmap[0].GetValue(x, y)!=-1);

					pltpfeat[0][tposmap[0].GetValue(x, y)]++;
					pltpfeat[1][tposmap[1].GetValue(x, y)]++;
					pltpfeat[2][tposmap[2].GetValue(x, y)]++;

					nltpfeat[0][tnegmap[0].GetValue(x, y)]++;
					nltpfeat[1][tnegmap[1].GetValue(x, y)]++;
					nltpfeat[2][tnegmap[2].GetValue(x, y)]++;

				}
			}
			{// normalization

				REAL *posptr = &features[offset], *negptr = &features[offset
						+ nbins];

				possum[0] = EPS;
				negsum[0] = EPS;

				if (norm == L2) {
					for (UINT biter = 0; biter < nbins; ++biter) { // add the three channels
						posptr[biter] = pltpfeat[0][biter] + pltpfeat[1][biter]
								+ pltpfeat[2][biter];
						negptr[biter] = nltpfeat[0][biter] + nltpfeat[1][biter]
								+ nltpfeat[2][biter];

						possum[0] += posptr[biter] * posptr[biter];
						negsum[0] += negptr[biter] * negptr[biter];
					}
					possum[0] = sqrt(possum[0]);
					negsum[0] = sqrt(negsum[0]);
					for (UINT biter = 0; biter < nbins; ++biter) {
						posptr[biter] = (posptr[biter] / possum[0]);
						negptr[biter] = (negptr[biter] / negsum[0]);
					}
				} else if (norm == NT_Average) {

					for (UINT biter = 0; biter < nbins; ++biter) { // add the three channels
						posptr[biter] = (pltpfeat[0][biter]
								+ pltpfeat[1][biter] + pltpfeat[2][biter])
								/ 3.0;
						negptr[biter] = (nltpfeat[0][biter]
								+ nltpfeat[1][biter] + nltpfeat[2][biter])
								/ 3.0;
					}

				} else {
					for (UINT biter = 0; biter < nbins; ++biter) { // add the three channels
						posptr[biter] = pltpfeat[0][biter] + pltpfeat[1][biter]
								+ pltpfeat[2][biter];
						negptr[biter] = nltpfeat[0][biter] + nltpfeat[1][biter]
								+ nltpfeat[2][biter];

						possum[0] += posptr[biter];
						negsum[0] += negptr[biter];
					}
					if (norm == L1) {
						for (UINT biter = 0; biter < nbins; ++biter) {
							posptr[biter] = (posptr[biter] / possum[0]);
							negptr[biter] = (negptr[biter] / negsum[0]);
						}
					} else {
						for (UINT biter = 0; biter < nbins; ++biter) {
							posptr[biter] = sqrt(posptr[biter] / possum[0]);
							negptr[biter] = sqrt(negptr[biter] / negsum[0]);
						}
					}

				}
			}
			for (UINT i = 0; i < 3; ++i) {// Zero out the temporary  histograms...
				fill(pltpfeat[i].begin(), pltpfeat[i].end(), 0);
				fill(nltpfeat[i].begin(), nltpfeat[i].end(), 0);
			}
			offset += 2 * nbins;
		}
	}
}

#ifdef MULTI_THRESH_LEVELS
void OnlyLTPFeaturesRGBCell::ComputeLBPFeatures(Image& imgref, UINT index,
		UINT winstride, UINT tlbpstride) {
	// for the speed modularization is sacrificed...

	//	LBPMap tposmap[3], tnegmap[3];
	//	lbpmap2d tposmap(nthrlevels, vector<LBPMap> (3)), tnegmap(nthrlevels,
	//			vector<LBPMap> (3));
	LBPMap *tposmap = new LBPMap[3 * nthrlevels], *tnegmap = new LBPMap[3
	* nthrlevels];
	//#ifdef FAST_FEATURES
	//TODO: Write LTPMap Fast function for this
	//	ComputeLTPMapFast(imgref, tlbpmap, tposmap, tnegmap);
	//#else
	ComputeLTPMap(imgref, tposmap, tnegmap);
	//#endif

	UINT x, y, xbmax = lbpfeatures.GetXBlock(index) * winstride, ybmax =
	lbpfeatures.GetYBlock(index) * winstride, offset = 0, endx, endy;

	vector<REAL> &features = lbpfeatures.GetFeatureRef(index);
	vector<vector<REAL> > pltpfeat(3, vector<REAL> (nbins)), nltpfeat(3,
			vector<REAL> (nbins));
	//	REAL *possum = new REAL[3 * nthrlevels], *negsum = new REAL[3 * nthrlevels]; // normalization constants...
	REAL possum[3], negsum[3];

	for (UINT rind = 0; rind < ybmax; rind += winstride) {
		y = rind + winstride/*foffset*/- radius;
		for (UINT colind = 0; colind < xbmax; colind += winstride) {
			x = colind + winstride /*foffset*/- radius;

			endx = MIN(x + tlbpstride, tposmap[0].nxpoints);
			endy = MIN(y + tlbpstride, tposmap[0].nypoints);
			for (int tl = 0; tl < 3 * nthrlevels; tl += 3) /*Remove Explicit Vectorization*/{
				for (UINT yiter = y; yiter < endy; ++yiter)
				for (UINT xiter = x; xiter < endx; ++xiter) {

					pltpfeat[0][tposmap[tl].GetValue(xiter, yiter)]++;
					pltpfeat[1][tposmap[tl + 1].GetValue(xiter, yiter)]++;
					pltpfeat[2][tposmap[tl + 2].GetValue(xiter, yiter)]++;

					nltpfeat[0][tnegmap[tl].GetValue(xiter, yiter)]++;
					nltpfeat[1][tnegmap[tl + 1].GetValue(xiter, yiter)]++;
					nltpfeat[2][tnegmap[tl + 2].GetValue(xiter, yiter)]++;
				}

				REAL *posptr = &features[offset], *negptr = &features[offset
				+ nbins];

				possum[0] = EPS;
				negsum[0] = EPS;
				if (norm == L2) {

					for (UINT biter = 0; biter < nbins; ++biter) { // add the three channels
						//							lbpptr[biter] = lbpfeat[0][biter]
						//									+ lbpfeat[1][biter] + lbpfeat[2][biter];
						posptr[biter] = pltpfeat[0][biter] + pltpfeat[1][biter]
						+ pltpfeat[2][biter];
						negptr[biter] = nltpfeat[0][biter] + nltpfeat[1][biter]
						+ nltpfeat[2][biter];

						possum[0] += posptr[biter] * posptr[biter];
						negsum[0] += negptr[biter] * negptr[biter];
					}
					possum[0] = sqrt(possum[0]);
					negsum[0] = sqrt(negsum[0]);
					for (UINT biter = 0; biter < nbins; ++biter) {
						posptr[biter] = posptr[biter] / possum[0];
						negptr[biter] = negptr[biter] / negsum[0];
					}

				} else {
					for (UINT biter = 0; biter < nbins; ++biter) { // add the three channels
						//							lbpptr[biter] = lbpfeat[0][biter]
						//									+ lbpfeat[1][biter] + lbpfeat[2][biter];
						posptr[biter] = pltpfeat[0][biter] + pltpfeat[1][biter]
						+ pltpfeat[2][biter];
						negptr[biter] = nltpfeat[0][biter] + nltpfeat[1][biter]
						+ nltpfeat[2][biter];

						//							lbpsum[0] += lbpptr[biter];
						possum[0] += posptr[biter];
						negsum[0] += negptr[biter];
					}

					for (UINT biter = 0; biter < nbins; ++biter) {
						//  							lbpptr[biter] = sqrt(lbpptr[biter] / lbpsum[0]);
						posptr[biter] = sqrt(posptr[biter] / possum[0]);
						negptr[biter] = sqrt(negptr[biter] / negsum[0]);
					}
				}
				for (UINT i = 0; i < 3; ++i) {// Zero out the temporary  histograms...
					fill(pltpfeat[i].begin(), pltpfeat[i].end(), 0);
					fill(nltpfeat[i].begin(), nltpfeat[i].end(), 0);
				}
				offset += 2 * nbins;
			}
		}

	}
	delete[] tposmap;
	delete[] tnegmap;
}

// With Multiple Threshold levels
void OnlyLTPFeaturesRGBCell::ComputeLBPFeatures(Image &image, UINT index,
		UINT sbin) {
	LBPMap tlbpmap[3], tposmap[3], tnegmap[3];
	ComputeLTPMap(imgref, tlbpmap, tposmap, tnegmap);
	vector<REAL> &tvec = lbpfeatures.GetFeatureRef(index);
	REAL *feat = &tvec[0];
	int dims[2];
	dims[0] = image.rows();
	dims[1] = image.columns();

	// memory for caching orientation histograms & their norms
	int blocks[2], histdim, bdim;
	blocks[0] = (int) round((double) dims[0] / (double) sbin);
	blocks[1] = (int) round((double) dims[1] / (double) sbin);
	histdim = blocks[0] * blocks[1] * nbins;
	bdim = blocks[0] * blocks[1];

	double *prhist = new double *[nthrlevels], *pghist =
	new double *[nthrlevels], *pbhist = new double *[nthrlevels],
	*phist = new double *[nthrlevels];
	double *nrhist = new double *[nthrlevels], *nghist =
	new double *[nthrlevels], *nbhist = new double *[nthrlevels],
	*nhist = new double *[nthrlevels];
	double *pnorm = new double*[bdim], *nnorm = new double*[bdim];

	for (UINT i = 0; i < nthrlevels; ++i) {
		prhist[i] = new double[(histdim)];
		pghist[i] = new double[(histdim)];
		pbhist[i] = new double[(histdim)];
		phist[i] = new double[(histdim)];

		nrhist[i] = new double[(histdim)];
		nghist[i] = new double[(histdim)];
		nbhist[i] = new double[(histdim)];
		nhist[i] = new double[(histdim)];

		pnorm[i] = new double[(bdim)];
		nnorm[i] = new double[(bdim)];

		fill(pnorm[i], pnorm[i] + bdim, 0);
		fill(nnorm[i], nnorm[i] + bdim, 0);
		fill(prhist[i], prhist[i] + histdim, 0);
		fill(pghist[i], pghist[i] + histdim, 0);
		fill(pbhist[i], pbhist[i] + histdim, 0);
		fill(phist[i], phist[i] + histdim, 0); // for pooling of RGB channels....
		fill(nrhist[i], nrhist[i] + histdim, 0);
		fill(nghist[i], nghist[i] + histdim, 0);
		fill(nbhist[i], nbhist[i] + histdim, 0);
		fill(nhist[i], nhist[i] + histdim, 0); // for pooling of RGB channels....
	}
	// memory for HOG features
	int out[2];
	out[0] = max(blocks[0] - 2, 0);
	out[1] = max(blocks[1] - 2, 0);

	int visible[2];
	visible[0] = blocks[0] * sbin;
	visible[1] = blocks[1] * sbin;
	UINT width = image.columns(), height = image.rows(), tx, ty;
	for (int x = radius; x < visible[1] - radius; x++) {
		for (int y = radius; y < visible[0] - radius; y++) {

			// first color channel
			/*because the value of map is offset by radius, value of x at x-radius */
			/*-1 because of < ; MIN(x, dims[1] - radius-1) for which lbp map is computed;
			 * final -radius the offset in lbpmap*/
			tx = MIN(x, dims[1] - radius-1) - radius;
			ty = MIN(y, dims[0] - radius-1) - radius;

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

			for (UINT i = 0; i < nthrlevels; ++i) // for each level compute histogram
			{
				int pro = tposmap[i][0].GetValue(tx, ty), // red, GREEN & BLUE
				pgo = tposmap[i][1].GetValue(tx, ty), //
				pbo = tposmap[i][2].GetValue(tx, ty); //
				int nro = tnegmap[i][0].GetValue(tx, ty), // red, GREEN & BLUE
				ngo = tnegmap[i][1].GetValue(tx, ty), //
				nbo = tnegmap[i][2].GetValue(tx, ty); //

				if (ixp >= 0 && iyp >= 0) {
					*(prhist[i] + ixp * blocks[0] + iyp + pro * bdim) += vx1
					* vy1 * v;
					*(pghist[i] + ixp * blocks[0] + iyp + pgo * bdim) += vx1
					* vy1 * v;
					*(pbhist[i] + ixp * blocks[0] + iyp + pbo * bdim) += vx1
					* vy1 * v;

					*(nrhist[i] + ixp * blocks[0] + iyp + nro * bdim) += vx1
					* vy1 * v;
					*(nghist[i] + ixp * blocks[0] + iyp + ngo * bdim) += vx1
					* vy1 * v;
					*(nbhist[i] + ixp * blocks[0] + iyp + nbo * bdim) += vx1
					* vy1 * v;

				}

				if (ixp + 1 < blocks[1] && iyp >= 0) {

					*(prhist[i] + (ixp + 1) * blocks[0] + iyp + pro * bdim)
					+= vx0 * vy1 * v;
					*(pghist[i] + (ixp + 1) * blocks[0] + iyp + pgo * bdim)
					+= vx0 * vy1 * v;
					*(pbhist[i] + (ixp + 1) * blocks[0] + iyp + pbo * bdim)
					+= vx0 * vy1 * v;

					*(nrhist[i] + (ixp + 1) * blocks[0] + iyp + nro * bdim)
					+= vx0 * vy1 * v;
					*(nghist[i] + (ixp + 1) * blocks[0] + iyp + ngo * bdim)
					+= vx0 * vy1 * v;
					*(nbhist[i] + (ixp + 1) * blocks[0] + iyp + nbo * bdim)
					+= vx0 * vy1 * v;
				}

				if (ixp >= 0 && iyp + 1 < blocks[0]) {

					*(prhist[i] + ixp * blocks[0] + (iyp + 1) + pro * bdim)
					+= vx1 * vy0 * v;
					*(pghist[i] + ixp * blocks[0] + (iyp + 1) + pgo * bdim)
					+= vx1 * vy0 * v;
					*(pbhist[i] + ixp * blocks[0] + (iyp + 1) + pbo * bdim)
					+= vx1 * vy0 * v;

					*(nrhist[i] + ixp * blocks[0] + (iyp + 1) + nro * bdim)
					+= vx1 * vy0 * v;
					*(nghist[i] + ixp * blocks[0] + (iyp + 1) + ngo * bdim)
					+= vx1 * vy0 * v;
					*(nbhist[i] + ixp * blocks[0] + (iyp + 1) + nbo * bdim)
					+= vx1 * vy0 * v;

				}

				if (ixp + 1 < blocks[1] && iyp + 1 < blocks[0]) {

					*(prhist[i] + (ixp + 1) * blocks[0] + (iyp + 1) + pro
							* bdim) += vx0 * vy0 * v;
					*(pghist[i] + (ixp + 1) * blocks[0] + (iyp + 1) + pgo
							* bdim) += vx0 * vy0 * v;
					*(pbhist[i] + (ixp + 1) * blocks[0] + (iyp + 1) + pbo
							* bdim) += vx0 * vy0 * v;

					*(nrhist[i] + (ixp + 1) * blocks[0] + (iyp + 1) + nro
							* bdim) += vx0 * vy0 * v;
					*(nghist[i] + (ixp + 1) * blocks[0] + (iyp + 1) + ngo
							* bdim) += vx0 * vy0 * v;
					*(nbhist[i] + (ixp + 1) * blocks[0] + (iyp + 1) + nbo
							* bdim) += vx0 * vy0 * v;

				}
			}
		}
	}
	// compute energy in each block by summing over orientations
	for (UINT i = 0; i < nthrlevels; ++i)
	for (int o = 0; o < nbins; o++) {

		double *psrcr = prhist[i] + o * bdim,
		*psrcg = pghist[i] + o * bdim, *psrcb = pbhist[i] + o
		* bdim, *psrc = phist[i] + o * bdim;

		double *nsrcr = nrhist[i] + o * bdim,
		*nsrcg = nghist[i] + o * bdim, *nsrcb = nbhist[i] + o
		* bdim, *nsrc = nhist[i] + o * bdim;

		double *pdst = pnorm[i], *ndst = *nnorm[i];
		double *end = pnorm[i] + blocks[1] * blocks[0];
		while (pdst < end) {
			*psrc = *psrcr + *psrcg + *psrcb; // pooling of rgb channels...
			*nsrc = *nsrcr + *nsrcg + *nsrcb; // pooling of rgb channels...
			if (norm == L2) {
				*(pdst++) += *psrc + *psrc;
				*(ndst++) += *nsrc + *nsrc;
			} else {
				*(pdst++) += *psrc;
				*(ndst++) += *nsrc;
			}

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
			REAL *pdst = feat + nbins * 2 * (x + y * out[1]); // write in column major order
			REAL *ndst = pdst + nbins;
			double *p, n1, n2, n3, n4;
			double *psrc = phist + (x + 1) * blocks[0] + (y + 1);
			double *nsrc = nhist + (x + 1) * blocks[0] + (y + 1);

			// single cell normalization...
			if (norm == L2) {
				REAL pnormval = 1 / sqrt(*(pnorm + (x + 1) * blocks[0] + y + 1)
						+ EPS);
				REAL nnormval = 1 / sqrt(*(nnorm + (x + 1) * blocks[0] + y + 1)
						+ EPS);
				for (int o = 0; o < nbins; o++) {
					*pdst = *psrc * pnormval;
					*ndst = *nsrc * nnormval;

					psrc += bdim;
					nsrc += bdim;

					++pdst;
					++ndst;
				}
			} else {
				REAL pnormval = *(pnorm + (x + 1) * blocks[0] + y + 1) + EPS;
				REAL nnormval = (nnorm + (x + 1) * blocks[0] + y + 1) + EPS;
				for (int o = 0; o < nbins; o++) {
					*pdst = sqrt(*psrc / pnormval);
					*ndst = sqrt(*nsrc / nnormval);
					psrc += bdim;
					nsrc += bdim;
					++pdst;
					++ndst;
				}
			}
		}
	}

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
/*********Multiple Threshold Levels****/

void OnlyLTPFeaturesRGBCell::InitalizeMaps(Image &image) {
	sspace = NoPyramid;
	UINT xbmax, ybmax;
	lbpfeatures.Initialize(1, 0, nbins * 2 * nthrlevels * dimmult);
	GetXYBlocks(image, xbmax, ybmax);
	lbpfeatures.SetFeature(xbmax, ybmax, 1.0);
	ComputeLBPFeatures(image, 0, cell, lbpstride);
}
void OnlyLTPFeaturesRGBCell::InitalizeMaps(Pyramid &pyobj, PyramidType sspace_) {
	UINT xbmax, ybmax, pycount = 0;
	sspace = sspace_;
	// Read the image....
	if (sspace == SingleRes) {
		pycount = pyobj.GetNLevels();
		lbpfeatures.Initialize(pycount, 0, nbins * 2 * nthrlevels * dimmult);
		for (UINT i = 0; i < pycount; i++) {
			Image& imgref = pyobj.GetImageRef(i);
			// Compute the Gradient Image
			GetXYBlocks(imgref, xbmax, ybmax);
			lbpfeatures.SetFeature(xbmax, ybmax, pyobj.GetScale(i));
			ComputeLBPFeatures(imgref, i, cell, lbpstride);
		}
	} else /*if (sspace == DoubleRes)*/{
		cout << "ERROR: Only LTP Not Implemented In Detail; ";
		exit(EXIT_FAILURE);

		/*		pycount = pyobj.GeneratePyramid(image);
		 UINT nlpyramid = pyobj.GetInitIndex();
		 lbpfeatures.Initialize(pycount + nlpyramid, nlpyramid, nbins * 2
		 * nthrlevels * dimmult);
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
		 lbpfeatures.SetFeature(i + nlpyramid, xbmax, ybmax, pyobj.GetScale(
		 i));
		 ComputeLBPFeatures(imgref, i + nlpyramid, cell, lbpstride);
		 }*/
	}
}
void OnlyLTPFeaturesRGBCell::InitalizeMaps(Image&image, PyramidType sspace_) {
	UINT xbmax, ybmax, pycount = 0;
	sspace = sspace_;
	// Read the image....
	if (sspace == SingleRes) {
		pycount = pyobj.GeneratePyramid(image);
		lbpfeatures.Initialize(pycount, 0, nbins * 2 * nthrlevels * dimmult);
		for (UINT i = 0; i < pycount; i++) {
			Image& imgref = pyobj.GetImageRef(i);
			// Compute the Gradient Image
			GetXYBlocks(imgref, xbmax, ybmax);
			lbpfeatures.SetFeature(xbmax, ybmax, pyobj.GetScale(i));
			ComputeLBPFeatures(imgref, i, cell, lbpstride);
		}
	} else if (sspace == DoubleRes) {
		cout << "ERROR: Only LTP Not Implemented In Detail; ";
		exit(EXIT_FAILURE);

		/*		pycount = pyobj.GeneratePyramid(image);
		 UINT nlpyramid = pyobj.GetInitIndex();
		 lbpfeatures.Initialize(pycount + nlpyramid, nlpyramid, nbins * 2
		 * nthrlevels * dimmult);
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
		 lbpfeatures.SetFeature(i + nlpyramid, xbmax, ybmax, pyobj.GetScale(
		 i));
		 ComputeLBPFeatures(imgref, i + nlpyramid, cell, lbpstride);
		 }*/
	} else {
		lbpfeatures.Initialize(1, 0, nbins * 2 * nthrlevels * dimmult);
		GetXYBlocks(image, xbmax, ybmax);
		lbpfeatures.SetFeature(xbmax, ybmax, 1.0);
		ComputeLBPFeatures(image, 0, cell, lbpstride);
	}
}
void OnlyLTPFeaturesRGBCell::ComputeLTPMap(Image &image, LBPMap *tpmap,
		LBPMap *tnmap) {

	UINT xdim = image.columns(), ydim = image.rows(), rcode, gcode, bcode;
	int nypoints = ydim - 2 * radius, reg = 1, fx, fy, cx, cy, nxpoints = xdim
	- 2 * radius, cenx, ceny, count, tval;
	//	int (*posvalue)[3] = new int[nthrlevels][3], ();
	int *posvalue = new int[3 * nthrlevels], *negvalue =
	new int[3 * nthrlevels];
	REAL tmpx, tmpy, fracx, fracy, pixval[3], cenval[3];
	REAL x, y, w1, w2, w3, w4;
	for (int i = 0; i < 3 * nthrlevels; ++i)// for 3 channels RGB
	{
		tpmap[i].Init(nxpoints, nypoints);
		tnmap[i].Init(nxpoints, nypoints);
	}
	vector<REAL> rimpix(xdim * ydim, 0), gimpix(xdim * ydim, 0), bimpix(xdim
			* ydim, 0);
	pim.Process(image, rimpix, gimpix, bimpix);
	for (int i = 0; i < nypoints; ++i) {
		for (int j = 0; j < nxpoints; ++j) {
			cenx = j + radius;
			ceny = i + radius;

			//			posvalue[0] = posvalue[1] = posvalue[2] = 0;
			//			negvalue[0] = negvalue[1] = negvalue[2] = 0;
			fill(posvalue, posvalue + nthrlevels * 3, 0);
			fill(negvalue, negvalue + nthrlevels * 3, 0);

			count = 0;
			rcode = gcode = bcode = 0;
			cenval[0] = rimpix[ceny * xdim + cenx];
			cenval[1] = gimpix[ceny * xdim + cenx];
			cenval[2] = bimpix[ceny * xdim + cenx];
			for (vector<Point<REAL> >::iterator iter = spoints.begin(); iter
					!= spoints.end(); ++iter, ++count) {
				iter->GetCoord(x, y);
				tmpy = ceny + y;
				tmpx = cenx + x;

				fx = floor(tmpx);
				fy = floor(tmpy);
				cx = ceil(tmpx);
				cy = ceil(tmpy);
				fracx = tmpx - fx;
				fracy = tmpy - fy;

				if (ABS(fracx) < 1e-6 && ABS(fracy) < 1e-6) {
					// no interpolation needed
					pixval[0] = rimpix[fy * xdim + fx];
					pixval[1] = gimpix[fy * xdim + fx];
					pixval[2] = bimpix[fy * xdim + fx];
				} else {

					w1 = (1 - fracx) * (1 - fracy);
					w2 = fracx * (1 - fracy);
					w3 = (1 - fracx) * fracy;
					w4 = fracx * fracy;

					pixval[0] = w1 * rimpix[fy * xdim + fx] + w2 * rimpix[fy
					* xdim + cx] + w3 * rimpix[cy * xdim + fx] + w4
					* rimpix[cy * xdim + cx];

					pixval[1] = w1 * gimpix[fy * xdim + fx] + w2 * gimpix[fy
					* xdim + cx] + w3 * gimpix[cy * xdim + fx] + w4
					* gimpix[cy * xdim + cx];

					pixval[2] = w1 * bimpix[fy * xdim + fx] + w2 * bimpix[fy
					* xdim + cx] + w3 * bimpix[cy * xdim + fx] + w4
					* bimpix[cy * xdim + cx];
				}

				tval = reg << count;

				// posvalue and negvalue contains the thresholded values for r,g & b
				// tcount is used to run over the values...
				for (int tl = 0, tcount = 0; tl < nthrlevels; ++tl)
				for (int iter = 0; iter < 3; ++iter, ++tcount) {

					if (pixval[iter] >= cenval[iter] + tollevels[tl]) {
						posvalue[tcount] += tval; // pos map for rgb

					} else if (pixval[iter] <= cenval[iter] - tollevels[tl]) {
						negvalue[tcount] += tval; // neg map for rgb

					}
				}
			}
			//			for (int tl = 0, count = 0; tl < nthrlevels; ++tl)
			for (int iter = 0; iter < 3 * nthrlevels; ++iter) {
				tpmap[iter](j, i) = map[posvalue[iter]];
				tnmap[iter](j, i) = map[negvalue[iter]];
			}
		}
	}
	delete[] posvalue;
	delete[] negvalue;
}
void OnlyLTPFeaturesRGBCell::ComputeLBPFeatures(Image& imgref, UINT index,
		UINT winstride, UINT tlbpstride) {
	// for the speed modularization is sacrificed...

	//	LBPMap tposmap[3], tnegmap[3];
	//	lbpmap2d tposmap(nthrlevels, vector<LBPMap> (3)), tnegmap(nthrlevels,
	//			vector<LBPMap> (3));
	LBPMap *tposmap = new LBPMap[3 * nthrlevels], *tnegmap = new LBPMap[3
	* nthrlevels];
	//#ifdef FAST_FEATURES
	//TODO: Write LTPMap Fast function for this
	//	ComputeLTPMapFast(imgref, tlbpmap, tposmap, tnegmap);
	//#else
	ComputeLTPMap(imgref, tposmap, tnegmap);
	//#endif

	UINT x, y, xbmax = lbpfeatures.GetXBlock(index) * winstride, ybmax =
	lbpfeatures.GetYBlock(index) * winstride, offset = 0, endx, endy;

	vector<REAL> &features = lbpfeatures.GetFeatureRef(index);
	vector<vector<REAL> > pltpfeat(3, vector<REAL> (nbins)), nltpfeat(3,
			vector<REAL> (nbins));
	//	REAL *possum = new REAL[3 * nthrlevels], *negsum = new REAL[3 * nthrlevels]; // normalization constants...
	REAL possum[3], negsum[3];

	for (UINT rind = 0; rind < ybmax; rind += winstride) {
		y = rind + winstride/*foffset*/- radius;
		for (UINT colind = 0; colind < xbmax; colind += winstride) {
			x = colind + winstride /*foffset*/- radius;

			endx = MIN(x + tlbpstride, tposmap[0].nxpoints);
			endy = MIN(y + tlbpstride, tposmap[0].nypoints);
			for (int tl = 0; tl < 3 * nthrlevels; tl += 3) /*Remove Explicit Vectorization*/{
				for (UINT yiter = y; yiter < endy; ++yiter)
				for (UINT xiter = x; xiter < endx; ++xiter) {

					pltpfeat[0][tposmap[tl].GetValue(xiter, yiter)]++;
					pltpfeat[1][tposmap[tl + 1].GetValue(xiter, yiter)]++;
					pltpfeat[2][tposmap[tl + 2].GetValue(xiter, yiter)]++;

					nltpfeat[0][tnegmap[tl].GetValue(xiter, yiter)]++;
					nltpfeat[1][tnegmap[tl + 1].GetValue(xiter, yiter)]++;
					nltpfeat[2][tnegmap[tl + 2].GetValue(xiter, yiter)]++;
				}

				REAL *posptr = &features[offset], *negptr = &features[offset
				+ nbins];

				possum[0] = EPS;
				negsum[0] = EPS;
				if (norm == L2) {

					for (UINT biter = 0; biter < nbins; ++biter) { // add the three channels
						//							lbpptr[biter] = lbpfeat[0][biter]
						//									+ lbpfeat[1][biter] + lbpfeat[2][biter];
						posptr[biter] = pltpfeat[0][biter] + pltpfeat[1][biter]
						+ pltpfeat[2][biter];
						negptr[biter] = nltpfeat[0][biter] + nltpfeat[1][biter]
						+ nltpfeat[2][biter];

						possum[0] += posptr[biter] * posptr[biter];
						negsum[0] += negptr[biter] * negptr[biter];
					}
					possum[0] = sqrt(possum[0]);
					negsum[0] = sqrt(negsum[0]);
					for (UINT biter = 0; biter < nbins; ++biter) {
						posptr[biter] = posptr[biter] / possum[0];
						negptr[biter] = negptr[biter] / negsum[0];
					}

				} else {
					for (UINT biter = 0; biter < nbins; ++biter) { // add the three channels
						//							lbpptr[biter] = lbpfeat[0][biter]
						//									+ lbpfeat[1][biter] + lbpfeat[2][biter];
						posptr[biter] = pltpfeat[0][biter] + pltpfeat[1][biter]
						+ pltpfeat[2][biter];
						negptr[biter] = nltpfeat[0][biter] + nltpfeat[1][biter]
						+ nltpfeat[2][biter];

						//							lbpsum[0] += lbpptr[biter];
						possum[0] += posptr[biter];
						negsum[0] += negptr[biter];
					}

					for (UINT biter = 0; biter < nbins; ++biter) {
						//  							lbpptr[biter] = sqrt(lbpptr[biter] / lbpsum[0]);
						posptr[biter] = sqrt(posptr[biter] / possum[0]);
						negptr[biter] = sqrt(negptr[biter] / negsum[0]);
					}
				}
				for (UINT i = 0; i < 3; ++i) {// Zero out the temporary  histograms...
					fill(pltpfeat[i].begin(), pltpfeat[i].end(), 0);
					fill(nltpfeat[i].begin(), nltpfeat[i].end(), 0);
				}
				offset += 2 * nbins;
			}
		}

	}
	delete[] tposmap;
	delete[] tnegmap;
}

void OnlyLTPFeaturesRGBCell::GetFeatures(UINT index, int x, int y, int width_,
		int height_, vector<REAL>&feat) {// There is offset of one cell so, input x=0 means x=8 by taking
//	into account offset and the removed border
	int xcell = x / cell, ycell = y / cell;
	width_ /= cell;
	height_ /= cell;
	lbpfeatures.GetSkippedHOGCells(index, xcell, ycell, width_, height_, skip,
			&feat[0]);
}
void OnlyLTPFeaturesRGBCell::GetFeatures(UINT index, int x, int y, int width_,
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
void OnlyLTPFeaturesRGBCell::GetFoldedFeatures(UINT index, int x, int y,
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
void OnlyLTPFeaturesRGBCell::GetFoldedFeatures(UINT index, int x, int y,
		int twidth, int theight, REAL*feat) {// There is offset of one cell so, input x=0 means x=8 by taking
//	into account offset and the removed border
	int xcell = x / cell, ycell = y / cell;
	twidth /= cell;
	theight /= cell;
	lbpfeatures.GetSkippedFoldedHOGCells(index, xcell, ycell, twidth, theight,
			skip, feat);
}

// Only Flipping hog images....
void OnlyLTPFeaturesRGBCell::GetFlippedFeatures(UINT index, int x, int y,
		int twidth, int theight, vector<REAL>&feat) {// There is offset of one cell so, input x=0 means x=8 by taking
//	into account offset and the removed border
	int xcell = x / cell, ycell = y / cell;
	twidth /= cell;
	theight /= cell;
	lbpfeatures.GetSkippedFlippedFeatures(index, xcell, ycell, twidth, theight,
			skip, &feat[0]);

}

void OnlyLTPFeaturesRGBCell::GetFlippedFeatures(int twidth, int theight,
		REAL *ifeat, REAL *ofeat) {// There is offset of one cell so, input x=0 means x=8 by taking //	into account offset and the removed border
	twidth /= cell;
	theight /= cell;
	HInfo::FlipFeatures(twidth / skip, theight / skip, nbins * 2 * nthrlevels,
			ifeat, ofeat);
}

void OnlyLTPFeaturesRGBCell::UnFoldWeights(REAL *input, UINT w, UINT h, vector<
		REAL>& weights) {

	UINT tcwidth = w / cell, tcheight = h / cell;
	HInfo hinfo((UINT) ceil((double) tcwidth / (2.0 * skip)),
			(tcheight / skip), 1, nbins * 2 * nthrlevels, input);
	hinfo.GetUnFoldedHOGCells(tcwidth / skip, &weights[0]);// Flip the features
}

void OnlyLTPFeaturesRGBCell::PadFeatureMap(UINT index, UINT padx, UINT pady) { // does padding in the pixel space
// first pad the hogFeatures map...
	padx /= cell;
	pady /= cell;
	lbpfeatures.PadFeatureMap(index, padx, pady);
}

void OnlyLTPFeaturesRGBCell::DotProduct(UINT index, UINT fwidth, UINT fheight,
		vector<REAL>&filter, Store&response) {
	fwidth /= cell;
	fheight /= cell;
	lbpfeatures.SkippedDotProduct(index, fwidth, fheight, skip, &filter[0],
			response);
}

void OnlyLTPFeaturesRGBCell::DotProduct(UINT index, UINT fwidth, UINT fheight,
		REAL *filter, Store&response) {
	fwidth /= cell;
	fheight /= cell;
	lbpfeatures.SkippedDotProduct(index, fwidth, fheight, skip, filter,
			response);
}
#endif

