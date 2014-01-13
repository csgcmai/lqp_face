/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/
#define EPS 0.0001
#include "lbpFeaturesRGBBlock.h"
void LBPFeaturesRGBBlock::InitalizeMaps(Pyramid &pyobj, PyramidType sspace_) {
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
			ComputeLBPFeatures(imgref, i, lbpstride);
		}
	} else /*if (sspace == DoubleRes) */{

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
				ComputeLBPFeatures(imgref, i, lbpstride / 2);
			}
			// Compute the features...
			//        imgref=pyobj.GetImageRef(i);
			GetXYBlocks(imgref, xbmax, ybmax);
			lbpfeatures.SetFeature(i + nlpyramid, xbmax, ybmax,
					pyobj.GetScale(i));
			ComputeLBPFeatures(imgref, i + nlpyramid, lbpstride);
		}
	}
}
void LBPFeaturesRGBBlock::InitalizeMaps(Image &image) {
	UINT xbmax, ybmax;
	sspace = NoPyramid;
	lbpfeatures.Initialize(1, 0, celldim * dimmult);
	GetXYBlocks(image, xbmax, ybmax);
	lbpfeatures.SetFeature(xbmax, ybmax, 1.0);
	ComputeLBPFeatures(image, 0, lbpstride);
}
void LBPFeaturesRGBBlock::InitalizeMaps(Image&image, PyramidType sspace_) {
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
			ComputeLBPFeatures(imgref, i, lbpstride);
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
				ComputeLBPFeatures(imgref, i, lbpstride / 2);
			}
			// Compute the features...
			//        imgref=pyobj.GetImageRef(i);
			GetXYBlocks(imgref, xbmax, ybmax);
			lbpfeatures.SetFeature(i + nlpyramid, xbmax, ybmax,
					pyobj.GetScale(i));
			ComputeLBPFeatures(imgref, i + nlpyramid, lbpstride);
		}
	} else {
		lbpfeatures.Initialize(1, 0, celldim * dimmult);
		GetXYBlocks(image, xbmax, ybmax);
		lbpfeatures.SetFeature(xbmax, ybmax, 1.0);
		ComputeLBPFeatures(image, 0, lbpstride);
	}
}
void LBPFeaturesRGBBlock::ComputeTextureFeatures(Image &imgref, UINT xcellsize,
		UINT ycellsize, vector<REAL>&features) {
	ComputeTextureFeatures(imgref, xcellsize, ycellsize, 0, 0,
			imgref.columns() - 1, imgref.rows() - 1, features);
}
void LBPFeaturesRGBBlock::ComputeTextureFeatures(Image &imgref, int sx, int sy,
		int ex, int ey, UINT xcellsize, UINT ycellsize, vector<REAL>&features) {
	LBPMap tlbpmap[3];
	ComputeLBPMap(imgref, tlbpmap);
	UINT width = ex - sx + 1, height = ey - sy + 1;
	int xoffset = -radius, yoffset = -radius;
	UINT xbmax = width / xcellsize, ybmax = height / ycellsize;
	UINT offset = 0;

	features.resize(xbmax * ybmax * celldim * dimmult, 0);
	vector<vector<REAL> > lbpfeat(3, vector<REAL> (celldim));
	REAL lbpsum[3]; // normalization constants...
	//	fill(features.begin(), features.end(), 0);

	for (UINT rind = 0; rind < ybmax; rind += ycellsize) {
		for (UINT colind = 0; colind < xbmax; colind += xcellsize) {

			for (int yiter = rind; yiter < rind + ycellsize; ++yiter) {
				int y = MIN(MAX(yiter+yoffset+sy,0),tlbpmap[0].nypoints-1);
				for (int xiter = colind; xiter < colind + xcellsize; ++xiter) {
					int x = MIN(MAX(xiter+xoffset+sx,0),tlbpmap[0].nxpoints-1);

					assert(tlbpmap[0].GetValue(x, y)!=-1);

					lbpfeat[0][tlbpmap[0].GetValue(x, y)]++;// red
					lbpfeat[1][tlbpmap[1].GetValue(x, y)]++;// green
					lbpfeat[2][tlbpmap[2].GetValue(x, y)]++;// blue

				}
			}
			{// normalization
				REAL *lbpptr = &features[offset];
				lbpsum[0] = EPS;
				if (norm == L2) {
					for (UINT biter = 0; biter < nbins; ++biter) { // add the three channels
						lbpptr[biter] = lbpfeat[0][biter] + lbpfeat[1][biter]
								+ lbpfeat[2][biter];
						lbpsum[0] += lbpptr[biter] * lbpptr[biter];
					}
					lbpsum[0] = sqrt(lbpsum[0]);
					for (UINT biter = 0; biter < nbins; ++biter) {
						lbpptr[biter] = (lbpptr[biter] / lbpsum[0]);
					}

				} else if (norm == NT_Average) {
					for (UINT biter = 0; biter < nbins; ++biter) { // add the three channels
						lbpptr[biter] = (lbpfeat[0][biter] + lbpfeat[1][biter]
								+ lbpfeat[2][biter]) / 3.0;
					}
				} else {
					for (UINT biter = 0; biter < nbins; ++biter) { // add the three channels
						lbpptr[biter] = lbpfeat[0][biter] + lbpfeat[1][biter]
								+ lbpfeat[2][biter];
						lbpsum[0] += lbpptr[biter];
					}
					if (norm == L1) {
						for (UINT biter = 0; biter < nbins; ++biter)
							lbpptr[biter] = (lbpptr[biter] / lbpsum[0]);
					} else {
						for (UINT biter = 0; biter < nbins; ++biter)
							lbpptr[biter] = sqrt(lbpptr[biter] / lbpsum[0]);
					}
				}
			}
			for (UINT i = 0; i < 3; ++i) {// Zero out the temporary  histograms...
				fill(lbpfeat[i].begin(), lbpfeat[i].end(), 0);
			}
			offset += nbins;
		}
	}
}
void LBPFeaturesRGBBlock::ComputeLBPFeatures(Image &image, UINT index,
		UINT sbin) {
	LBPMap tlbpmap[3];
	ComputeLBPMap(image, tlbpmap);
	if (hmethod == HM_Bilinear) {
		vector<REAL> &tvec = lbpfeatures.GetFeatureRef(index);
		ComputeHistBilinear(image, sbin, tlbpmap, &tvec[0]);
	} else
		ComputeHistDiscrete(image, index, sbin, tlbpmap);
}
void LBPFeaturesRGBBlock::ComputeHistDiscrete(Image& imgref, UINT index,
		UINT winstride, LBPMap *tlbpmap) {

	UINT x, y, xbmax = lbpfeatures.GetXBlock(index) * winstride, ybmax =
			lbpfeatures.GetYBlock(index) * winstride, offset = 0, endx, endy;

	vector<REAL> &features = lbpfeatures.GetFeatureRef(index);
	vector<vector<REAL> > lbpfeat(3, vector<REAL> (nbins));
	REAL lbpsum[3]; // normalization constants...
	//	fill(features.begin(), features.end(), 0);

	for (UINT rind = 0; rind < ybmax; rind += winstride) {
		y = rind + winstride/*foffset*/- radius;
		for (UINT colind = 0; colind < xbmax; colind += winstride) {
			x = colind + winstride/*foffset*/- radius;

			endx = MIN(x+winstride,tlbpmap[0].nxpoints);
			endy = MIN(y+winstride,tlbpmap[0].nypoints);

			for (UINT yiter = y; yiter < endy; ++yiter)
				for (UINT xiter = x; xiter < endx; ++xiter) {
					lbpfeat[0][tlbpmap[0].GetValue(xiter, yiter)]++;// red
					lbpfeat[1][tlbpmap[1].GetValue(xiter, yiter)]++;// green
					lbpfeat[2][tlbpmap[2].GetValue(xiter, yiter)]++;// blue

				}

			if (nseparate) {

				lbpsum[0] = lbpsum[1] = lbpsum[2] = EPS;
				// Currently only doing L1-Sqrt Normalization further to be added...
				for (UINT biter = 0; biter < nbins; ++biter) {
					lbpsum[0] += lbpfeat[0][biter];
					lbpsum[1] += lbpfeat[1][biter];
					lbpsum[2] += lbpfeat[2][biter];

				}

				for (UINT biter = 0; biter < nbins; ++biter) {
					lbpfeat[0][biter] = sqrt(lbpfeat[0][biter] / lbpsum[0]);
					lbpfeat[1][biter] = sqrt(lbpfeat[1][biter] / lbpsum[1]);
					lbpfeat[2][biter] = sqrt(lbpfeat[2][biter] / lbpsum[2]);
				}

				for (UINT biter = 0; biter < nbins; ++biter) {
					features[offset + biter] = lbpfeat[0][biter]
							+ lbpfeat[1][biter] + lbpfeat[2][biter];
				}
			} else {

				REAL *lbpptr = &features[offset];

				lbpsum[0] = EPS;

				for (UINT biter = 0; biter < nbins; ++biter) { // add the three channels
					lbpptr[biter] = lbpfeat[0][biter] + lbpfeat[1][biter]
							+ lbpfeat[2][biter];
					lbpsum[0] += lbpptr[biter];
				}

				for (UINT biter = 0; biter < nbins; ++biter) {
					lbpptr[biter] = sqrt(lbpptr[biter] / lbpsum[0]);
				}

			}

			for (UINT i = 0; i < 3; ++i) {// Zero out the temporary  histograms...
				fill(lbpfeat[i].begin(), lbpfeat[i].end(), 0);
			}
			offset += nbins;
		}

	}

}
void LBPFeaturesRGBBlock::ComputeHistBilinear(Image &image, UINT sbin,
		LBPMap *tlbpmap, REAL *feat) {

	int dims[2];
	dims[0] = image.rows();
	dims[1] = image.columns();

	// memory for caching histograms & their norms
	int blocks[2];
	blocks[0] = (int) round((double) dims[0] / (double) sbin);
	blocks[1] = (int) round((double) dims[1] / (double) sbin);
	UINT bdim = blocks[0] * blocks[1];
	double *rhist = new double[(bdim * celldim)], *ghist = new double[(bdim
			* celldim)], *bhist = new double[(bdim * celldim)], *hist =
			new double[(bdim * celldim)];
	double *lbpnorm = new double[bdim * dimmult];
	fill(lbpnorm, lbpnorm + bdim * dimmult, 0);
	fill(rhist, rhist + bdim * celldim, 0);
	fill(ghist, ghist + bdim * celldim, 0);
	fill(bhist, bhist + bdim * celldim, 0);
	fill(hist, hist + bdim * celldim, 0); // for pooling of RGB channels....

	// memory for HOG features
	int out[2];
	out[0] = max(blocks[0] - 2, 0);
	out[1] = max(blocks[1] - 2, 0);

	int visible[2];
	visible[0] = blocks[0] * sbin;
	visible[1] = blocks[1] * sbin;
	UINT tx, ty;

	UINT width = image.columns(), height = image.rows(), mx, my;
	const PixelPacket *pix = image.getConstPixels(0, 0, width, height);
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

			int ro = tlbpmap[0].GetValue(tx, ty), // red, GREEN & BLUE
					go = tlbpmap[1].GetValue(tx, ty), //
					bo = tlbpmap[2].GetValue(tx, ty); //
			// add to 4 histograms around pixel using linear interpolation
			double xp = ((double) x + 0.5) / (double) sbin - 0.5;
			double yp = ((double) y + 0.5) / (double) sbin - 0.5;
			int ixp = (int) floor(xp);
			int iyp = (int) floor(yp);
			double vx0 = xp - ixp;
			double vy0 = yp - iyp;
			double vx1 = 1.0 - vx0;
			double vy1 = 1.0 - vy0;
			// finding gradient values only for testing purposes...
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

			if (ixp >= 0 && iyp >= 0) {
				*(rhist + ixp * blocks[0] + iyp + ro * bdim) += vx1 * vy1 * v;
				*(ghist + ixp * blocks[0] + iyp + go * bdim) += vx1 * vy1 * v;
				*(bhist + ixp * blocks[0] + iyp + bo * bdim) += vx1 * vy1 * v;
			}

			if (ixp + 1 < blocks[1] && iyp >= 0) {
				*(rhist + (ixp + 1) * blocks[0] + iyp + ro * bdim) += vx0 * vy1
						* v;
				*(ghist + (ixp + 1) * blocks[0] + iyp + go * bdim) += vx0 * vy1
						* v;
				*(bhist + (ixp + 1) * blocks[0] + iyp + bo * bdim) += vx0 * vy1
						* v;
			}

			if (ixp >= 0 && iyp + 1 < blocks[0]) {
				*(rhist + ixp * blocks[0] + (iyp + 1) + ro * bdim) += vx1 * vy0
						* v;
				*(ghist + ixp * blocks[0] + (iyp + 1) + go * bdim) += vx1 * vy0
						* v;
				*(bhist + ixp * blocks[0] + (iyp + 1) + bo * bdim) += vx1 * vy0
						* v;
			}

			if (ixp + 1 < blocks[1] && iyp + 1 < blocks[0]) {
				*(rhist + (ixp + 1) * blocks[0] + (iyp + 1) + ro * bdim) += vx0
						* vy0 * v;
				*(ghist + (ixp + 1) * blocks[0] + (iyp + 1) + go * bdim) += vx0
						* vy0 * v;
				*(bhist + (ixp + 1) * blocks[0] + (iyp + 1) + bo * bdim) += vx0
						* vy0 * v;
			}
		}
	}
	// compute energy in each block by summing over orientations
	if (add) {
		for (int o = 0; o < celldim; o++) {
			double *srcr = rhist + o * bdim, *srcg = ghist + o * bdim, *srcb =
					bhist + o * bdim, *src = hist + o * bdim;
			double *dst = lbpnorm;
			double *end = lbpnorm + bdim;
			while (dst < end) {
				*src = *srcr + *srcg + *srcb; // pooling of rgb channels...
				if (norm == L2)
					*(dst++) += *src * *src;
				else
					*(dst++) += *src;
				srcr++;
				srcg++;
				srcb++;
				src++;
			}
		}
		// compute features
		REAL sc1, sc2, sc3, sc4;
		for (int x = 0; x < out[1]; x++) {
			for (int y = 0; y < out[0]; y++) {
				REAL *dst = feat + celldim * dimmult * (x + y * out[1]); // write in column major order
				double *p, n1, n2, n3, n4;
				double *src = hist + (x + 1) * blocks[0] + (y + 1);
				if (useblocknorm) {/*normlaize by lbp energy of four neighbouring blocks*/
					p = lbpnorm + (x + 1) * blocks[0] + y + 1;
					n1 = 1.0 / (*p + *(p + 1) + *(p + blocks[0]) + *(p
							+ blocks[0] + 1) + EPS);
					p = lbpnorm + (x + 1) * blocks[0] + y;
					n2 = 1.0 / (*p + *(p + 1) + *(p + blocks[0]) + *(p
							+ blocks[0] + 1) + EPS);
					p = lbpnorm + x * blocks[0] + y + 1;
					n3 = 1.0 / (*p + *(p + 1) + *(p + blocks[0]) + *(p
							+ blocks[0] + 1) + EPS);
					p = lbpnorm + x * blocks[0] + y;
					n4 = 1.0 / (*p + *(p + 1) + *(p + blocks[0]) + *(p
							+ blocks[0] + 1) + EPS);

					for (int o = 0; o < celldim; o++) {
						double sum = *src;
						double h1 = sqrt(sum * n1);
						double h2 = sqrt(sum * n2);
						double h3 = sqrt(sum * n3);
						double h4 = sqrt(sum * n4);
						*(dst + o) = h4;
						*(dst + celldim + o) = h2;
						*(dst + celldim * 2 + o) = h3;
						*(dst + celldim * 3 + o) = h1;
						src += bdim;
					}
				} else {
					// single cell normalization...
					if (norm == L2) {
						REAL normval = 1 / sqrt(
								*(lbpnorm + (x + 1) * blocks[0] + y + 1) + EPS);
						for (int o = 0; o < celldim; o++) {
							*dst = *src * normval;
							src += bdim;
							++dst;
						}
					} else if (norm == LOG2) {
						for (int o = 0; o < celldim; o++) {
							*dst = log2(*src + 1);
							src += bdim;
							++dst;
						}

					} else {
						REAL normval = *(lbpnorm + (x + 1) * blocks[0] + y + 1)
								+ EPS;
						for (int o = 0; o < celldim; o++) {
							*dst = sqrt(*src / normval);
							src += bdim;
							++dst;
						}
					}
				}
			}
		}
	} else {// separate normalization blocks for each color channel
		for (int o = 0; o < celldim; o++) {
			double *srcr = rhist + o * bdim, *srcg = ghist + o * bdim, *srcb =
					bhist + o * bdim;
			double *rdst = lbpnorm, *gdst = lbpnorm + bdim, *bdst = gdst + bdim;
			double *end = lbpnorm + bdim;
			while (rdst < end) {
				if (norm == L2) {
					*(rdst++) += *srcr * *srcr;
					*(gdst++) += *srcg * *srcg;
					*(bdst++) += *srcb * *srcb;
				} else {
					*(rdst++) += *srcr;
					*(gdst++) += *srcg;
					*(bdst++) += *srcb;
				}
				srcr++;
				srcg++;
				srcb++;
			}
		}

		// compute features
		REAL sc1, sc2, sc3, sc4;
		UINT offset;
		for (int x = 0; x < out[1]; x++) {
			for (int y = 0; y < out[0]; y++) {
				REAL *rdst = feat + celldim * dimmult * (x + y * out[1]),
						*gdst = rdst + celldim, *bdst = gdst + celldim; // write in column major order
				offset = (x + 1) * blocks[0] + (y + 1);
				double *rsrc = rhist + offset, *gsrc = ghist + offset, *bsrc =
						bhist + offset;
				double *rnptr = lbpnorm + offset, *gnptr = rnptr + bdim,
						*bnptr = gnptr + bdim;
				// single cell normalization...
				// single cell normalization...
				if (norm == L2) {

					REAL rnormval = 1 / sqrt(*rnptr + EPS);
					REAL gnormval = 1 / sqrt(*gnptr + EPS);
					REAL bnormval = 1 / sqrt(*bnptr + EPS);
					for (int o = 0; o < celldim; o++) {
						*rdst = *rsrc * rnormval;
						*gdst = *gsrc * gnormval;
						*bdst = *bsrc * bnormval;
						rsrc += bdim;
						gsrc += bdim;
						bsrc += bdim;

						++rdst;
						++gdst;
						++bdst;
					}
				} else {

					REAL rnormval = 1 / (*rnptr + EPS);
					REAL gnormval = 1 / (*gnptr + EPS);
					REAL bnormval = 1 / (*bnptr + EPS);
					for (int o = 0; o < celldim; o++) {
						*rdst = sqrt(*rsrc * rnormval);
						*gdst = sqrt(*gsrc * gnormval);
						*bdst = sqrt(*bsrc * bnormval);
						rsrc += bdim;
						gsrc += bdim;
						bsrc += bdim;

						++rdst;
						++gdst;
						++bdst;

					}
				}
			}
		}
	}

	delete[] hist;
	delete[] rhist;
	delete[] ghist;
	delete[] bhist;
	delete[] lbpnorm;
}
void LBPFeaturesRGBBlock::GetFeatures(UINT index, int x, int y, int width_,
		int height_, vector<REAL>&feat) {// There is offset of one cell so, input x=0 means x=8 by taking
	//    into account offset and the removed border
	int xcell = x / cell, ycell = y / cell;
	width_ /= cell;
	height_ /= cell;
	lbpfeatures.GetSkippedHOGCells(index, xcell, ycell, width_, height_, skip,
			&feat[0]);
}
void LBPFeaturesRGBBlock::GetFeatures(UINT index, int x, int y, int width_,
		int height_, REAL*feat) {// There is offset of one cell so, input x=0 means x=8 by taking
	//    into account offset and the removed border
	int xcell = x / cell, ycell = y / cell;
	width_ /= cell;
	height_ /= cell;
	lbpfeatures.GetSkippedHOGCells(index, xcell, ycell, width_, height_, skip,
			feat);
}
/*
 * Contains only folded HOG's....*/
void LBPFeaturesRGBBlock::GetFoldedFeatures(UINT index, int x, int y,
		int twidth, int theight, vector<REAL>&feat) {// There is offset of one cell so, input x=0 means x=8 by taking
	//    into account offset and the removed border
	int xcell = x / cell, ycell = y / cell;
	twidth /= cell;
	theight /= cell;
	lbpfeatures.GetSkippedFoldedHOGCells(index, xcell, ycell, twidth, theight,
			skip, &feat[0]);
}

/*
 * Contains only folded HOG's....*/
void LBPFeaturesRGBBlock::GetFoldedFeatures(UINT index, int x, int y,
		int twidth, int theight, REAL*feat) {// There is offset of one cell so, input x=0 means x=8 by taking
	//    into account offset and the removed border
	int xcell = x / cell, ycell = y / cell;
	twidth /= cell;
	theight /= cell;
	lbpfeatures.GetSkippedFoldedHOGCells(index, xcell, ycell, twidth, theight,
			skip, feat);
}

// Only Flipping hog images....
void LBPFeaturesRGBBlock::GetFlippedFeatures(UINT index, int x, int y,
		int twidth, int theight, vector<REAL>&feat) {// There is offset of one cell so, input x=0 means x=8 by taking
	//    into account offset and the removed border
	int xcell = x / cell, ycell = y / cell;
	twidth /= cell;
	theight /= cell;
	lbpfeatures.GetSkippedFlippedFeatures(index, xcell, ycell, twidth, theight,
			skip, &feat[0]);

}

void LBPFeaturesRGBBlock::GetFlippedFeatures(int twidth, int theight,
		REAL *ifeat, REAL *ofeat) {// There is offset of one cell so, input x=0 means x=8 by taking
	//    into account offset and the removed border
	twidth /= cell;
	theight /= cell;
	HInfo::FlipFeatures(twidth / skip, theight / skip, celldim * dimmult,
			ifeat, ofeat);
}

void LBPFeaturesRGBBlock::UnFoldWeights(REAL *input, UINT w, UINT h,
		vector<REAL>& weights) {

	UINT tcwidth = w / cell, tcheight = h / cell;
	HInfo hinfo((UINT) ceil((double) tcwidth / (2.0 * skip)),
			(tcheight / skip), 1, celldim * dimmult, input);
	hinfo.GetUnFoldedHOGCells(tcwidth / skip, &weights[0]);// Flip the features
}

void LBPFeaturesRGBBlock::PadFeatureMap(UINT index, UINT padx, UINT pady) { // does padding in the pixel space
	// first pad the hogFeatures map...
	padx /= cell;
	pady /= cell;
	lbpfeatures.PadFeatureMap(index, padx, pady);
}

void LBPFeaturesRGBBlock::DotProduct(UINT index, UINT fwidth, UINT fheight,
		vector<REAL>&filter, Store&response) {
	fwidth /= cell;
	fheight /= cell;
	lbpfeatures.SkippedDotProduct(index, fwidth, fheight, skip, &filter[0],
			response);
}
void LBPFeaturesRGBBlock::DotProduct(UINT index, UINT fwidth, UINT fheight,
		REAL *filter, Store&response) {
	fwidth /= cell;
	fheight /= cell;
	lbpfeatures.SkippedDotProduct(index, fwidth, fheight, skip, filter,
			response);
}

void CSLBPFeaturesRGBCell::ComputeLBPMap(Image &image, LBPMap lbpmap[]) {
	ColorGray mycolor, c1, c2, c3, c4;
	UINT xdim = image.columns(), ydim = image.rows();
	int nypoints = ydim - 2 * radius, reg = 1, fx, fy, cx, cy, nxpoints = xdim
			- 2 * radius, cenx, ceny, count, rvalue, gvalue, bvalue;
	REAL tmpx, tmpy, fracx, fracy, rcenval, gcenval, bcenval;
	REAL x, y, w1, w2, w3, w4;
	vector<REAL> pixval(npoints * 3, 0);
	REAL *rpixval = &pixval[0], *gpixval = rpixval + npoints, *bpixval =
			gpixval + npoints;
	lbpmap[0].Init(nxpoints, nypoints); // for RGB channels
	lbpmap[1].Init(nxpoints, nypoints);
	lbpmap[2].Init(nxpoints, nypoints);
	vector<REAL> rimpix(xdim * ydim, 0), bimpix(xdim * ydim, 0), gimpix(
			xdim * ydim, 0);
	//	image.flop();
	pim.Process(image, rimpix, gimpix, bimpix);
	//	int fbins[]={4,3,2,1,0,7,6,5};

	//		ofstream ofile("test_map.lst");
	for (int i = 0; i < nypoints; ++i) {
		for (int j = 0; j < nxpoints; ++j) {
			cenx = j + radius;
			ceny = i + radius;
			count = 0;
			rcenval = rimpix[ceny * xdim + cenx];
			gcenval = gimpix[ceny * xdim + cenx];
			bcenval = bimpix[ceny * xdim + cenx];
			//ofile<<endl<< i<<" " << j<< " "<< rcenval <<" " ;
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

					if (ABS(fracx) < 1e-6 && ABS(fracy) < 1e-6)
					// no interpolation needed
					{
						rpixval[count] = rimpix[fy * xdim + fx];
						gpixval[count] = gimpix[fy * xdim + fx];
						bpixval[count] = bimpix[fy * xdim + fx];
					} else {
						w1 = (1 - fracx) * (1 - fracy);
						w2 = fracx * (1 - fracy);
						w3 = (1 - fracx) * fracy;
						w4 = fracx * fracy;
						rpixval[count] = w1 * rimpix[fy * xdim + fx] + w2
								* rimpix[fy * xdim + cx] + w3 * rimpix[cy
								* xdim + fx] + w4 * rimpix[cy * xdim + cx];

						gpixval[count] = w1 * gimpix[fy * xdim + fx] + w2
								* gimpix[fy * xdim + cx] + w3 * gimpix[cy
								* xdim + fx] + w4 * gimpix[cy * xdim + cx];

						bpixval[count] = w1 * bimpix[fy * xdim + fx] + w2
								* bimpix[fy * xdim + cx] + w3 * bimpix[cy
								* xdim + fx] + w4 * bimpix[cy * xdim + cx];

					}
				} else {
					rpixval[count] = rimpix[tmpy * xdim + tmpx];
					gpixval[count] = gimpix[tmpy * xdim + tmpx];
					bpixval[count] = bimpix[tmpy * xdim + tmpx];
				}

			}
			rvalue = 0;
			gvalue = 0;
			bvalue = 0;
			int wt = 0;
			if (issigned) { /*by default use the unsigned for cslbp,ltp*/
				for (UINT iter = 0; iter < midpoint; ++iter) {
					wt = reg << iter;
					rvalue += ((rpixval[iter] - rpixval[iter + midpoint])
							> tolerance ? wt : 0);
					gvalue += ((gpixval[iter] - gpixval[iter + midpoint])
							> tolerance ? wt : 0);
					bvalue += ((bpixval[iter] - bpixval[iter + midpoint])
							> tolerance ? wt : 0);
				}
			} else {
				for (UINT iter = 0; iter < midpoint; ++iter) {
					wt = reg << iter;
					rvalue += (ABS(rpixval[iter] - rpixval[iter + midpoint])
							> tolerance ? wt : 0);
					gvalue += (ABS(gpixval[iter] - gpixval[iter + midpoint])
							> tolerance ? wt : 0);
					bvalue += (ABS(bpixval[iter] - bpixval[iter + midpoint])
							> tolerance ? wt : 0);
				}

			}

			lbpmap[0](j, i) = rvalue;
			lbpmap[1](j, i) = gvalue;
			lbpmap[2](j, i) = bvalue;
		}
	}
}
