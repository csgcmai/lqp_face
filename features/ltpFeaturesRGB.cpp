/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/
#include "ltpFeaturesRGB.h"
void LTPFeaturesRGB::ComputeLTPMap(Image &image, vector<LBPMap> &lbpmap,
		vector<LBPMap>& tpmap, vector<LBPMap> &tnmap) {
	UINT xdim = image.columns(), ydim = image.rows();
	int nypoints = ydim - 2 * radius, reg = 1, fx, fy, cx, cy, nxpoints = xdim
			- 2 * radius, cenx, ceny, count, tval, value[3], posvalue[3],
			negvalue[3];
	REAL tmpx, tmpy, fracx, fracy, pixval[3], cenval[3];
	REAL x, y, w1, w2, w3, w4;

	for (int i = 0; i < 3; ++i)// for 3 channels RGB
	{
		lbpmap[i].Init(nxpoints, nypoints);
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
					value[iter] += pixval[iter] >= cenval[iter] ? tval : 0;
					posvalue[iter]
							+= pixval[iter] >= cenval[iter] + tolerance ? tval
									: 0;
					negvalue[iter]
							+= pixval[iter] <= cenval[iter] - tolerance ? tval
									: 0;
				}
			}

			for (int iter = 0; iter < 3; ++iter) {
				lbpmap[iter](j, i) = map[value[iter]];
				tpmap[iter](j, i) = map[posvalue[iter]];
				tnmap[iter](j, i) = map[negvalue[iter]];
			}
		}
	}
}
void LTPFeaturesRGB::ComputeLTPMap(Image &image, LBPMap lbpmap[],
		LBPMap tpmap[], LBPMap tnmap[]) {
	UINT xdim = image.columns(), ydim = image.rows();
	int nypoints = ydim - 2 * radius, reg = 1, fx, fy, cx, cy, nxpoints = xdim
			- 2 * radius, cenx, ceny, count, tval, value[3], posvalue[3],
			negvalue[3];
	REAL tmpx, tmpy, fracx, fracy, pixval[3], cenval[3];
	REAL x, y, w1, w2, w3, w4;

	for (int i = 0; i < 3; ++i)// for 3 channels RGB
	{
		lbpmap[i].Init(nxpoints, nypoints);
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
					if (pixval[iter] >= cenval[iter]) {
						value[iter] += tval;
					}
					if (pixval[iter] >= (cenval[iter] + tolerance))
						posvalue[iter] += tval;
					if (pixval[iter] <= (cenval[iter] - tolerance))
						negvalue[iter] += tval;
				}
			}

			for (int iter = 0; iter < 3; ++iter) {
				lbpmap[iter](j, i) = map[value[iter]];
				tpmap[iter](j, i) = map[posvalue[iter]];
				tnmap[iter](j, i) = map[negvalue[iter]];
			}
		}
	}
}
void LTPFeaturesRGB::ComputeLTPMapFast(Image &image, LBPMap lbpmap[],
		LBPMap tpmap[], LBPMap tnmap[]) {
	UINT xdim = image.columns(), ydim = image.rows();
	int nypoints = ydim - 2 * radius, reg = 1, fx, fy, cx, cy, nxpoints = xdim
			- 2 * radius, cenx, ceny, count, tval, rvalue, gvalue, bvalue,
			rposvalue, gposvalue, bposvalue, rnegvalue, gnegvalue, bnegvalue;
	REAL tmpx, tmpy, fracx, fracy, rpixval, gpixval, bpixval, rcenval, gcenval,
			bcenval;
	REAL x, y, w1, w2, w3, w4;

	for (int i = 0; i < 3; ++i)// for 3 channels RGB
	{
		lbpmap[i].Init(nxpoints, nypoints);
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

			rvalue = gvalue = bvalue = 0;
			rposvalue = gposvalue = bposvalue = 0;
			rnegvalue = gnegvalue = bnegvalue = 0;

			count = 0;

			rcenval = rimpix[ceny * xdim + cenx];
			gcenval = gimpix[ceny * xdim + cenx];
			bcenval = bimpix[ceny * xdim + cenx];
			for (UINT count = 0, pcount = 0; count < npoints; ++count, pcount
					+= 10) {
				fx = cenx + biinfo[pcount];
				fy = ceny + biinfo[pcount + 1];
				cx = cenx + biinfo[pcount + 2];
				cy = ceny + biinfo[pcount + 3];
				UINT m1 = fy * xdim, p1 = m1 + fx;
				if (dointerpol[count]) {
					w1 = biinfo[pcount + 6];
					w2 = biinfo[pcount + 7];
					w3 = biinfo[pcount + 8];
					w4 = biinfo[pcount + 9];

					UINT p2 = m1 + cx, m2 = cy * xdim, p3 = m2 + fx, p4 = m2
							+ cx;

					rpixval = w1 * rimpix[p1] + w2 * rimpix[p2] + w3
							* rimpix[p3] + w4 * rimpix[p4];

					gpixval = w1 * gimpix[p1] + w2 * gimpix[p2] + w3
							* gimpix[p3] + w4 * gimpix[p4];

					bpixval = w1 * bimpix[p1] + w2 * bimpix[p2] + w3
							* bimpix[p3] + w4 * bimpix[p4];
				} else {
					rpixval = rimpix[p1];
					gpixval = gimpix[p1];
					bpixval = bimpix[p1];
				}

				tval = reg << count;
				/*				for (int iter = 0; iter < 3; ++iter) {
				 if (pixval[iter] >= cenval[iter])
				 value[iter] += tval;
				 if (pixval[iter] >= cenval[iter] + tolerance)
				 posvalue[iter] += tval;
				 if (pixval[iter] <= cenval[iter] - tolerance)
				 negvalue[iter] += tval;
				 }*/
				if (rpixval >= rcenval)
					rvalue += tval;
				if (gpixval >= gcenval)
					gvalue += tval;
				if (bpixval >= bcenval)
					bvalue += tval;

				if (rpixval >= rcenval + tolerance)
					rposvalue += tval;
				if (gpixval >= gcenval + tolerance)
					gposvalue += tval;
				if (bpixval >= bcenval + tolerance)
					bposvalue += tval;

				if (rpixval <= rcenval - tolerance)
					rnegvalue += tval;
				if (gpixval <= gcenval - tolerance)
					gnegvalue += tval;
				if (bpixval <= bcenval - tolerance)
					bnegvalue += tval;
			}

			//			for (int iter = 0; iter < 3; ++iter) {
			lbpmap[0](j, i) = map[rvalue];
			lbpmap[1](j, i) = map[gvalue];
			lbpmap[2](j, i) = map[bvalue];

			tpmap[0](j, i) = map[rposvalue];
			tpmap[1](j, i) = map[gposvalue];
			tpmap[2](j, i) = map[bposvalue];

			tnmap[0](j, i) = map[rnegvalue];
			tnmap[1](j, i) = map[gnegvalue];
			tnmap[2](j, i) = map[bnegvalue];
		}
	}
}
// to compute the LTP Features on the modulus of edge image...sqrt(gradx^2+grady^2)
void LTPFeaturesRGB::ComputeLTPMap(vector<REAL>&impix, Point<UINT>& dim,
		LBPMap &lbpmap, LBPMap& tpmap, LBPMap &tnmap) {
	UINT xdim, ydim;
	dim.GetCoord(xdim, ydim);
	int nypoints = ydim - 2 * radius, reg = 1, fx, fy, cx, cy, nxpoints = xdim
			- 2 * radius, cenx, ceny, count, value, posvalue, negvalue;
	REAL tmpx, tmpy, fracx, fracy, pixval, cenval;
	REAL x, y;
	lbpmap.Init(nxpoints, nypoints);
	tpmap.Init(nxpoints, nypoints);
	tnmap.Init(nxpoints, nypoints);
	//	vector<REAL> impix(xdim * ydim, 0);
	//	pim.Process(image, impix);
	for (int i = 0; i < nypoints; ++i) {
		for (int j = 0; j < nxpoints; ++j) {
			cenx = j + radius;
			ceny = i + radius;
			value = posvalue = negvalue = 0;
			count = 0;
			cenval = impix[ceny * xdim + cenx];//column major order
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
				//
				if (ABS(fracx) < 1e-6 && ABS(fracy) < 1e-6)
					// no interpolation needed
					pixval = impix[fy * xdim + fx];
				else
					pixval = (1 - fracx) * (1 - fracy) * impix[fy * xdim + fx]
							+ fracx * (1 - fracy) * impix[fy * xdim + cx] + (1
							- fracx) * fracy * impix[cy * xdim + fx] + fracx
							* fracy * impix[cy * xdim + cx];

				value += pixval >= cenval ? reg << count : 0;
				posvalue += pixval >= cenval + tolerance ? reg << count : 0;
				negvalue += pixval <= cenval - tolerance ? reg << count : 0;

			}
			lbpmap(j, i) = map[value];
			tpmap(j, i) = map[posvalue];
			tnmap(j, i) = map[negvalue];
		}
	}
}
void LTPFeaturesRGB::InitalizeMaps(Image &image) {
	lbpmaps.resize(1, vector<LBPMap> (3));
	posmap.resize(1, vector<LBPMap> (3));
	negmap.resize(1, vector<LBPMap> (3));
	ComputeLTPMap(image, lbpmaps[0], posmap[0], negmap[0]);
}
void LTPFeaturesRGB::InitalizeMaps(Pyramid&pyobj, PyramidType sspace_) {
	UINT pycount = 0;
	// Read the image....
	sspace = sspace_;
	if (sspace_ == SingleRes) {
		if (!IsFolded()) {
			pycount = pyobj.GetNLevels();
			lbpmaps.resize(pycount, vector<LBPMap> (3));
			posmap.resize(pycount, vector<LBPMap> (3));
			negmap.resize(pycount, vector<LBPMap> (3));
			for (UINT i = 0; i < pycount; i++) {
				Image& imgref = pyobj.GetImageRef(i);
				// Compute the Gradient Image
				ComputeLTPMap(imgref, lbpmaps[i], posmap[i], negmap[i]);
			}
		}
	} else {

	}
}
void LTPFeaturesRGB::InitalizeMaps(Image& image, PyramidType sspace_) {
	UINT pycount = 0;
	// Read the image....
	sspace = sspace_;
	if (sspace_ == SingleRes) {
		if (!IsFolded()) {
			pycount = pyobj.GeneratePyramid(image);
			lbpmaps.resize(pycount, vector<LBPMap> (3));
			posmap.resize(pycount, vector<LBPMap> (3));
			negmap.resize(pycount, vector<LBPMap> (3));
			for (UINT i = 0; i < pycount; i++) {
				Image& imgref = pyobj.GetImageRef(i);
				// Compute the Gradient Image
				ComputeLTPMap(imgref, lbpmaps[i], posmap[i], negmap[i]);
			}
		} else {
		}
	} else {
		if (!IsFolded()) {
			lbpmaps.resize(1, vector<LBPMap> (3));
			posmap.resize(1, vector<LBPMap> (3));
			negmap.resize(1, vector<LBPMap> (3));
			ComputeLTPMap(image, lbpmaps[0], posmap[0], negmap[0]);
		} else {
		}
	}
}
#ifndef OLD_CODE
void LTPFeaturesRGB::GetFeatures(UINT index, int x, int y, vector<REAL>&feat) {

	UINT featsize = nxcells * nycells * nbins;
	fill(feat.begin(), feat.end(), 0);// Initialize with 0's;
	vector<REAL> rfeat(featsize, 0), gfeat(featsize, 0), bfeat(featsize, 0);
	vector<REAL> prfeat(featsize, 0), pgfeat(featsize, 0), pbfeat(featsize, 0);
	vector<REAL> nrfeat(featsize, 0), ngfeat(featsize, 0), nbfeat(featsize, 0);
	ComputeHistogram(x, y, lbpmaps[index][0], &rfeat[0]);
	ComputeHistogram(x, y, lbpmaps[index][1], &gfeat[0]);
	ComputeHistogram(x, y, lbpmaps[index][2], &bfeat[0]);

	ComputeHistogram(x, y, posmap[index][0], &prfeat[0]);
	ComputeHistogram(x, y, posmap[index][1], &pgfeat[0]);
	ComputeHistogram(x, y, posmap[index][2], &pbfeat[0]);

	ComputeHistogram(x, y, negmap[index][0], &nrfeat[0]);
	ComputeHistogram(x, y, negmap[index][1], &ngfeat[0]);
	ComputeHistogram(x, y, negmap[index][2], &nbfeat[0]);

	if (add) {
		if (nseparate) {// normalize separately each channel and then add them together...
			NormalizeFeatures(rfeat);
			NormalizeFeatures(gfeat);
			NormalizeFeatures(bfeat);
			NormalizeFeatures(prfeat);
			NormalizeFeatures(pgfeat);
			NormalizeFeatures(pbfeat);
			NormalizeFeatures(nrfeat);
			NormalizeFeatures(ngfeat);
			NormalizeFeatures(nbfeat);
			// add red, blue and green channels
			for (UINT k = 0; k < featsize; ++k)
				feat[k] = rfeat[k] + gfeat[k] + bfeat[k];
			for (UINT k = 0; k < featsize; ++k)
				feat[k + featsize] = prfeat[k] + pgfeat[k] + pbfeat[k];
			for (UINT k = 0; k < featsize; ++k)
				feat[k + 2 * featsize] = nrfeat[k] + ngfeat[k] + nbfeat[k];
		} else {
			for (UINT k = 0; k < featsize; ++k)
				feat[k] = rfeat[k] + gfeat[k] + bfeat[k];
			NormalizeFeatures(feat);
			for (UINT k = 0; k < featsize; ++k)
				feat[k + featsize] = prfeat[k] + pgfeat[k] + pbfeat[k];
			NormalizeFeatures(&feat[featsize]);
			for (UINT k = 0; k < featsize; ++k)
				feat[k + 2 * featsize] = nrfeat[k] + ngfeat[k] + nbfeat[k];
			NormalizeFeatures(&feat[2 * featsize]);
		}
	} else {
		cerr << endl
				<< "This Option is Not defined Concat of Normalized Features LTPFeaturesRGB "
				<< endl;
		exit(EXIT_FAILURE);
		/*NormalizeFeatures(rfeat);
		 NormalizeFeatures(gfeat);
		 NormalizeFeatures(bfeat);
		 // combine the three channels
		 UINT size = nxcells * nycells * nbins;
		 copy(rfeat.begin(), rfeat.end(), feat.begin());
		 copy(gfeat.begin(), gfeat.end(), feat.begin() + size);
		 copy(bfeat.begin(), bfeat.end(), feat.begin() + 2 * size);*/
	}
}
#else
void LTPFeaturesRGB::GetFeatures(UINT index, int x, int y, vector<REAL>&feat) {

	UINT offset = 0;
	fill(feat.begin(), feat.end(), 0);// Initialize with 0's;

	x = MAX(0,x - radius); // to cater for offset
	y = MAX(0,y - radius);
	UINT twidth = x + width < lbpmaps[index][0].nxpoints ? x + width
	: lbpmaps[index][0].nxpoints, theight = y + height
	< lbpmaps[index][0].nypoints ? y + height
	: lbpmaps[index][0].nypoints;
	vector<REAL> gfeat(feat.size(), 0), bfeat(feat.size(), 0);
	if (add) {
		if (nseparate)// if normalize separately...
		{
			for (UINT i = y; i < theight; i += lbpstride)
			for (UINT j = x; j < twidth; j += lbpstride) {
				lbpmaps[index][0].ComputeHist(j, i, cellsize, cellsize,
						&feat[offset]);
				NormalizeFeatures(&feat[offset]);
				lbpmaps[index][1].ComputeHist(j, i, cellsize, cellsize,
						&gfeat[offset]);
				NormalizeFeatures(&gfeat[offset]);
				lbpmaps[index][2].ComputeHist(j, i, cellsize, cellsize,
						&bfeat[offset]);
				NormalizeFeatures(&bfeat[offset]);
				offset += nbins;

				posmap[index][0].ComputeHist(j, i, cellsize, cellsize,
						&feat[offset]);
				NormalizeFeatures(&feat[offset]);
				posmap[index][1].ComputeHist(j, i, cellsize, cellsize,
						&gfeat[offset]);
				NormalizeFeatures(&gfeat[offset]);
				posmap[index][2].ComputeHist(j, i, cellsize, cellsize,
						&bfeat[offset]);
				NormalizeFeatures(&bfeat[offset]);
				offset += nbins;

				negmap[index][0].ComputeHist(j, i, cellsize, cellsize,
						&feat[offset]);
				NormalizeFeatures(&feat[offset]);
				negmap[index][1].ComputeHist(j, i, cellsize, cellsize,
						&gfeat[offset]);
				NormalizeFeatures(&gfeat[offset]);
				negmap[index][2].ComputeHist(j, i, cellsize, cellsize,
						&bfeat[offset]);
				NormalizeFeatures(&bfeat[offset]);
				offset += nbins;
			}

			// add blue and green channels
			for (UINT k = 0; k < feat.size(); ++k)
			feat[k] += gfeat[k] + bfeat[k];

		} else {
			for (UINT i = y; i < theight; i += lbpstride)
			for (UINT j = x; j < twidth; j += lbpstride) {
				lbpmaps[index][0].ComputeHist(j, i, cellsize, cellsize,
						&feat[offset]);
				lbpmaps[index][1].ComputeHist(j, i, cellsize, cellsize,
						&feat[offset]);
				lbpmaps[index][2].ComputeHist(j, i, cellsize, cellsize,
						&feat[offset]);
				NormalizeFeatures(&feat[offset]);
				offset += nbins;

				posmap[index][0].ComputeHist(j, i, cellsize, cellsize,
						&feat[offset]);
				posmap[index][1].ComputeHist(j, i, cellsize, cellsize,
						&feat[offset]);
				posmap[index][2].ComputeHist(j, i, cellsize, cellsize,
						&feat[offset]);
				NormalizeFeatures(&feat[offset]);
				offset += nbins;

				negmap[index][0].ComputeHist(j, i, cellsize, cellsize,
						&feat[offset]);
				negmap[index][1].ComputeHist(j, i, cellsize, cellsize,
						&feat[offset]);
				negmap[index][2].ComputeHist(j, i, cellsize, cellsize,
						&feat[offset]);
				NormalizeFeatures(&feat[offset]);
				offset += nbins;
			}
		}
	} else {
		for (UINT i = y; i < theight; i += cellsize)
		for (UINT j = x; j < twidth; j += cellsize) {
			// For LBP
			lbpmaps[index][0].ComputeHist(j, i, cellsize, cellsize,
					&feat[offset]);
			NormalizeFeatures(&feat[offset]);
			offset += nbins;
			lbpmaps[index][1].ComputeHist(j, i, cellsize, cellsize,
					&feat[offset]);
			NormalizeFeatures(&feat[offset]);
			offset += nbins;
			lbpmaps[index][2].ComputeHist(j, i, cellsize, cellsize,
					&feat[offset]);
			NormalizeFeatures(&feat[offset]);
			offset += nbins;

			// FOR LTP +ve map
			posmap[index][0].ComputeHist(j, i, cellsize, cellsize,
					&feat[offset]);
			NormalizeFeatures(&feat[offset]);
			offset += nbins;
			posmap[index][1].ComputeHist(j, i, cellsize, cellsize,
					&feat[offset]);
			NormalizeFeatures(&feat[offset]);
			offset += nbins;
			posmap[index][2].ComputeHist(j, i, cellsize, cellsize,
					&feat[offset]);
			NormalizeFeatures(&feat[offset]);
			offset += nbins;

			//FOR LTP -ve map
			negmap[index][0].ComputeHist(j, i, cellsize, cellsize,
					&feat[offset]);
			NormalizeFeatures(&feat[offset]);
			offset += nbins;
			negmap[index][1].ComputeHist(j, i, cellsize, cellsize,
					&feat[offset]);
			NormalizeFeatures(&feat[offset]);
			offset += nbins;
			negmap[index][2].ComputeHist(j, i, cellsize, cellsize,
					&feat[offset]);
			NormalizeFeatures(&feat[offset]);
			offset += nbins;
		}
	}
}
#endif
void LTPFeaturesRGB::GetFoldedFeatures(UINT index, int x, int y,
		vector<REAL>&feat) {

	//	UINT offset=0;
	//	x = MAX(0,x - radius); // to cater for offset
	//	y = MAX(0,y - radius);
	//	UINT	twidth = x+width < lbpmaps[index].nxpoints ?
	//			x+width : lbpmaps[index].nxpoints,
	//			theight = y+height < lbpmaps[index].nypoints ?
	//					y+height : lbpmaps[index].nypoints;
	//	fill(feat.begin(),feat.end(),0);// Initialize with 0's;
	//	twidth /= 2; //half of the width
	//	int fx =  x + width - lbpstride; // flipped x index
	//	vector<REAL> tfeat(nbins,0),tpfeat(nbins,0),tnfeat(nbins,0);
	//	for(UINT i=y ; i < theight; i+=lbpstride)
	//	{
	//		fx =  x + width - lbpstride;
	//		for(UINT j=x ; j < twidth; j+=lbpstride)
	//		{
	//			// LBP Features
	//			lbpmaps[index].ComputeHist(j,i,cellsize,cellsize,&feat[offset]);
	//			NormalizeFeatures(&feat[offset]);
	//
	//			lbpmaps[index].ComputeHist(fx,i,cellsize,cellsize,&tfeat[0]);
	//			NormalizeFeatures(&tfeat[0]);
	//
	//			// LTP +Ve features
	//			posmap[index].ComputeHist(j,i,cellsize,cellsize,&feat[offset+nbins]);
	//			NormalizeFeatures(&feat[offset+nbins]);
	//
	//			posmap[index].ComputeHist(fx,i,cellsize,cellsize,&tpfeat[0]);
	//			NormalizeFeatures(&tpfeat[0]);
	//
	//
	//			// LTP +Ve features
	//			negmap[index].ComputeHist(j,i,cellsize,cellsize,&feat[offset+2*nbins]);
	//			NormalizeFeatures(&feat[offset+2*nbins]);
	//
	//			negmap[index].ComputeHist(fx,i,cellsize,cellsize,&tnfeat[0]);
	//			NormalizeFeatures(&tnfeat[0]);
	//
	//
	//
	//			for(UINT iter=0; iter < nbins;++iter)
	//			{
	//				feat[offset+iter] += tfeat[iter];
	//				feat[offset+nbins+iter] += tpfeat[iter];
	//				feat[offset+2*nbins+iter] += tnfeat[iter];
	//			}
	//
	//			fill(tfeat.begin(),tfeat.end(),0);// Initialize with 0's;
	//			fill(tpfeat.begin(),tpfeat.end(),0);// Initialize with 0's;
	//			fill(tnfeat.begin(),tnfeat.end(),0);// Initialize with 0's;
	//
	//			offset += 3*nbins;
	//			fx -= lbpstride;
	//		}
	//	}
}

void LTPFeaturesRGB::UnFoldWeights(REAL *input, vector<REAL>&weights) {
	UINT twidth = width / 2, offset = 0, count = 0, foffset = 0, boffset = 0,
			inoffset = 0;// forward & backward offsets

	for (UINT i = 0; i < height; i += lbpstride, ++count) {
		boffset = ((width - lbpstride) / lbpstride) * nbins * 3;
		foffset = 0;
		offset = count * (width / lbpstride) * nbins * 3;
		for (UINT j = 0; j < twidth; j += lbpstride) {
			for (UINT iter = 0; iter < 3 * nbins; ++iter) {
				weights[iter + offset + foffset] = input[inoffset + iter];
				weights[iter + offset + boffset] = input[inoffset + iter];
			}

			foffset += nbins * 3;
			boffset -= nbins * 3;
			inoffset += nbins * 3;
		}

	}

}

// for fast computations....
//void LTPFeaturesRGB::ComputeHist(UINT index,UINT x, UINT y,
//		UINT width,UINT height,REAL *feat)
//{
//
//	height = MIN(height+y,lbpmaps[index].nypoints);
//	width  = MIN(width+x,lbpmaps[index].nxpoints);
//	UINT nxpoints=lbpmaps[index].nxpoints;
//	REAL *lbpmap=&lbpmaps[index].lbpmap[0],
//	*tpmap = &posmap[index].lbpmap[0],
//	*tnmap = &negmap[index].lbpmap[0];
//	for(UINT i=y; i < height;++i)
//		for(UINT j=x; j < width;++j)
//		{
//			hist[lbpmap[i*nxpoints+j]]++;
//			hist[tpmap[i*nxpoints+j]+nbins]++;
//			hist[tnmap[i*nxpoints+j]+2*nbins]++;
//		}
//}
