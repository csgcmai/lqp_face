/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/
#include "lbpFeaturesRGB.h"

void LBPFeaturesRGB::GenerateMap() {

	UINT niter = unsigned(1 << npoints);
	int index = 0;
	map.resize(niter, 0);
	if (rinvariant) {
		int rimapping[] = { 0, 1, 1, 2, 1, 9, 2, 3, 1, 9, 9, 9, 2, 9, 3, 4, 1,
				9, 9, 9, 9, 9, 9, 9, 2, 9, 9, 9, 3, 9, 4, 5, 1, 9, 9, 9, 9, 9,
				9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 2, 9, 9, 9, 9, 9, 9, 9, 3, 9, 9,
				9, 4, 9, 5, 6, 1, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
				9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 2, 9, 9, 9, 9,
				9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 3, 9, 9, 9, 9, 9, 9, 9, 4, 9,
				9, 9, 5, 9, 6, 7, 1, 2, 9, 3, 9, 9, 9, 4, 9, 9, 9, 9, 9, 9, 9,
				5, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 6, 9, 9, 9, 9,
				9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
				9, 9, 9, 9, 9, 9, 7, 2, 3, 9, 4, 9, 9, 9, 5, 9, 9, 9, 9, 9, 9,
				9, 6, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 7, 3, 4, 9,
				5, 9, 9, 9, 6, 9, 9, 9, 9, 9, 9, 9, 7, 4, 5, 9, 6, 9, 9, 9, 7,
				5, 6, 9, 7, 6, 7, 7, 8 };
		copy(rimapping, rimapping + niter, map.begin());
		return;
	}
	if (nbins != pow(2.0, (double) npoints)) {/*Generate Uniform Patterns*/
		for (UINT i = 0; i < niter; i++) {
			if (Transitions(i, npoints) <= 2) //uniform
				map[i] = index++;
			else
				map[i] = nbins - 1;
		}

	} else {
		for (UINT i = 0; i < niter; i++)
			map[i] = i; // map each bin to itself...
	}
}

int LBPFeaturesRGB::Transitions(UINT c, int nbits) {
	int base = 1;
	int current = c & base, current2, changes = 0;
	for (int i = 1; i < nbits; i++) {
		base <<= 1;
		current2 = (c & base) >> i;
		if (current ^ current2)
			changes++;
		current = current2;
	}
	return changes; //(changes <= 2)? 1 : 0;
}

void LBPFeaturesRGB::ComputeSamplingPoints() {
	REAL step = (M_PI * 2) / npoints;
	REAL x, y;
	int count = 0;
	if (gridtype == Circular && ptype == PT_Circular) {
		UINT j = 0;
		for (UINT i = 0; i < npoints; ++i, j += 10) {
			x = radius * cos((double) (i * step));
			y = -radius * sin((double) (i * step));

			x = ABS(x) <= 1e-6 ? 0 : x;
			y = ABS(y) <= 1e-6 ? 0 : y;

			spoints[i].SetCoord(x, y);

			biinfo[j] = floor(x);
			biinfo[j + 1] = floor(y);
			biinfo[j + 2] = ceil(x);
			biinfo[j + 3] = ceil(y);
			biinfo[j + 4] = x - biinfo[j];
			biinfo[j + 5] = y - biinfo[j + 1];
			biinfo[j + 6] = (1 - biinfo[j + 4]) * (1 - biinfo[j + 5]);//w1
			biinfo[j + 7] = biinfo[j + 4] * (1 - biinfo[j + 5]);//w1
			biinfo[j + 8] = (1 - biinfo[j + 4]) * biinfo[j + 5];//w2
			biinfo[j + 9] = biinfo[j + 4] * biinfo[j + 5];
			dointerpol[i] = !(ABS(x-round(x)) <= 1e-6 && ABS(y-round(y))
					<= 1e-6);
			cout << endl << x << " " << y << " " << biinfo[j] << " "
					<< biinfo[j + 1] << " " << biinfo[j + 2] << " " << biinfo[j
					+ 3] << endl;
		}
		count = spoints.size();
	} else {
		for (int x = -radius; x <= radius; ++x)
			switch (ptype) {
			case PT_Circular:
				for (int y = -radius; y <= radius; ++y)
					if (x != 0 || y != 0)
						spoints[count++].SetCoord(x, y);
				break;
			case PT_HStrip:
				spoints[count++].SetCoord(x, 0);
				break;
			case PT_VStrip:
				spoints[count++].SetCoord(0, x);
				break;
			case PT_Diag:
				spoints[count++].SetCoord(x, x);
				break;
			case PT_ADiag:
				spoints[count++].SetCoord(-x, x);
				break;
			}
	}
	cout << "\n Sample Points = " << count << "\n ";
	for (int i = 0; i < count; ++i) {
		cout << "Point # " << i + 1 << " X : " << spoints[i].GetX() << " Y : "
				<< spoints[i].GetY() << endl;
	}
}
void LBPFeaturesRGB::ComputeGradientMap(Image &image, GenericMap<REAL> &gradmap) {

	UINT width = image.columns(), height = image.rows();
	gradmap.Init(width, height); // for RGB channels
	const PixelPacket *pix = image.getConstPixels(0, 0, width, height);
	REAL factor = 255.0 / (pow(2.0, 16.0) - 1);// scaling factor to have real values[0,1]
	for (int tx = 1; tx < width - 1; tx++) {
		for (int ty = 1; ty < height - 1; ty++) {
			// first color channel
			double dy = (pix[tx + (ty + 1) * width].red - pix[tx + (ty - 1)
					* width].red) * factor, dy2 =
					(pix[tx + (ty + 1) * width].green - pix[tx + (ty - 1)
							* width].green) * factor, dy3 = (pix[tx + (ty + 1)
					* width].blue - pix[tx + (ty - 1) * width].blue) * factor,
					dx = (pix[(tx + 1) + ty * width].red - pix[(tx - 1) + ty
							* width].red) * factor, dx2 = (pix[(tx + 1) + ty
							* width].green - pix[(tx - 1) + ty * width].green)
							* factor, dx3 = (pix[(tx + 1) + ty * width].blue
							- pix[(tx - 1) + ty * width].blue) * factor, v = dx
							* dx + dy * dy, v2 = dx2 * dx2 + dy2 * dy2, v3 =
							dx3 * dx3 + dy3 * dy3;
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
			gradmap(tx, ty) = sqrt(v);
		}
	}
}
void LBPFeaturesRGB::ComputeLBPMap(Image &image, vector<LBPMap>& lbpmap) {
	ColorGray mycolor, c1, c2, c3, c4;
	UINT xdim = image.columns(), ydim = image.rows();
	int nypoints = ydim - 2 * radius, reg = 1, fx, fy, cx, cy, nxpoints = xdim
			- 2 * radius, cenx, ceny, count, rvalue, gvalue, bvalue;
	REAL tmpx, tmpy, fracx, fracy, rpixval, gpixval, bpixval, rcenval, gcenval,
			bcenval;
	REAL x, y, w1, w2, w3, w4;
	lbpmap[0].Init(nxpoints, nypoints); // for RGB channels
	lbpmap[1].Init(nxpoints, nypoints);
	lbpmap[2].Init(nxpoints, nypoints);
	vector<REAL> rimpix(xdim * ydim, 0), bimpix(xdim * ydim, 0), gimpix(
			xdim * ydim, 0);
	//	image.flop();
	pim.Process(image, rimpix, gimpix, bimpix);
	//	int fbins[]={4,3,2,1,0,7,6,5};
	//	ofstream ofile("lbpmap_slow");
	//	ofile << nypoints << "  " << nxpoints;
	for (int i = 0; i < nypoints; ++i) {
		for (int j = 0; j < nxpoints; ++j) {
			cenx = j + radius;
			ceny = i + radius;
			rvalue = gvalue = bvalue = 0;
			count = 0;
			rcenval = rimpix[ceny * xdim + cenx];
			gcenval = gimpix[ceny * xdim + cenx];
			bcenval = bimpix[ceny * xdim + cenx];
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
						rpixval = rimpix[fy * xdim + fx];
						gpixval = gimpix[fy * xdim + fx];
						bpixval = bimpix[fy * xdim + fx];
					} else {
						w1 = (1 - fracx) * (1 - fracy);
						w2 = fracx * (1 - fracy);
						w3 = (1 - fracx) * fracy;
						w4 = fracx * fracy;
						rpixval = w1 * rimpix[fy * xdim + fx] + w2 * rimpix[fy
								* xdim + cx] + w3 * rimpix[cy * xdim + fx] + w4
								* rimpix[cy * xdim + cx];

						gpixval = w1 * gimpix[fy * xdim + fx] + w2 * gimpix[fy
								* xdim + cx] + w3 * gimpix[cy * xdim + fx] + w4
								* gimpix[cy * xdim + cx];

						bpixval = w1 * bimpix[fy * xdim + fx] + w2 * bimpix[fy
								* xdim + cx] + w3 * bimpix[cy * xdim + fx] + w4
								* bimpix[cy * xdim + cx];

					}
				} else {
					rpixval = rimpix[tmpy * xdim + tmpx];
					gpixval = gimpix[tmpy * xdim + tmpx];
					bpixval = bimpix[tmpy * xdim + tmpx];
				}
				int tval = reg << count;
				rvalue += rpixval >= rcenval + lbptolerance ? tval : 0;
				gvalue += gpixval >= gcenval + lbptolerance ? tval : 0;
				bvalue += bpixval >= bcenval + lbptolerance ? tval : 0;

			}
			lbpmap[0](j, i) = map[rvalue];
			lbpmap[1](j, i) = map[gvalue];
			lbpmap[2](j, i) = map[bvalue];
		}
	}
	//	ofile.close();
}

#ifdef CROPPEDWINDOWS
void LBPFeaturesRGB::ComputeLBPMap(UINT index, int xoffset, int yoffset) {
	ColorGray mycolor, c1, c2, c3, c4;
	int nypoints = height - 2 * radius, reg = 1, fx, fy, cx, cy, nxpoints =
			width - 2 * radius, cenx, ceny, count, rvalue, gvalue, bvalue;
	REAL tmpx, tmpy, fracx, fracy, rpixval, gpixval, bpixval, rcenval, gcenval,
			bcenval;
	REAL x, y, w1, w2, w3, w4;
	UINT xdim = xdims[index];
	for (int i = 0; i < nypoints; ++i) {
		for (int j = 0; j < nxpoints; ++j) {
			cenx = j + radius + xoffset;
			ceny = i + radius + yoffset;
			rvalue = gvalue = bvalue = 0;
			count = 0;
			rcenval = rimpix[index][ceny * xdim + cenx];
			gcenval = gimpix[index][ceny * xdim + cenx];
			bcenval = bimpix[index][ceny * xdim + cenx];
			//ofile<<endl<< i<<" " << j<< " "<< rcenval <<" " ;
			for (vector<Point<REAL> >::iterator iter = spoints.begin(); iter
					!= spoints.end(); ++iter, ++count) {
				iter->GetCoord(x, y);

				tmpx = cenx + x;
				tmpy = ceny + y;
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
						rpixval = rimpix[index][fy * xdim + fx];
						gpixval = gimpix[index][fy * xdim + fx];
						bpixval = bimpix[index][fy * xdim + fx];
					} else {
						w1 = (1 - fracx) * (1 - fracy);
						w2 = fracx * (1 - fracy);
						w3 = (1 - fracx) * fracy;
						w4 = fracx * fracy;
						rpixval = w1 * rimpix[index][fy * xdim + fx] + w2
								* rimpix[index][fy * xdim + cx] + w3
								* rimpix[index][cy * xdim + fx] + w4
								* rimpix[index][cy * xdim + cx];

						gpixval = w1 * gimpix[index][fy * xdim + fx] + w2
								* gimpix[index][fy * xdim + cx] + w3
								* gimpix[index][cy * xdim + fx] + w4
								* gimpix[index][cy * xdim + cx];

						bpixval = w1 * bimpix[index][fy * xdim + fx] + w2
								* bimpix[index][fy * xdim + cx] + w3
								* bimpix[index][cy * xdim + fx] + w4
								* bimpix[index][cy * xdim + cx];

					}
				} else {
					rpixval = rimpix[index][tmpy * xdim + tmpx];
					gpixval = gimpix[index][tmpy * xdim + tmpx];
					bpixval = bimpix[index][tmpy * xdim + tmpx];
				}
				int tval = reg << count;
				rvalue += rpixval >= rcenval ? tval : 0;
				gvalue += gpixval >= gcenval ? tval : 0;
				bvalue += bpixval >= bcenval ? tval : 0;

				//ofile<<x<<" "<< y<< " "<<rpixval << " "<<rvalue<<" ";
				///*for flipped case*/						reg << fbins[count] : 0;
			}

			cwlbpmaps[0](j, i) = map[rvalue];
			cwlbpmaps[1](j, i) = map[gvalue];
			cwlbpmaps[2](j, i) = map[bvalue];
		}
	}
}
#endif
void LBPFeaturesRGB::ComputeMeanLBPMap(Image &image, LBPMap lbpmap[]) {
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
	REAL rmeanval, gmeanval, bmeanval;
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
			rmeanval = gmeanval = bmeanval = 0;
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
				rmeanval += rpixval[count];
				gmeanval += gpixval[count];
				bmeanval += bpixval[count];
			}
			rvalue = 0;
			gvalue = 0;
			bvalue = 0;
			int wt = 0;
			rmeanval /= npoints;
			gmeanval /= npoints;
			bmeanval /= npoints;

			for (UINT iter = 0; iter < npoints; ++iter) {
				wt = reg << iter;
				rvalue += ((rpixval[iter] - rmeanval - lbptolerance) >= 0 ? wt
						: 0);
				gvalue += ((gpixval[iter] - gmeanval - lbptolerance) >= 0 ? wt
						: 0);
				bvalue += ((bpixval[iter] - bmeanval - lbptolerance) >= 0 ? wt
						: 0);
			}

			lbpmap[0](j, i) = map[rvalue];
			lbpmap[1](j, i) = map[gvalue];
			lbpmap[2](j, i) = map[bvalue];
		}
	}
}
void LBPFeaturesRGB::ComputeLBPMap(Image &image, LBPMap lbpmap[]) {
	ColorGray mycolor, c1, c2, c3, c4;
	UINT xdim = image.columns(), ydim = image.rows();
	int nypoints = ydim - 2 * radius, reg = 1, fx, fy, cx, cy, nxpoints = xdim
			- 2 * radius, cenx, ceny, count, rvalue, gvalue, bvalue;
	REAL tmpx, tmpy, fracx, fracy, rpixval, gpixval, bpixval, rcenval, gcenval,
			bcenval;
	REAL x, y, w1, w2, w3, w4;
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
			rvalue = gvalue = bvalue = 0;
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
						rpixval = rimpix[fy * xdim + fx];
						gpixval = gimpix[fy * xdim + fx];
						bpixval = bimpix[fy * xdim + fx];
					} else {
						w1 = (1 - fracx) * (1 - fracy);
						w2 = fracx * (1 - fracy);
						w3 = (1 - fracx) * fracy;
						w4 = fracx * fracy;

						rpixval = w1 * rimpix[fy * xdim + fx] + w2 * rimpix[fy
								* xdim + cx] + w3 * rimpix[cy * xdim + fx] + w4
								* rimpix[cy * xdim + cx];

						gpixval = w1 * gimpix[fy * xdim + fx] + w2 * gimpix[fy
								* xdim + cx] + w3 * gimpix[cy * xdim + fx] + w4
								* gimpix[cy * xdim + cx];

						bpixval = w1 * bimpix[fy * xdim + fx] + w2 * bimpix[fy
								* xdim + cx] + w3 * bimpix[cy * xdim + fx] + w4
								* bimpix[cy * xdim + cx];

					}
				} else {
					rpixval = rimpix[tmpy * xdim + tmpx];
					gpixval = gimpix[tmpy * xdim + tmpx];
					bpixval = bimpix[tmpy * xdim + tmpx];
				}
				int tval = reg << count;
				rvalue += rpixval >= rcenval + lbptolerance ? tval : 0;
				gvalue += gpixval >= gcenval + lbptolerance ? tval : 0;
				bvalue += bpixval >= bcenval + lbptolerance ? tval : 0;

				//ofile<<x<<" "<< y<< " "<<rpixval << " "<<rvalue<<" ";
				///*for flipped case*/						reg << fbins[count] : 0;
			}

			lbpmap[0](j, i) = map[rvalue];
			lbpmap[1](j, i) = map[gvalue];
			lbpmap[2](j, i) = map[bvalue];
		}
	}
}
void LBPFeaturesRGB::ComputeLBPMapFast(Image &image, vector<LBPMap> &lbpmap) {
	ColorGray mycolor, c1, c2, c3, c4;
	UINT xdim = image.columns(), ydim = image.rows();
	int nypoints = ydim - 2 * radius, reg = 1, fx, fy, cx, cy, nxpoints = xdim
			- 2 * radius, cenx, ceny, count, rvalue, gvalue, bvalue;
	REAL tmpx, tmpy, fracx, fracy, rpixval, gpixval, bpixval, rcenval, gcenval,
			bcenval;
	REAL x, y, w1, w2, w3, w4;
	lbpmap[0].Init(nxpoints, nypoints); // for RGB channels
	lbpmap[1].Init(nxpoints, nypoints);
	lbpmap[2].Init(nxpoints, nypoints);
	vector<REAL> rimpix(xdim * ydim, 0), bimpix(xdim * ydim, 0), gimpix(
			xdim * ydim, 0);
	pim.Process(image, rimpix, gimpix, bimpix);
	//	ofstream ofile("lbpmap_fast");
	//	ofile << nypoints << "  " << nxpoints;
	//		ofstream ofile("test_map.lst");
	for (int i = 0; i < nypoints; ++i) {
		for (int j = 0; j < nxpoints; ++j) {
			cenx = j + radius;
			ceny = i + radius;
			rvalue = gvalue = bvalue = 0;
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
				int tval = reg << count;
				if (rpixval >= rcenval + lbptolerance)
					rvalue += tval;
				if (gpixval >= gcenval + lbptolerance)
					gvalue += tval;
				if (bpixval >= bcenval + lbptolerance)
					bvalue += tval;

			}
			//			ofile << " " << map[rvalue] << " ";
			lbpmap[0](j, i) = map[rvalue];
			lbpmap[1](j, i) = map[gvalue];
			lbpmap[2](j, i) = map[bvalue];
		}
	}
}
void LBPFeaturesRGB::ComputeLBPMapFast(Image &image, LBPMap lbpmap[]) {
	ColorGray mycolor, c1, c2, c3, c4;
	UINT xdim = image.columns(), ydim = image.rows();
	int nypoints = ydim - 2 * radius, reg = 1, fx, fy, cx, cy, nxpoints = xdim
			- 2 * radius, cenx, ceny, count, rvalue, gvalue, bvalue;
	REAL tmpx, tmpy, fracx, fracy, rpixval, gpixval, bpixval, rcenval, gcenval,
			bcenval;
	REAL x, y, w1, w2, w3, w4;

	lbpmap[0].Init(nxpoints, nypoints); // for RGB channels
	lbpmap[1].Init(nxpoints, nypoints);
	lbpmap[2].Init(nxpoints, nypoints);
	vector<REAL> rimpix(xdim * ydim, 0), bimpix(xdim * ydim, 0), gimpix(
			xdim * ydim, 0);
	pim.Process(image, rimpix, gimpix, bimpix);

	//		ofstream ofile("test_map.lst");
	for (int i = 0; i < nypoints; ++i) {
		for (int j = 0; j < nxpoints; ++j) {
			cenx = j + radius;
			ceny = i + radius;
			rvalue = gvalue = bvalue = 0;
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

				if (ABS(biinfo[pcount+4]) < 1e-6 && ABS(biinfo[pcount+5])
						< 1e-6)
				// no interpolation needed
				{
					rpixval = rimpix[fy * xdim + fx];
					gpixval = gimpix[fy * xdim + fx];
					bpixval = bimpix[fy * xdim + fx];
				} else {
					w1 = biinfo[pcount + 6];
					w2 = biinfo[pcount + 7];
					w3 = biinfo[pcount + 8];
					w4 = biinfo[pcount + 9];
					rpixval = w1 * rimpix[fy * xdim + fx] + w2 * rimpix[fy
							* xdim + cx] + w3 * rimpix[cy * xdim + fx] + w4
							* rimpix[cy * xdim + cx];

					gpixval = w1 * gimpix[fy * xdim + fx] + w2 * gimpix[fy
							* xdim + cx] + w3 * gimpix[cy * xdim + fx] + w4
							* gimpix[cy * xdim + cx];

					bpixval = w1 * bimpix[fy * xdim + fx] + w2 * bimpix[fy
							* xdim + cx] + w3 * bimpix[cy * xdim + fx] + w4
							* bimpix[cy * xdim + cx];

				}
				int tval = reg << count;
				if (rpixval >= rcenval + lbptolerance)
					rvalue += tval;
				if (gpixval >= gcenval + lbptolerance)
					gvalue += tval;
				if (bpixval >= bcenval + lbptolerance)
					bvalue += tval;
			}
			lbpmap[0](j, i) = map[rvalue];
			lbpmap[1](j, i) = map[gvalue];
			lbpmap[2](j, i) = map[bvalue];
		}
	}
}

void LBPFeaturesRGB::InitalizeMapsFast(Image& image, PyramidType sspace_) {
	UINT pycount = 0;
	// Read the image....
	sspace = sspace_;
	if (sspace_ == SingleRes) {
		if (!IsFolded()) {
			pycount = pyobj.GeneratePyramid(image);// for 3 Channels RGB
			lbpmaps.resize(pycount, vector<LBPMap> (3));
			for (UINT i = 0; i < pycount; i++) {
				Image& imgref = pyobj.GetImageRef(i);
				// Compute the Gradient Image
				ComputeLBPMapFast(imgref, lbpmaps[i]);
			}
		} else {
		}
	} else {
		if (!IsFolded()) {
			lbpmaps.resize(1, vector<LBPMap> (3));
			ComputeLBPMapFast(image, lbpmaps[0]);
		} else {
		}
	}

}
void LBPFeaturesRGB::InitalizeMaps(Image &image) {
	lbpmaps.resize(1, vector<LBPMap> (3));
	ComputeLBPMap(image, lbpmaps[0]);
	if (usegradmag) {
		gradientmaps.resize(1);
		ComputeGradientMap(image, gradientmaps[0]);
	}
}
void LBPFeaturesRGB::InitalizeMaps(Pyramid &pyobj, PyramidType sspace_) {
	UINT pycount = 0;
	// Read the image....
	sspace = sspace_;
	if (sspace_ == SingleRes) {
		if (!IsFolded()) {
			pycount = pyobj.GetNLevels();
			lbpmaps.resize(pycount, vector<LBPMap> (3));
			gradientmaps.resize(pycount);
			for (UINT i = 0; i < pycount; i++) {
				Image& imgref = pyobj.GetImageRef(i);
				// Compute the Gradient Image
				ComputeLBPMap(imgref, lbpmaps[i]);
				if (usegradmag)
					ComputeGradientMap(imgref, gradientmaps[i]);
			}
		} else {
		}
	}
}
void LBPFeaturesRGB::InitalizeMaps(Image& image, PyramidType sspace_) {
	UINT pycount = 0;
	// Read the image....
	sspace = sspace_;
	if (sspace_ == SingleRes) {
		if (!IsFolded()) {
			pycount = pyobj.GeneratePyramid(image);// for 3 Channels RGB
#ifndef CROPPEDWINDOWS
			lbpmaps.resize(pycount, vector<LBPMap> (3));
#else
			xdims.resize(pycount);
			rimpix.resize(pycount);
			gimpix.resize(pycount);
			bimpix.resize(pycount);
#endif
			gradientmaps.resize(pycount);
			for (UINT i = 0; i < pycount; i++) {
				Image& imgref = pyobj.GetImageRef(i);
				// Compute the Gradient Image
#ifndef CROPPEDWINDOWS
				ComputeLBPMap(imgref, lbpmaps[i]);
#else
				UINT dim = imgref.rows() * imgref.columns();
				xdims[i] = imgref.columns();
				rimpix[i].resize(dim);
				gimpix[i].resize(dim);
				bimpix[i].resize(dim);
				pim .Process(imgref, rimpix[i], gimpix[i], bimpix[i]);
#endif
				if (usegradmag)
					ComputeGradientMap(imgref, gradientmaps[i]);
			}
		} else {
		}
	} else {
		if (!IsFolded()) {

#ifndef CROPPEDWINDOWS
			lbpmaps.resize(1, vector<LBPMap> (3));
			ComputeLBPMap(image, lbpmaps[0]);
#else
			xdims.resize(1);
			rimpix.resize(1);
			gimpix.resize(1);
			bimpix.resize(1);
			UINT dim = image.rows() * image.columns();
			xdims[0] = image.columns();
			rimpix[0].resize(dim);
			gimpix[0].resize(dim);
			bimpix[0].resize(dim);
			pim.Process(image, rimpix[0], gimpix[0], bimpix[0]);
#endif
			if (usegradmag) {
				gradientmaps.resize(1);
				ComputeGradientMap(image, gradientmaps[0]);
			}
		} else {
		}
	}
}

#ifndef OLD_CODE
void LBPFeaturesRGB::ComputeHistogram(int initx, int inity, LBPMap &lbpmap,
		REAL *hist) /*Height & Width refer to window width & height, initx, inity are the histogram offsets...*/{

	if (hmethod == HM_Discrete) {
#if 0
		/*Fastest
		 * but if initx < radius,
		 * initx is ignored and the grid is placed over the image in range [radius,width+radius]
		 * In other words feature computed for the range initx < radius are same as initx=radius, irrespective of initx value
		 * */
		// Note: Normalization routine has been updated so
		UINT offset = 0;
		initx = MAX(0,initx - radius); // to cater for offset
		inity = MAX(0,inity - radius);
		UINT twidth = MIN(initx+width, lbpmap.nxpoints), theight =
		MIN(inity+height , lbpmap.nypoints );
		for (UINT i = inity; i < theight; i += cellsize)
		for (UINT j = initx; j < twidth; j += cellsize) {
			lbpmap.ComputeHist(i, j, cellsize, cellsize, hist + offset);
			offset += nbins;
		}
#endif
		int tx, ty;
		for (int y = 0; y < height; ++y) {
			/*Boundary Points are included 2 times (due to boundary extension) in the boundary cell-histograms*/
			/*LBP Map is translated by radius*/
			ty = MIN(MAX(inity + y - radius,0),lbpmap.nypoints-1);
			int iyp = y / nycells * nxcells;
			for (int x = 0; x < width; ++x) {
				tx = MIN(MAX(initx + x - radius,0),lbpmap.nxpoints-1);
				int ixp = x / nxcells;
				int o = lbpmap.GetValue(tx, ty);
				*(hist + ixp + iyp + o * nycells * nxcells) += 1;
			}
		}
	} else {
		ComputeHistBilinear(initx, inity, lbpmap, hist);
	}
}
void LBPFeaturesRGB::ComputeHistBilinear(int initx, int inity, LBPMap &lbpmap,
		REAL *hist) {
	int tx, ty;
	for (int y = 0; y < height; ++y) {
		/*Boundary Points are included 2 times (due to boundary extension) in the boundary cell-histograms*/
		ty = MIN(MAX(inity + y - radius,0),lbpmap.nypoints-1);
		for (int x = 0; x < width; ++x) {
			tx = MIN(MAX(initx + x - radius,0),lbpmap.nxpoints-1);
			/*because the value of map is offset by radius, value of x at x-radius */
			int o = lbpmap.GetValue(tx, ty);
			// add to 4 histograms around pixel using linear interpolation
			double xp = ((double) x + 0.5) / (double) cellsize - 0.5; // cell interpolation grid
			double yp = ((double) y + 0.5) / (double) cellsize - 0.5;
			int ixp = (int) floor(xp);
			int iyp = (int) floor(yp);
			double vx0 = xp - ixp;
			double vy0 = yp - iyp;
			double vx1 = 1.0 - vx0;
			double vy1 = 1.0 - vy0;
			double v = 1; // pooling of gradient magnitude can be an option

			if (ixp >= 0 && iyp >= 0) {
				*(hist + ixp + nxcells * iyp + o * nycells * nxcells) += vx1
						* vy1 * v;
			}

			if (ixp + 1 < nxcells && iyp >= 0) {
				*(hist + (ixp + 1) + nxcells * iyp + o * nycells * nxcells)
						+= vx0 * vy1 * v;
			}

			if (ixp >= 0 && iyp + 1 < nycells) {
				*(hist + ixp + nxcells * (iyp + 1) + o * nycells * nxcells)
						+= vx1 * vy0 * v;
			}

			if (ixp + 1 < nxcells && iyp + 1 < nycells) {
				*(hist + (ixp + 1) + nxcells * (iyp + 1) + o * nycells
						* nxcells) += vx0 * vy0 * v;
			}
		}
	}
}

void LBPFeaturesRGB::ComputeHistBilinear(int initx, int inity, LBPMap &lbpmap,
		GenericMap<REAL> &gradmap, REAL *hist) {
	int tx, ty;
	for (int y = 0; y < height; ++y) {
		/*Boundary Points are included 2 times (due to boundary extension) in the boundary cell-histograms*/
		ty = MIN(MAX(inity + y - radius,0),lbpmap.nypoints-1);
		for (int x = 0; x < width; ++x) {
			tx = MIN(MAX(initx + x - radius,0),lbpmap.nxpoints-1);
			/*because the value of map is offset by radius, value of x at x-radius */
			int o = lbpmap.GetValue(tx, ty);
			// add to 4 histograms around pixel using linear interpolation
			double xp = ((double) x + 0.5) / (double) cellsize - 0.5; // cell interpolation grid
			double yp = ((double) y + 0.5) / (double) cellsize - 0.5;
			int ixp = (int) floor(xp);
			int iyp = (int) floor(yp);
			double vx0 = xp - ixp;
			double vy0 = yp - iyp;
			double vx1 = 1.0 - vx0;
			double vy1 = 1.0 - vy0;
			double v = gradmap.GetValue(MIN(initx+x,gradmap.nxpoints-1),
					MIN(inity+y,gradmap.nypoints-1)); // pooling gradient magnitude
			if (ixp >= 0 && iyp >= 0) {
				*(hist + ixp + nxcells * iyp + o * nycells * nxcells) += vx1
						* vy1 * v;
			}

			if (ixp + 1 < nxcells && iyp >= 0) {
				*(hist + (ixp + 1) + nxcells * iyp + o * nycells * nxcells)
						+= vx0 * vy1 * v;
			}

			if (ixp >= 0 && iyp + 1 < nycells) {
				*(hist + ixp + nxcells * (iyp + 1) + o * nycells * nxcells)
						+= vx1 * vy0 * v;
			}

			if (ixp + 1 < nxcells && iyp + 1 < nycells) {
				*(hist + (ixp + 1) + nxcells * (iyp + 1) + o * nycells
						* nxcells) += vx0 * vy0 * v;
			}
		}
	}

}
void LBPFeaturesRGB::GetFeatures(UINT index, int x, int y, vector<REAL>&feat) {
	// offset in ComputeHist
	UINT featsize = nxcells * nycells * nbins;
	fill(feat.begin(), feat.end(), 0);// Initialize with 0's;
	vector<REAL> rfeat(featsize, 0), gfeat(featsize, 0), bfeat(featsize, 0);
	/*RGB channels based code*/
	if (!usegradmag) {
#ifdef CROPPEDWINDOWS
		ComputeLBPMap(index, 0, 0);
		ComputeHistogram(x, y, cwlbpmaps[0], &rfeat[0]);
		ComputeHistogram(x, y, cwlbpmaps[1], &gfeat[0]);
		ComputeHistogram(x, y, cwlbpmaps[2], &bfeat[0]);
#else
		ComputeHistogram(x, y, lbpmaps[index][0], &rfeat[0]);
		ComputeHistogram(x, y, lbpmaps[index][1], &gfeat[0]);
		ComputeHistogram(x, y, lbpmaps[index][2], &bfeat[0]);

#endif
	} else {
		ComputeHistBilinear(x, y, lbpmaps[index][0], gradientmaps[index],
				&rfeat[0]);
		ComputeHistBilinear(x, y, lbpmaps[index][1], gradientmaps[index],
				&gfeat[0]);
		ComputeHistBilinear(x, y, lbpmaps[index][2], gradientmaps[index],
				&bfeat[0]);
	}
	if (add) {
		if (nseparate) {// normalize separately each channel and then add them together...
			NormalizeFeatures(rfeat);
			NormalizeFeatures(gfeat);
			NormalizeFeatures(bfeat);
			// add reg,blue and green channels
			for (UINT k = 0; k < featsize; ++k)
				feat[k] = rfeat[k] + gfeat[k] + bfeat[k];

		} else {
			for (UINT k = 0; k < featsize; ++k)
				feat[k] = rfeat[k] + gfeat[k] + bfeat[k];
			NormalizeFeatures(feat);
		}
	} else {
		NormalizeFeatures(rfeat);
		NormalizeFeatures(gfeat);
		NormalizeFeatures(bfeat);
		// combine the three channels

		copy(rfeat.begin(), rfeat.end(), feat.begin());
		copy(gfeat.begin(), gfeat.end(), feat.begin() + featsize);
		copy(bfeat.begin(), bfeat.end(), feat.begin() + 2 * featsize);

	}
}
void LBPFeaturesRGB::NormalizeFeatures(REAL *feat) {
	vector<REAL> normvec(nxcells * nycells, 0);
	REAL *normptr = &normvec[0];

	for (UINT o = 0; o < nbins; o++) {
		REAL *src = feat + o * nycells * nxcells;
		REAL *dst = normptr;
		REAL *end = normptr + nxcells * nycells;
		while (dst < end) {
			if (norm == L1 || norm == L1sqrt || norm == LOG2Norm) // for optimization
				*(dst++) += *src;
			else
				*(dst++) += *src * *src;
			src++;
		}
	}
	vector<REAL> hist(nxcells * nycells * nbins, 0);
	REAL *histptr = &hist[0];
	REAL epsilon = 1e-3;
	for (UINT x = 0; x < nxcells; x++) {
		for (UINT y = 0; y < nycells; y++) {
			REAL *src = feat + x + y * nxcells;
			REAL *dst = histptr + nbins * (x + y * nxcells); //
			// single cell normalization...
			REAL normval = *(normptr + x + y * nxcells);
			if (norm == L2Hystersis || norm == L2)
				normval = sqrt(normval + epsilon * epsilon);
			for (UINT o = 0; o < nbins; o++) {
				switch (norm) {
				case L1:
					*dst = (*src / (normval + epsilon));
					break;
				case L2:
					*dst = *src / normval;
					break;
				case L2Hystersis:
					*dst = MIN(*src /normval ,0.2);
					break;
				case LOG2Norm:
					*dst = log2(*src / (normval + epsilon) + 1);
					break;
				case LOG2:
					*dst = log2(*src + 1);
					break;
				case L1sqrt:
				default:
					*dst = sqrt(*src / (normval + epsilon));
				}
				src += nxcells * nycells;
				++dst;
			}
		}
	}
	//TODO:
	copy(hist.begin(), hist.end(), feat);
}
#else
void LBPFeaturesRGB::GetFeatures(UINT index, int x, int y, vector<REAL>&feat) {
	UINT offset = 0;
	x = MAX(0,x - radius); // to cater for offset
	y = MAX(0,y - radius);
	UINT twidth = x + width < lbpmaps[index][0].nxpoints ? x + width
	: lbpmaps[index][0].nxpoints, theight = y + height
	< lbpmaps[index][0].nypoints ? y + height
	: lbpmaps[index][0].nypoints;
	fill(feat.begin(), feat.end(), 0);// Initialize with 0's;
	vector<REAL> gfeat(feat.size(), 0), bfeat(feat.size(), 0);
	if (add) {
		if (nseparate) {// normalize separately each channel and then add them together...
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
			}
		}
	} else {
		for (UINT i = y; i < theight; i += lbpstride)
		for (UINT j = x; j < twidth; j += lbpstride) {
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
		}
	}
}
void LBPFeaturesRGB::NormalizeFeatures(REAL *feat) {
	REAL normval = 0, EPS = 1e-3;
	if (norm == L1 || norm == L1sqrt || norm == LOG2Norm) {
		for (UINT k = 0; k < nbins; ++k)
		normval += feat[k];
		normval += EPS;
		if (norm == L1)
		for (UINT k = 0; k < nbins; ++k)
		feat[k] /= normval;
		else if (norm == LOG2Norm)
		for (UINT k = 0; k < nbins; ++k)
		feat[k] = log2(feat[k] / normval + 1);
		else
		for (UINT k = 0; k < nbins; ++k)
		feat[k] = sqrt(feat[k] / normval);
	} else if (norm == LOG2) {
		for (UINT k = 0; k < nbins; ++k)
		feat[k] = log2(feat[k] + 1);
	} else {
		EPS = 1e-4;
		for (UINT k = 0; k < nbins; ++k)
		normval += feat[k] * feat[k];
		normval = sqrt(normval) + EPS;

		for (UINT k = 0; k < nbins; ++k)
		feat[k] /= normval;
	}
}
#endif

void LBPFeaturesRGB::GetFoldedFeatures(UINT index, int x, int y,
		vector<REAL>&feat) {

}
void LBPFeaturesRGB::UnFoldWeights(REAL *input, vector<REAL>&weights) {
	UINT count = 0;// forward & backward offsets
	UINT xblocks = width / cellsize, yblocks = height / cellsize, fxblocks =
			xblocks / 2.0, cxblocks = ceil((REAL) xblocks / 2.0);
	int findex;
	for (UINT i = 0; i < yblocks; ++i, ++count) {
		findex = xblocks - 1;
		for (UINT j = 0; j < fxblocks; ++j, --findex) {
			for (UINT iter = 0; iter < nbins; ++iter) {
				weights[(i * xblocks + j) * nbins + iter] = input[(i * cxblocks
						+ j) * nbins + iter];
				weights[(i * xblocks + findex) * nbins + iter] = input[(i
						* cxblocks + j) * nbins + iter];
			}

		}
		if (fxblocks != cxblocks) // for odd width window
		{
			for (UINT iter = 0; iter < nbins; ++iter)
				weights[(i * xblocks + fxblocks) * nbins + iter] = input[(i
						* cxblocks + fxblocks) * nbins + iter];
		}
	}
}

// to be  used fo the LTP's and HOG-LTP's...
void LBPFeaturesRGB::UnFoldWeights(REAL *input, REAL *weights) {
	UINT count = 0;// forward & backward offsets
	UINT xblocks = width / cellsize, yblocks = height / cellsize, fxblocks =
			xblocks / 2.0, cxblocks = ceil((REAL) xblocks / 2.0);
	int findex;
	for (UINT i = 0; i < yblocks; ++i, ++count) {
		findex = xblocks - 1;
		for (UINT j = 0; j < fxblocks; ++j, --findex) {
			for (UINT iter = 0; iter < nbins; ++iter) {
				weights[(i * xblocks + j) * nbins + iter] = input[(i * cxblocks
						+ j) * nbins + iter];
				weights[(i * xblocks + findex) * nbins + iter] = input[(i
						* cxblocks + j) * nbins + iter];
			}
		}
		if (fxblocks != cxblocks) // for odd width window
		{
			for (UINT iter = 0; iter < nbins; ++iter)
				weights[(i * xblocks + fxblocks) * nbins + iter] = input[(i
						* cxblocks + fxblocks) * nbins + iter];
		}
	}
}
