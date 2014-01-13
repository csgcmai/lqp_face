/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/
#include "lbpFeatures.h"

void LBPFeatures::GenerateMap() {

	UINT niter = unsigned(1 << npoints);
	int index = 0;
	map.resize(niter, 0);
	if (nbins != 256) {
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
int LBPFeatures::Transitions(UINT c, int nbits) {
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

void LBPFeatures::ComputeSamplingPoints() {
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
void LBPFeatures::ComputeGradientMap(Image &image, GenericMap<REAL> &gradmap) {

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

void LBPFeatures::ComputeLBPMap(Image &image, LBPMap &lbpmap) {
	UINT xdim = image.columns(), ydim = image.rows();
	int nypoints = ydim - 2 * radius, reg = 1, fx, fy, cx, cy, nxpoints = xdim
			- 2 * radius, cenx, ceny, count, value;
	REAL tmpx, tmpy, fracx, fracy, pixval, cenval;
	REAL x, y;
	lbpmap.Init(nxpoints, nypoints);
	vector<REAL> impix(xdim * ydim, 0);
	//	image.flop();
	pim.Process(image, impix);
	//	int fbins[]={4,3,2,1,0,7,6,5};

	for (int i = 0; i < nypoints; ++i) {
		for (int j = 0; j < nxpoints; ++j) {
			cenx = j + radius;
			ceny = i + radius;
			value = 0;
			count = 0;
			cenval = impix[ceny * xdim + cenx];
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
						pixval = impix[fy * xdim + fx];
					else
						pixval = (1 - fracx) * (1 - fracy) * impix[fy * xdim
								+ fx] + fracx * (1 - fracy) * impix[fy * xdim
								+ cx] + (1 - fracx) * fracy * impix[cy * xdim
								+ fx] + fracx * fracy * impix[cy * xdim + cx];
				} else {
					pixval = impix[tmpy * xdim + tmpx];
					//					pixval[1] = gimpix[tmpy * xdim + tmpx];
					//					pixval[2] = bimpix[tmpy * xdim + tmpx];
				}
				value += pixval >= cenval ? reg << count : 0;

			}
			lbpmap(j, i) = map[value];
		}
	}
}
/*
 void LBPFeatures::ComputeLBPMap(Image &image, LBPMap &lbpmap, LBPMap &flbpmap) {
 ColorGray mycolor, c1, c2, c3, c4;
 UINT xdim = image.columns(), ydim = image.rows();
 int nypoints = ydim - 2 * radius, reg = 1, fx, fy, cx, cy, nxpoints = xdim
 - 2 * radius, cenx, ceny, count, value, fvalue;
 REAL tmpx, tmpy, fracx, fracy, pixval, cenval;
 REAL x, y;
 lbpmap.Init(nxpoints, nypoints);
 flbpmap.Init(nxpoints, nypoints);
 vector<REAL> impix(xdim * ydim, 0);
 pim.Process(image, impix);
 int fbins[] = { 4, 3, 2, 1, 0, 7, 6, 5 };
 for (int i = 0; i < nypoints; ++i) {
 for (int j = 0; j < nxpoints; ++j) {
 cenx = j + radius;
 ceny = i + radius;
 value = fvalue = 0;
 count = 0;
 cenval = impix[ceny * xdim + cenx];
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

 if (ABS(fracx) < 1e-6 && ABS(fracy) < 1e-6)
 // no interpolation needed
 pixval = impix[fy * xdim + fx];
 else
 pixval = (1 - fracx) * (1 - fracy) * impix[fy * xdim + fx]
 + fracx * (1 - fracy) * impix[fy * xdim + cx] + (1
 - fracx) * fracy * impix[cy * xdim + fx] + fracx
 * fracy * impix[cy * xdim + cx];

 if (pixval >= cenval) {
 value += reg << count;
 fvalue += reg << fbins[count];
 }
 }
 lbpmap(j, i) = map[value];
 flbpmap(j, i) = map[fvalue];
 }
 }
 }
 */
void LBPFeatures::InitalizeMaps(Image& image) {

	UINT pycount = 0;
	// Read the image....
	sspace = NoPyramid;
	lbpmaps.resize(1);
	gradientmaps.resize(1);
	ComputeLBPMap(image, lbpmaps[0]);
	if (usegradmag)
		ComputeGradientMap(image, gradientmaps[0]);
}
void LBPFeatures::InitalizeMaps(Image& image, PyramidType sspace_) {
	UINT pycount = 0;
	// Read the image....
	sspace = sspace_;
	if (sspace_ == SingleRes) {
		pycount = pyobj.GeneratePyramid(image);
		lbpmaps.resize(pycount);
		gradientmaps.resize(pycount);
		for (UINT i = 0; i < pycount; i++) {
			Image& imgref = pyobj.GetImageRef(i);
			// Compute the Gradient Image
			ComputeLBPMap(imgref, lbpmaps[i]);
			if (usegradmag)
				ComputeGradientMap(imgref, gradientmaps[i]);
		}
	} else if (sspace == NoPyramid) {
		InitalizeMaps(image);
	}
}
#ifndef OLDCODE
void LBPFeatures::ComputeHistogram(int initx, int inity, LBPMap &lbpmap,
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
void LBPFeatures::ComputeHistBilinear(int initx, int inity, LBPMap &lbpmap,
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
void LBPFeatures::ComputeHistBilinear(int initx, int inity, LBPMap &lbpmap,
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
void LBPFeatures::GetFeatures(UINT index, int x, int y, vector<REAL>&feat) {
	// offset in ComputeHist
	UINT featsize = nxcells * nycells * nbins;
	fill(feat.begin(), feat.end(), 0);// Initialize with 0's;
	/*RGB channels based code*/
	if (!usegradmag) {
		ComputeHistogram(x, y, lbpmaps[index], &feat[0]);
	} else {
		ComputeHistBilinear(x, y, lbpmaps[index], gradientmaps[index], &feat[0]);
	}
	NormalizeFeatures(&feat[0]);
}
void LBPFeatures::NormalizeFeatures(REAL *feat) {
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
			REAL *dst = histptr + nbins * (x + y * nxcells); // write
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
void LBPFeatures::GetFeatures(UINT index, int x, int y, vector<REAL>&feat) {

	UINT offset = 0;
	x = MAX(0,x - radius); // to cater for offset
	y = MAX(0,y - radius);
	UINT twidth = x + width < lbpmaps[index].nxpoints ? x + width
	: lbpmaps[index].nxpoints, theight = y + height
	< lbpmaps[index].nypoints ? y + height : lbpmaps[index].nypoints;
	fill(feat.begin(), feat.end(), 0);// Initialize with 0's;
	for (UINT i = y; i < theight; i += lbpstride)
	for (UINT j = x; j < twidth; j += lbpstride) {
		lbpmaps[index].ComputeHist(j, i, cellsize, cellsize, &feat[offset]);
		NormalizeFeatures(&feat[offset]);
		offset += nbins;
	}
	if (nflevels > 1) // multi-level features
	{
		for (UINT i = y; i < theight; i += 2 * lbpstride)
		for (UINT j = x; j < twidth; j += 2 * lbpstride) {
			lbpmaps[index].ComputeHist(j, i, 2 * cellsize, 2 * cellsize,
					&feat[offset]);
			NormalizeFeatures(&feat[offset]);
			offset += nbins;
		}
	}
}
void LBPFeatures::NormalizeFeatures(REAL *feat) {
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

void LBPFeatures::GetFoldedFeatures(UINT index, int x, int y, vector<REAL>&feat) {

	UINT offset = 0;
	x = MAX(0,x - radius); // to cater for offset
	y = MAX(0,y - radius);
	UINT twidth = x + width < lbpmaps[index].nxpoints ? x + width
	: lbpmaps[index].nxpoints, theight = y + height
	< lbpmaps[index].nypoints ? y + height : lbpmaps[index].nypoints;
	fill(feat.begin(), feat.end(), 0);// Initialize with 0's;
	twidth = ceil((REAL) twidth / 2.0); //half of the width
	int fx = x + width - lbpstride; // flipped x index
	vector<REAL> tfeat(nbins, 0);
	//	vector<UINT> tmap(cellsize*cellsize,0);
	for (UINT i = y; i < theight; i += lbpstride) {
		fx = x + width - lbpstride;
		for (UINT j = x; j < twidth; j += lbpstride) {
			lbpmaps[index].ComputeHist(j, i, cellsize, cellsize, &feat[offset]);
			NormalizeFeatures(&feat[offset]);

			flbpmaps[index].ComputeHist(fx, i, cellsize, cellsize, &tfeat[0]);
			NormalizeFeatures(&tfeat[0]);

			for (UINT iter = 0; iter < nbins; ++iter)
			feat[offset + iter] += tfeat[iter];
			fill(tfeat.begin(), tfeat.end(), 0);// Initialize with 0's;

			offset += nbins;
			fx -= lbpstride;
		}
	}
}

// temporary
//			lbpmaps[index].GetMap(j,i,cellsize,cellsize,&tmap[0]);
//			for(UINT iter=0; iter < tmap.size();++iter)
//				cout<< " "<<tmap[iter];
//			flbpmaps[index].GetMap(fx,i,cellsize,cellsize,&tmap[0]);
//			for(UINT iter=0; iter < tmap.size();++iter)
//				cout<< " "<<tmap[iter];
//			cout<<endl;

//void LBPFeatures::GetFoldedFeatures(UINT index,int x, int y ,vector<REAL>&feat)
//{
//
//	vector<REAL> tfeat(GetDim(width,height),0);
//	GetFeatures(index,x,y,tfeat);
//	UINT twidth = width/2;
//	int fx;
//	fill(feat.begin(),feat.end(),0);
//	for(UINT i=0 ; i < height; i+=lbpstride)
//	{
//		bx =  (width - lbpstride)/lbpstride * nbins;
//		for(UINT j=0 ; j < twidth; j+=lbpstride)
//		{
//			for(UINT iter=0; iter < nbins;++iter)
//				feat[offset+iter] += tfeat[offset+bx];
//
//			fill(tfeat.begin(),tfeat.end(),0);// Initialize with 0's;
//			offset += nbins;
//			bx -= nbins;
//		}
//	}
//}


/*void LBPFeatures::UnFoldWeights(REAL *input,vector<REAL>&weights)
 {
 UINT twidth = width/2,offset=0,count=0,
 foffset=0,boffset=0,inoffset=0;// forward & backward offsets

 for(UINT i=0 ; i < height; i+=lbpstride,++count)
 {
 boffset = ((width - lbpstride)/lbpstride)*nbins;
 foffset=0;
 offset = count * (width/lbpstride) * nbins;
 for(UINT j=0 ; j < twidth; j+=lbpstride)
 {
 for(UINT iter=0; iter < nbins;++iter)
 {
 weights[iter+offset+foffset] = input[inoffset + iter];
 weights[iter+offset+boffset] = input[inoffset + iter];
 }

 foffset += nbins;
 boffset -= nbins;
 inoffset += nbins;
 }
 }
 }

 void LBPFeatures::UnFoldWeights(REAL *input, vector<REAL>&weights) {
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
 void LBPFeatures::UnFoldWeights(REAL *input, REAL *weights) {
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

 */
#endif
