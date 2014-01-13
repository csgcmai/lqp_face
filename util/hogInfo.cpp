/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/
// Data to be shared among class objects
#include "hogInfo.h"
//
//UINT HInfo::norient =9;
HInfo::HInfo(UINT xbmax_, UINT ybmax_, REAL scale_, UINT norient_) :
	xbmax(xbmax_), ybmax(ybmax_), scale(scale_), norient(norient_) {
	if (xbmax != 0 && ybmax != 0)
		hogvalue.resize(xbmax * ybmax * norient, 0);
	permvec.resize(norient, 0);
	GetPermutationFeatures(&permvec[0], norient);
}
HInfo::HInfo(UINT xbmax_, UINT ybmax_, UINT norient_, REAL *weights) :
	xbmax(xbmax_), ybmax(ybmax_), scale(1), norient(norient_) {
	if (xbmax != 0 && ybmax != 0) {
		hogvalue.resize(xbmax * ybmax * norient, 0);
		for (UINT i = 0; i < hogvalue.size(); ++i)
			hogvalue[i] = weights[i];
	} else
		cout << "\n Error in Initialization of the HInfo "
				<< "Either xbmax ==0 or ybmax==0" << xbmax << "  " << ybmax
				<< endl;
	permvec.resize(norient, 0);
	GetPermutationFeatures(&permvec[0], norient);
}
HInfo::HInfo(UINT xbmax_, UINT ybmax_, REAL scale_, UINT norient_,
		REAL *weights) :
	xbmax(xbmax_), ybmax(ybmax_), scale(scale_), norient(norient_) {
	if (xbmax != 0 && ybmax != 0) {
		hogvalue.resize(xbmax * ybmax * norient, 0);
		for (UINT i = 0; i < hogvalue.size(); ++i)
			hogvalue[i] = weights[i];
	} else
		cout << "\n Error in Initialization of the HInfo "
				<< "Either xbmax ==0 or ybmax==0" << xbmax << "  " << ybmax
				<< endl;
	permvec.resize(norient, 0);
	GetPermutationFeatures(&permvec[0], norient);
}
void HInfo::Initialize(UINT w, UINT h, REAL scale_, UINT norient_) {
	xbmax = w;
	ybmax = h;
	scale = scale_;
	norient = norient_;
	// total_blocks * number of cells in one block * number of orientations...
	hogvalue.resize(xbmax * ybmax * norient, 0);
	permvec.resize(norient, 0);
	GetPermutationFeatures(&permvec[0], norient);
}
void HInfo::Initialize(UINT w, UINT h, REAL scale_, UINT norient_,
		const vector<REAL>&feat) {
	xbmax = w; // number of xcells;
	ybmax = h; // number of ycell;
	scale = scale_;
	norient = norient_;
	hogvalue = feat;
	permvec.resize(norient, 0);
	GetPermutationFeatures(&permvec[0], norient);
}
void HInfo::Initialize(UINT w, UINT h, REAL scale_, UINT norient_, REAL*feat) {
	xbmax = w;
	ybmax = h;
	scale = scale_;
	norient = norient_;
	hogvalue.resize(xbmax * ybmax * norient, 0);
	copy(feat, feat + hogvalue.size(), hogvalue.begin());
	permvec.resize(norient, 0);
	GetPermutationFeatures(&permvec[0], norient);
}
HInfo::HInfo(const HInfo & obj) {
	*this = obj;
}

void HInfo::operator =(const HInfo &obj) {
	xbmax = obj.xbmax;
	ybmax = obj.ybmax;
	hogvalue = obj.hogvalue;
	scale = obj.scale;
	norient = obj.norient;
	permvec = obj.permvec;
}
void HInfo::GetPermutationFeatures(int *pfeat, int norient) {

	if (norient == 9) {
		int p[] = { 0, 8, 7, 6, 5, 4, 3, 2, 1 };
		copy(p, p + norient, pfeat);

	} else if (norient == 36) {
		int p[] = { 0, 8, 7, 6, 5, 4, 3, 2, 1, 9, 17, 16, 15, 14, 13, 12, 11,
				10, 18, 26, 25, 24, 23, 22, 21, 20, 19, 27, 35, 34, 33, 32, 31,
				30, 29, 28 };
		copy(p, p + norient, pfeat);

	} else if (norient == 31) {
		int p[] = { 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 17, 16, 15, 14, 13, 12, 11,
				10, 18, 26, 25, 24, 23, 22, 21, 20, 19, 29, 30, 27, 28 };
		copy(p, p + norient, pfeat);

	}/*FOR CSLBP's*/else if (norient == 16) {
		//		int p[] = { 4, 3, 2, 1, 0, 7, 6, 5, 12, 11, 10, 9, 8, 15, 14, 13 };/*LBPHIST*/
		int p[] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };
		copy(p, p + norient, pfeat);

	} else if (norient == 24) {
		int p[] = { 4, 3, 2, 1, 0, 7, 6, 5, 12, 11, 10, 9, 8, 15, 14, 13, 20,
				19, 18, 17, 16, 23, 22, 21 };/*LBPHIST*/
		copy(p, p + norient, pfeat);

	} else if (norient % 256 == 0) {
		int p[] = { 0, 16, 8, 24, 4, 20, 12, 28, 2, 18, 10, 26, 6, 22, 14, 30,
				1, 17, 9, 25, 5, 21, 13, 29, 3, 19, 11, 27, 7, 23, 15, 31, 128,
				144, 136, 152, 132, 148, 140, 156, 130, 146, 138, 154, 134,
				150, 142, 158, 129, 145, 137, 153, 133, 149, 141, 157, 131,
				147, 139, 155, 135, 151, 143, 159, 64, 80, 72, 88, 68, 84, 76,
				92, 66, 82, 74, 90, 70, 86, 78, 94, 65, 81, 73, 89, 69, 85, 77,
				93, 67, 83, 75, 91, 71, 87, 79, 95, 192, 208, 200, 216, 196,
				212, 204, 220, 194, 210, 202, 218, 198, 214, 206, 222, 193,
				209, 201, 217, 197, 213, 205, 221, 195, 211, 203, 219, 199,
				215, 207, 223, 32, 48, 40, 56, 36, 52, 44, 60, 34, 50, 42, 58,
				38, 54, 46, 62, 33, 49, 41, 57, 37, 53, 45, 61, 35, 51, 43, 59,
				39, 55, 47, 63, 160, 176, 168, 184, 164, 180, 172, 188, 162,
				178, 170, 186, 166, 182, 174, 190, 161, 177, 169, 185, 165,
				181, 173, 189, 163, 179, 171, 187, 167, 183, 175, 191, 96, 112,
				104, 120, 100, 116, 108, 124, 98, 114, 106, 122, 102, 118, 110,
				126, 97, 113, 105, 121, 101, 117, 109, 125, 99, 115, 107, 123,
				103, 119, 111, 127, 224, 240, 232, 248, 228, 244, 236, 252,
				226, 242, 234, 250, 230, 246, 238, 254, 225, 241, 233, 249,
				229, 245, 237, 253, 227, 243, 235, 251, 231, 247, 239, 255 };
		//		copy(p, p + norient, pfeat);
		int nbins = 256, nlbps = norient / nbins;
		for (int i = 0; i < nlbps; ++i) {
			int offset = i * nbins, *tptr = pfeat + offset;
			for (int j = 0; j < nbins; ++j)
				tptr[j] = p[j] + offset;
		}

	} else if (norient % 59 == 0) { // for all different types of LBP's and LTP's...
		int nlbps = norient / 59;
		int p[] = { 0, 11, 7, 12, 4, 8, 13, 2, 5, 9, 14, 1, 3, 6, 10, 15, 29,
				30, 31, 32, 33, 34, 22, 36, 37, 38, 39, 40, 41, 16, 17, 18, 19,
				20, 21, 35, 23, 24, 25, 26, 27, 28, 42, 47, 51, 54, 56, 43, 48,
				52, 55, 44, 49, 53, 45, 50, 46, 57, 58, /*LBP*/};
		//		copy(p, p + norient, pfeat);
		for (int i = 0; i < nlbps; ++i) {
			int offset = i * 59, *tptr = pfeat + offset;
			for (int j = 0; j < 59; ++j)
				tptr[j] = p[j] + offset;
		}
	}
}
void HInfo::Get2XHOGCells(UINT xstart, UINT ystart, UINT ncols, UINT nrows,
		REAL *outhog) const
// Performs the binlinear interpolation on the given  HInfo image...
{
	UINT x0, y0, x1, y1;
	//	outhog.resize(nrows*ncols*norient,0);//
	REAL U, V, u, v, xfactor = (REAL) (ncols / 2 - 1) / (ncols - 1), yfactor =
			(REAL) (nrows / 2 - 1) / (nrows - 1);

	nrows += ystart;
	ncols += xstart;
	for (UINT k = 0; k < norient; ++k) {
		for (UINT i = ystart; i < nrows; ++i) {
			u = i * yfactor;
			U = u - floor(u);
			//			++u ;
			y0 = (UINT) floor(u);
			y1 = (UINT) ceil(u);
			// so copy the original pixels at boundary...
			for (UINT j = xstart; j < ncols; ++j) {
				v = j * xfactor;
				V = v - floor(v);
				//				++v;
				x0 = (UINT) floor(v);
				x1 = (UINT) ceil(v);
				outhog[norient * (j - xstart) + k + (i - ystart) * ncols
						* norient] = (V - 1) * ((U - 1) * hogvalue[norient * x0
						+ k + y0 * xbmax * norient] - U * hogvalue[norient * x0
						+ k + y1 * xbmax * norient]) - V * ((U - 1)
						* hogvalue[norient * x1 + k + y0 * xbmax * norient] - U
						* hogvalue[norient * x1 + k + y1 * xbmax * norient]);
			}
		}
	}
}
void HInfo::FlipFeatures() {
	vector<REAL> ffeat(hogvalue.size(), 0);
	GetFlippedFeatures(0, 0, xbmax, ybmax, &ffeat[0]);
	hogvalue = ffeat;
}
void HInfo::GetFlippedFeatures(HInfo& hinfo) {
	hinfo.Initialize(xbmax, ybmax, scale, norient);
	GetFlippedFeatures(0, 0, xbmax, ybmax, &(hinfo.hogvalue[0]));
}
void HInfo::FlipFeatures(UINT xbmax_, UINT ybmax_, UINT norient_,
		REAL* ihogvalue, REAL*hvalue) {
	vector<int> permvec(norient_, 0);
	GetPermutationFeatures(&permvec[0], norient_);
	UINT loffset, roffset;
	for (UINT i = 0; i < ybmax_; ++i)
		for (int j = xbmax_ - 1, l = 0; j >= 0; --j, ++l) {
			loffset = (l + i * xbmax_) * norient_;
			roffset = (j + i * xbmax_) * norient_;
			for (UINT k = 0; k < norient_; ++k)
				hvalue[loffset + k] = ihogvalue[roffset + permvec[k]];
		}
}

void HInfo::GetFlippedFeatures(UINT xstart, UINT ystart, UINT xsize,
		UINT ysize, REAL *hvalue) const {
	UINT roffset, counter = 0;
	for (UINT i = ystart; i < ystart + ysize; ++i)
		for (int j = xstart + xsize - 1; j >= (int) xstart; --j) {
			roffset = (j + i * xbmax) * norient;
			for (UINT k = 0; k < norient; ++k)
				hvalue[counter++] = hogvalue[roffset + permvec[k]];
		}
}

void HInfo::GetUnFoldedHOGCells(UINT xsize, REAL *outhog) const
// All the inputs are in cell coordinates...
{
	UINT width1 = (UINT) ceil((double) xsize / 2.0), width2 = (UINT) floor(
			(double) xsize / 2.0), counter = 0, roffset, coffset;

	vector<REAL> feat2(width2 * ybmax * norient);
	GetFlippedFeatures(0, 0, width2, ybmax, &feat2[0]);

	ofstream ofile("unfold.txt");
	// Now copy the original and un flipped weights...
	for (UINT i = 0; i < ybmax; ++i) {

		roffset = i * width1 * norient;
		for (UINT j = 0; j < width1; ++j) {
			coffset = roffset + j * norient;
			for (UINT k = 0; k < norient; ++k)
				outhog[counter++] = hogvalue[coffset + k];
		}
		roffset = i * width2 * norient;
		for (UINT j = 0; j < width2; ++j) {
			coffset = roffset + j * norient;
			for (UINT k = 0; k < norient; ++k)
				outhog[counter++] = feat2[coffset + k];
		}
	}

	for (UINT i = 0; i < feat2.size(); ++i)
		ofile << feat2[i] << " ";
	ofile.close();
	//
	//
	//	ofile.open("temp_f.txt");
	//	for(int i=0; i < feat2.size();++i)
	//		ofile<<feat2[i]<<" ";
	//	ofile.close();

	/*

	 ofile.open("temp_o.txt");
	 feat2.resize(ysize*xsize*norient,0);
	 GetHOGCells(xstart,ystart,xsize,ysize,&feat2[0]);
	 for(int i=0; i < ysize*xsize*norient;++i)
	 ofile<<feat2[i]<<" ";
	 ofile.close();*/

}

void HInfo::GetFoldedHOGCells(UINT xstart, UINT ystart, UINT xsize, UINT ysize,
		REAL *outhog) const
// All the inputs are in cell coordinates...
{
	UINT width1 = ceil((double) xsize / 2.0), width2 = floor((double) xsize
			/ 2.0), roffset1, coffset1, roffset2, coffset2;
	if (xstart + width1 > xbmax) {
		cout << "\n Invalid Size for Folded HoG cells" << endl;
		return;
	}

	vector<REAL> feat2(width2 * ysize * norient);
	GetFlippedFeatures(xstart + width1, ystart, width2, ysize, &feat2[0]);
	GetHOGCells(xstart, ystart, width1, ysize, outhog);
	//	ofstream ofile("temp.txt");
	for (UINT i = 0; i < ysize; ++i) {
		roffset1 = i * width1 * norient;
		roffset2 = i * width2 * norient;
		for (UINT j = 0; j < width2; ++j) {
			coffset1 = roffset1 + j * norient;
			coffset2 = roffset2 + j * norient;
			for (UINT k = 0; k < norient; ++k)
				outhog[coffset1 + k] = outhog[coffset1 + k] + feat2[coffset2
						+ k];
		}
	}

	/*for(int i=0; i < ysize*width1*norient;++i)
	 ofile<<outhog[i]<<" ";
	 ofile.close();*/

	/*
	 ofile.open("temp_f.txt");
	 for(int i=0; i < feat2.size();++i)
	 ofile<<feat2[i]<<" ";
	 ofile.close();



	 ofile.open("temp_o.txt");
	 feat2.resize(ysize*xsize*norient,0);
	 GetHOGCells(xstart,ystart,xsize,ysize,&feat2[0]);
	 for(int i=0; i < ysize*xsize*norient;++i)
	 ofile<<feat2[i]<<" ";
	 ofile.close();*/

}
void HInfo::GetSkippedFlippedFeatures(UINT xstart, UINT ystart, UINT xsize,
		UINT ysize, UINT skip, REAL *hvalue) const {

	UINT roffset, counter = 0;
	for (UINT i = ystart; i < ystart + ysize; i += skip)
		for (int j = xstart + xsize - skip; j >= (int) xstart; j -= skip) {
			roffset = (j + i * xbmax) * norient;
			for (UINT k = 0; k < norient; ++k)
				hvalue[counter++] = hogvalue[roffset + permvec[k]];
		}

}

void HInfo::GetSkippedFoldedHOGCells(UINT xstart, UINT ystart, UINT xsize,
		UINT ysize, UINT skip, REAL *outhog) const {
	// All the inputs are in cell coordinates...
	UINT width1 = ceil((double) xsize / 2.0), width2 = floor((double) xsize
			/ 2.0), roffset1, coffset1, roffset2, coffset2;
	if (xstart + width1 > xbmax) {
		cout << "\n Invalid Size for Folded HoG cells" << endl;
		return;
	}

	vector<REAL> feat2(width2 * ysize * norient, 0);
	GetSkippedFlippedFeatures(xstart + width1, ystart, width2, ysize, skip,
			&feat2[0]);
	GetSkippedHOGCells(xstart, ystart, width1, ysize, skip, outhog);
	//	ofstream ofile("temp.txt");
	ysize /= skip;
	width2 /= skip;
	width1 /= skip;
	for (UINT i = 0; i < ysize; ++i) {
		roffset1 = i * width1 * norient;
		roffset2 = i * width2 * norient;
		for (UINT j = 0; j < width2; ++j) {
			coffset1 = roffset1 + j * norient;
			coffset2 = roffset2 + j * norient;
			for (UINT k = 0; k < norient; ++k)
				outhog[coffset1 + k] = outhog[coffset1 + k] + feat2[coffset2
						+ k];
		}
	}
}

void HInfo::GetHOGCells(UINT xstart, UINT ystart, UINT xsize, UINT ysize,
		REAL *outhog) const
// All the inputs are in cell coordinates...
{
	//
	//	UINT oxsize = xsize,oysize=ysize;// original xsize
	//	xsize = xstart+xsize > xbmax ? xsize - ((xstart+xsize) - xbmax) -1 : xsize;
	//	ysize = ystart+ysize > ybmax ? ysize - ((ystart+ysize) - ybmax) -1 : ysize;
	//
	//	// Pad by Constant = Zeros for the extended boundaries....
	//	if(oxsize!=xsize || oysize != ysize)
	//	{
	//		for(UINT i=0; i < oxsize*oysize*norient; ++i)
	//			outhog[i]=0;
	//	}

	UINT counter = 0, offset = xstart * norient + ystart * xbmax * norient,
			coffset, roffset;
	for (UINT i = 0; i < ysize; ++i) {
		roffset = offset + i * xbmax * norient;
		for (UINT j = 0; j < xsize; ++j) {
			coffset = roffset + j * norient;
			for (UINT k = 0; k < norient; ++k)
				outhog[counter++] = hogvalue[coffset + k];
		}
	}
}
void HInfo::GetSkippedHOGCells(UINT xstart, UINT ystart, UINT xsize,
		UINT ysize, UINT skip, REAL *outhog) const
// All the inputs are in cell coordinates...
// Skip: contains number of cells to skip after reading a cell..
{

	//	UINT oxsize = xsize,oysize=ysize;// original xsize
	//	xsize = xstart+xsize > xbmax ? xsize - ((xstart+xsize) - xbmax) : xsize;
	//	ysize = ystart+ysize > ybmax ? ysize - ((ystart+ysize) - ybmax) : ysize;
	//
	//	// Pad by Constant = Zeros for the extended boundaries....
	//	if(oxsize!=xsize || oysize != ysize)
	//	{
	//		for(UINT i=0; i < oxsize/skip*oysize/skip*norient; ++i)
	//			outhog[i]=0;
	//	}

	UINT counter = 0, offset = xstart * norient + ystart * xbmax * norient,
			coffset, roffset;
	for (UINT i = 0; i < ysize; i += skip) {
		roffset = offset + i * xbmax * norient;
		for (UINT j = 0; j < xsize; j += skip) {
			coffset = roffset + j * norient;
			for (UINT k = 0; k < norient; ++k)
				outhog[counter++] = hogvalue[coffset + k];
		}
	}
}

void HInfo::SetHOGCells(UINT xstart, UINT ystart, UINT xsize, UINT ysize,
		REAL ivalue)
// All the inputs are in cell coordinates...
{

	xsize = xstart + xsize > xbmax ? xsize - ((xstart + xsize) - xbmax) : xsize;
	ysize = ystart + ysize > ybmax ? ysize - ((ystart + ysize) - ybmax) : ysize;

	UINT offset = xstart * norient + ystart * xbmax * norient, coffset, roffset;
	for (UINT i = 0; i < ysize; ++i) {
		roffset = offset + i * xbmax * norient;
		for (UINT j = 0; j < xsize; ++j) {
			coffset = roffset + j * norient;
			for (UINT k = 0; k < norient; ++k)
				hogvalue[coffset + k] = ivalue;
		}
	}
}
void HInfo::TrimWeights(const UINT &offset)
// Trim the weights by considering only the center information of the feature vector...
{
	if (offset <= 0 || offset >= xbmax / 2)
		return;
	vector<REAL> tvec((xbmax - 2 * offset) * (ybmax - 2 * offset) * norient, 0);
	GetHOGCells(offset, offset, xbmax - 2 * offset, ybmax - 2 * offset,
			&tvec[0]);
	xbmax -= 2 * offset;
	ybmax -= 2 * offset;
	hogvalue = tvec;
}
REAL HInfo::GetSumHOGCells(UINT xstart, UINT ystart, UINT xsize, UINT ysize,
		bool pweights) const
// All the inputs are in cell coordinates...
{

	xsize = xstart + xsize > xbmax ? xsize - ((xstart + xsize) - xbmax) : xsize;
	ysize = ystart + ysize > ybmax ? ysize - ((ystart + ysize) - ybmax) : ysize;

	REAL sum = 0;
	UINT offset = xstart * norient + ystart * xbmax * norient, coffset, roffset;
	for (UINT i = 0; i < ysize; ++i) {
		roffset = offset + i * xbmax * norient;
		for (UINT j = 0; j < xsize; ++j) {
			coffset = roffset + j * norient;
			for (UINT k = 0; k < norient; ++k)
				if (pweights)
					sum += (hogvalue[coffset + k] >= 0 ? hogvalue[coffset + k]
							: 0);
				else
					sum += hogvalue[coffset + k];
		}
	}
	return sum;
}

REAL HInfo::DotProduct(UINT xstart, UINT ystart, UINT xsize, UINT ysize,
		const vector<REAL>&inhog) const
// All the inputs are cell coordinates...
{
	xsize = xstart + xsize > xbmax ? xsize - ((xstart + xsize) - xbmax) : xsize;
	ysize = ystart + ysize > ybmax ? ysize - ((ystart + ysize) - ybmax) : ysize;
	REAL sum = 0;
	UINT offset = xstart * norient + ystart * xbmax * norient, coffset,
			roffset, counter = 0;
	for (UINT i = 0; i < ysize; ++i) {
		roffset = offset + i * xbmax * norient;
		for (UINT j = 0; j < xsize; ++j) {
			coffset = roffset + j * norient;
			for (UINT k = 0; k < norient; ++k)
				sum += hogvalue[coffset + k] * inhog[counter++];
		}
	}
	return sum;
}
void HInfo::DotProduct(UINT fxsize, UINT fysize, const vector<REAL>&filter,
		vector<REAL>&outhog) const
//Computes the dot product over complete HOG-Image in Column-Major order
{
	int nxwindows = xbmax - fxsize + 1, nywindows = ybmax - fysize + 1, count =
			0;
	outhog.resize(nxwindows * nywindows, 0);
	for (int ystart = 0; ystart < nywindows; ++ystart)
		for (int xstart = 0; xstart < nxwindows; ++xstart) {
			REAL sum = 0;
			UINT offset = xstart * norient + ystart * xbmax * norient, coffset,
					roffset, counter = 0;
			for (UINT i = 0; i < fysize; ++i) {
				roffset = offset + i * xbmax * norient;
				for (UINT j = 0; j < fxsize; ++j) {
					coffset = roffset + j * norient;
					for (UINT k = 0; k < norient; ++k)
						sum += hogvalue[coffset + k] * filter[counter++];
				}
			}
			outhog[count++] = sum;
		}

}
void HInfo::DotProduct(UINT fxsize, UINT fysize, const vector<REAL>&filter,
		Store&outhog) const
//Computes the dot product over complete HOG-Image in Column-Major order
{
	DotProduct(fxsize, fysize, &filter[0], outhog);
}
//TODO For LBP_LTP Features to have the Dot Product with Quadratic coefficient..
void HInfo::DotProductWithQuadratic(UINT fxsize, UINT fysize,
		const REAL &meanssd, const REAL *mpvector, const REAL *filter,
		Store&outhog) const
//Computes the dot product over complete HOG-Image in Column-Major order
{
	int count = 0;
	for (int ystart = 0; ystart < outhog.height; ++ystart)
		for (int xstart = 0; xstart < outhog.width; ++xstart) {
			REAL sum = 0, ssd = 0;
			UINT offset = xstart * norient + ystart * xbmax * norient, coffset,
					roffset, counter = 0;
			for (UINT i = 0; i < fysize; ++i) {
				roffset = offset + i * xbmax * norient;
				for (UINT j = 0; j < fxsize; ++j) {
					coffset = roffset + j * norient;
					for (UINT k = 0; k < norient; ++k, ++counter) {
						sum += hogvalue[coffset + k] * filter[counter];
						ssd += ((hogvalue[coffset + k] - mpvector[counter])
								* (hogvalue[coffset + k] - mpvector[counter]));
					}
				}
			}
			sum += sqrt(ssd) / meanssd * filter[counter];
			outhog.SetValue(xstart, ystart, sum);
		}
}
void HInfo::DotProductWithProjection(UINT fxsize, UINT fysize,
		const VectorXf &sigma, const MatrixXf &projMat, const REAL *filter,
		Store&outhog) const
//Computes the dot product over complete HOG-Image in Column-Major order
{
	int count = 0;
	VectorXf feat(projMat.rows()), projfeat;
	for (int ystart = 0; ystart < outhog.height; ++ystart)
		for (int xstart = 0; xstart < outhog.width; ++xstart) {
			REAL sum = 0, ssd = 0;
			UINT offset = xstart * norient + ystart * xbmax * norient, coffset,
					roffset, counter = 0;
			for (UINT i = 0; i < fysize; ++i) {
				roffset = offset + i * xbmax * norient;
				for (UINT j = 0; j < fxsize; ++j) {
					coffset = roffset + j * norient;
					for (UINT k = 0; k < norient; ++k, ++counter) {
						sum += hogvalue[coffset + k] * filter[counter];
						feat[counter] = hogvalue[coffset + k];
					}
				}
			}
			projfeat = feat.transpose() * projMat;
//			projfeat = projfeat.cwise() / sigma;
			projfeat = projfeat.array() / sigma.array();

			for (UINT k = 0; k < sigma.size(); ++k) {
				for (UINT tvar = k; tvar < sigma.size(); ++tvar)
					sum += projfeat[tvar] * projfeat[k] * filter[counter++];
			}

			outhog.SetValue(xstart, ystart, sum);
		}

}
void HInfo::DotProduct(UINT fxsize, UINT fysize, const REAL *filter,
		Store&outhog) const
//Computes the dot product over complete HOG-Image in Column-Major order
{
	int count = 0;
	for (int ystart = 0; ystart < outhog.height; ++ystart)
		for (int xstart = 0; xstart < outhog.width; ++xstart) {
			REAL sum = 0;
			UINT offset = xstart * norient + ystart * xbmax * norient, coffset,
					roffset, counter = 0;
			for (UINT i = 0; i < fysize; ++i) {
				roffset = offset + i * xbmax * norient;
				for (UINT j = 0; j < fxsize; ++j) {
					coffset = roffset + j * norient;
					for (UINT k = 0; k < norient; ++k)
						sum += hogvalue[coffset + k] * filter[counter++];
				}
			}
			outhog.SetValue(xstart, ystart, sum);
		}

}

void HInfo::SkippedDotProduct(UINT fxsize, UINT fysize, UINT skip,
		const REAL *filter, Store&outhog) const
//Computes the dot product over complete HOG-Image in Column-Major order
{
	int count = 0;
	for (int ystart = 0; ystart < outhog.height; ++ystart)
		for (int xstart = 0; xstart < outhog.width; ++xstart) {
			REAL sum = 0;
			UINT offset = xstart * norient + ystart * xbmax * norient, coffset,
					roffset, counter = 0;
			for (UINT i = 0; i < fysize; i += skip) {
				roffset = offset + i * xbmax * norient;
				for (UINT j = 0; j < fxsize; j += skip) {
					coffset = roffset + j * norient;
					for (UINT k = 0; k < norient; ++k)
						sum += hogvalue[coffset + k] * filter[counter++];
				}
			}
			outhog.SetValue(xstart, ystart, sum);
		}

}

REAL HInfo::DotProduct(UINT xstart, UINT ystart, UINT xsize, UINT ysize,
		REAL *inhog) const
// All the inputs are cell coordinates...
{

	xsize = xstart + xsize > xbmax ? xsize - ((xstart + xsize) - xbmax) : xsize;
	ysize = ystart + ysize > ybmax ? ysize - ((ystart + ysize) - ybmax) : ysize;
	REAL sum = 0;
	UINT coffset, roffset, counter = 0;
	const REAL *thogvalue = &hogvalue[xstart * norient + ystart * xbmax
			* norient];
	for (UINT i = 0; i < ysize; ++i) {
		roffset = i * xbmax * norient;
		for (UINT j = 0; j < xsize; ++j) {
			coffset = roffset + j * norient;
			for (UINT k = 0; k < norient; ++k)
				sum += thogvalue[coffset + k] * inhog[counter++];
		}
	}
	return sum;
}

void HInfo::PadFeatureMap(UINT padx, UINT pady) {
	UINT nxbmax = xbmax + 2 * padx, nybmax = ybmax + 2 * pady;
	vector<REAL> feat(nxbmax * nybmax * norient, 0);
	UINT counter = 0, offset = (padx + pady * nxbmax) * norient, coffset,
			roffset;
	for (UINT i = 0; i < ybmax; ++i) {
		roffset = offset + i * nxbmax * norient;
		for (UINT j = 0; j < xbmax; ++j) {
			coffset = roffset + j * norient;
			for (UINT k = 0; k < norient; ++k)
				feat[coffset + k] = hogvalue[counter++];
		}
	}
	xbmax = nxbmax;
	ybmax = nybmax;
	hogvalue = feat;
}
void HInfo::PruneNegativeWeights() {// to be used for the selection of parts
	// Square the value and theshold it...
	for (UINT i = 0; i < hogvalue.size(); ++i)
		hogvalue[i] = hogvalue[i] > 0 ? hogvalue[i] * hogvalue[i] : 0;

	vector<REAL> feat2(xbmax * ybmax * norient);
	GetFlippedFeatures(0, 0, xbmax, ybmax, &feat2[0]);

	for (UINT i = 0; i < hogvalue.size(); ++i)
		hogvalue[i] = hogvalue[i] + feat2[i];
}
void HInfo::Write(fstream &ofile) {
	ofile.write((char*) &scale, sizeof(REAL));
	ofile.write((char*) &xbmax, sizeof(UINT));
	ofile.write((char*) &ybmax, sizeof(UINT));
	UINT size = hogvalue.size();
	ofile.write((char*) &size, sizeof(UINT));
	ofile.write((char*) &hogvalue[0], sizeof(REAL) * size);
}
void HInfo::WriteText(fstream &ofile) {
	ofile << " Scale = " << scale << " Xbmax = " << xbmax << " Ybmax = "
			<< ybmax << " Dim of Hogfeature = " << hogvalue.size() << endl;

	for (UINT i = 0; i < hogvalue.size(); ++i)
		ofile << hogvalue[i] << " ";
	ofile << endl;
}
void HInfo::Read(fstream &ifile) {
	ifile.read((char*) &scale, sizeof(REAL));
	ifile.read((char*) &xbmax, sizeof(UINT));
	ifile.read((char*) &ybmax, sizeof(UINT));
	UINT size;
	ifile.read((char*) &size, sizeof(UINT));
	hogvalue.resize(size, 0);
	ifile.read((char*) &hogvalue[0], sizeof(REAL) * size);
}

void HInfo::BilinearInterpolation(REAL sfactor, HInfo& output)
// Performs the binlinear interpolation on the given  HInfo image...
{
	UINT nrows = (UINT) (ybmax * sfactor), ncols = (UINT) (xbmax * sfactor),
			x0, y0, x1, y1;
	output.hogvalue.resize(nrows * ncols * norient, 0);//
	REAL U, V, u, v, xfactor = (REAL) (xbmax - 1) / (ncols - 1), yfactor =
			(REAL) (ybmax - 1) / (nrows - 1);

	for (UINT k = 0; k < norient; ++k) {
		for (UINT i = 0; i < nrows; ++i) {
			u = i * yfactor;
			U = u - floor(u);
			//			++u ;
			y0 = floor(u);
			y1 = ceil(u);
			// so copy the original pixels at boundary...
			for (UINT j = 0; j < ncols; ++j) {
				v = j * xfactor;
				V = v - floor(v);
				//				++v;
				x0 = floor(v);
				x1 = ceil(v);
				output.hogvalue[norient * j + k + i * ncols * norient]
						= (V - 1) * ((U - 1) * hogvalue[norient * x0 + k + y0
								* xbmax * norient] - U * hogvalue[norient * x0
								+ k + y1 * xbmax * norient]) - V * ((U - 1)
								* hogvalue[norient * x1 + k + y0 * xbmax
										* norient] - U * hogvalue[norient * x1
								+ k + y1 * xbmax * norient]);
			}
		}
	}
}
void HInfo::BlackOut(UINT x1, UINT y1, UINT x2, UINT y2) {
	x1 = MAX(floor((REAL)x1/scale),0);
	y1 = MAX(floor((REAL)y1/scale),0);
	if (x1 < xbmax && y1 < ybmax) {
		x2 = MIN(ceil((REAL)x2/scale),xbmax);
		y2 = MIN(ceil((REAL)y2/scale),ybmax);
		SetHOGCells(x1, y1, (x2 - x1) + 1, y2 - y1 + 1, 0);
	}
}

/*void HInfo::BilinearInterpolation(REAL sfactor,HInfo& output)
 // Performs the binlinear interpolation on the given  HInfo image...
 {
 UINT nrows =(UINT)( ybmax*sfactor),
 ncols = (UINT)(xbmax*sfactor),x0,y0,x1,y1;
 output.hogvalue.resize(nrows*ncols*norient,0);//
 REAL ii,ij,dx,dy;

 for(UINT k=0; k < norient;++k)
 {
 for(UINT i=0; i < nrows; ++i)
 {
 ii = i/sfactor;
 dy = ii - floor(ii);
 y0 =(UINT) floor(ii); y1 =(UINT) ceil(ii); // y0 = ii, y1 = y0+1;
 // so copy the original pixels at boundary...
 y1 = y1 > ybmax-1 ? ybmax-1 : y1;

 for(UINT j=0; j < ncols;++j)
 {
 ij = j/sfactor;
 dx = ij - floor(ij);
 x0 =(UINT) floor(ij); x1 = (UINT)ceil(ij);// x0 = ij, x1 = x0+1;
 x1 = x1 > xbmax-1 ? xbmax-1 : x1;
 output.hogvalue[norient*j+k+i*ncols*norient] = (1-dx)*(1-dy)*hogvalue[norient*x0+k+y0*xbmax*norient] +
 (1-dx)*dy*hogvalue[norient*x0+k+y1*xbmax*norient]+ (1-dy)*dx*hogvalue[norient*x1+k+y0*xbmax*norient]+
 dx*dy*hogvalue[norient*x1+k+y1*xbmax];
 }
 }
 }
 }

 */

//------------------------ HOGInfo Class definition-------------------------

HOGInfo::HOGInfo(UINT nimages, UINT nalevel_) :
	nalevel(nalevel_), counter(0) {
	imfeatures.resize(nimages);
}

HOGInfo::HOGInfo(const HOGInfo & obj) {
	*this = obj;
}

void HOGInfo::operator =(const HOGInfo &obj) {
	imfeatures = obj.imfeatures;
	nalevel = obj.nalevel;
	counter = obj.counter;
	norient = obj.norient;
}
void HOGInfo::SetInfo(UINT nalevel_, REAL sratio_) {
	nalevel = nalevel_;
}
void HOGInfo::SetFeature(UINT xmax, UINT ymax, REAL scale,
		const vector<REAL>&feat) {

	if (counter < imfeatures.size()) {
		imfeatures[counter++].Initialize(xmax, ymax, scale, norient, feat);
		//		imfeatures[counter++].SetFeature(feat);
	}
}
void HOGInfo::SetFeature(UINT xmax, UINT ymax, REAL scale, REAL* feat) {
	//	imfeatures.push_back(HInfo(xmax,ymax,scale));
	//	(imfeatures.back()).SetFeature(feat);
	if (counter < imfeatures.size()) {
		imfeatures[counter].Initialize(xmax, ymax, scale, norient, feat);
		//		imfeatures[counter++].SetFeature(feat);
	}
}
void HOGInfo::SetFeature(UINT index, UINT xmax, UINT ymax, REAL scale) {
	//	imfeatures.push_back(Info(xmax,ymax,scale,feat));
	//	imfeatures.push_back(HInfo(xmax,ymax,scale));
	if (index < imfeatures.size()) {
		imfeatures[index].Initialize(xmax, ymax, scale, norient);
		++counter;
	}
}
void HOGInfo::SetFeature(UINT xmax, UINT ymax, REAL scale) {
	//	imfeatures.push_back(Info(xmax,ymax,scale,feat));
	//	imfeatures.push_back(HInfo(xmax,ymax,scale));
	if (counter < imfeatures.size())
		imfeatures[counter++].Initialize(xmax, ymax, scale, norient);
}
void HOGInfo::Initialize(UINT nimages, UINT nalevel_, UINT norient_) {

	nalevel = nalevel_;
	if (!imfeatures.empty())
		imfeatures.clear();
	counter = 0;
	norient = norient_;
	imfeatures.resize(nimages);
}

void HOGInfo::WriteInfo(fstream &ofile) {
	ofile.write((char*) &nalevel, sizeof(UINT));
	UINT size = imfeatures.size();
	ofile.write((char*) &size, sizeof(UINT));
	for (UINT i = 0; i < imfeatures.size(); ++i)
		imfeatures[i].Write(ofile);
}
void HOGInfo::WriteTextInfo(fstream &ofile) {
	ofile << " Number of Above Levels =" << nalevel << endl
			<< " Total Number of Levels = " << imfeatures.size() << endl;
	//	ofile.write((char*)&sratio,sizeof(REAL));
	//	ofile.write((char*)&nalevel,sizeof(UINT));
	//	UINT size=imfeatures.size();
	for (UINT i = 0; i < imfeatures.size(); ++i)
		imfeatures[i].WriteText(ofile);
}
bool HOGInfo::ReadInfo(fstream &ifile) {
	if (!ifile)
		return false;

	ifile.read((char*) &nalevel, sizeof(UINT));
	UINT size;
	ifile.read((char*) &size, sizeof(UINT));
	imfeatures.resize(size);
	counter = size;
	for (UINT i = 0; i < size; ++i)
		imfeatures[i].Read(ifile);

	return true;
}

void HOGInfo::DisplayHOG(UINT index) {
	for (UINT i = 0; i < imfeatures[index].hogvalue.size(); ++i)
		cout << imfeatures[index].hogvalue[i] << " ";
	cout << endl;
}

UINT HOGInfo::GetNearestScaleIndex(REAL scale) const {
	REAL diff = ABS(imfeatures[0].scale-scale);
	UINT index = 0;
	for (UINT i = 1; i < imfeatures.size(); ++i)
		if (ABS(imfeatures[i].scale-scale) < diff) {
			diff = ABS(imfeatures[i].scale-scale);
			index = i;
		}
	return index;
}
void HOGInfo::FlipFeatures() {
	for (UINT i = 0; i < imfeatures.size(); ++i)
		imfeatures[i].FlipFeatures();
}
void HOGInfo::BlackOut(UINT x1, UINT y1, UINT x2, UINT y2) {
	for (UINT i = 0; i < imfeatures.size(); ++i)
		imfeatures[i].BlackOut(x1, y1, x2, y2);
}
