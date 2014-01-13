/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/
#include "lbpMap.h"
void LBPMap::ComputeHist(UINT x, UINT y, UINT width, UINT height,
		vector<REAL>&hist) {
	height = MIN(height+y,nypoints);
	width = MIN(width+x,nxpoints);
	int binloc;
	for (UINT i = y; i < height; ++i)
		for (UINT j = x; j < width; ++j) {
			binloc = lbpmap[i * nxpoints + j];
			if (binloc >= 0)
				hist[binloc]++;
		}
}
void LBPMap::ComputeHist(UINT x, UINT y, UINT width, UINT height, REAL *hist) {
	height = MIN(height+y,nypoints);
	width = MIN(width+x,nxpoints);
	int binloc;
	for (UINT i = y; i < height; ++i)
		for (UINT j = x; j < width; ++j) {
			binloc = lbpmap[i * nxpoints + j];
			if (binloc >= 0)
				hist[binloc]++;
		}
}
// Linear interpolation based histogram computation ...
// build cell based histogram...

void LBPMap::GetMap(UINT x, UINT y, UINT width, UINT height, UINT *hist) const {
	height = MIN(height+y,nypoints);
	width = MIN(width+x,nxpoints);
	UINT count = 0;
	for (UINT i = y; i < height; ++i)
		for (UINT j = x; j < width; ++j)
			hist[count++] = lbpmap[i * nxpoints + j];
}

// histograms are concatenated in the column major order
// and also y*x*nbins+x*nbins is used to access the particular histogram...
//void LBPMap::ComputeHist(UINT xstart, UINT ystart, UINT width, UINT height,UINT nbins, REAL *hist) {
//
//	int binloc;
//	int xind=0,yind=0,foffset=0,hs = log2(height), ws =log2(width);
//	for (UINT i = ystart; i < nypoints; ++i)
//	{
//		yind = i >> hs;
//		foffset  = yind *
//		for (UINT j = xstart; j < nxpoints; ++j) {
////			foffset = (j >> ws) * nbins;
//			++hist[ (j >> ws) * nbins + lbpmap[i * nxpoints + j]];
//		}
//	}
//}

void LBPMap::PrintMap() {
	for (UINT i = 0, count = 0; i < nypoints; ++i) {
		for (UINT j = 0; j < nxpoints; ++j, ++count)
			cout << lbpmap[count] << " ";
		cout << endl;
	}
	cout << endl;
}
void LBPMap::WriteMap(const string &fname) const {
	//	ofile<<" NxPoints  "<< nxpoints << " NyPoints "
	//	<<nypoints<< endl;
	ofstream ofile(fname.c_str(), ios::out);
	for (UINT i = 0, count = 0; i < nypoints; ++i) {
		for (UINT j = 0; j < nxpoints; ++j, ++count)
			ofile << lbpmap[count] << " ";
		ofile << endl;
	}
	ofile << endl;
}
void LBPMap::WriteMap(ofstream &ofile) const {
	//	ofile<<" NxPoints  "<< nxpoints << " NyPoints "
	//	<<nypoints<< endl;
	for (UINT i = 0, count = 0; i < nypoints; ++i) {
		for (UINT j = 0; j < nxpoints; ++j, ++count)
			ofile << lbpmap[count] << " ";
		ofile << endl;
	}
	ofile << endl;
}

/*
 * UINT nxbmax = xbmax+2*padx,nybmax=ybmax+2*pady;
 vector<REAL> feat(nxbmax*nybmax*norient,0);
 UINT counter=0,
 offset = (padx + pady * nxbmax) * norient,coffset,roffset;
 for(UINT i= 0; i < ybmax; ++i)
 {
 roffset = offset +  i * nxbmax* norient;
 for(UINT j=0; j < xbmax;++j)
 {
 coffset = roffset + j*norient;
 for(UINT k=0; k < norient;++k)
 feat[coffset+k]=hogvalue[counter++];
 }
 }
 xbmax = nxbmax; ybmax=nybmax;
 hogvalue=feat;*/
void LBPMap::PadFeatureMap(UINT padx, UINT pady) {// to pad the boundary values with -1
	UINT xpoints = nxpoints + 2 * padx, ypoints = nypoints + 2 * pady;
	vector<int> feat(xpoints * ypoints, -1);
	UINT counter = 0, roffset, offset = (padx + pady * xpoints);
	for (UINT i = 0; i < nypoints; ++i) {
		roffset = offset + i * xpoints;
		for (UINT j = 0; j < nxpoints; ++j) {
			feat[roffset + j] = lbpmap[counter++];
		}
	}
	nxpoints = xpoints;
	nypoints = ypoints;
	lbpmap = feat;
}

