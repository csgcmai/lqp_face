/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/
#include "onlyLTPFeatures.h"
void OnlyLTPFeatures::ComputeLTPMap(Image &image, LBPMap &lbpmap, LBPMap &tnmap) {
	ColorGray mycolor, c1, c2, c3, c4;
	UINT xdim = image.columns(), ydim = image.rows();
	int nypoints = ydim - 2 * radius, reg = 1, fx, fy, cx, cy, nxpoints = xdim
			- 2 * radius, cenx, ceny, count, value, posvalue, negvalue;
	REAL tmpx, tmpy, fracx, fracy, pixval, cenval;
	REAL x, y;
	lbpmap.Init(nxpoints, nypoints);
	tnmap.Init(nxpoints, nypoints);
	vector<REAL> impix(xdim * ydim, 0);
	pim.Process(image, impix);
	for (int i = 0; i < nypoints; ++i) {
		for (int j = 0; j < nxpoints; ++j) {
			cenx = j + radius;
			ceny = i + radius;
			value = posvalue = negvalue = 0;
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

				value += pixval >= cenval + tolerance ? reg << count : 0;
				negvalue += pixval <= cenval - tolerance ? reg << count : 0;

			}
			lbpmap(j, i) = map[value];
			tnmap(j, i) = map[negvalue];
		}
	}
}
void OnlyLTPFeatures::InitalizeMaps(Image& image, PyramidType sspace_) {
	UINT pycount = 0;
	// Read the image....
	sspace = sspace_;
	if (sspace_ == SingleRes) {
		pycount = pyobj.GeneratePyramid(image);
		lbpmaps.resize(pycount);
		negmap.resize(pycount);
		for (UINT i = 0; i < pycount; i++) {
			Image& imgref = pyobj.GetImageRef(i);
			// Compute the Gradient Image
			ComputeLTPMap(imgref, lbpmaps[i], negmap[i]);
		}
	} else {
		lbpmaps.resize(1);
		negmap.resize(1);
		ComputeLTPMap(image, lbpmaps[0], negmap[0]);
	}
}
void OnlyLTPFeatures::GetFeatures(UINT index, int x, int y, vector<REAL>&feat) {

	if (!usegradmag) {
		UINT offset = 0;
		fill(feat.begin(), feat.end(), 0);// Initialize with 0's;
		ComputeHistogram(x, y, lbpmaps[index], &feat[offset]);
		NormalizeFeatures(&feat[offset]);
		offset += nxcells * nycells * nbins;
		ComputeHistogram(x, y, negmap[index], &feat[offset]);
		NormalizeFeatures(&feat[offset]);
	} else {
		cerr
				<< " \n Gradient Magnitude is only used during LBP Feature Computations"
				<< endl;
		exit(EXIT_FAILURE);
	}

}
/*

 void OnlyLTPFeatures::GetFoldedFeatures(UINT index,int x, int y ,vector<REAL>&feat)
 {

 UINT offset=0;
 x = MAX(0,x - radius); // to cater for offset
 y = MAX(0,y - radius);
 UINT	twidth = x+width < lbpmaps[index].nxpoints ?
 x+width : lbpmaps[index].nxpoints,
 theight = y+height < lbpmaps[index].nypoints ?
 y+height : lbpmaps[index].nypoints;
 fill(feat.begin(),feat.end(),0);// Initialize with 0's;
 twidth /= 2; //half of the width
 int fx =  x + width - lbpstride; // flipped x index
 vector<REAL> tfeat(nbins,0),tpfeat(nbins,0),tnfeat(nbins,0);
 for(UINT i=y ; i < theight; i+=lbpstride)
 {
 fx =  x + width - lbpstride;
 for(UINT j=x ; j < twidth; j+=lbpstride)
 {
 // LBP Features
 lbpmaps[index].ComputeHist(j,i,cellsize,cellsize,&feat[offset]);
 NormalizeFeatures(&feat[offset]);

 lbpmaps[index].ComputeHist(fx,i,cellsize,cellsize,&tfeat[0]);
 NormalizeFeatures(&tfeat[0]);


 // LTP +Ve features
 negmap[index].ComputeHist(j,i,cellsize,cellsize,&feat[offset+2*nbins]);
 NormalizeFeatures(&feat[offset+2*nbins]);

 negmap[index].ComputeHist(fx,i,cellsize,cellsize,&tnfeat[0]);
 NormalizeFeatures(&tnfeat[0]);



 for(UINT iter=0; iter < nbins;++iter)
 {
 feat[offset+iter] += tfeat[iter];
 feat[offset+nbins+iter] += tpfeat[iter];
 feat[offset+2*nbins+iter] += tnfeat[iter];
 }

 fill(tfeat.begin(),tfeat.end(),0);// Initialize with 0's;
 fill(tpfeat.begin(),tpfeat.end(),0);// Initialize with 0's;
 fill(tnfeat.begin(),tnfeat.end(),0);// Initialize with 0's;

 offset += 3*nbins;
 fx -= lbpstride;
 }
 }
 }

 /*void OnlyLTPFeatures::UnFoldWeights(REAL *input,vector<REAL>&weights)
 {
 UINT twidth = width/2,offset=0,count=0,
 foffset=0,boffset=0,inoffset=0;// forward & backward offsets

 for(UINT i=0 ; i < height; i+=lbpstride,++count)
 {
 boffset = ((width - lbpstride)/lbpstride)*nbins*3;
 foffset=0;
 offset = count * (width/lbpstride) * nbins*3;
 for(UINT j=0 ; j < twidth; j+=lbpstride)
 {
 for(UINT iter=0; iter < 3*nbins;++iter)
 {
 weights[iter+offset+foffset] = input[inoffset + iter];
 weights[iter+offset+boffset] = input[inoffset + iter];
 }

 foffset += nbins*3;
 boffset -= nbins*3;
 inoffset += nbins*3;
 }

 }

 }


 // for fast computations....
 //void OnlyLTPFeatures::ComputeHist(UINT index,UINT x, UINT y,
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
 */
