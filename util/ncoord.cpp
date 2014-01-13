/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/
#include "ncoord.h"
void Coord::Write(ofstream &ofile) {
	ofile.write((char*) &xmin, sizeof(int));
	ofile.write((char*) &ymin, sizeof(int));
	ofile.write((char*) &xmax, sizeof(int));
	ofile.write((char*) &ymax, sizeof(int));
}
void Coord::WriteText(ofstream &ofile) {
	ofile << " xmin = " << xmin << ", ymin = " << ymin << " ,xmax = " << xmax
			<< " ,ymax = " << ymax << endl;
}
void Coord::Display() {
	cout << " xmin = " << xmin << ", ymin = " << ymin << " ,xmax = " << xmax
			<< " ,ymax = " << ymax << endl;
}
void Coord::Read(ifstream &ifile) {
	ifile.read((char*) &xmin, sizeof(int));
	ifile.read((char*) &ymin, sizeof(int));
	ifile.read((char*) &xmax, sizeof(int));
	ifile.read((char*) &ymax, sizeof(int));
}
REAL Coord::FindOverlap(const Coord &winfo, REAL scale) const
// Overlap = Area(A n B)/Area(A u B)
{
	//
	int w = MIN(xmax/scale, winfo.xmax) - MAX(xmin/scale,winfo.xmin), h =
			MIN(ymax/scale, winfo.ymax) - MAX(ymin/scale,winfo.ymin), anb,
			aub;
	if (w > 0 && h > 0) {
		anb = w * h;
		aub = GetHeight() / scale * GetWidth() / scale + winfo.GetWidth()
				* winfo.GetHeight() - anb;

		return (REAL) anb / aub;
	} else
		return 0;
}

REAL Coord::FindOverlap(int xmin_, int ymin_, int xmax_, int ymax_, REAL scale) const {
	return FindOverlap(Coord(xmin_, ymin_, xmax_, ymax_), scale);
}

//REAL Coord::FindOverlap(UINT xmini,UINT ymini,UINT width,UINT height)const
//{
//	return FindOverlap(Coord(xmini,ymini,xmini+width,ymini+height));
//}
REAL Coord::FindOverlap(const Coord &winfo) const
// Overlap = Area(A n B)/Area(A u B)
{
	//
	
	int w = MIN(xmax, winfo.xmax) - MAX(xmin,winfo.xmin) + 1, h =
			MIN(ymax, winfo.ymax) - MAX(ymin,winfo.ymin) + 1, anb, aub;
	if (w > 0 && h > 0) {
		anb = w * h;
		aub = GetHeight() * GetWidth() + winfo.GetWidth() * winfo.GetHeight()
				- anb;

		return (REAL) anb / aub;
	} else
		return 0;
}
REAL Coord::FindNMSOverlap(int ixmin, int iymin, int ixmax, int iymax) const {
	int w = MIN(xmax, ixmax) - MAX(xmin,ixmin) + 1, h = MIN(ymax, iymax)
			- MAX(ymin,iymin) + 1, anb;
	if (w > 0 && h > 0) {
		anb = w * h;
		return (REAL) anb / ((ixmax - ixmin + 1) * (iymax - iymin + 1));
	} else
		return 0;
}
REAL Coord::FindOverlap(int ixmin, int iymin, int ixmax, int iymax) const {
	return FindOverlap(Coord(ixmin, iymin, ixmax, iymax));
}
REAL Coord::FindOverlap(const Coord &c1, const Coord &c2) {
	return c1.FindOverlap(c2);
}
