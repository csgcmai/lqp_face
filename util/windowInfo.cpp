/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/
#include "windowInfo.h"
// takes the upper coordinates along with window width & height....
WindowInfo::WindowInfo(int xmin_, int ymin_, int width, int height,
		UINT sindex_, REAL res, REAL scale, int cat, REAL ol) {
	SetValues(xmin_, ymin_, width, height, sindex_, res, scale, cat, ol);
}

WindowInfo::WindowInfo(const WindowInfo &winfo) :
	Coord(winfo), response(winfo.response), overlap(winfo.overlap), cid(
			winfo.cid), sindex(winfo.sindex), scale(winfo.scale) {
	//	features = winfo.features;
	//	xmin = winfo.xmin;
	//	ymin = winfo.ymin;
	//	xmax = winfo.xmax;
	//	ymax = winfo.ymax;
}
bool WindowInfo::CompareResponse(const WindowInfo & obj1,
		const WindowInfo & obj2) {
	if (obj1.response <= obj2.response)
		return true;
	return false;
}
bool WindowInfo::CompareDescResponse(const WindowInfo & obj1,
		const WindowInfo & obj2) {
	if (obj1.response >= obj2.response)
		return true;
	return false;
}
bool WindowInfo::CompareOverlap(const WindowInfo & obj1,
		const WindowInfo & obj2) {
	if (obj1.overlap <= obj2.overlap)
		return true;
	return false;
}
void WindowInfo::Write(ofstream &ofile) {
	Coord::Write(ofile);
	ofile.write((char*) &sindex, sizeof(UINT));
	ofile.write((char*) &response, sizeof(REAL));
	ofile.write((char*) &cid, sizeof(int));
	ofile.write((char*) &overlap, sizeof(REAL));
	ofile.write((char*) &scale, sizeof(REAL));
}
void WindowInfo::Read(ifstream &ifile) {
	Coord::Read(ifile);
	ifile.read((char*) &sindex, sizeof(UINT));
	ifile.read((char*) &response, sizeof(REAL));
	ifile.read((char*) &cid, sizeof(int));
	ifile.read((char*) &overlap, sizeof(REAL));
	ifile.read((char*) &scale, sizeof(REAL));
}
void WindowInfo::WriteText(ofstream &ofile) {
	//	Coord::WriteText(ofile);
	ofile << xmin << " " << ymin << " " << xmax << " " << ymax << " " << sindex
			<< " " << response << " " << overlap << " " << cid << " " << scale
			<< endl;
}
void WindowInfo::WriteTextPascal(ofstream &ofile) {
	ofile << " " << response << " " << xmin << " " << ymin << " " << xmax
			<< " " << ymax << endl;
}
void WindowInfo::operator =(const WindowInfo & winfo) {
	Coord::operator=(winfo);
	overlap = winfo.overlap;
	sindex = winfo.sindex;
	response = winfo.response;
	cid = winfo.cid;
	scale = winfo.scale;
	//	features = winfo.features;


}
void WindowInfo::SetValues(int xmin_, int ymin_, int xmax_, int ymax_,
		UINT sindex_, REAL res, REAL scale_, int cat, REAL overlap_) {

	Coord::SetCoord(xmin_, ymin_, xmax_, ymax_);
	sindex = sindex_;
	response = res;
	cid = cat;
	overlap = overlap_;
	scale = scale_;
}

REAL WindowInfo::FindMaximumResponse(const vector<WindowInfo>& lwinfo,
		UINT &index) {
	REAL maxres = lwinfo[0].response;
	index = 0;
	UINT k = 1;
	for (vector<WindowInfo>::const_iterator iter = lwinfo.begin() + 1; iter
			!= lwinfo.end(); ++iter, ++k)
		if (iter->response > maxres) {
			index = k;
			maxres = iter->response;
		}

	return maxres;
}
REAL WindowInfo::FindMaximumResponse(const vector<WindowInfo>& lwinfo,
		UINT olindex, REAL percen, UINT &index) {
	REAL maxres = lwinfo[olindex].response, olthres = lwinfo[olindex].overlap
			* percen;
	index = olindex;
	UINT k = 1;
	for (vector<WindowInfo>::const_iterator iter = lwinfo.begin(); iter
			!= lwinfo.end(); ++iter, ++k)
		if (iter->response > maxres && iter->overlap > olthres) {
			index = k;
			maxres = iter->response;
		}

	return maxres;
}

REAL WindowInfo::FindMaximumOverlap(const vector<WindowInfo>& lwinfo,
		UINT &index) {
	REAL maxol = lwinfo[0].overlap;
	index = 0;
	UINT k = 1;
	for (vector<WindowInfo>::const_iterator iter = lwinfo.begin() + 1; iter
			!= lwinfo.end(); ++iter, ++k)
		if (iter->overlap > maxol) {
			index = k;
			maxol = iter->overlap;
		}

	return maxol;
}
// Scaling for the  Non-Max-Suppression....
void WindowInfo::ScaleCoord(REAL scale) {
	xmin = xmin * scale;
	ymin = ymin * scale;
	xmax = xmax * scale;
	ymax = ymax * scale;
}

void WindowInfo::Shrink(UINT svalue, UINT rwidth, UINT rheight) {
	UINT sloc = (REAL) svalue / 2.0 * scale;
	xmin += sloc;
	ymin += sloc;
	xmax = xmin + (rwidth - svalue) * scale;
	ymax = ymin + (rheight - svalue) * scale;
}
