/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/
#ifndef WINDOWINFO_H_
#define WINDOWINFO_H_
#include "ncoord.h"

class WindowInfo: public Coord {
public:
	WindowInfo() {
		response = overlap = 0;
		cid = sindex = 0;
		scale = 1;
	}
	WindowInfo(const WindowInfo&);
	// xmin,ymin, width & height
//	WindowInfo(int, int, int, int, UINT, REAL = -10e3, REAL = 1, int = 1, REAL =
//			0);
	WindowInfo(int xmin_, int ymin_, int width, int height,
				UINT sindex_=1, REAL res=10e-3, REAL scale=1, int cat=1, REAL ol=0);
	void WriteText(ofstream &);
	void Write(ofstream &);
	void Read(ifstream &);

	void SetValues(int, int, int, int, UINT = 1, REAL = 10e-3, REAL = 1, int =
			1, REAL = 0);
	void SetScale(REAL scale_) {
		scale = scale_;
	}
	void SetOverlap(REAL ol) {
		overlap = ol;
	}
	void SetResponse(REAL res) {
		response = res;
	}
	REAL GetScale() const {
		return scale;
	}
	UINT GetScaleIndex() const {
		return sindex;
	} // used Because of hogInfo respective coordinates
	REAL GetResponse() const {
		return response;
	}
	REAL GetOverlap() const {
		return overlap;
	}
	void operator=(const WindowInfo&);
	double FindNMSOverlap(const WindowInfo& winfo) {
		return Coord::FindNMSOverlap(winfo.xmin, winfo.ymin, winfo.xmax,
				winfo.ymax);
	}
	void ScaleCoord(REAL);
	//	void Shrink(UINT); // Shrink the detections window final location by given size
	void Shrink(UINT svalue, UINT rwidth, UINT rheight);
	static bool CompareResponse(const WindowInfo &, const WindowInfo &);
	static bool CompareDescResponse(const WindowInfo &, const WindowInfo &);
	static bool CompareOverlap(const WindowInfo &, const WindowInfo &);
	static REAL FindMaximumResponse(const vector<WindowInfo>&, UINT&);
	static REAL FindMaximumResponse(const vector<WindowInfo>&, UINT, REAL,
			UINT&);
	static REAL FindMaximumOverlap(const vector<WindowInfo>&, UINT&);
	void WriteTextPascal(ofstream &ofile);
	int GetCID() {
		return cid;
	}
	void SetCID(const int &value) {
		cid = value;
	}
	void SetSIndex(const UINT &sindex_) {
		sindex = sindex_;
	}
	UINT GetSIndex() {
		return sindex;
	}
protected:
	REAL response, overlap;
	int cid; // component id
	UINT sindex;// scale index;
	REAL scale;
};

#endif /*WINDOWINFO_H_*/
