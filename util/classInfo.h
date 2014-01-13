/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/

#ifndef CLASSINFO_H_
#define CLASSINFO_H_
//#include "util.h"
#include "ncoord.h"
struct ImageInfo {
	ImageInfo(const string &name) :
		iname(name) {
	}
	vector<Coord> coord;
	string iname;
};
class ClassInfo {
public:
	ClassInfo() {
	}
	ClassInfo(const string &fname) :
		width(0), height(0), counter(0), ocounter(0) {
		ReadInfo(fname);
	}
	void Initialize(const string &fname) {
		width = 0;
		height = 0;
		counter = ocounter = 0;
		ReadInfo(fname);
	}
	bool GetNextImageAnnotation(string&, vector<Coord>&);
	bool GetImageAnnotation(string&, vector<Coord>&);
	void ReadInfo(const string&);
	UINT GetWidth() {
		return width;
	}
	UINT GetHeight() {
		return height;
	}
	void ResetCounter() {
		counter = 0;
	}
	UINT GetNObjects() {
		return ocounter;
	}
	void GetImageList(vector<string>& ilist) {
		ilist.resize(iinfo.size());
		for (UINT i = 0; i < iinfo.size(); ++i)
			ilist[i] = iinfo[i].iname;
	}
	UINT GetNImages() {
		return iinfo.size();
	}
private:
	UINT width, height;
	vector<ImageInfo> iinfo;
	UINT counter, ocounter;// image & object counters.
};

#endif /* CLASSINFO_H_ */
