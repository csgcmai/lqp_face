/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/

#ifndef WINDOWINFOFEAT_H_
#define WINDOWINFOFEAT_H_
#include "windowInfo.h"

class WindowInfoFeat: public WindowInfo {
public:
	WindowInfoFeat() {
	}
	WindowInfoFeat(UINT size) {
		features.resize(size, 0);
	}
	WindowInfoFeat(int cindex_, UINT sindex_, vector<REAL>& feat) {
		cid = cindex_;
		sindex = sindex_;
		features = feat;
		response = overlap = 0;
		scale = 1;
	}
	//	WindowInfoFeat(int xmin_,int ymin_,int width,
	//			int height,UINT sindex_,REAL res,REAL scale,UINT cid,REAL ol)
	//	:WindowInfo(xmin_,ymin_,width,height,sindex_,res,scale,cid,ol)
	//	{
	//
	//	}
	WindowInfoFeat(int xmin_, int ymin_, int width, int height, UINT sindex_,
			REAL res, REAL scale, UINT cid, REAL ol, vector<REAL>&feat) :
		WindowInfo(xmin_, ymin_, width, height, sindex_, res, scale, cid, ol) {
		features = feat;
	}
	WindowInfoFeat(int xmin_, int ymin_, int width, int height, UINT sindex_,
			REAL res, REAL scale, UINT cid, REAL ol, UINT size) :
		WindowInfo(xmin_, ymin_, width, height, sindex_, res, scale, cid, ol) {
		if (size != 0)
			features.resize(size, 0);
	}

	WindowInfoFeat(const WindowInfoFeat& winfo) :
		WindowInfo(winfo) {
		features = winfo.features;
	}
	WindowInfoFeat(const WindowInfo& winfo, UINT featsize) :
		WindowInfo(winfo) {
		features.resize(featsize);
	}
	void SetFeatures(vector<REAL>&feat) {
		features = feat;
	}
	vector<REAL>& GetFeatureRef() {
		return features;
	}
	REAL *GetFeaturePtr() {
		return &features[0];
	}
	void GetFeat(vector<REAL>& feat) {
		feat = features;
	}

	void operator=(const WindowInfoFeat & winfo) {
		WindowInfo::operator=(winfo);
		features = winfo.features;
	}

	void SetValues(int xmin_, int ymin_, int width, int height, UINT sindex_,
			REAL res, REAL scale, UINT cid, REAL ol, vector<REAL>&feat) {
		WindowInfo::SetValues(xmin_, ymin_, width, height, sindex_, res, scale,
				cid, ol);
		features = feat;
	}
	void SetValues(int xmin_, int ymin_, int width, int height, UINT sindex_,
			REAL res, REAL scale, UINT cid, REAL ol, UINT size) {
		WindowInfo::SetValues(xmin_, ymin_, width, height, sindex_, res, scale,
				cid, ol);
		if (features.size() != size)
			features.resize(size, 0);
	}
	void SetValues(WindowInfo winfo, UINT size) {
		WindowInfo::operator=(winfo);
		features.resize(size, 0);
	}
	void SetValues(int cindex_, UINT sindex_, REAL *feat) {
		//		features.resize(size, 0);
		cid = cindex_;
		sindex = sindex_;
		copy(feat, feat + features.size(), features.begin());
	}
	void WriteText(ofstream& ofile) {
		WindowInfo::WriteText(ofile);
	}
	static REAL FindMaximumResponse(const vector<WindowInfoFeat>& lwinfo,
			UINT &index) {
		REAL maxres = lwinfo[0].response;
		index = 0;
		UINT k = 1;
		for (vector<WindowInfoFeat>::const_iterator iter = lwinfo.begin() + 1; iter
				!= lwinfo.end(); ++iter, ++k)
			if (iter->response > maxres) {
				index = k;
				maxres = iter->response;
			}

		return maxres;
	}

private:
	vector<REAL> features;
};

#endif /* WINDOWINFOFEAT_H_ */
