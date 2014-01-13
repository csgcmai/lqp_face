/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/
/*version 2: Replaces the UINT with integer to cater for
 * padding as -1 bin location means don't bin it....*/
#ifndef LBPMAP_H_
#define LBPMAP_H_
#include "util.hpp"
struct LBPMap {
	void Init(UINT nx, UINT ny, int val = 0) {
		nxpoints = nx;
		nypoints = ny;
		lbpmap.resize(nx * ny, val);
	}
	void ComputeHist(UINT, UINT, UINT, UINT, vector<REAL>&);
	void ComputeHist(UINT, UINT, UINT, UINT, REAL*);

	void GetMap(UINT, UINT, UINT, UINT, UINT*) const;
	int & operator()(UINT x, UINT y) {
		return lbpmap[y * nxpoints + x];
	}
	int operator()(UINT x, UINT y) const {
		return lbpmap[y * nxpoints + x];
	}
	void PadFeatureMap(UINT padx, UINT pady);
	void PrintMap();
	void SetValue(int val) {
		fill(lbpmap.begin(), lbpmap.end(), val);
	}
	void WriteMap(ofstream &) const;

	int GetValue(const UINT &x, const UINT &y) {
		return lbpmap[y * nxpoints + x];
	}
	int GetValueInBounds(const int &x, const int &y) {
		int tx = MIN(MAX(x,0),nxpoints-1), ty = MIN(MAX(y,0),nypoints-1);
		return lbpmap[ty * nxpoints + tx];
	}
	void WriteMap(UINT x, UINT y, UINT width, UINT height, const string &fname) {

		ofstream ofile(fname.c_str(), ios::out);
		for (UINT j = 0; j < height - 2; ++j) {
			for (UINT i = 0; i < width - 2; ++i) {
				ofile << GetValue(i + x, j + y) << " ";
			}
			ofile << endl;
		}
		ofile.close();
	}
	void WriteMap(const string &fname) const;
	UINT nxpoints, nypoints;
	vector<int> lbpmap;
};
template<class T>
struct GenericMap {
	void Init(UINT nx, UINT ny, int val = 0) {
		nxpoints = nx;
		nypoints = ny;
		map.resize(nx * ny, val);
	}
	T & operator()(UINT x, UINT y) {
		return map[y * nxpoints + x];
	}
	T operator()(UINT x, UINT y) const {
		return map[y * nxpoints + x];
	}
	void SetValue(int val) {
		fill(map.begin(), map.end(), val);
	}
	T GetValue(const UINT &x, const UINT &y) {
		return map[y * nxpoints + x];
	}
	UINT nxpoints, nypoints;
	vector<T> map;
};
#endif /* LBPMAP_H_ */
