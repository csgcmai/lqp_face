/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/

#ifndef CACHE_H_
#define CACHE_H_
//#include "../features/features.h"
#include <cmath>
struct Store {
	Store() :
		width(0), height(0) {
	}
	Store(const Store& st) {
		*this = st;
	}
	Store(UINT w, UINT h, UINT val = 0) {
		Initialize(w, h, val);
	}
	void clear() {
		width = 0;
		height = 0;
		values.clear();
	}
	void operator=(const Store &st) {
		width = st.width;
		height = st.height;
		values = st.values;
	}
	void Initialize(UINT w, UINT h, UINT val = 0) {
		width = w;
		height = h;
		values.resize(w * h, val/*NAN*/);
		fill(values.begin(),values.end(),val);
	}
	REAL GetValue(UINT x, UINT y) {
		return values[y * width + x];
	}
	REAL GetValue(UINT x, UINT y, UINT xoffset, UINT yoffset) {
		return values[(y + yoffset) * width + (x + xoffset)];
	}
	void SetValue(UINT x, UINT y, REAL val) {
		values[y * width + x] = val;
	}
	void AddValue(UINT x, UINT y, REAL val) {
		values[y * width + x] += val;
	}
	void PrintVector() {
		cout << endl;
		for (vector<REAL>::iterator lhiter = values.begin(); lhiter
				!= values.end(); ++lhiter)
			cout << *lhiter << " ";
		cout << flush;

	}
	void operator+=(const Store & st) {

		vector<REAL>::const_iterator rhiter = st.values.begin();
		for (vector<REAL>::iterator lhiter = values.begin(); lhiter
				!= values.end(); ++lhiter, ++rhiter)
			*lhiter += *rhiter;

	}

	UINT width, height;
	vector<REAL> values;
};
/*
 class Cache2D
 {
 public:
 Cache2D(Features & fobj)
 {
 offset = fobj.GetInitIndex();
 cache.resize(fobj.GetNLevels()-offset);
 for(UINT i=offset; i < fobj.GetNLevels();++i)
 cache[i-offset].Initialize(fobj.GetXDim(i),fobj.GetYDim(i));
 }

 REAL GetValue(Features&,UINT,UINT,UINT);
 REAL GetValue(UINT index,UINT x,UINT y )
 {
 return cache[index-offset].GetValue(x,y);
 }
 void  SetValue(UINT index,UINT x,UINT y,REAL val )
 {
 cache[index-offset].SetValue(x,y,val);
 }
 private:
 vector<Store> cache;
 UINT offset;
 };



 class Cache3D
 {
 public:

 Cache3D(Features & fobj,UINT nparts)
 {
 cache.resize(nparts,vector<Store>(fobj.GetNLevels()));

 for(UINT j=0; j < nparts;++j)
 {
 for(UINT i=0; i < fobj.GetNLevels();++i)
 //#ifndef CELL_ORDER
 cache[j][i].Initialize(fobj.GetXDim(i),fobj.GetYDim(i));
 //#else
 //			cache[j][i].Initialize(fobj.GetXDim(i)/fobj.xcell,fobj.GetYDim(i)/fobj.ycell);
 //#endif
 }
 }
 Cache3D(HOGInfo & hoginfo,UINT nparts)
 {
 cache.resize(nparts,vector<Store>(hoginfo.GetNLevels()));

 for(int j=0; j < nparts;++j)
 {
 for(UINT i=0; i < hoginfo.GetNLevels();++i)
 cache[j][i].Initialize(hoginfo.GetXBlock(i),hoginfo.GetYBlock(i));
 }
 }


 REAL GetValue(Features&,UINT,UINT,UINT);
 REAL GetValue(UINT pindex,UINT index,UINT x,UINT y )
 {
 return cache[pindex][index].GetValue(x,y);
 }
 void  SetValue(UINT pindex,UINT index,UINT x,UINT y,REAL val )
 {
 cache[pindex][index].SetValue(x,y,val);
 }
 private:
 vector<vector<Store> > cache;
 };
 */

#endif /* CACHE_H_ */
