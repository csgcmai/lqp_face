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

class WindowInfoFeat:public WindowInfo
{
public:
	WindowInfoFeat()
	{}
	WindowInfoFeat(UINT size)
	{
		features.resize(size,0);
	}
	WindowInfoFeat(int xmin_,int ymin_,int width,
			int height,UINT sindex_,REAL res,REAL scale,REAL ol)
	:WindowInfo(xmin_,ymin_,width,height,sindex_,res,scale,0,ol)
	{

	}
	WindowInfoFeat(int xmin_,int ymin_,int width,
			int height,UINT sindex_,REAL res,REAL scale,REAL ol,vector<REAL>&feat)
	:WindowInfo(xmin_,ymin_,width,height,sindex_,res,scale,0,ol)
	{
		features=feat;
	}
	WindowInfoFeat(int xmin_,int ymin_,int width,
			int height,UINT sindex_,REAL res,REAL scale,REAL ol,UINT size)
	:WindowInfo(xmin_,ymin_,width,height,sindex_,res,scale,0,ol)
	{
		features.resize(size,0);
	}


	WindowInfoFeat(const WindowInfoFeat& winfo)
	:WindowInfo(winfo)
	{
		features = winfo.features;
	}
	void SetFeatures(vector<REAL>&feat)
	{
		features=feat;
	}
	vector<REAL>& GetFeatureRef(){return features;}
	void GetFeat(vector<REAL>& feat){feat=features;}

	void operator=(const WindowInfoFeat & winfo)
	{
		WindowInfo::operator=(winfo);
		features = winfo.features;
	}

	void SetValues(int xmin_,int ymin_,int width,
			int height,UINT sindex_,REAL res,REAL scale,REAL ol,vector<REAL>&feat)
	{
		WindowInfo::SetValues(xmin_,ymin_,width,height,sindex_,res,scale
				,0,ol);
		features = feat;
	}
	void SetValues(int xmin_,int ymin_,int width,
			int height,UINT sindex_,REAL res,REAL scale,REAL ol,UINT size)
	{
		WindowInfo::SetValues(xmin_,ymin_,width,height,sindex_,res,scale
				,0,ol);
		if(features.size()!=size)
		features.resize(size,0);
	}
	void SetValues(WindowInfo winfo,UINT size)
		{
			WindowInfo::operator=(winfo);
			features.resize(size,0);
		}
	void WriteText(ofstream& ofile)
	{
		WindowInfo::WriteText(ofile);
//		ofile<<endl;
//		for(int i=0; i < features.size();++i)
//			ofile<<features[i]<<" ";
	}
private:
	vector<REAL> features;
};


#endif /* WINDOWINFOFEAT_H_ */
