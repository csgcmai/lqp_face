/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/
#ifndef ONLYLTPFEATURES_H_
#define ONLYLTPFEATURES_H_
#include "lbpFeatures.h"
class OnlyLTPFeatures:public LBPFeatures
{
public:
	OnlyLTPFeatures(UINT cellsize_,UINT npoints_,UINT radius_,
			UINT width_,UINT height_,UINT nplevels,NORM norm_,
			UINT nflevels_,UINT lbpstride_,ProcessImage &pim_,bool folded_,REAL tol)
	:LBPFeatures(cellsize_,npoints_,radius_,width_,
			height_,nplevels,norm_,nflevels_,lbpstride_,pim_,folded_){
		tolerance = tol;
		cout<<" Dim Feature = "<< GetDim(width,height)<< endl;
		cout<<" Tolerance = " << tol <<endl;
	}
	virtual	UINT GetDim(UINT width_,UINT height_)const{
		return nflevels==1 ? (width_/lbpstride) * (height_/lbpstride) * nbins*2
				: ( ( (width_/lbpstride) * (height_/lbpstride) ) +
						( (width_/(2*lbpstride)) * (height_/(2*lbpstride))  ) ) * nbins*2;
	}
	virtual	UINT GetFoldedDim(UINT width_,UINT height_)const{ // for horizontal symmetry
		return nflevels==1 ? (width_/(2*lbpstride)) * (height_/lbpstride) * 2 *nbins
				: ( ( (width_/(2*lbpstride)) * (height_/lbpstride) ) +
						( width_/(4*lbpstride) * height_/(2*lbpstride)  ) ) * 2 *nbins;
	}


	virtual	 void InitalizeMaps( Image &,PyramidType);
	void ComputeLTPMap( Image &image,LBPMap &lbpmap,LBPMap &tnmap);
	virtual ~OnlyLTPFeatures(){}
	virtual	void GetFeatures(UINT,int,int,vector<REAL>&);
//	virtual	void GetFoldedFeatures(UINT,int,int,vector<REAL>&);
//	virtual void UnFoldWeights(REAL*,vector<REAL>&);
	virtual	UINT GetXSize(UINT index){
		return lbpmaps[index].nxpoints+radius;
	}
	virtual UINT GetYSize(UINT index){
		return lbpmaps[index].nypoints+radius;
	}
protected:
	vector<LBPMap> negmap;
	REAL tolerance;
};

#endif /* HOGLBPFEATURES_H_ */
