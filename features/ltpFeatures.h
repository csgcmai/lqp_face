/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/
#ifndef LTPFEATURES_H_
#define LTPFEATURES_H_
#include "lbpFeatures.h"
class LTPFeatures:public LBPFeatures
{
public:
	LTPFeatures(UINT cellsize_,UINT npoints_,UINT radius_,
			UINT width_,UINT height_,UINT nplevels,NORM norm_,
			UINT nflevels_,UINT lbpstride_,ProcessImage &pim_,bool folded_,REAL tol)
	:LBPFeatures(cellsize_,npoints_,radius_,width_,
			height_,nplevels,norm_,nflevels_,lbpstride_,pim_,folded_){
		tolerance = tol;
		cout<<" Dim Feature = "<< GetDim(width,height)<< endl;
		cout<< " Folded Dim = "<< GetFoldedDim(width,height)<< endl;
	}
	virtual	UINT GetDim(UINT width_,UINT height_)const{
		return nflevels==1 ? (width_/lbpstride) * (height_/lbpstride) * nbins*3
				: ( ( (width_/lbpstride) * (height_/lbpstride) ) +
						( (width_/(2*lbpstride)) * (height_/(2*lbpstride))  ) ) * nbins*3;
	}
	virtual	UINT GetFoldedDim(UINT width_,UINT height_)const{ // for horizontal symmetry
		return nflevels==1 ? ceil((REAL)width_/(2*lbpstride)) * (height_/lbpstride) * 3 *nbins
				: ( ( (width_/(2*lbpstride)) * (height_/lbpstride) ) +
						( width_/(4*lbpstride) * height_/(2*lbpstride)  ) ) * 3 *nbins;
	}


	virtual	 void InitalizeMaps( Image &,PyramidType);
	void ComputeLTPMap( Image &image,LBPMap &lbpmap,LBPMap& tpmap,LBPMap &tnmap);
	void ComputeLTPMap( Image &image,
			LBPMap &lbpmap,LBPMap & flbmap, LBPMap& tpmap,
			LBPMap& ftpmap,
			LBPMap &tnmap,LBPMap &ftnmap);
	virtual ~LTPFeatures(){}
	virtual	void GetFeatures(UINT,int,int,vector<REAL>&);
	virtual	void GetFoldedFeatures(UINT,int,int,vector<REAL>&);
	virtual void UnFoldWeights(REAL*,vector<REAL>&);
	virtual	UINT GetXSize(UINT index){
		return lbpmaps[index].nxpoints+radius;
	}
	virtual UINT GetYSize(UINT index){
		return lbpmaps[index].nypoints+radius;
	}
protected:
	vector<LBPMap> posmap;
	vector<LBPMap> fposmap;
	vector<LBPMap> negmap;
	vector<LBPMap> fnegmap;
	REAL tolerance;
};

#endif /* HOGLBPFEATURES_H_ */
