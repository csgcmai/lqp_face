/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/
#ifndef LTPFEATURESRGB_H_
#define LTPFEATURESRGB_H_
#include "lbpFeaturesRGB.h"
class LTPFeaturesRGB: public LBPFeaturesRGB {
public:
	LTPFeaturesRGB(const LBPParams &lbpparam, UINT width_, UINT height_,
			UINT nplevels, ProcessImage &pim_) :
		LBPFeaturesRGB(lbpparam, width_, height_, nplevels, pim_) {
		tolerance = lbpparam.tolerance;
		cout << " LTPDim = " << (IsFolded() == true ? GetFoldedDim(width_,
				height_) : GetDim(width_, height_)) << endl;
		punicodes = pnunicodes = 0;
		nunicodes = nnunicodes = 0;
	}
	LTPFeaturesRGB(UINT cellsize_, UINT npoints_, UINT radius_, UINT width_,
			UINT height_, UINT nplevels, NORM norm_, UINT nflevels_,
			UINT lbpstride_, ProcessImage &pim_, bool folded_, bool add_,
			UINT nbins_, bool nseparate_, REAL tol, GridType gridtype_) :
		LBPFeaturesRGB(cellsize_, npoints_, radius_, width_, height_, nplevels,
				norm_, nflevels_, lbpstride_, pim_, folded_, add_, nbins_,
				nseparate_, gridtype_) {
		tolerance = tol;
		cout << " LTPDim = " << (IsFolded() == true ? GetFoldedDim(width_,
				height_) : GetDim(width_, height_)) << endl;
		punicodes = pnunicodes = 0;
		nunicodes = nnunicodes = 0;
	}
	virtual UINT GetDim(UINT width_, UINT height_) const {
		return nflevels == 1 ? (width_ / lbpstride) * (height_ / lbpstride)
				* nbins * 3 * dimmult : (((width_ / lbpstride) * (height_
				/ lbpstride)) + ((width_ / (2 * lbpstride)) * (height_ / (2
				* lbpstride)))) * nbins * 3 * dimmult;
	}
	virtual UINT GetFoldedDim(UINT width_, UINT height_) const { // for horizontal symmetry
		return nflevels == 1 ? ceil((REAL) width_ / (2 * lbpstride)) * (height_
				/ lbpstride) * 3 * nbins * dimmult : (((width_
				/ (2 * lbpstride)) * (height_ / lbpstride)) + (width_ / (4
				* lbpstride) * height_ / (2 * lbpstride))) * 3 * nbins
				* dimmult;
	}
	virtual void InitalizeMaps(Image &image);
	virtual void InitalizeMaps(Pyramid &pyobj_, PyramidType sspace);
	virtual void InitalizeMaps(Image &, PyramidType);
	void ComputeLTPMap(Image &image, LBPMap &lbpmap, LBPMap& tpmap,
			LBPMap &tnmap);
	void ComputeLTPMap(Image &image, LBPMap &lbpmap, LBPMap & flbmap,
			LBPMap& tpmap, LBPMap& ftpmap, LBPMap &tnmap, LBPMap &ftnmap);
	void ComputeLTPMap(Image &image, vector<LBPMap> &, vector<LBPMap> &,
			vector<LBPMap> &);
	void ComputeLTPMap(Image &image, LBPMap[], LBPMap[], LBPMap[]);
	void ComputeLTPMapFast(Image &image, LBPMap[], LBPMap[], LBPMap[]);
	void ComputeLTPMap(vector<REAL>&impix, Point<UINT>&, LBPMap &lbpmap,
			LBPMap& tpmap, LBPMap &tnmap);// for lbp and ltp over the edge images....
	virtual ~LTPFeaturesRGB() {
	}
	virtual void GetFeatures(UINT, int, int, vector<REAL>&);
	virtual void GetFoldedFeatures(UINT, int, int, vector<REAL>&);
	virtual void UnFoldWeights(REAL*, vector<REAL>&);
	virtual UINT GetXSize(UINT index) {
		return lbpmaps[index][0].nxpoints + radius;
	}
	virtual UINT GetYSize(UINT index) {
		return lbpmaps[index][0].nypoints + radius;
	}
	virtual UINT GetNumberBins() {
		return nbins * 3;
	}
	virtual UINT PrintCodesInfo() {
		cout << endl << " Positive Codes: uniform =" << punicodes
				<< " , non-uniform =" << pnunicodes << endl
				<< " Negative Codes: uniform =" << nunicodes
				<< " , non-uniform =" << nnunicodes << endl;
	}
protected:
	vector<vector<LBPMap> > posmap;
	vector<vector<LBPMap> > negmap;
	vector<LBPMap> fnegmap;
	vector<LBPMap> fposmap;
	REAL tolerance;
	HOGInfo ltppos, ltpneg;
	UINT punicodes, nunicodes; // counts of uniform codes in the +ve and -ve maps...
	UINT pnunicodes, nnunicodes; // counts of uniform codes in the +ve and -ve maps...
};

#endif /* HOGLBPFEATURES_H_ */
