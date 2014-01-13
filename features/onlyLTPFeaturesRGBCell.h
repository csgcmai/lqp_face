/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/

#ifndef ONLYLTPFEATURESRGBCELL_H_
#define ONLYLTPFEATURESRGBCELL_H_
#include "ltpFeaturesRGBCell.h"
#include "ltpFeaturesRGB.h"

//#include "encodeFeatures.h"
#include "hogInfo.h"
//#include "util.hpp"
// Its symmetric for y coordinates...
// Modifications to include the multiLevel LBP,LTP and Multi-Level Thresholds....
// Because as the threshold is varied the detail encoded by LTP varies...
// Threshold Level defines number of threshold used for the computation of LTPs....
// Currently Using Power of 2  as base power of 2,4,8,16,32
// 2 is same as having local texture features, 32 as having the most coarsest information...
// Important Note: If number o threshold levels =1 then tol value is used as it is,
// otherwise the tol value is used as base for the computation of threshold rings...

class OnlyLTPFeaturesRGBCell: public LTPFeaturesRGBCell {
public:
	OnlyLTPFeaturesRGBCell(const LBPParams &lbpparm, UINT width_, UINT height_,
			UINT nplevels, ProcessImage &pim_) :
		LTPFeaturesRGBCell(lbpparm, width_, height_, nplevels, pim_) {
		useblocknorm = lbpparm.useblocknorm;
		//		dimmult = useblocknorm == true ? dimmult * 4 : dimmult * 1;
		clbpdim = nbins * (useblocknorm == true ? 4 : 1);
		celldim = 2 * clbpdim;

		lbpdim = GetLBPDim(width_, height_);
		UINT flbpdim = GetFoldedLBPDim(width_, height_);
		cout << "\n LBP Dim = " << lbpdim << " Folded LBP Dim = " << flbpdim
				<< endl << " Original Feature Dim = "
				<< GetDim(width_, height_) << " Folded Feature Dim = " << endl
				<< GetFoldedDim(width_, height_) << " UseBnorm ="
				<< useblocknorm << " CellDim =" << celldim << " Dimmult ="
				<< dimmult;
	}
	OnlyLTPFeaturesRGBCell(UINT cellsize_, UINT npoints_, UINT radius_,
			UINT width_, UINT height_, UINT nplevels, NORM norm_,
			UINT nflevels_, UINT lbpstride_, ProcessImage &pim_, bool folded_,
			bool add_, UINT nbins_, bool nseparate_, REAL tol,
			GridType gridtype_, bool usebnorm = false) :
				LTPFeaturesRGBCell(cellsize_, npoints_, radius_, width_,
						height_, nplevels, norm_, nflevels_, lbpstride_, pim_,
						folded_, add_, nbins_, nseparate_, tol, gridtype_) {
		useblocknorm = usebnorm;
		//		dimmult = useblocknorm == true ? dimmult * 4 : dimmult * 1;
		clbpdim = nbins * (useblocknorm == true ? 4 : 1);
		celldim = 2 * clbpdim;
		lbpdim = GetLBPDim(width_, height_);
		UINT flbpdim = GetFoldedLBPDim(width_, height_);
		cout << "\n LBP Dim = " << lbpdim << " Folded LBP Dim = " << flbpdim
				<< endl << " Original Feature Dim = "
				<< GetDim(width_, height_) << " Folded Feature Dim = "
				<< GetFoldedDim(width_, height_);
	}
	// GetDim returns the size of feature in a window without considering offset...
	virtual UINT GetDim(UINT width_, UINT height_) const {
		return (GetLBPDim(width_, height_));
	}

	virtual UINT GetHOGDim(UINT width_, UINT height_) const {
		return 0;
	}
	virtual UINT GetLBPDim(UINT width_, UINT height_) const {
		return (width_ / cellsize) * (height_ / cellsize) * celldim * dimmult;
	}
	virtual UINT GetFoldedLBPDim(UINT width_, UINT height_) const {
		return ceil((REAL) width_ / (2 * cellsize)) * (height_ / cellsize)
				* celldim * dimmult;
	}
	virtual UINT GetFoldedHOGDim(UINT width_, UINT height_) const {
		return 0;
	}
	virtual UINT GetFoldedDim(UINT width_, UINT height_) const {
		return GetFoldedLBPDim(width_, height_);
	}
	virtual ~OnlyLTPFeaturesRGBCell() {
	}
	virtual void ComputeLTPMap(Image &image, LBPMap tpmap[], LBPMap tnmap[]);
	void ComputeLBPFeatures(Image& imgref, UINT index, UINT winstride,
			UINT tlbpstride);
	void ComputeLBPFeatures(Image &image, UINT index, UINT sbin);

	virtual void GetFlippedFeatures(int twidth, int theight, REAL *ifeat,
			REAL *ofeat);
	virtual void UnFoldWeights(REAL *input, UINT w, UINT h,
			vector<REAL>& weights);
	virtual void ComputeTextureFeatures(Image &imgref, UINT xcellsize,
			UINT ycellsize, vector<REAL>&features);
	virtual void ComputeTextureFeatures(Image &imgref, int sx, int sy, int ex,
			int ey, UINT xcellsize, UINT ycellsize, vector<REAL>&features);
	virtual void InitalizeMaps(Image &, PyramidType);
	virtual void InitalizeMaps(Pyramid & pyobj_, PyramidType);
	virtual void InitalizeMaps(Image &image);
	virtual UINT GetNumberBins() {
		return 2 * clbpdim;
	}
	virtual UINT PrintCodesInfo() {
		cout << endl << " Positive Codes: uniform =" << punicodes
				<< " , non-uniform =" << pnunicodes << endl
				<< " Negative Codes: uniform =" << nunicodes
				<< " , non-uniform =" << nnunicodes << endl;
	}
protected:
	UINT celldim, clbpdim;
	bool useblocknorm; // use blocks normalization.
	void ComputeHistDiscrete(Image& imgref, UINT index, UINT winstride,
			LBPMap *tposmap, LBPMap *tnegmap);
	void ComputeHistBilinear(Image &image, UINT sbin, LBPMap *tposmap,
			LBPMap *tnegmap, REAL *feat);
};

#endif
