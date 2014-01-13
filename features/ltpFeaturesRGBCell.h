/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/

#ifndef LTPFEATURESRGBCELL_H_
#define LTPFEATURESRGBCELL_H_
#include "ltpFeaturesRGB.h"
#include "hogInfo.h"
//#include "util.hpp"
// xcell size of a cell;
// xblock size of block in cell coordinates 2 means each block contain 2cell
// xstride difference between two block in number of pixels....
// e.g for non-overlapping HOGs xstride = xcell*xblock;
// Its symmetric for y coordinates...
class LTPFeaturesRGBCell: public LTPFeaturesRGB {
public:
	LTPFeaturesRGBCell(const LBPParams &lbpparam, UINT width_, UINT height_,
			UINT nplevels, ProcessImage &pim_) :
		LTPFeaturesRGB(lbpparam, width_, height_, nplevels, pim_) {
		cell = lbpparam.cellsize;
		cellsize = cell;
		cwidth = width / cell;
		cheight = height / cell;
		lbpdim = GetLBPDim(width_, height_);
		UINT flbpdim = GetFoldedLBPDim(width_, height_);
		foffset = cell; // offset by One cell
		cout << "Cell Size = " << cell << " LBP Dim = " << lbpdim
				<< " Folded LBP Dim = " << flbpdim << endl
				<< " Original Feature Dim = " << GetDim(width_, height_)
				<< " Folded Feature Dim = " << GetFoldedDim(width_, height_)
				<< " Feature Offset =" << foffset;

		skip = lbpstride / cell; // how many cells to skip
	}
	LTPFeaturesRGBCell(UINT cellsize_, UINT npoints_, UINT radius_,
			UINT width_, UINT height_, UINT nplevels, NORM norm_,
			UINT nflevels_, UINT lbpstride_, ProcessImage &pim_, bool folded_,
			bool add_, UINT nbins_, bool nseparate_, REAL tol,
			GridType gridtype_) :
				LTPFeaturesRGB(cellsize_, npoints_, radius_, width_, height_,
						nplevels, norm_, nflevels_, lbpstride_, pim_, folded_,
						add_, nbins_, nseparate_, tol, gridtype_) {
		cell = 8;
		cwidth = width / cell;
		cheight = height / cell;
		lbpdim = GetLBPDim(width_, height_);
		UINT flbpdim = GetFoldedLBPDim(width_, height_);
		cout << " LBP Dim = " << lbpdim << " Folded LBP Dim = " << flbpdim
				<< endl << " Original Feature Dim = "
				<< GetDim(width_, height_) << " Folded Feature Dim = "
				<< GetFoldedDim(width_, height_);
		foffset = 8; // offset by One cell

		skip = lbpstride / cell; // how many cells to skip
	}
	// GetDim returns the size of feature in a window without considering offset...
	virtual UINT GetDim(UINT width_, UINT height_) const {
		return (GetLBPDim(width_, height_));
	}

	virtual UINT GetHOGDim(UINT width_, UINT height_) const {
		return 0;
	}
	virtual UINT GetLBPDim(UINT width_, UINT height_) const {
		return (width_ / cellsize) * (height_ / cellsize) * nbins * 3 * dimmult;
	}
	virtual UINT GetFoldedLBPDim(UINT width_, UINT height_) const {
		return ceil((REAL) width_ / (2 * cellsize)) * (height_ / cellsize)
				* nbins * 3 * dimmult;
	}
	virtual UINT GetFoldedHOGDim(UINT width_, UINT height_) const {
		return 0;
	}
	virtual UINT GetFoldedDim(UINT width_, UINT height_) const {
		return GetFoldedLBPDim(width_, height_);
	}
	virtual ~LTPFeaturesRGBCell() {
	}
	virtual void ComputeTextureFeatures(Image &imgref, UINT xcellsize,
			UINT ycellsize, vector<REAL>&features);
	virtual void ComputeTextureFeatures(Image &imgref, int sx, int sy, int ex,
			int ey, UINT xcellsize, UINT ycellsize, vector<REAL>&features);
	virtual void InitalizeMaps(Image &, PyramidType);
	virtual void InitalizeMaps(Pyramid & pyobj_, PyramidType);
	virtual void InitalizeMaps(Image &image);

	virtual void GetFeatures(UINT, int, int, int, int, vector<REAL>&);
	virtual void GetFeatures(UINT, int, int, int, int, REAL*);
	virtual void GetFoldedFeatures(UINT, int, int, int, int, vector<REAL>&);
	virtual void GetFoldedFeatures(UINT, int, int, int, int, REAL*);
	virtual void
	GetFlippedFeatures(UINT, int, int, int, int, vector<REAL>&);
	virtual void UnFoldWeights(REAL*, UINT, UINT, vector<REAL>&);

	virtual UINT GetXSize(UINT index) {
		return lbpfeatures.GetXBlock(index) * cell;
	}
	virtual UINT GetYSize(UINT index) {
		return lbpfeatures.GetYBlock(index) * cell;
	}
	virtual REAL GetScale(UINT index) {
		return lbpfeatures.GetScale(index);
	}
	virtual UINT GetInitIndex() {
		return sspace == DoubleRes ? pyobj.GetInitIndex() : 0;
	}
	virtual void PadFeatureMap(UINT index, UINT padx, UINT pady);
	void WriteMap() {
		ofstream ofile("maps.txt");
		for (UINT i = 0; i < lbpmaps.size(); ++i) {
			for (int j = 0; j < 3; ++j) {
				lbpmaps[i][j].WriteMap(ofile);
				posmap[i][j].WriteMap(ofile);
				negmap[i][j].WriteMap(ofile);
			}
		}
	}
	virtual void DotProduct(UINT index, UINT width, UINT height, vector<REAL>&,
			Store&response);
	virtual void DotProduct(UINT index, UINT width, UINT height, REAL *filter,
			Store&);
	virtual void DotProductWithQuadratic(UINT index, UINT width, UINT height,
			REAL, vector<REAL>&, vector<REAL>&, Store&) {
	}
	virtual void DotProductWithProjection(UINT index, UINT xsize, UINT ysize,
			const VectorXf &sigma, const MatrixXf &projMat, vector<REAL>&,
			Store&inhog) const {

	}
	void ComputeLBPFeatures(Image& imgref, UINT index, UINT winstride,
			UINT tlbpstride);
	/*****Test Code ****/
	// All the inputs are in Pixels .....
	virtual void GetHOGFeatures(UINT index, int x, int y, int width_,
			int height_, REAL*feat) { // There is offset of one cell so, input x=0 means x=8 by taking
		//	into account offset and the removed border
		//		hoginfo.GetHOGCells(index, x / cell, y / cell, width_ / cell, height_
		//				/ cell, feat);
	}

	virtual void GetFoldedHOGFeatures(UINT index, int x, int y, int twidth,
			int theight, REAL*feat) { // There is offset of one cell so, input x=0 means x=8 by taking
		//	into account offset and the removed border
		//		hoginfo.GetFoldedHOGCells(index, x / cell, y / cell, twidth / cell,
		//				theight / cell, feat);
	}

	virtual void GetFlippedHOGFeatures(UINT index, int x, int y, int twidth,
			int theight, vector<REAL>&feat) { // There is offset of one cell so, input x=0 means x=8 by taking
		//	into account offset and the removed border
		//		hoginfo.GetFlippedFeatures(index, x / cell, y / cell, twidth / cell,
		//				theight / cell, &feat[0]);
	}
	virtual void HOGDotProduct(UINT index, UINT fwidth, UINT fheight,
			vector<REAL>&filter, Store&response) {
		//		hoginfo.DotProduct(index, fwidth / cell, fheight / cell, &filter[0],
		//				response);
	}
	virtual void LBPDotProduct(UINT index, UINT fwidth, UINT fheight,
			vector<REAL>&filter, Store&response) {
		lbpfeatures.SkippedDotProduct(index, fwidth / cell, fheight / cell,
				skip, &filter[0], response);
	}
	virtual void GetLBPFeatures(UINT index, int x, int y, int w, int h,
			REAL* feat) {
		lbpfeatures.GetSkippedHOGCells(index, x / cell, y / cell, w / cell,
				h / cell, skip, feat);
	}
	void GetLBPRGB(Image&image, vector<REAL>&feat, vector<REAL>&ffeat) {

		UINT fdim = LBPFeaturesRGB::GetDim(64, 160);
		feat.resize(fdim, 0);
		ffeat.resize(fdim, 0);
		LBPFeaturesRGB::InitalizeMaps(image, NoPyramid);
		LBPFeaturesRGB::GetFeatures(0, 8, 8, feat);
		image.flop();
		LBPFeaturesRGB::InitalizeMaps(image, NoPyramid);
		LBPFeaturesRGB::GetFeatures(0, 8, 8, ffeat);

	}
	//*****/
	virtual void GetFlippedFeatures(int twidth, int theight, REAL *ifeat,
			REAL *ofeat);
	void ComputeLBPFeatures(Image &image, UINT index, UINT sbin);
	virtual UINT PrintCodesInfo() {
		cout << endl << " Positive Codes: uniform =" << punicodes
				<< " , non-uniform =" << pnunicodes << endl
				<< " Negative Codes: uniform =" << nunicodes
				<< " , non-uniform =" << nnunicodes << endl;
	}

protected:
	UINT cell, cwidth, cheight, lbpdim, skip;

	void Compute(const Image&, UINT, REAL*);
	void GetXYBlocks(const Image&image, UINT&xbmax, UINT&ybmax) { // look in the code of Compute as the last two boundary blocks are not used
		xbmax = (UINT) (round((REAL) image.columns() / (REAL) cell) - 2);
		ybmax = (UINT(round((REAL) image.rows() / (REAL) cell) - 2));
	}

	void Get2XYBlocks(const Image&image, UINT&xbmax, UINT&ybmax) { // look in the code of Compute as the last two boundary blocks are not used

		xbmax = (UINT) (round((REAL) (image.columns() * 2) / (REAL) cell) - 2);
		ybmax = (UINT) (round((REAL) (image.rows() * 2) / (REAL) cell) - 2);
	}

	void GetLBPFeatures(UINT, UINT, UINT, UINT, UINT, REAL*); // scale index with starting location and size also size

	void ComputeHistDiscrete(Image& imgref, UINT index, UINT winstride,
			LBPMap *tlbpmap, LBPMap *tposmap, LBPMap *tnegmap);
	void ComputeHistBilinear(Image &image, UINT sbin, LBPMap *tlbpmap,
			LBPMap *tposmap, LBPMap *tnegmap, REAL *feat);
};

#endif
