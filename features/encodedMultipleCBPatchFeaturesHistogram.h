/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/

#ifndef ENCODEDMULTIPLECBPATCHFEATURESHISTOGRAM_H_
#define ENCODEDMULTIPLECBPATCHFEATURESHISTOGRAM_H_
#include"encodedMultipleCBPatchFeatures.h"
// Histogram of all the features at Window Level, Twice Normalization once at cell level at 2nd at Window Level..
// Current Implementation only Supports Fix Dimension Windows...
// Two Possibilties..
// First To Encode each cell against all the Cell CodeBooks and Negative CodeBooks and then to build histogram
// over codebooks for the detection window.
// Second To Encode each cell against its respective codeBook just like Pyramid style approach, and -ve
// code book.. Now the possiblties are whether to keep the neg-cell code book as cell level or histogram all
// the negative codebook contributions into a histogram at window level features...
#define LIVE_COMPUTE // set int make file this flag
#ifndef LIVE_COMPUTE // works on 64bit systems with more than 4gb of ram
class EncodedMultipleCBPatchFeaturesHistogram: public EncodedMultipleCBPatchFeatures {
public:
	EncodedMultipleCBPatchFeaturesHistogram(const LBPParams &lbpparam,
			UINT width_, UINT height_, UINT nplevels, ProcessImage &pim_,
			const PatchParams& pparam) :
				EncodedMultipleCBPatchFeatures(lbpparam, width_, height_,
						nplevels, pim_, pparam) {
		ncbs = cwidth * cheight + 1; // code books for each cell...
		featdim = ncenters * ncbs; // featdim to be used internally
		// windim = dimension of a window with size cwidth, cheight...
		if (cbtype == CBT_PosCellNegSep) {
			windim = (ncenters * 2 * (ncbs - 1)/**/);
		} else
			windim = featdim; // to quantize each cell against its own code book and negative code book
		hnormtype = pparam.hnormtype;
		tmpfeat.resize(featdim * cwidth * cheight, 0);
	}
	virtual ~EncodedMultipleCBPatchFeaturesHistogram() {

	}
	void BuildWindowHistogram(REAL *feat);
	void ComputeLBPFeatures(Image& imgref, UINT index, UINT winstride,
			UINT tlbpstride);
	virtual void GetFeatures(UINT index, int x, int y, int width_, int height_,
			vector<REAL>&feat);
	virtual void GetFeatures(UINT index, int x, int y, int width_, int height_,
			REAL *feat);
	virtual void GetFoldedFeatures(UINT index, int x, int y, int width_,
			int height_, vector<REAL>&feat);
	virtual void GetFoldedFeatures(UINT index, int x, int y, int width_,
			int height_, REAL *feat);
	virtual void DotProduct(UINT index, UINT fwidth, UINT fheight,
			vector<REAL>&filter, Store&response);
	virtual void DotProduct(UINT index, UINT fwidth, UINT fheight,
			REAL *filter, Store&response);

	virtual UINT GetDim(UINT width_, UINT height_) const {
		return GetFeatDim(width_, height_);
	}
	virtual UINT GetHOGDim(UINT width_, UINT height_) const {
		return 0;
	}
	virtual UINT GetLBPDim(UINT width_, UINT height_) const {
		return GetFeatDim(width_, height_);;
	}
	virtual UINT GetFoldedLBPDim(UINT width_, UINT height_) const {
		return GetFeatDim(width_, height_);;
	}
	virtual UINT GetFoldedHOGDim(UINT width_, UINT height_) const {
		return 0;
	}
	virtual UINT GetFoldedDim(UINT width_, UINT height_) const {
		return GetFeatDim(width_, height_);
	}
private:

	UINT GetFeatDim(UINT width_, UINT height_) const {
		return /*(width_ / cellsize) * (height_ / cellsize) * celldim*/windim;
	}
protected:
	UINT windim; // used for external feature dimension calculations...
	EncodedHistogramNormType hnormtype;
	vector<REAL> tmpfeat; // used during computation of features...
};
#else

class EncodedMultipleCBPatchFeaturesHistogram: public EncodedMultipleCBPatchFeatures {
public:
	EncodedMultipleCBPatchFeaturesHistogram(const LBPParams &lbpparam,
			UINT width_, UINT height_, UINT nplevels, ProcessImage &pim_,
			const PatchParams& pparam) :
	EncodedMultipleCBPatchFeatures(lbpparam, width_, height_, nplevels,
			pim_, pparam) {
		ncbs = cwidth * cheight + 1; // code books for each cell...
		featdim = ncenters * ncbs; // featdim to be used internally
		if (cbtype == CBT_PosCellNegSep) {
			windim = (ncenters * 2 * (ncbs - 1)/**/);
		} else
		windim = featdim;
		hnormtype = pparam.hnormtype;
		tmpfeat.resize(featdim * cwidth * cheight, 0);
	}
	virtual ~EncodedMultipleCBPatchFeaturesHistogram() {

	}
	void BuildWindowHistogram(REAL *feat);
	void ComputeLBPFeatures(Image& imgref, UINT index, UINT winstride,
			UINT tlbpstride);
	virtual void GetFeatures(UINT index, int x, int y, int width_, int height_,
			vector<REAL>&feat);
	virtual void GetFeatures(UINT index, int x, int y, int width_, int height_,
			REAL *feat);
	virtual void GetFoldedFeatures(UINT index, int x, int y, int width_,
			int height_, vector<REAL>&feat);
	virtual void GetFoldedFeatures(UINT index, int x, int y, int width_,
			int height_, REAL *feat);
	virtual void DotProduct(UINT index, UINT fwidth, UINT fheight,
			vector<REAL>&filter, Store&response);
	virtual void DotProduct(UINT index, UINT fwidth, UINT fheight,
			REAL *filter, Store&response);

	virtual void PadFeatureMap(UINT index, UINT padx, UINT pady);

	virtual UINT GetDim(UINT width_, UINT height_) const {
		return GetFeatDim(width_, height_);
	}
	virtual UINT GetHOGDim(UINT width_, UINT height_) const {
		return 0;
	}
	virtual UINT GetLBPDim(UINT width_, UINT height_) const {
		return GetFeatDim(width_, height_);;
	}
	virtual UINT GetFoldedLBPDim(UINT width_, UINT height_) const {
		return GetFeatDim(width_, height_);;
	}
	virtual UINT GetFoldedHOGDim(UINT width_, UINT height_) const {
		return 0;
	}
	virtual UINT GetFoldedDim(UINT width_, UINT height_) const {
		return GetFeatDim(width_, height_);
	}
	virtual void InitalizeMaps(Image &, PyramidType);
	virtual void InitalizeMaps(Image &image);
	virtual void InitalizeMaps(Pyramid &pyobj, PyramidType sspace_);
	virtual UINT GetXSize(UINT index) {
		return blockinfo[index].GetX() * cell;
	}
	virtual UINT GetYSize(UINT index) {
		return blockinfo[index].GetY() * cell;
	}
	virtual REAL GetScale(UINT index) {
		return scaleinfo[index];
	}
	void CheckIsComputed(const UINT &index) {
		if (computedindex != index && index < pyobj.GetNLevels()) {
			Image & imgref = pyobj.GetImageRef(index);
			lbpfeatures.Initialize(1, 0, featdim);
			UINT xbmax, ybmax, padx, pady;
			GetXYBlocks(imgref, xbmax, ybmax);
			lbpfeatures.SetFeature(xbmax, ybmax, pyobj.GetScale(index));
			ComputeLBPFeatures(imgref, 0, cell, lbpstride);
			computedindex = index;
			padinfo[index].GetCoord(padx, pady);
			if (padx != 0 && pady != 0) {
				lbpfeatures.PadFeatureMap(0, padx, pady);
				xbmax = lbpfeatures.GetXBlock(0);
				ybmax = lbpfeatures.GetYBlock(0);
				blockinfo[index].SetCoord(xbmax, ybmax);
			}
		}
	}
private:

	UINT GetFeatDim(UINT width_, UINT height_) const {
		return /*(width_ / cellsize) * (height_ / cellsize) * celldim*/windim;
	}
protected:
	UINT windim; // used for external feature dimension calculations...
	EncodedHistogramNormType hnormtype;
	vector<REAL> tmpfeat; // used during computation of features...
	vector<Point<UINT> > padinfo;
	vector<Point<UINT> > blockinfo;
	vector<REAL> scaleinfo;
	int computedindex;
};
#endif
#endif /* ENCODEDMULTIPLECBPATCHFEATURESHISTOGRAM_H_ */
