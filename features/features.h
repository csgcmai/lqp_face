/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/

#ifndef FEATURES_H_
#define FEATURES_H_

#include "../util/npyramid.h"
#include "../util/point.hpp"
#include "../util/util.hpp"
#include "../util/cache.h"
//#include "../Kmeans/cluster.hpp"
#include "../util/classInfo.h"
#include "../util/windowInfoFeat.h"
#include "lbpMap.h"
#include "definitions.h"

class Features {
public:
	Features(UINT w, UINT h, UINT nlevels, bool folded_ = false) :
		width(w), height(h), pyobj(width, height, nlevels, Bilinear) {
		sspace = NoPyramid;

		folded = folded_;
		foffset = 0;
	}
	virtual ~Features() {
	}
	virtual UINT GetDim(UINT width_, UINT height_) const=0;
	virtual UINT GetFoldedDim(UINT width_, UINT height_) const =0;

	//	virtual void InitalizeMaps(Image &image) {
	//		PrintSpecializationError(
	//				"(InitalizeMaps) Must define a specialization for this object");
	//	}


	//	virtual void InitalizeMaps(Pyramid &pyobj_, PyramidType sspace) { // to use with multi - Levels LBP's LTP's features...
	//		PrintSpecializationError(
	//				" Error: (InitalizeMaps) Must define a specialization for this object");
	//	}
	virtual void InitalizeMaps(vector<Image> &image)/*Required for Gabor-Filtered Images...*/{
		PrintSpecializationError(
				" Error: (InitalizeMaps(vector<Image> &image )) Must define a specialization for this object");
	}
	virtual void InitalizeMaps(Image &image, PyramidType sspace)=0;

	virtual void GetFeatures(UINT index, int x, int y, vector<REAL>& feat) {
		PrintSpecializationError(
				" GetFeatures Specialization needs to be defined");
	}
	virtual void GetFoldedFeatures(UINT index, int x, int y, vector<REAL>& feat) {
		PrintSpecializationError(
				" GetFoldedFeatures Specialization needs to be defined");
	}
	void SetPyramidInterval(UINT nlevels) {
		pyobj.SetOctaveSize(nlevels);
	}
	UINT GetOctaveSize() {
		return pyobj.GetOctaveSize();
	}
	UINT GetNLevels() {
		return sspace == NoPyramid ? 1
				: sspace == SingleRes ? pyobj.GetNLevels() : pyobj.GetNLevels()
						+ pyobj.GetInitIndex();
	}

	virtual UINT GetXSize(UINT index)=0;
	virtual UINT GetYSize(UINT index)=0;
	bool IsFolded() {
		return folded;
	}
	UINT GetFeaturesOffset() {
		return foffset;
	}

#ifdef WINDOW_CLASSIFIER
	virtual REAL GetScale(UINT index) {
		return sspace == NoPyramid ? 1 : sspace == SingleRes ? pyobj.GetScale(
				index) : 1;
	}
	virtual void UnFoldWeights(REAL*iweight, vector<REAL>&oweight) {
		cout << " Error: Specific Definition Missing " << endl;
	}

#else
	virtual UINT GetInitIndex()=0;
	virtual REAL GetScale(UINT index)=0;
	virtual void PadFeatureMap(UINT index, UINT padx, UINT pady)=0; //	{ // does padding in the pixel space	}
	virtual void DotProduct(UINT index, UINT width, UINT height, vector<REAL>&,
			Store&)=0;//{}
	virtual void DotProduct(UINT index, UINT width, UINT height, REAL *filter,
			Store&) {
		PrintSpecializationError(
				" SpecializationError: DotProduct(Filter*) Define the Specialization of the Function ");
	}
	//	virtual void DotProductWithQuadratic(UINT index, UINT width, UINT height,
	//			REAL, vector<REAL>&, vector<REAL>&, Store&)=0;
	//	virtual void DotProductWithProjection(UINT index, UINT xsize, UINT ysize,
	//			const VectorXf &sigma, const MatrixXf &projMat, vector<REAL>&,
	//			Store&inhog) const=0;

	virtual UINT GetHOGDim(UINT width, UINT height) const=0;//{	return 0;	}
	virtual UINT GetFoldedHOGDim(UINT width, UINT height) const=0;//{	return 0;}
	virtual UINT GetLBPDim(UINT width, UINT height) const=0;//{return 0;}
	virtual UINT GetFoldedLBPDim(UINT width, UINT height) const=0;//{return 0;}

	virtual void GetFeatures(UINT index, int x, int y, int w, int h,
			vector<REAL>& feat)=0; //{}
	virtual void
	GetFeatures(UINT index, int x, int y, int w, int h, REAL* feat)=0;
	virtual void GetFoldedFeatures(UINT index, int x, int y, int w, int h,
			vector<REAL>& feat) =0;
	virtual void GetFoldedFeatures(UINT index, int x, int y, int w, int h,
			REAL* feat) =0;//{}

	virtual void GetFlippedFeatures(UINT index, int x, int y, int w, int h,
			vector<REAL>& feat)=0;// {}
	virtual void GetFlippedFeatures(int twidth, int theight, REAL *ifeat,
			REAL *ofeat)=0;
#endif

	/* Generally features are computed once by calling the function InitializeMaps(...) first and then calling the GetFeatures()
	 * However, ComputeTextureFeatures is a special case and can be used to compute features over complete image
	 * (e.g by setting  xcellsize=image.width() and ycellsize=image.height() ) as well
	 * as grid-based cell level features. It is almost a stand-alone function.
	 */
	virtual void ComputeTextureFeatures(Image &imgref, UINT xcellsize,
			UINT ycellsize, vector<REAL>&features) {
		PrintSpecializationError("ComputeTextureFeatures_NoSpecialization");
	}
	virtual void ComputeTextureFeatures(Image &imgref, UINT xcellsize,
			UINT ycellsize, int sx, int sy, int ex, int ey,
			vector<REAL>&features) {
		PrintSpecializationError("ComputeTextureFeatures_NoSpecialization");
	}

	void PrintSpecializationError(const string& err) const {
		cerr << endl << err << endl;
		exit(EXIT_FAILURE);
	}

	virtual UINT GetMaxFeatureOffset() { /*To be used for initial Cropping of Windows....*/
		return foffset;
	}
	// Generate the list of all windows in a scale space
	UINT GenerateWindows(Image &image, vector<WindowInfoFeat>& winfo) {
		sspace = SingleRes;
		UINT pycount = pyobj.GeneratePyramid(image);
		UINT wCount = 0;
		UINT ybmax, xbmax, ystride = ycell, xstride = xcell;
		REAL response = 0, overlap = 0, scale;
		for (UINT i = 0; i < pycount; i++) {

			wCount = 0;
			ybmax = GetYSize(i);
			xbmax = GetXSize(i);
			if (ybmax >= height && xbmax >= width) {
				int trwindows = (int) floor((ybmax - height) / ystride) + 1,
						tcwindows = (int) floor((xbmax - width) / xstride) + 1;
				scale = GetScale(i);
				for (UINT rcount = 0; rcount < trwindows * ystride; rcount
						+= ystride)
					for (UINT ccount = 0; ccount < tcwindows * xstride; ccount
							+= xstride, ++wCount) {
						winfo.push_back(
								WindowInfoFeat(ccount, rcount, width, height,
										i, response, scale, 0, overlap, 0));
					}
			}
		}
		return wCount;
	}
	virtual UINT GetNumberBins() {
		return nbins;
	}
	virtual UINT PrintCodesInfo() {
		PrintSpecializationError("GenerateCodeBook_NoSpecialization");
	}
	virtual void GenerateCodeBook() {
		PrintSpecializationError("GenerateCodeBook_NoSpecialization");
	}

	static UINT foffset, xcell, ycell, norient, cellsize;
protected:
	UINT width, height;
	UINT nbins;
	bool folded;
	PyramidType sspace;
	Pyramid pyobj;
};
#endif /* FEATURES_H_ */
