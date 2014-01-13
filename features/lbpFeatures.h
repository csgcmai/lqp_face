/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/
#ifndef LBPFEATURES_H
#define LBPFEATURES_H
#include "features.h"
#include "processImage.h"
// xcell size of a cell;
class LBPFeatures: public Features {
public:
	LBPFeatures(const LBPParams &lbpparam, UINT width_, UINT height_,
			UINT nplevels, ProcessImage &pim_) :
		Features(width_, height_, nplevels, false), pim(pim_) {
		//		 npoints(npoints_), unif(2),
		//						radius(radius_), lbpstride(lbpstride_), pim(pim_)
		npoints = lbpparam.npoints;
		unif = 2;
		radius = lbpparam.radius;
		lbpstride = lbpparam.stride;
		cellsize = lbpparam.cellsize;
		gridtype = lbpparam.gridtype;
		ptype = lbpparam.ptype;
		cout << "\n LBP Grid Type = " << (gridtype == Circular ? "Circular"
				: "Rectangular") << flush << endl;
		nbins = lbpparam.nbins != 256 ? npoints * (npoints - 1) + 3
				: lbpparam.nbins;
		GenerateMap();
		cout << endl;
		spoints.resize(npoints);
		norm = lbpparam.featnorm;
		nflevels = lbpparam.nflevels;
		width = width_;
		height = height_;
		ComputeSamplingPoints();
		nycells = height / cellsize;
		nxcells = width / cellsize;
		usegradmag = lbpparam.usegradmag;
		hmethod = lbpparam.hmethod;
		lbpparam.PrintParamInfo();
		cout << endl << " ------------------------------------" << endl
				<< " Width = " << width << " Height = " << height
				<< " LBPStride= " << lbpstride << endl
				<< " Feature Normalization Type = " << norm << endl
				<< " LBPbins = " << nbins << endl
				<< " Geometrical Organization =  "
				<< PatchParams::GetPatchType(ptype) << "\t  Histogram Method= "
				<< LBPParams::GetHistogramMethod(hmethod) << endl << flush;
		cout << " Nxcells = " << nxcells << " Nycells = " << nycells << flush;
		cout << " LBPDim = " << (IsFolded() == true ? GetFoldedDim(width_,
				height_) : GetDim(width_, height_)) << endl;

	}
	LBPFeatures(UINT cellsize_, UINT npoints_, UINT radius_, UINT width_,
			UINT height_, UINT nplevels, NORM norm_, UINT nflevels_,
			UINT lbpstride_, ProcessImage &pim_, bool folded_ = false) :
		Features(width_, height_, nplevels, folded_), npoints(npoints_),
				unif(2), radius(radius_), lbpstride(lbpstride_), pim(pim_) {
		nbins = npoints * (npoints - 1) + 3;
		cellsize = cellsize_;
		GenerateMap();
		spoints.resize(npoints);
		ComputeSamplingPoints();
		norm = norm_;
		width = width_;
		height = height_;
		nflevels = nflevels_;
		cout << "\n Feature Normalization = " << norm << endl;
		cout << "\n Gamma Normalization = " << pim.GetGamma() << endl;
		cout << "Dim of Feature = " << (IsFolded() == true ? GetFoldedDim(
				width_, height_) : GetDim(width_, height_)) << endl;
	}
	virtual ~LBPFeatures() {
	}
	virtual UINT GetDim(UINT width_, UINT height_) const {
		return nflevels == 1 ? (width_ / lbpstride) * (height_ / lbpstride)
				* nbins : (((width_ / lbpstride) * (height_ / lbpstride))
				+ (width_ / (2 * lbpstride) * height_ / (2 * lbpstride)))
				* nbins;
	}
	virtual UINT GetFoldedDim(UINT width_, UINT height_) const { // for horizontal symmetry
		return nflevels == 1 ? ceil((REAL) width_ / (2 * lbpstride)) * (height_
				/ lbpstride) * nbins : ((ceil((REAL) width_ / (2 * lbpstride))
				* (height_ / lbpstride)) + (ceil(width_ / (4 * lbpstride))
				* height_ / (2 * lbpstride))) * nbins;
	}
	virtual void InitalizeMaps(Image &);
	virtual void InitalizeMaps(Image &, PyramidType);
	virtual void GetFeatures(UINT, int, int, vector<REAL>&);
	//	virtual void GetFoldedFeatures(UINT, int, int, vector<REAL>&);
	//	virtual void UnFoldWeights(REAL*, vector<REAL>&);
	//	void UnFoldWeights(REAL*, REAL*);
	void GenerateMap();
	int Transitions(UINT, int);
	virtual void ComputeLBPMap(Image&, LBPMap&);
	void ComputeLBPMap(Image&, LBPMap&, LBPMap&);
	void ComputeSamplingPoints();
	void NormalizeFeatures(REAL*);
	void ComputeHistogram(int initx, int inity, LBPMap &lbpmap, REAL *hist);
	void ComputeHistBilinear(int initx, int inity, LBPMap &lbpmap, REAL *hist);
	void ComputeHist(int initx, int inity, LBPMap &lbpmap, REAL *hist);

	void ComputeGradientMap(Image &image, GenericMap<REAL> &gradmap);
	void ComputeHistBilinear(int initx, int inity, LBPMap &lbpmap,
			GenericMap<REAL> &gradmap, REAL *hist); // Histogram using Gradient Maps
	int GetFeatureOffset() {
		return foffset;
	}
	UINT GetNBins() {
		return nbins;
	}
	const vector<int>& GetMap(UINT index) {
		return lbpmaps[index].lbpmap;
	}
	const vector<int>& GetFoldedMap(UINT index) {
		return flbpmaps[index].lbpmap;
	}
	const LBPMap& GetMapObj(UINT index) {
		return lbpmaps[index];
	}
	const LBPMap& GetFoldedMapObj(UINT index) {
		return flbpmaps[index];
	}
	virtual UINT GetXSize(UINT index) {
		return lbpmaps[index].nxpoints + 2 * radius;
	}
	virtual UINT GetYSize(UINT index) {
		return lbpmaps[index].nypoints + 2 * radius;
	}

protected:
	vector<LBPMap> lbpmaps, flbpmaps;// flipped lbpmaps for folded features...
	vector<GenericMap<REAL> > gradientmaps;// flipped lbpmaps for folded features...
	UINT npoints, unif, nflevels;
	vector<Point<REAL> > spoints;// sampling points on the circle...
	int radius;
	vector<UINT> map;
	NORM norm;
	UINT lbpstride;
	ProcessImage &pim;
	GridType gridtype;
	int nxcells, nycells;
	PatchType ptype;
	bool usegradmag; // whether to quantize gradient magnitude or not...
	HistogramMethod hmethod;
};

#endif
