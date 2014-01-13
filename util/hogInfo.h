/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/
#ifndef HOGINFO_H_
#define HOGINFO_H_
#include "util.hpp"
#include<iostream>
#include<fstream>
#include "cache.h"

using namespace std;
/*HInfo is the major data structure used to store and access features and perform dot products of features with a given filter...
 * Features are stored in a linear row-major order i.e: first features belonging to first cell of first row are stored and then
 * the next cell of first row and so.
 *
 * Function: GetPermutationFeatures defines permutations (for the HOG,LBP,LTP and their different combinations) for flipping
 * or folding the features..
 *
 * Caution: There are some redundant interfaces to maintain backward compatability...
 * */
class HInfo {
	friend class HOGInfo;
public:
	HInfo(UINT = 0, UINT = 0, REAL = 1, UINT = 36);
	HInfo(UINT, UINT, UINT, REAL*);
	HInfo(UINT, UINT, REAL, UINT, REAL*);

	HInfo(const HInfo &);
	void Initialize(UINT, UINT, REAL, UINT); // number of orientations...
	void Initialize(UINT, UINT, REAL, UINT, const vector<REAL>&);
	void Initialize(UINT, UINT, REAL, UINT, REAL*);
	void SetOrient(UINT orient) {
		norient = orient;
	}

	void GetHOGCells(UINT, UINT, UINT, UINT, REAL*) const;
	void GetSkippedHOGCells(UINT, UINT, UINT, UINT, UINT, REAL*) const;// skips the cell in between...
	void Get2XHOGCells(UINT, UINT, UINT, UINT, REAL*) const;
	void GetFoldedHOGCells(UINT, UINT, UINT, UINT, REAL*) const;
	void GetUnFoldedHOGCells(UINT, REAL *) const;
	// Get as parameter the folding & unfolding indeces...
	void GetSkippedFoldedHOGCells(UINT, UINT, UINT, UINT, UINT, REAL*) const;
	//	void GetSkippedFlippedFeatures(UINT,UINT,UINT,UINT,UINT,vector<REAL>&)const;
	void GetSkippedFlippedFeatures(UINT, UINT, UINT, UINT, UINT, REAL*) const;

	REAL GetSumHOGCells(UINT, UINT, UINT, UINT, bool = true) const;
	void SetHOGCells(UINT, UINT, UINT, UINT, REAL = 0);
	void TrimWeights(const UINT&);
	void FlipFeatures();// Reflect the features around vertical axis....
	void GetFlippedFeatures(HInfo&);
	//	void GetFlippedFeatures(UINT,UINT,UINT,UINT,vector<REAL>&)const;
	void GetFlippedFeatures(UINT, UINT, UINT, UINT, REAL*) const;
	static void FlipFeatures(UINT, UINT, UINT, REAL*, REAL*);

	void PadFeatureMap(UINT, UINT);
	REAL DotProduct(UINT, UINT, UINT, UINT, const vector<REAL>&) const;
	REAL DotProduct(UINT, UINT, UINT, UINT, REAL*) const;
	void DotProduct(UINT, UINT, const vector<REAL>&, vector<REAL>&) const;
	void DotProduct(UINT, UINT, const vector<REAL>&, Store&) const;
	void DotProduct(UINT, UINT, const REAL*, Store&) const;
	void DotProductWithQuadratic(UINT, UINT, const REAL&, const REAL*,
			const REAL*, Store&) const;
	void DotProductWithProjection(UINT fxsize, UINT fysize,
			const VectorXf &sigma, const MatrixXf &projMat, const REAL *filter,
			Store&outhog) const;
	void SkippedDotProduct(UINT, UINT, UINT, const REAL*, Store&) const;
	REAL SymmetricDotProduct(UINT, UINT, UINT, UINT, const REAL*) const;
	void PruneNegativeWeights();
	void BlackOut(UINT, UINT, UINT, UINT);
	UINT GetWidth() {
		return xbmax;
	}
	UINT GetHeight() {
		return ybmax;
	}
	REAL GetScale() {
		return scale;
	}
	void operator =(const HInfo &);
	vector<REAL>& GetFeatureReference() {
		return hogvalue;
	}
	REAL *GetFeaturePtr() {
		return (&hogvalue[0]);
	}
	void BilinearInterpolation(REAL, HInfo &);
	void Write(fstream&);
	void WriteText(fstream&);
	void Read(fstream&);
	static void GetPermutationFeatures(int *pfeat, int norient);
private:
	UINT xbmax, ybmax; // number of maximum cell per dimension...
	REAL scale;
	UINT norient;
	vector<REAL> hogvalue; // hogvalue
	vector<int> permvec;// permutation vector used for flipping and folding features...
};

class HOGInfo {
public:
	HOGInfo() {
		nalevel = 0;
		counter = 0;
	}
	HOGInfo(UINT, UINT = 0);
	HOGInfo(const HOGInfo&);
	void Initialize(UINT, UINT, UINT);
	void SetFeature(UINT, UINT, REAL, const vector<REAL>&);
	void SetFeature(UINT, UINT, REAL, REAL*);
	void SetFeature(UINT, UINT, REAL);
	void SetFeature(UINT, UINT, UINT, REAL);
	void SetInfo(UINT, REAL);
	bool ReadInfo(fstream&);
	void WriteInfo(fstream&);
	void WriteTextInfo(fstream&);
	vector<REAL>* GetFeature(UINT index) {
		return &imfeatures[index].hogvalue;
	}
	vector<REAL>& GetFeatureRef(UINT index) {
		return imfeatures[index].hogvalue;
	}
	REAL GetScale(UINT index) const {
		return imfeatures[index].scale;
	}
	void operator=(const HOGInfo&);
	void ReInitialize(UINT, UINT, REAL);
	UINT GetNLevels() const {
		return imfeatures.size();
	}
	UINT GetInitIndex() const {
		return nalevel;
	}
	UINT GetXBlock(UINT index) const {
		return imfeatures[index].xbmax;
	}
	UINT GetYBlock(UINT index) const {
		return imfeatures[index].ybmax;
	}
	HInfo& GetHInfo(UINT index) {
		return imfeatures[index];
	}
	UINT Get2XScaleIndex(UINT ulimit) const {
		return (ulimit - nalevel < 0 ? 0 : ulimit - nalevel);
	}
	UINT GetNearestScaleIndex(REAL) const;

	void GetHOGCells(UINT index, UINT xmin, UINT ymin, UINT width, UINT height,
			REAL *outhog) const {
		imfeatures[index].GetHOGCells(xmin, ymin, width, height, outhog);
	}
	void GetSkippedHOGCells(UINT index, UINT xmin, UINT ymin, UINT width,
			UINT height, UINT skip, REAL *outhog) const {
		imfeatures[index].GetSkippedHOGCells(xmin, ymin, width, height, skip,
				outhog);
	}
	void GetFoldedHOGCells(UINT index, UINT xmin, UINT ymin, UINT width,
			UINT height, REAL *outhog) const {
		imfeatures[index].GetFoldedHOGCells(xmin, ymin, width, height, outhog);
	}
	void GetSkippedFoldedHOGCells(UINT index, UINT xmin, UINT ymin, UINT width,
			UINT height, UINT skip, REAL *outhog) const {
		imfeatures[index].GetSkippedFoldedHOGCells(xmin, ymin, width, height,
				skip, outhog);
	}
	void DisplayHOG(UINT);
	REAL DotProduct(UINT index, UINT xstart, UINT ystart, UINT xsize,
			UINT ysize, const vector<REAL>&inhog) const {
		return imfeatures[index].DotProduct(xstart, ystart, xsize, ysize, inhog);
	}
	REAL DotProduct(UINT index, UINT xstart, UINT ystart, UINT xsize,
			UINT ysize, REAL*inhog) const {
		return imfeatures[index].DotProduct(xstart, ystart, xsize, ysize, inhog);
	}
	void DotProduct(UINT index, UINT xsize, UINT ysize, vector<REAL>& filter,
			Store&inhog) const {
		imfeatures[index].DotProduct(xsize, ysize, filter, inhog);
	}
	void DotProduct(UINT index, UINT xsize, UINT ysize, const REAL *filter,
			Store&inhog) const {
		imfeatures[index].DotProduct(xsize, ysize, filter, inhog);
	}
	void DotProductv1(UINT index, UINT xsize, UINT ysize, const REAL *filter,
			Store&inhog) const {
		imfeatures[index].DotProduct(xsize, ysize, filter, inhog);
	}
	void DotProductWithQuadratic(UINT index, UINT xsize, UINT ysize,
			const REAL &meanssd, const REAL *mpvector, const REAL *filter,
			Store&inhog) const {
		imfeatures[index].DotProductWithQuadratic(xsize, ysize, meanssd,
				mpvector, filter, inhog);
	}
	void DotProductWithProjection(UINT index, UINT xsize, UINT ysize,
			const VectorXf &sigma, const MatrixXf &projMat, const REAL *filter,
			Store&inhog) const {
		imfeatures[index].DotProductWithProjection(xsize, ysize, sigma,
				projMat, filter, inhog);
	}

	void SkippedDotProduct(UINT index, UINT xsize, UINT ysize, UINT skip,
			const REAL *filter, Store&inhog) const {
		imfeatures[index].SkippedDotProduct(xsize, ysize, skip, filter, inhog);
	}
	REAL SymmetricDotProduct(UINT index, UINT xstart, UINT ystart, UINT xsize,
			UINT ysize, REAL*inhog) const {
		return imfeatures[index].SymmetricDotProduct(xstart, ystart, xsize,
				ysize, inhog);
	}
	REAL SymmetricDotProduct(UINT index, UINT xstart, UINT ystart, UINT xsize,
			UINT ysize, const vector<REAL>&inhog) const {
		return imfeatures[index].SymmetricDotProduct(xstart, ystart, xsize,
				ysize, &inhog[0]);
	}
	void FlipFeatures(UINT index) {
		imfeatures[index].FlipFeatures();
	}
	void FlipFeatures();
	void GetFlippedFeatures(UINT index, HInfo&hinfo) {
		imfeatures[index].GetFlippedFeatures(hinfo);
	}
	void GetFlippedFeatures(UINT index, UINT xmin, UINT ymin, UINT width,
			UINT height, REAL*outhog) {
		imfeatures[index].GetFlippedFeatures(xmin, ymin, width, height, outhog);
	}
	void GetSkippedFlippedFeatures(UINT index, UINT xmin, UINT ymin,
			UINT width, UINT height, UINT skip, REAL*outhog) {
		imfeatures[index].GetSkippedFlippedFeatures(xmin, ymin, width, height,
				skip, outhog);
	}
	void PadFeatureMap(UINT index, UINT padx, UINT pady) {
		imfeatures[index].PadFeatureMap(padx, pady);
	}
	void BlackOut(UINT, UINT, UINT, UINT);
private:
	vector<HInfo> imfeatures;
	UINT nalevel; // number of levels above zero.
	UINT counter;
	UINT norient;
};

#endif /*HOGINFO_H_*/
