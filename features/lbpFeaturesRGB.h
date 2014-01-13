/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/
#ifndef LBPFEATURESRGB_H
#define LBPFEATURESRGB_H
#include "features.h"
#include "../util/processImage.h"
#include "../util/hogInfo.h"
/*LBP Features for Color-Scale Images: Note that features are computed at window level and thus can be used
 * with any step size or can be computed for any given location*/
class LBPFeaturesRGB: public Features {
public:
	LBPFeaturesRGB(const LBPParams &lbpparam, UINT width_, UINT height_,
			UINT nplevels, ProcessImage &pim_) :
		Features(width_, height_, nplevels, false), pim(pim_) {
		//		 npoints(npoints_), unif(2),
		//						radius(radius_), lbpstride(lbpstride_), pim(pim_)
		npoints = lbpparam.npoints;
		ptype = lbpparam.ptype;
		unif = 2;
		lbpstride = lbpparam.stride;
		add = lbpparam.addfeat;
		cellsize = lbpparam.cellsize;

		hmethod = lbpparam.hmethod;
		nseparate = lbpparam.normfeatsep;
		gridtype = lbpparam.gridtype;
		radius = lbpparam.radius;
		norm = lbpparam.featnorm;
		nflevels = lbpparam.nflevels;
		usemeanval = lbpparam.usemeanval;
		lbptolerance = lbpparam.lbptolerance;
		if (ptype == PT_HStrip || ptype == PT_VStrip || ptype == PT_Diag
				|| ptype == PT_ADiag)
			radius = npoints / 2;

		cout << "\n LBP Grid Type = " << (gridtype == Circular ? "Circular"
				: "Rectangular") << flush << endl;
		dimmult = add == true ? 1 : 3; // if !add concatenate features on RGB channel...
		rinvariant = lbpparam.rinvariant;
		nbins = rinvariant == true ? 10 : (lbpparam.nbins != pow(2.0,
				(double) npoints) ? npoints * (npoints - 1) + 3
				: lbpparam.nbins);
		GenerateMap();
		cout << endl;
		spoints.resize(npoints);
		biinfo = new double[npoints * 10]; // 11 variables per sampling point
		dointerpol = new bool[npoints]; // Interpolation flag for

		width = width_;
		height = height_;

		ComputeSamplingPoints();
		nycells = height / cellsize;
		nxcells = width / cellsize;
		usegradmag = lbpparam.usegradmag;

		lbpparam.PrintParamInfo();
		cout << endl << " ------------------------------------" << endl
				<< " Width = " << width << " Height = " << height
				<< " LBPStride= " << lbpstride << endl
				<< " Feature Normalization Type = " << norm << endl
				<< " LBPbins = " << nbins
				<< (add == true ? "  & Adding RGB Channels at Cell levels"
						: "  Concatenating RGB Channels at Cell levels ")
				<< endl << "Patch Geometrical Organization =  "
				<< PatchParams::GetPatchType(ptype) << endl << flush;
		cout << " Nxcells = " << nxcells << " Nycells = " << nycells << flush;
		cout << " LBPDim = " << (IsFolded() == true ? GetFoldedDim(width_,
				height_) : GetDim(width_, height_)) << endl;
		if (hmethod == HM_Bilinear)
			cout << "~~~~~~~~~~~~~~~~ Interpolated Version ~~~~~~~~~~~~~~~~~~~";
		else
			cout << "~~~~~~~~~~~~~~~~ Discrete Version ~~~~~~~~~~~~~~~~~~~";

#ifdef TMPDEBUG
		fcount = 'a';
		fname = "debug";
#endif
#ifdef CROPPEDWINDOWS
		cout << "~~~~~~~~~~~~~~~~CROPPED WINDOWS ~~~~~~~~~~~~~~~~~" << endl;
		cwlbpmaps[0].Init(width - 2 * radius, height - 2 * radius);
		cwlbpmaps[1].Init(width - 2 * radius, height - 2 * radius);
		cwlbpmaps[2].Init(width - 2 * radius, height - 2 * radius);
#endif
	}
	void ComputeLBPMap(UINT index, int xoffset, int yoffset);
	virtual UINT PrintCodesInfo() {

	}
	LBPFeaturesRGB(UINT cellsize_, UINT npoints_, UINT radius_, UINT width_,
			UINT height_, UINT nplevels, NORM norm_, UINT nflevels_,
			UINT lbpstride_, ProcessImage &pim_, bool folded_ = false,
			bool add_ = true, UINT nbins_ = 59, bool nseparate_ = true,
			GridType gridtype_ = Circular) :
		Features(width_, height_, nplevels, folded_), npoints(npoints_),
				unif(2), radius(radius_), lbpstride(lbpstride_), pim(pim_) {
		add = add_;
		cellsize = cellsize_;
		nseparate = nseparate_;
		gridtype = gridtype_;
		usemeanval = false;
		lbptolerance = 0;
		cout << "\n LBP Grid Type = " << (gridtype == Circular ? "Circular"
				: "Rectangular") << flush << endl;
		dimmult = add == true ? 1 : 3; // if !add concatenate features on RGB channel...
		nbins = nbins_ != 256 ? npoints * (npoints - 1) + 3 : nbins_;
		GenerateMap();
		cout << endl;
		spoints.resize(npoints);
		biinfo = new double[npoints * 10]; // 11 variables per sampling point
		dointerpol = new bool[npoints]; // Interpolation flag for
		ComputeSamplingPoints();
		norm = norm_;
		width = width_;
		nflevels = nflevels_;
		height = height_;
		ptype = PT_Circular; // here
		nycells = height / cellsize;
		nxcells = width / cellsize;
		cout << endl << " ------------------------------------" << endl
				<< " Width = " << width << " Height = " << height
				<< " LBPStride= " << lbpstride << endl
				<< " Feature Normalization Type = " << norm << endl
				<< " LBPbins = " << nbins
				<< (add == true ? "  & Adding RGB Channels at Cell levels"
						: "  Concatenating RGB Channels at Cell levels ")
				<< endl << flush;
		cout << " Nxcells = " << nxcells << " Nycells = " << nycells << flush;
		cout << " LBPDim = " << (IsFolded() == true ? GetFoldedDim(width_,
				height_) : GetDim(width_, height_)) << endl;

	}
	virtual ~LBPFeaturesRGB() {
		delete[] biinfo;
		delete[] dointerpol;
	}
	virtual UINT GetDim(UINT width_, UINT height_) const {
		return nflevels == 1 ? (width_ / lbpstride) * (height_ / lbpstride)
				* nbins * dimmult : (((width_ / lbpstride) * (height_
				/ lbpstride)) + (width_ / (2 * lbpstride) * height_ / (2
				* lbpstride))) * nbins * dimmult;
	}
	virtual UINT GetFoldedDim(UINT width_, UINT height_) const { // for horizontal symmetry
		return nflevels == 1 ? ceil((REAL) width_ / (2 * lbpstride)) * (height_
				/ lbpstride) * nbins * dimmult : ((ceil(
				(REAL) width_ / (2 * lbpstride)) * (height_ / lbpstride))
				+ (ceil(width_ / (4 * lbpstride)) * height_ / (2 * lbpstride)))
				* nbins * dimmult;
	}
	virtual void InitalizeMaps(Image &image);
	virtual void InitalizeMaps(Pyramid &pyobj_, PyramidType sspace);
	virtual void InitalizeMaps(Image &, PyramidType);
	virtual void InitalizeMapsFast(Image &, PyramidType);
	virtual void GetFeatures(UINT, int, int, vector<REAL>&);
	virtual void GetFoldedFeatures(UINT, int, int, vector<REAL>&);
	virtual void UnFoldWeights(REAL*, vector<REAL>&);
	void UnFoldWeights(REAL*, REAL*);
	void GenerateMap();
	int Transitions(UINT, int);
	//	void ComputeLBPMap(Image&, LBPMap&);
	//	void ComputeLBPMap(Image&, LBPMap&, LBPMap&);
	virtual void ComputeLBPMap(Image &image, vector<LBPMap>&);
	void ComputeLBPMapFast(Image &image, vector<LBPMap>&);
	void ComputeGradientMap(Image &image, GenericMap<REAL> &gradmap);
	virtual void ComputeLBPMap(Image &image, LBPMap[]);
	void ComputeMeanLBPMap(Image &image, LBPMap lbpmap[]);
	void ComputeLBPMapFast(Image &image, LBPMap[]);
	void ComputeSamplingPoints();
	void NormalizeFeatures(REAL*);
	void NormalizeFeatures(vector<REAL>&feat) {
		NormalizeFeatures(&feat[0]);
	}
	void ComputeHistogram(int initx, int inity, LBPMap &lbpmap, REAL *hist);
	void ComputeHistBilinear(int initx, int inity, LBPMap &lbpmap, REAL *hist);
	void ComputeHistBilinear(int initx, int inity, LBPMap &lbpmap,
			GenericMap<REAL> &gradmap, REAL *hist);
	int GetFeatureOffset() {
		return foffset;
	}
	UINT GetNBins() {
		return nbins;
	}
	const vector<int>& GetMap(UINT index) {
		return lbpmaps[index][0].lbpmap;
	}
	const vector<int>& GetFoldedMap(UINT index) {
		return flbpmaps[index].lbpmap;
	}
	const LBPMap& GetMapObj(UINT index) {
		return lbpmaps[index][0];
	}
	const LBPMap& GetFoldedMapObj(UINT index) {
		return flbpmaps[index];
	}
#ifndef CROPPEDWINDOWS
	virtual UINT GetXSize(UINT index) {
		return lbpmaps[index][0].nxpoints + 2 * radius;
	}
	virtual UINT GetYSize(UINT index) {
		return lbpmaps[index][0].nypoints + 2 * radius;
	}
#else
	virtual UINT GetXSize(UINT index) {
		Image & ref = pyobj.GetImageRef(index);
		return ref.columns();
	}
	virtual UINT GetYSize(UINT index) {
		Image & ref = pyobj.GetImageRef(index);
		return ref.rows();
	}
#endif
	vector<UINT>& GetInitMap() {
		return map;
	}
	vector<Point<REAL> >& GetSamplingPoints() {
		return spoints;
	}
protected:
#ifdef CROPPEDWINDOWS
	LBPMap cwlbpmaps[3];
	vector<vector<REAL> > rimpix;
	vector<vector<REAL> > gimpix;
	vector<vector<REAL> > bimpix;
	vector<UINT> xdims;
#endif
	vector<vector<LBPMap> > lbpmaps;// flipped lbpmaps for folded features...
	vector<GenericMap<REAL> > gradientmaps;// flipped lbpmaps for folded features...
	vector<LBPMap> flbpmaps;
	UINT npoints, unif, nflevels;
	vector<Point<REAL> > spoints;// sampling points on the circle...
	int radius;
	vector<UINT> map;
	REAL lbptolerance;
	UINT fmap[59];
	NORM norm;
	UINT lbpstride;
	ProcessImage &pim;
	bool add, nseparate; // if add add all the color-channels else concatenates
	PatchType ptype; // by default circular , Geometrical  organization of the features...
	UINT dimmult;// dimensional multipliear
	HOGInfo lbpfeatures;
	GridType gridtype;
	double *biinfo; // contains the
	bool *dointerpol;
	HistogramMethod hmethod;
	int nxcells, nycells;
	bool usegradmag; // whether to quantize gradient magnitude or not...
	bool usemeanval; // use mean of the circular points for the comparison...
	bool rinvariant; /// rotation invariant coding
#ifdef TMPDEBUG
	char fcount;
	string fname;
#endif
};

class MeanLBPFeaturesRGB: public LBPFeaturesRGB {
public:
	MeanLBPFeaturesRGB(const LBPParams &lbpparam, UINT width_, UINT height_,
			UINT nplevels, ProcessImage &pim_) :
		LBPFeaturesRGB(lbpparam, width_, height_, nplevels, pim_) {

	}
	virtual void ComputeLBPMap(Image &image, vector<LBPMap> &lbpmap) {
		UINT xdim = image.columns(), ydim = image.rows();
		int nypoints = ydim - 2 * radius, reg = 1, fx, fy, cx, cy, nxpoints =
				xdim - 2 * radius, cenx, ceny, count, rvalue, gvalue, bvalue;
		REAL tmpx, tmpy, fracx, fracy, rcenval, gcenval, bcenval;
		REAL x, y, w1, w2, w3, w4;
		vector<REAL> pixval(npoints * 3, 0);
		REAL *rpixval = &pixval[0], *gpixval = rpixval + npoints, *bpixval =
				gpixval + npoints;
		lbpmap[0].Init(nxpoints, nypoints); // for RGB channels
		lbpmap[1].Init(nxpoints, nypoints);
		lbpmap[2].Init(nxpoints, nypoints);
		vector<REAL> rimpix(xdim * ydim, 0), bimpix(xdim * ydim, 0), gimpix(
				xdim * ydim, 0);
		//	image.flop();
		pim.Process(image, rimpix, gimpix, bimpix);
		REAL rmeanval, gmeanval, bmeanval;
		//	int fbins[]={4,3,2,1,0,7,6,5};

		//		ofstream ofile("test_map.lst");
		for (int i = 0; i < nypoints; ++i) {
			for (int j = 0; j < nxpoints; ++j) {
				cenx = j + radius;
				ceny = i + radius;
				count = 0;
				rcenval = rimpix[ceny * xdim + cenx];
				gcenval = gimpix[ceny * xdim + cenx];
				bcenval = bimpix[ceny * xdim + cenx];
				rmeanval = gmeanval = bmeanval = 0;
				//ofile<<endl<< i<<" " << j<< " "<< rcenval <<" " ;
				for (vector<Point<REAL> >::iterator iter = spoints.begin(); iter
						!= spoints.end(); ++iter, ++count) {
					iter->GetCoord(x, y);
					tmpy = ceny + y;
					tmpx = cenx + x;
					if (ptype == PT_Circular && gridtype == Circular) { // do bilinear interpolation
						fx = floor(tmpx);
						fy = floor(tmpy);
						cx = ceil(tmpx);
						cy = ceil(tmpy);
						fracx = tmpx - fx;
						fracy = tmpy - fy;

						if (ABS(fracx) < 1e-6 && ABS(fracy) < 1e-6)
						// no interpolation needed
						{
							rpixval[count] = rimpix[fy * xdim + fx];
							gpixval[count] = gimpix[fy * xdim + fx];
							bpixval[count] = bimpix[fy * xdim + fx];
						} else {
							w1 = (1 - fracx) * (1 - fracy);
							w2 = fracx * (1 - fracy);
							w3 = (1 - fracx) * fracy;
							w4 = fracx * fracy;
							rpixval[count] = w1 * rimpix[fy * xdim + fx] + w2
									* rimpix[fy * xdim + cx] + w3 * rimpix[cy
									* xdim + fx] + w4 * rimpix[cy * xdim + cx];

							gpixval[count] = w1 * gimpix[fy * xdim + fx] + w2
									* gimpix[fy * xdim + cx] + w3 * gimpix[cy
									* xdim + fx] + w4 * gimpix[cy * xdim + cx];

							bpixval[count] = w1 * bimpix[fy * xdim + fx] + w2
									* bimpix[fy * xdim + cx] + w3 * bimpix[cy
									* xdim + fx] + w4 * bimpix[cy * xdim + cx];

						}
					} else {
						rpixval[count] = rimpix[tmpy * xdim + tmpx];
						gpixval[count] = gimpix[tmpy * xdim + tmpx];
						bpixval[count] = bimpix[tmpy * xdim + tmpx];
					}
					rmeanval += rpixval[count];
					gmeanval += gpixval[count];
					bmeanval += bpixval[count];
				}
				rvalue = 0;
				gvalue = 0;
				bvalue = 0;
				int wt = 0;
				rmeanval /= npoints;
				gmeanval /= npoints;
				bmeanval /= npoints;

				for (UINT iter = 0; iter < npoints; ++iter) {
					wt = reg << iter;
					rvalue += ((rpixval[iter] - rmeanval) >= 0 ? wt : 0);
					gvalue += ((gpixval[iter] - gmeanval) >= 0 ? wt : 0);
					bvalue += ((bpixval[iter] - bmeanval) >= 0 ? wt : 0);
				}

				lbpmap[0](j, i) = map[rvalue];
				lbpmap[1](j, i) = map[gvalue];
				lbpmap[2](j, i) = map[bvalue];
			}
		}
	}
};

#endif
