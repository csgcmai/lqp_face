/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/
#ifndef _PYRAMID_H
#define _PYRAMID_H
#ifdef WITH_BLITZ
#include "../numericutil/gauss.h"
#include "../blitz/ext/convolve.h"
#include "../flow/rescale.h"
#include "def.h"
#endif
#include "util.hpp"
//typedef unsigned int UINT;
//typedef float REAL;
// Working on the color Images...

#define SRATIO 1.2
// For the isotropic scaling of the images
// Though octave is referred as 8, but here octave specifies number of intervals
class Pyramid {
public:
	Pyramid(UINT, UINT, UINT, FilterType = Bilinear, UINT minsize_ = 40);
	Pyramid(REAL, UINT, FilterType);
	void
	Initialize(UINT, UINT, UINT, FilterType = Bilinear, UINT minsize_ = 40);
	static void SmoothGaussian(Image&, UINT, UINT, bool = true);
	static void SmoothBilinear(Image&, UINT, UINT, bool = true);
	void SubtractGaussian(REAL factor, UINT winsize = 1, REAL sigma = 1);
	void SubtractImage(Image& image1, const Image& image2, REAL factor);

#ifdef WITH_BLITZ
	UINT ComputePyramidInfo(const Image &iImage) {
		REAL sigma = 1;
		int kernelsize = 7;
		kernel.reference(LJK::numericutil::discreteGauss(sigma, kernelsize));

		cout << " Rows = " << iImage.rows() << " Cols = " << iImage.columns()
				<< endl;
		UINT pycount = 0;

#ifdef __DEBUG
		nlevels = round(log((REAL) MIN( iImage.rows()/(REAL) minsize ,
								iImage.columns()/(REAL) minsize)) / log(sratio)) + 1;
#else
		nlevels = floor(log((REAL) MIN( iImage.rows()/(REAL) minsize ,
				iImage.columns()/(REAL) minsize)) / log(sratio)) + 1;
#endif
		cout << " Number of Levels = " << nlevels << endl;
		scale.resize(nlevels, 0);
		for (; pycount < nlevels; ++pycount) {
			scale[pycount] = pow((double) sratio, (double) pycount);
		}
		cout << " pycount = " << pycount << endl;
		return nlevels;
	}
	mImageType ScaleImage(Image&);
	mImageType ScaleAndSmoothImage(Image&, REAL factor);
	mImageType ScaleImage(UINT index, const Image&);
	mImageType ScaleAndSmoothImage(UINT index, const Image&, REAL factor);
	mImageType ScaleAndSmoothImage2(UINT index, const Image&, REAL factor);

	mImageType& GetSmoothImage(UINT index) {
			return simgobj[index];
	}
	void WriteImage(const std::string& iname, blitz::Array<PixelType, 2> &img) {
		Magick::Image outimg(Magick::Geometry(img.columns(), img.rows()),
				"Black");
		blitz::TinyVector<float, 3> tvec;
		Magick::ColorRGB cvar;
		for (int i = 0; i < img.rows(); ++i)
			for (int j = 0; j < img.columns(); ++j) {
				tvec = img(i, j);
				cvar.red(tvec[0]);
				cvar.green(tvec[1]);
				cvar.blue(tvec[2]);
				outimg.pixelColor(j, i, cvar);
			}
		outimg.write(iname);

	}
	void ReadImage(const Image & image, mImageType & imgout) {
		int m = image.rows(), n = image.columns();
		imgout.resize(m, n);
		Magick::ColorRGB mycolor;
		blitz::TinyVector<REAL, 3> tvec;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				mycolor = image.pixelColor(j, i);
				tvec[0] = mycolor.red();
				tvec[1] = mycolor.green();
				tvec[2] = mycolor.blue();
				imgout(i, j) = tvec;
			}
		}
	}
#endif
	inline UINT GetNLevels() {
		return nlevels;
	}
	inline UINT GetInitIndex() {
		return nloctave;
	}
	float GetScale(UINT index) {
		return scale[index];
	}
	UINT GeneratePyramid(const Image &);
	UINT GeneratePyramidSimple(const Image &);

	REAL GetScaleRatio() {
		return sratio;
	}
	Image & GetImageRef(UINT index) {
		return imgobj[index];
	}
	void GetImage(UINT index, Image &image) {
		image = imgobj[index];
	}
	void SetOctaveSize(UINT nloctave_) {
		nloctave = nloctave_;
		sratio = pow(2.0, 1.0 / nloctave);
	}
	UINT GetOctaveSize() {
		return nloctave;
	}
	UINT GetXSize(UINT index) {
		return imgobj[index].columns();
	}
	UINT GetYSize(UINT index) {
		return imgobj[index].rows();

	}
	void Clear() {
		imgobj.clear();
		scale.clear();
	}
private:
	//	vector<Image> imgObj;
	REAL sratio;
	UINT nlevels; // total number of levels;
	UINT nloctave; // number of level in an octave...
	vector<REAL> scale;
	vector<Image> imgobj;

	UINT rows, cols;
	UINT minsize;
	FilterType ftype;
#ifdef WITH_BLITZ
	vector<mImageType> simgobj; // smooth Image Obj
	blitz::Array<REAL, 1> kernel;
#endif
};
#endif
