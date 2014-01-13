/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/
#include "npyramid.h"

const static double gkernel[] =
		{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 7.22562316040198e-06, 3.23829963443551e-05,
				5.33905348819428e-05, 3.23829963443551e-05,
				7.22562316040198e-06, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				1.96412801362135e-05, 0.000239279776639906, 0.00107237755972000,
				0.00176805169293186, 0.00107237755972000, 0.000239279776639906,
				1.96412801362135e-05, 0, 0, 0, 0, 0, 0, 0, 7.22562316040198e-06,
				0.000239279776639906, 0.00291502443383413, 0.0130642331448828,
				0.0215392790713540, 0.0130642331448828, 0.00291502443383413,
				0.000239279776639906, 7.22562316040198e-06, 0, 0, 0, 0, 0, 0,
				3.23829963443551e-05, 0.00107237755972000, 0.0130642331448828,
				0.0585498308977697, 0.0965323515970485, 0.0585498308977697,
				0.0130642331448828, 0.00107237755972000, 3.23829963443551e-05,
				0, 0, 0, 0, 0, 0, 5.33905348819428e-05, 0.00176805169293186,
				0.0215392790713540, 0.0965323515970485, 0.159154941388757,
				0.0965323515970485, 0.0215392790713540, 0.00176805169293186,
				5.33905348819428e-05, 0, 0, 0, 0, 0, 0, 3.23829963443551e-05,
				0.00107237755972000, 0.0130642331448828, 0.0585498308977697,
				0.0965323515970485, 0.0585498308977697, 0.0130642331448828,
				0.00107237755972000, 3.23829963443551e-05, 0, 0, 0, 0, 0, 0,
				7.22562316040198e-06, 0.000239279776639906, 0.00291502443383413,
				0.0130642331448828, 0.0215392790713540, 0.0130642331448828,
				0.00291502443383413, 0.000239279776639906, 7.22562316040198e-06,
				0, 0, 0, 0, 0, 0, 0, 1.96412801362135e-05, 0.000239279776639906,
				0.00107237755972000, 0.00176805169293186, 0.00107237755972000,
				0.000239279776639906, 1.96412801362135e-05, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 7.22562316040198e-06, 3.23829963443551e-05,
				5.33905348819428e-05, 3.23829963443551e-05,
				7.22562316040198e-06, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, };

Pyramid::Pyramid(UINT cols_, UINT rows_, UINT nlevels_, FilterType ftype_,
		UINT minsize_) { // nolevels_ is the number of level in one octave...
	UINT nolevels_ = nlevels_;
	Initialize(cols_, rows_, nolevels_, ftype_, minsize_);
}
Pyramid::Pyramid(REAL scaleratio_, UINT minsize_, FilterType ftype_) {
	sratio = scaleratio_;
	minsize = MAX(16, minsize_);
	ftype = ftype_;
	nloctave = round(log(2.0) / log(sratio));
	cout << " Number of Levels Per Octave = " << nloctave << endl;
}
void Pyramid::Initialize(UINT cols_, UINT rows_, UINT nolevels_,
		FilterType ftype_, UINT minsize_) {
	nlevels = 0;
	cols = cols_;
	rows = rows_;
	minsize = MAX(16,minsize_);
	nloctave = nolevels_;
	ftype = ftype_;
	sratio = pow(2.0, (double) 1.0 / (nolevels_)); // log_2(2) = nolevels * log_2(sRatio);
	if (nolevels_ != 1)
		cout << "Image Pyramid Scale ratio = " << sratio << endl;
}
// for more detail on the scale & resize methods look convert documentation in Magick++
void Pyramid::SmoothGaussian(Image &simage, UINT rows, UINT cols, bool aspect) {
	if (simage.rows() == rows && simage.columns() == cols)
		return;
	int w = (UINT) (2 * round(cols / simage.columns()) + 1), // width of gaussian mask
			h = (UINT) (2 * round(rows / simage.rows()) + 1); // height of gaussian mask
	int size = w < h ? w : h;
	float sigma = (size / 2 - 1) * 0.3 + 0.8;
	Geometry gobj = Geometry(cols, rows, 0, 0);
	gobj.aspect(aspect);
	simage.gaussianBlur(size, sigma);
	simage.scale(gobj); // should use a sample function instead of scale
}
// Scaling the image using bilinear interpolation
void Pyramid::SmoothBilinear(Image &simage, UINT rows, UINT cols, bool aspect) {
	// Image to be scale, output rows and columns, aspect=true keep the aspect ratio...
	if (simage.rows() == rows && simage.columns() == cols)
		return;
	Geometry gobj = Geometry(cols, rows, 0, 0);
	gobj.aspect(aspect); // if aspect true then resize without preserving aspect ratio...
	simage.filterType(TriangleFilter);
	simage.resize(gobj);
}

UINT Pyramid::GeneratePyramid(const Image &iImage)
// Octave Based Pyramid Implementation....

		{

	cout << " Rows = " << iImage.rows() << " Cols = " << iImage.columns()
			<< endl;
	//	float sScale = 1, eScale, trows = ((float) iImage.rows() / rows), tcols =
	//			((float) iImage.columns() / cols);
	//	eScale = (trows < tcols ? trows : tcols);//end_scale
	UINT pycount = 0;

	//	nlevels= floor( log(eScale/sScale)/ log(sratio) ) +1 ; old criteria
	// New criteria (5*xcell=40 )
	nlevels = floor(
			log((REAL) MIN(iImage.rows(), iImage.columns()) / (REAL) minsize)
					/ log(sratio)) + 1;
	cout << " Number of Levels = " << nlevels << endl;
	//nlevels = nlevels < nloctave ? nloctave:nlevels;
	nloctave = nlevels < nloctave ? nlevels : nloctave;
	try {
		imgobj.resize(nlevels);
		scale.resize(nlevels, 0);
		for (; pycount < nloctave; ++pycount) {
			scale[pycount] = pow((double) sratio, (double) pycount);
			imgobj[pycount] = iImage;

			if (ftype == Gaussian)
				SmoothGaussian(imgobj[pycount],
						(UINT) (iImage.rows() / scale[pycount]),
						(UINT) (iImage.columns() / scale[pycount]), true);
			else
				SmoothBilinear(imgobj[pycount],
						(UINT) (iImage.rows() / scale[pycount]),
						(UINT) (iImage.columns() / scale[pycount]), true);
			// Generate next image in the pyrmaid
			for (UINT k = nloctave + pycount; k < nlevels; k += nloctave) {
				scale[k] = scale[k - nloctave] / 0.5;
				imgobj[k] = imgobj[k - nloctave];
				if (ftype == Gaussian)
					SmoothGaussian(imgobj[k], (UINT) (iImage.rows() / scale[k]),
							(UINT) (iImage.columns() / scale[k]), true);
				else
					SmoothBilinear(imgobj[k], (UINT) (iImage.rows() / scale[k]),
							(UINT) (iImage.columns() / scale[k]), true);
			}
		}
	} catch (Exception &exp) {
		cout << "Exception Caught:";
		exp.what();
	}
	cout << " Number of Levels = " << pycount << endl;
	return nlevels;
}

UINT Pyramid::GeneratePyramidSimple(const Image &iImage)
// Simple Sequential Image Pyramid Implementation....

		{
	cout << " Rows = " << iImage.rows() << " Cols = " << iImage.columns()
			<< endl;
	UINT pycount = 0;
	nlevels = floor(log((REAL) MIN( iImage.rows()/(REAL) minsize,
			iImage.columns()/(REAL) minsize)) / log(sratio)) + 1;
	cout << " Number of Levels = " << nlevels << endl;
	//nlevels = nlevels < nloctave ? nloctave:nlevels;
	try {
		imgobj.resize(nlevels);
		scale.resize(nlevels, 0);
		for (; pycount < nlevels; ++pycount) {
			scale[pycount] = pow((double) sratio, (double) pycount);
			imgobj[pycount] = iImage;

			if (ftype == Gaussian)
				SmoothGaussian(imgobj[pycount],
						(UINT) (iImage.rows() / scale[pycount]),
						(UINT) (iImage.columns() / scale[pycount]), true);
			else
				SmoothBilinear(imgobj[pycount],
						(UINT) (iImage.rows() / scale[pycount]),
						(UINT) (iImage.columns() / scale[pycount]), true);
		}
	} catch (Exception &exp) {
		cout << "Exception Caught:";
		exp.what();
	}
	cout << " Number of Levels = " << pycount << endl;
	return nlevels;
}
#ifdef WITH_BLITZ
void Pyramid::SubtractGaussian(REAL factor, UINT winsize, REAL sigma) {
	if (winsize % 2 != 1) {
		cout << " Kernel Size Must be Odd as 2*kernel+1 is total size ";
		exit(EXIT_FAILURE);
	}

	string tiname = "Image_ .png";
	simgobj.resize(nlevels);
	for (UINT i = 0; i < nlevels; ++i) {
		Image timage = imgobj[i];
		timage.convolve(15, gkernel); // not happy with Gaussian Blur...
		tiname[tiname.length() - 5] = i + 'a';
		mImageType image2;
		ReadImage(imgobj[i], simgobj[i]);
		ReadImage(timage, image2);
		simgobj[i] -= factor * image2;
#ifdef __DEBUG
		//		simgobj[i] = image1 - factor * image2;
		cout << simgobj[i].rows() << " " << simgobj[i].columns() << endl
		<< flush;
		tiname[tiname.length() - 5] = i + 'A';
		//		WriteImage(tiname, simgobj[i]);
		WriteImage(tiname, simgobj[i]);
#endif
	}
}
// Image - = factor * Image2;
void Pyramid::SubtractImage(Image& image1, const Image& image2, REAL factor) {
	UINT m = image1.rows(), n = image1.columns();
	if (m != image2.rows() && n != image2.columns()) {
		cout << " Subtract Image: Image Dimensions doesn't match " << endl;
		return;
	}
	const PixelPacket *pix1 = image1.getConstPixels(0, 0, n, m), *pix2 =
	image2.getConstPixels(0, 0, n, m);
	PixelPacket tpix1, tpix2;
	UINT count = 0, tvar;
	for (UINT i = 0; i < m; ++i) {
		tvar = i * n;
		for (UINT j = 0; j < n; ++j, ++count) {
			//			tpix1 = pix1[tvar + j];
			//			tpix2 = pix2[tvar + j];
			//			REAL r = pix1[tvar + j].red - factor * pix2[tvar + j].red, g =
			//					pix1[tvar + j].green - factor * pix2[tvar + j].green, b =
			//					pix1[tvar + j].blue - factor * pix2[tvar + j].blue;
			ColorRGB pix2 = image2.pixelColor(j, i), pix1 = image1.pixelColor(
					j, i);
			pix1.red(pix1.red() - pix2.red() * factor);
			pix1.green(pix1.green() - pix2.green() * factor);
			pix1.blue(pix1.blue() - pix2.blue() * factor);
			image1.pixelColor(j, i, pix1);
			//			pix1 = pix1 - pix2;
			//			image1.pixelColor(j, i, image1.pixelColor(j, i) - pix2);
			//			ColorRGB diff = image1.pixelColor(j, i) - factor * image2.pixelColor(j,
			//					i);
			//			image1.pixelColor(j, i, ColorRGB(r / 65535, g / 65535, b / 65535));
		}
	}
	image1.syncPixels();
}

mImageType Pyramid::ScaleImage(Image &iImage) {

	if (ftype == Gaussian)
	SmoothGaussian(iImage, (UINT) (iImage.rows() / sratio),
			(UINT) (iImage.columns() / sratio), true);
	else
	SmoothBilinear(iImage, (UINT) (iImage.rows() / sratio),
			(UINT) (iImage.columns() / sratio), true);
	mImageType timage;
	ReadImage(iImage, timage);
	return timage;
}
mImageType Pyramid::ScaleAndSmoothImage(Image &iImage, REAL factor) {

	if (ftype == Gaussian)
	SmoothGaussian(iImage, (UINT) (iImage.rows() / sratio),
			(UINT) (iImage.columns() / sratio), true);
	else
	SmoothBilinear(iImage, (UINT) (iImage.rows() / sratio),
			(UINT) (iImage.columns() / sratio), true);
	Image bimage = iImage;
	mImageType timage, timage2;
	ReadImage(iImage, timage);
	bimage.convolve(15, gkernel); // convolve with Gaussain
	ReadImage(bimage, timage2);
	timage -= factor * timage2;
	return timage;
}

mImageType Pyramid::ScaleImage(UINT index, const Image &tmpImage) {
	Image iImage = tmpImage;
	if (ftype == Gaussian)
	SmoothGaussian(iImage, (UINT) (iImage.rows() / scale[index]),
			(UINT) (iImage.columns() / scale[index]), false);
	else
	SmoothBilinear(iImage, (UINT) (iImage.rows() / scale[index]),
			(UINT) (iImage.columns() / scale[index]), false);
	mImageType timage;
	ReadImage(iImage, timage);
	return timage;
}
mImageType Pyramid::ScaleAndSmoothImage(UINT index, const Image &tmpimage,
		REAL factor) {
	Image iImage = tmpimage;
#ifdef __DEBUG
	cout << " \n: Input Image  Rows = " << iImage.rows() << " Columns = "
	<< iImage.columns() << "\n : Scale = " << scale[index] << flush;
	if (ftype == Gaussian)
	SmoothGaussian(iImage, round(iImage.rows() / scale[index]), round(
					iImage.columns() / scale[index]), true);
	else
	SmoothBilinear(iImage, round(iImage.rows() / scale[index]),
			round(iImage.columns() / scale[index]), true);
#else
	if (ftype == Gaussian)
	SmoothGaussian(iImage, (UINT) (iImage.rows() / scale[index]),
			(UINT) (iImage.columns() / scale[index]), true);
	else
	SmoothBilinear(iImage, (UINT) (iImage.rows() / scale[index]),
			(UINT) (iImage.columns() / scale[index]), true);
#endif
	Image bimage = iImage;
	mImageType timage, timage2;
	ReadImage(iImage, timage);
	//	Write2File("bilin.txt", timage[0]);
#ifndef __WITHOUT_GAUSSIAN
	bimage.convolve(15, gkernel); // convolve with Gaussain
	ReadImage(bimage, timage2);
	timage -= factor * timage2;
#endif
	return timage;
}

mImageType Pyramid::ScaleAndSmoothImage2(UINT index, const Image &iimage,
		REAL factor) {
	/*	blitz::TinyVector<int, 2> extent(round(iimage.rows() / scale[index]),
	 round(iimage.columns() / scale[index]));
	 Image iImage = iimage;
	 SmoothBilinear(iImage, (UINT) (iImage.rows() / scale[index]),
	 (UINT) (iImage.columns() / scale[index]), true);
	 mImageType image, tmpimage(extent), snimage(extent);

	 ReadImage(iImage, image);
	 //	mImageType nimage = LJK::rescale(image, extent, LJK::BILINEAR_FILTER,
	 //			detail::ValueUnbounded()), tmpimage(extent), snimage(extent);
	 blitz::ConvolveExact<blitz::CPolicy_Stretch> conv;
	 conv.dim2(image, kernel, tmpimage);
	 conv.dim1(tmpimage, kernel, snimage);
	 //	Write2File("conv_res", res);
	 image -= factor * snimage;
	 return image;*/

	blitz::TinyVector<int, 2> extent(round(iimage.rows() / scale[index]),
			round(iimage.columns() / scale[index]));
	mImageType image;
	ReadImage(iimage, image);
	//	mImageType nimage = LJK::rescale(image, extent, LJK::BILINEAR_FILTER,
	//			detail::ValueUnbounded());
	mImageType nimage = LJK::fastrescale(image, 1 / scale[index]);
#ifndef __WITHOUT_GAUSSIAN
	mImageType tmpimage(extent), snimage(extent);
	blitz::ConvolveExact<blitz::CPolicy_Stretch> conv;
	conv.dim2(nimage, kernel, tmpimage);
	conv.dim1(tmpimage, kernel, snimage);
	nimage -= factor * snimage;
#endif
	//	Write2File("conv_res", res);

	return nimage;
}
#endif
