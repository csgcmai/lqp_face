/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/
#include "processImage.h"

void ProcessImage::ProcessSaveRes(Image &image, vector<REAL>&r,
		vector<REAL>& g, vector<REAL>& b, const string&ofile) {// Work on each color channel RGB

	if (gamma == 1 && !dog)
		ReadImage(image, r, g, b);

	if (gamma != 1)// Perform the gamma normalization of the image...
	// image magick gamma is x=y^(1/gamma); so sending
	// 1/gamma leads to x=y^gamma;
	// to have tip style normalization...
	{
		image.gamma(1 / gamma);
		if (!dog)
			ReadImage(image, r, g, b);
	}
	if (dog) {
		// Now using ImageMagick Implementation for padding
		UINT k2size = ceil(k2 * 3);
		Image image2 = image;
		// 3 to include 3 sigma bell
		image.blur(ceil(k1 * 3), k1);// applies separable filters...
		image2.blur(k2size, k2);
		SubtractImage(image, image2, r, g, b);
	}
	WriteImageFile(ofile + "_dog", r, g, b, image.rows(), image.columns());
	if (constr || coneq) {
		REAL imgstat[6];
		ExtractChannelStats(r, g, b, imgstat);
		if (constr)
			DoContrastStretching(r, g, b, imgstat);
		if (coneq)
			DoContrastEqualization(r, g, b, imgstat);
	}

}
/*Returns same pixel values for red, green & blue channels if a single channel information is demanded.
 * Thus same computational methods (atleast for Avg, L1 & L1Sqrt Norms)can be used for the computation of features*/
void ProcessImage::Process(Image &image, vector<REAL>&r, vector<REAL>& g,
		vector<REAL>& b) {
	switch (channel) {
	case WandellOCS:
	case SandeOCS:
		ProcessOCS(image, r, g, b);
		break;
	case HSL:
		ProcessHSL(image, r, g, b);
		break;
	case RGB:
		ProcessRGB(image, r, g, b);
		break;
	case Gray:
		ProcessGray(image, r);
		copy(r.begin(), r.end(), g.begin());
		copy(r.begin(), r.end(), b.begin());
		break;
	case Red:
	case Green:
	case Blue:
	case ABSGradientMag:
		Process(image, r);
		copy(r.begin(), r.end(), g.begin());
		copy(r.begin(), r.end(), b.begin());
		break;
	}
}
void ProcessImage::ProcessRGB(Image &image, vector<REAL>&r, vector<REAL>& g,
		vector<REAL>& b) {// Work on each color channel RGB

	if (gamma == 1 && !dog) {
		ReadImage(image, r, g, b);
	}
	if (gamma != 1)// Perform the gamma normalization of the image...
	// image magick gamma is x=y^(1/gamma); so sending
	// 1/gamma leads to x=y^gamma;
	// to have tip style normalization...
	{
		image.gamma(1 / gamma);
		if (!dog)
			ReadImage(image, r, g, b);
	}
	if (dog) {
		// Now using ImageMagick Implementation for padding
		UINT k2size = ceil(k2 * 3);
		Image image2 = image;
		// 3 to include 3 sigma bell
		image.blur(ceil(k1 * 3), k1);// applies separable filters...
		image2.blur(k2size, k2);
		SubtractImage(image, image2, r, g, b);
	}
	if (constr || coneq) {
		REAL imgstat[6];
		ExtractChannelStats(r, g, b, imgstat);
		if (constr)
			DoContrastStretching(r, g, b, imgstat);
		if (coneq)
			DoContrastEqualization(r, g, b, imgstat);
	}

}
void ProcessImage::ProcessOCS(Image &image, vector<REAL>&rpixval,
		vector<REAL>& gpixval, vector<REAL>& bpixval) {// Work on each color channel HSL
	// Currently only reading the OCS 1 & OCS 2 channels
	REAL wtm[3][3] = { { 0.2661, 0.6019, 0.0006 },
			{ -0.1250, 0.0377, -0.1332 }, { -0.0803, -0.3315, 0.4490 } },
			stm[3][3] = { { 0.7071, -0.7071, 0 }, { 0.4082, 0.4082, -0.8165 },
					{ 0.5774, 0.5774, 0.5774 } };
	REAL (*tm)[3];
	if (channel == WandellOCS)
		tm = wtm;
	else
		tm = stm;
	UINT n = image.columns(), m = image.rows();
	const PixelPacket *pix1 = image.getConstPixels(0, 0, n, m);
	PixelPacket tpix;
	UINT count = 0, tvar;
	REAL factor = 1.0 / MaxRGB;
	for (UINT i = 0; i < m; ++i) {
		tvar = i * n;
		for (UINT j = 0; j < n; ++j, ++count) {
			tpix = pix1[tvar + j];//0.2814 0.6938 0.0638

			rpixval[count] = (tpix.red * tm[0][0] + tpix.green * tm[0][1]
					+ tpix.blue * tm[0][2]) * factor;
			gpixval[count] = (tpix.red * tm[1][0] + tpix.green * tm[1][1]
					+ tpix.blue * tm[1][2]) * factor;
			bpixval[count] = (tpix.red * tm[2][0] + tpix.green * tm[2][1]
					+ tpix.blue * tm[2][2]) * factor;
		}
	}
}
// for specific Channel..
void ProcessImage::Process(Image &image, vector<REAL>&des1) {

	if (gamma == 1 && !dog)
		ReadImage(image, des1);

	if (gamma != 1)// Perform the gamma normalization of the image...
	{
		image.gamma(1 / gamma);
		if (!dog) {
			ReadImage(image, des1);
		}
	}
	if (dog) {
		// Now using ImageMagick Implementation for padding
		UINT k2size = ceil(k2 * 3);
		Image image2 = image;
		// 2 to include 2 sigma bell
		image.blur(ceil(k1 * 3), k1);// applies separable filters...
		image2.blur(k2size, k2);
		SubtractImage(image, image2, des1);
	}
	if (constr || coneq) {
		REAL imgstat[2];
		ExtractChannelStats(des1, imgstat);
		if (constr)
			DoContrastStretching(des1, imgstat);
		if (coneq)
			DoContrastEqualization(des1, imgstat);
	}

}
void ProcessImage::ProcessHSL(Image &image, vector<REAL>&h, vector<REAL>& s,
		vector<REAL>& l) {// Work on each color channel RGB

	if (gamma == 1 && !dog)
		ReadHSLImage(image, h, s, l);

	if (gamma != 1)// Perform the gamma normalization of the image...
	// image magick gamma is x=y^1/gamma; so sending
	// 1/gamma leads to x=y^gamma;
	// to have tip(Bill) style normalization...
	{
		image.gamma(1 / gamma);
		if (!dog)
			ReadHSLImage(image, h, s, l);
	}
	if (dog) {
		// Now using ImageMagick Implementation for padding
		UINT k2size = ceil(k2 * 3);
		Image image2 = image;
		// 2 to include 2 sigma bell
		image.blur(ceil(k1 * 3), k1);// applies separable filters...
		image2.blur(k2size, k2);
		SubtractImageHSL(image, image2, h, s, l);
	}
	if (constr || coneq) {
		REAL imgstat[6];
		ExtractChannelStats(h, s, l, imgstat);
		if (constr)
			DoContrastStretching(h, s, l, imgstat);
		if (coneq)
			DoContrastEqualization(h, s, l, imgstat);
	}
}

void ProcessImage::ProcessGray(Image &image, vector<REAL>&des1) {
	REAL imgstat[2];
	if (gamma == 1 && !dog)
		ReadImage(image, des1);

	if (gamma != 1)// Perform the gamma normalization of the image...
	{
		image.gamma(1 / gamma);
		if (!dog) {
			ReadImage(image, des1);
		}
	}
	if (dog) {
		// Now using ImageMagick Implementation for padding
		UINT k2size = ceil(k2 * 3);
		Image image2 = image;
		// 2 to include 2 sigma bell
		image.blur(ceil(k1 * 3), k1);// applies separable filters...
		image2.blur(k2size, k2);
		SubtractImage(image, image2, des1);
	}
	if (constr || coneq) {
		REAL imgstat[2];
		ExtractChannelStats(des1, imgstat);
		if (constr)
			DoContrastStretching(des1, imgstat);
		if (coneq)
			DoContrastEqualization(des1, imgstat);
	}
}
void ProcessImage::Normalize(vector<REAL>& ivector, REAL factor) {
	/*Multiply the Given input Image with the factor value*/
	for (vector<REAL>::iterator iter = ivector.begin(); iter != ivector.end(); ++iter)
		*iter *= factor;
}
void ProcessImage::Normalize(vector<REAL>& rch, vector<REAL>& gch,
		vector<REAL>& bch, REAL factor) {

	vector<REAL>::iterator riter = rch.begin(), giter = gch.begin(), biter =
			bch.begin();
	for (; riter != rch.end(); ++riter, ++giter, ++biter) {
		*riter *= factor;
		*giter *= factor;
		*biter *= factor;
	}
}
void ProcessImage::ReadImage(const Image& image, vector<REAL>&r,
		vector<REAL>&g, vector<REAL>&b) {
	UINT m = image.rows(), n = image.columns();
	const PixelPacket *pix1 = image.getConstPixels(0, 0, n, m);
	PixelPacket tpix;
	UINT count = 0, tvar;
	REAL factor = 1.0 / MaxRGB;
	for (UINT i = 0; i < m; ++i) {
		tvar = i * n;
		for (UINT j = 0; j < n; ++j, ++count) {
			tpix = pix1[tvar + j];
			r[count] = tpix.red * factor;
			g[count] = tpix.green * factor;
			b[count] = tpix.blue * factor;
		}
	}
}
/*void ProcessImage::ReadGrayImage(const Image& image, vector<REAL>&pixval) {
 UINT m = image.rows(), n = image.columns();
 ColorGray mycolor;
 UINT count = 0;
 for (UINT i = 0; i < m; ++i) {
 for (UINT j = 0; j < n; ++j) {
 ColorGray mycolor = image.pixelColor(j, i);
 pixval[count++] = mycolor.shade();
 }
 }
 }*/
void ProcessImage::ReadHSLImage(const Image &image, vector<REAL>&h,
		vector<REAL>&s, vector<REAL>&l) {
	int count = 0;
	for (unsigned int row = 0; row < image.rows(); row++) {
		for (unsigned int col = 0; col < image.columns(); col++, ++count) {
			Magick::Color pixCol = image.pixelColor(col, row);
			Magick::ColorHSL hslPixCol = Magick::ColorHSL(pixCol); // white: 0.0h 0.0s 1.0l, black: 0.0h 0.0s 0.0h
			h[count] = hslPixCol.hue();
			s[count] = hslPixCol.saturation();
			l[count] = hslPixCol.luminosity();
		}
	}
}

void ProcessImage::ReadImage(Image& image, vector<REAL>&pixval) {
	UINT m = image.rows(), n = image.columns();
	if (channel == Gray) {
		image.type(GrayscaleType);
		ColorGray mycolor;
		UINT count = 0;
		for (UINT i = 0; i < m; ++i) {
			for (UINT j = 0; j < n; ++j) {
				ColorGray mycolor = image.pixelColor(j, i);
				pixval[count++] = mycolor.shade();
			}
		}
	}
	// look into folder ~/ocs/
	// sRGB2XYZ ==> Standard RGB values to XYZ color Space, then we convert from the XYZ to opponent Color Space
	// Where O1== luminance component...
	// O2 is the red-green channel
	// O3 is the blue-yellow channel....
	else if (channel == O1 || channel == O2 || channel == O3) {
		//	Older Transformation matrix...
		//		REAL tm[3][3]={{0.2814, 0.6938, 0.0638},
		//		 {-0.0971 0.1458 -0.0250},
		//		{-0.0930 -0.2529 0.4665}};
		//
		// RGB2OCS rgb 2 opponent Color Space...
		REAL tm[3][3] = { { 0.2661, 0.6019, 0.0006 }, { -0.1250, 0.0377,
				-0.1332 }, { -0.0803, -0.3315, 0.4490 } };

		const PixelPacket *pix1 = image.getConstPixels(0, 0, n, m);
		PixelPacket tpix;
		UINT count = 0, tvar;
		REAL factor = 1.0 / MaxRGB;
		UINT row = channel == O1 ? 0 : channel == O2 ? 1 : 2;
		for (UINT i = 0; i < m; ++i) {
			tvar = i * n;
			for (UINT j = 0; j < n; ++j) {
				tpix = pix1[tvar + j];//0.2814 0.6938 0.0638

				pixval[count++] = (tpix.red * tm[row][0] + tpix.green
						* tm[row][1] + tpix.blue * tm[row][2]) * factor;
			}
		}

	} else if (channel == ABSGradientMag) {

		UINT width = image.columns(), height = image.rows();
		const PixelPacket *pix = image.getConstPixels(0, 0, width, height);
		REAL factor = 255.0 / (pow(2.0, 16.0) - 1);// scaling factor to have real values[0,1]
		UINT count = 0;
		for (int ty = 1; ty < height - 1; ty++) {
			for (int tx = 1; tx < width - 1; tx++) {
				//abs(g)=abs(gx)+abs(gy)
				// first color channel
				double dy = (pix[tx + (ty + 1) * width].red - pix[tx + (ty - 1)
						* width].red) * factor, dy2 = (pix[tx + (ty + 1)
						* width].green - pix[tx + (ty - 1) * width].green)
						* factor, dy3 = (pix[tx + (ty + 1) * width].blue
						- pix[tx + (ty - 1) * width].blue) * factor, dx =
						(pix[(tx + 1) + ty * width].red - pix[(tx - 1) + ty
								* width].red) * factor, dx2 = (pix[(tx + 1)
						+ ty * width].green - pix[(tx - 1) + ty * width].green)
						* factor, dx3 = (pix[(tx + 1) + ty * width].blue
						- pix[(tx - 1) + ty * width].blue) * factor, v =
						ABS(dx) + ABS(dy), v2 = ABS(dx2) + ABS(dy2), v3 =
						ABS(
								dx3) + ABS(dy3);
				// pick channel with strongest ABSOLUTE gradient
				if (v2 > v) {
					v = v2;
				}
				if (v3 > v) {
					v = v3;
				}
				pixval[count++] = v;
			}
		}

	} else {
		const PixelPacket *pix1 = image.getConstPixels(0, 0, n, m);
		UINT count = 0, tvar;
		REAL factor = 1.0 / MaxRGB;
		if (channel == Red) {
			for (UINT i = 0; i < m; ++i) {
				tvar = i * n;
				for (UINT j = 0; j < n; ++j)
					pixval[count++] = pix1[tvar + j].red * factor;
			}
		} else if (channel == Blue) {
			for (UINT i = 0; i < m; ++i) {
				tvar = i * n;
				for (UINT j = 0; j < n; ++j)
					pixval[count++] = pix1[tvar + j].blue * factor;
			}
		} else {
			for (UINT i = 0; i < m; ++i) {
				tvar = i * n;
				for (UINT j = 0; j < n; ++j)
					pixval[count++] = pix1[tvar + j].green * factor;
			}
		}
	}

}
void ProcessImage::SubtractImage(const Image& image1, const Image& image2,
		vector<REAL>&r, vector<REAL>&g, vector<REAL>&b) {
	UINT m = image1.rows(), n = image1.columns();
	if (m != image2.rows() && n != image2.columns()) {
		cout << " Subtract Image: Image Dimensions doesn't match " << endl;
		return;
	}
	ColorRGB mycolor1, mycolor2;
	const PixelPacket *pix1 = image1.getConstPixels(0, 0, n, m), *pix2 =
			image2.getConstPixels(0, 0, n, m);
	PixelPacket tpix1, tpix2;
	UINT count = 0, tvar;
	for (UINT i = 0; i < m; ++i) {
		tvar = i * n;
		for (UINT j = 0; j < n; ++j, ++count) {
			tpix1 = pix1[tvar + j];
			tpix2 = pix2[tvar + j];
			r[count] = ((REAL) (tpix1.red - tpix2.red)) / MaxRGB;
			g[count] = ((REAL) (tpix1.green - tpix2.green)) / MaxRGB;
			b[count] = ((REAL) (tpix1.blue - tpix2.blue)) / MaxRGB;
		}
	}
}
void ProcessImage::SubtractImageHSL(const Image& image1, const Image& image2,
		vector<REAL>&h, vector<REAL>&s, vector<REAL>&l) {
	int count = 0;
	Magick::ColorHSL hslPixCol1, hslPixCol2;
	for (unsigned int row = 0; row < image1.rows(); row++) {
		for (unsigned int col = 0; col < image1.columns(); col++, ++count) {
			hslPixCol1 = Magick::ColorHSL(image1.pixelColor(col, row));
			hslPixCol2 = Magick::ColorHSL(image2.pixelColor(col, row)); // white: 0.0h 0.0s 1.0l, black: 0.0h 0.0s 0.0h
			h[count] = hslPixCol1.hue() - hslPixCol2.hue();
			s[count] = hslPixCol1.saturation() - hslPixCol2.saturation();
			l[count] = hslPixCol1.luminosity() - hslPixCol2.luminosity();
		}
	}
}
void ProcessImage::ExtractChannelStats(vector<REAL>&gray, REAL imgstat[]) {
	imgstat[0] = -30000;
	imgstat[1] = +30000;
	for (vector<REAL>::iterator giter = gray.begin(); giter != gray.end(); ++giter) {
		if (*giter > imgstat[0])
			imgstat[0] = *giter;
		if (*giter < imgstat[1])
			imgstat[1] = *giter;
	}
}
void ProcessImage::ExtractChannelStats(vector<REAL>&r, vector<REAL>&g,
		vector<REAL>&b, REAL imgstat[]) {
	imgstat[0] = imgstat[2] = imgstat[4] = -30000;
	imgstat[1] = imgstat[3] = imgstat[5] = +30000;
	for (vector<REAL>::iterator riter = r.begin(), giter = g.begin(), biter =
			b.begin(); riter != r.end(); ++riter, ++giter, ++biter) {

		if (*riter > imgstat[0])
			imgstat[0] = *riter;
		if (*riter < imgstat[1])
			imgstat[1] = *riter;
		//green
		if (*giter > imgstat[2])
			imgstat[2] = *giter;
		if (*giter < imgstat[3])
			imgstat[3] = *giter;
		//blue
		if (*biter > imgstat[4])
			imgstat[4] = *biter;
		if (*biter < imgstat[5])
			imgstat[5] = *biter;

	}
}
void ProcessImage::SubtractImage(const Image& image1, const Image& image2,
		vector<REAL>&pixval) {
	UINT m = image1.rows(), n = image1.columns();
	if (m != image2.rows() && n != image2.columns()) {
		cout << " Subtract Image: Image Dimensions doesn't match " << endl;
		return;
	}

	if (channel == Gray) {
		ColorGray mycolor1, mycolor2;
		const PixelPacket *pix1 = image1.getConstPixels(0, 0, n, m), *pix2 =
				image2.getConstPixels(0, 0, n, m);
		UINT count = 0;
		for (UINT i = 0; i < m; ++i) {
			for (UINT j = 0; j < n; ++j) {
				mycolor1 = *(pix1 + i * n + j);
				mycolor2 = *(pix2 + i * n + j);
				pixval[count++] = mycolor1.shade() - mycolor2.shade();
			}
		}
	} else {
		ColorRGB mycolor1, mycolor2;
		const PixelPacket *pix1 = image1.getConstPixels(0, 0, n, m), *pix2 =
				image2.getConstPixels(0, 0, n, m);
		UINT count = 0, tvar;
		if (channel == Red) {
			for (UINT i = 0; i < m; ++i) {
				tvar = i * n;
				for (UINT j = 0; j < n; ++j)
					pixval[count++] = pix1[tvar + j].red - pix2[tvar + j].red;
			}
		} else if (channel == Green) {
			for (UINT i = 0; i < m; ++i) {
				tvar = i * n;
				for (UINT j = 0; j < n; ++j)
					pixval[count++] = pix1[tvar + j].green
							- pix2[tvar + j].green;
			}

		} else {
			for (UINT i = 0; i < m; ++i) {
				tvar = i * n;
				for (UINT j = 0; j < n; ++j)
					pixval[count++] = pix1[tvar + j].blue - pix2[tvar + j].blue;
			}
		}

	}
}
void ProcessImage::DoContrastStretching(vector<REAL>&gray, REAL imgstat[]) {

	REAL mult = ((maxr - minr) / (imgstat[0] - imgstat[1]));
	for (vector<REAL>::iterator giter = gray.begin(); giter != gray.end(); ++giter)
		*giter = (*giter - imgstat[1]) * mult + minr;

}
void ProcessImage::DoContrastStretching(vector<REAL>&r, vector<REAL>&g,
		vector<REAL>&b, REAL imgstat[]) {
	// TODO: Should Use the histogram base statistics...
	if (uchstat) // Use Channel Based Statistics for contrast Stretching
	{
		for (vector<REAL>::iterator riter = r.begin(), giter = g.begin(),
				biter = b.begin(); riter != r.end(); ++riter, ++giter, ++biter) {
			*riter = (*riter - imgstat[1]) * ((maxr - minr) / (imgstat[0]
					- imgstat[1])) + minr;
			*giter = (*giter - imgstat[3]) * ((maxr - minr) / (imgstat[2]
					- imgstat[3])) + minr;
			*biter = (*biter - imgstat[5]) * ((maxr - minr) / (imgstat[4]
					- imgstat[5])) + minr;
		}
	} else // Use Image base Statistic
	{
		REAL maxi, mini;
		maxi = *max_element(imgstat, imgstat + 6);
		mini = *min_element(imgstat, imgstat + 6);
		REAL mult = ((maxr - minr) / (maxi - mini));
		for (vector<REAL>::iterator riter = r.begin(), giter = g.begin(),
				biter = b.begin(); riter != r.end(); ++riter, ++giter, ++biter) {

			*riter = (*riter - mini) * mult + minr;
			*giter = (*giter - mini) * mult + minr;
			*biter = (*biter - mini) * mult + minr;
		}
	}
}
void ProcessImage::DoContrastEqualization(vector<REAL>&gray, REAL imgstat[]) {

	REAL mag = 0;
	for (vector<REAL>::iterator giter = gray.begin(); giter != gray.end(); ++giter) {
		mag += pow(ABS(*giter), alpha);
	}
	UINT len = gray.size();
	mag /= len;

	for (vector<REAL>::iterator giter = gray.begin(); giter != gray.end(); ++giter) {
		*giter /= mag;
	}
	mag = 0;
	for (vector<REAL>::iterator giter = gray.begin(); giter != gray.end(); ++giter) {
		mag += pow(MIN(tou, ABS(*giter)), alpha);
	}

	mag /= len;
	mag = pow(mag, 1 / alpha);

	for (vector<REAL>::iterator giter = gray.begin(); giter != gray.end(); ++giter) {
		*giter /= mag;
	}

	// Tanh transformation...% Optional...
	for (vector<REAL>::iterator giter = gray.begin(); giter != gray.end(); ++giter) {
		*giter = tou * tanh(*giter / tou);
	}
}
// Bill's Method for Contrast Stretching See TIP paper
void ProcessImage::DoContrastEqualization(vector<REAL>&r, vector<REAL>&g,
		vector<REAL>&b, REAL imgstat[]) {

	REAL mar = 0, mag = 0, mab = 0;
	for (vector<REAL>::iterator riter = r.begin(), giter = g.begin(), biter =
			b.begin(); riter != r.end(); ++riter, ++giter, ++biter) {
		mar += pow(ABS(*riter), alpha);
		mag += pow(ABS(*giter), alpha);
		mab += pow(ABS(*biter), alpha);
	}
	UINT len = r.size();
	mar /= len;
	mag /= len;
	mab /= len;

	for (vector<REAL>::iterator riter = r.begin(), giter = g.begin(), biter =
			b.begin(); riter != r.end(); ++riter, ++giter, ++biter) {
		*riter /= mar;
		*giter /= mag;
		*biter /= mar;
	}

	mar = mag = mab = 0;

	for (vector<REAL>::iterator riter = r.begin(), giter = g.begin(), biter =
			b.begin(); riter != r.end(); ++riter, ++giter, ++biter) {
		mar += pow(MIN(tou, ABS(*riter)), alpha);
		mag += pow(MIN(tou, ABS(*giter)), alpha);
		mab += pow(MIN(tou, ABS(*biter)), alpha);
	}

	mar /= len;
	mag /= len;
	mab /= len;

	mar = pow(mar, 1 / alpha);
	mag = pow(mag, 1 / alpha);
	mab = pow(mab, 1 / alpha);

	for (vector<REAL>::iterator riter = r.begin(), giter = g.begin(), biter =
			b.begin(); riter != r.end(); ++riter, ++giter, ++biter) {
		*riter /= mar;
		*giter /= mag;
		*biter /= mar;
	}

	// Tanh transformation...
	for (vector<REAL>::iterator riter = r.begin(), giter = g.begin(), biter =
			b.begin(); riter != r.end(); ++riter, ++giter, ++biter) {
		*riter = tou * tanh(*riter / tou);
		*giter = tou * tanh(*giter / tou);
		*biter = tou * tanh(*biter / tou);
	}
}
/*void ProcessImage::Gaussian(vector<REAL>& input,vector<REAL>&ker,
 UINT imwidth,UINT imheight,
 ,vector<REAL>&output)// already padded image as input..
 {
 UINT fsize = ker.size()/2;
 REAL sum =0;
 REAL *fil;
 // apply in the X-direction
 for(UINT i=0; i < imheight-2*fsize;++i)
 for(UINT j=0;j < imwidth-2*fsize;++j)
 {
 sum=0;
 fil = &input[i*imwidth+j];
 for(UINT k=0; k < ker.size();++k)
 sum+= ker[k] * fil[k];
 output[i*(imwidth-2*fsize)+j]=sum;
 }
 //apply in the Y-direction
 for(UINT i=0; i < imheight-2*fsize;++i)
 for(UINT j=0;j < imwidth-2*fsize;++j)
 {
 sum=0;
 fil = &input[i*imwidth+j];
 for(UINT k=0; k < ker.size();++k)
 sum+	= ker[k] * fil[k];
 output[i*(imwidth-2*fsize)+j]=sum;
 }
 }
 */
