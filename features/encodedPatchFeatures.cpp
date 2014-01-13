/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/
#define EPS 0.0001
#include "encodedPatchFeatures.h"
#ifdef STATIC_CODES
vector<vector<CodeVecType> > ComputeCode::codevec;
vector<UINT> ComputeCode::codeid;
bool ComputeCode::initcode = false; // by default no initalization
#endif
void EncodedPatchFeatures::InitalizeMaps(Image &image) {
	sspace = NoPyramid;
	UINT xbmax, ybmax;
	lbpfeatures.Initialize(1, 0, celldim);
	GetXYBlocks(image, xbmax, ybmax);
	lbpfeatures.SetFeature(xbmax, ybmax, 1.0);
	ComputeLBPFeatures(image, 0, cell, lbpstride);
}
void EncodedPatchFeatures::InitalizeMaps(Pyramid &pyobj, PyramidType sspace_) {
	UINT xbmax, ybmax, pycount = 0;
	sspace = sspace_;
	// Read the image....
	if (sspace == SingleRes) {
		pycount = pyobj.GetNLevels();
		lbpfeatures.Initialize(pycount, 0, celldim);
		for (UINT i = 0; i < pycount; i++) {
			Image& imgref = pyobj.GetImageRef(i);
			// Compute the Gradient Image
			GetXYBlocks(imgref, xbmax, ybmax);
			lbpfeatures.SetFeature(xbmax, ybmax, pyobj.GetScale(i));
			ComputeLBPFeatures(imgref, i, cell, lbpstride);
		}
	} else /*if (sspace == DoubleRes)*/{
		cout << "ERROR: Only LTP Not Implemented In Detail; ";
		exit(EXIT_FAILURE);

		/*		pycount = pyobj.GeneratePyramid(image);
		 UINT nlpyramid = pyobj.GetInitIndex();
		 lbpfeatures.Initialize(pycount + nlpyramid, nlpyramid, nbins * 2
		 * nthrlevels );
		 // Number of level in 2x...
		 for (UINT i = 0; i < pycount; i++) {
		 Image &imgref = pyobj.GetImageRef(i);

		 if (i < nlpyramid) {
		 Get2XYBlocks(imgref, xbmax, ybmax);
		 lbpfeatures.SetFeature(i, xbmax, ybmax, pyobj.GetScale(i) / 2);
		 ComputeLBPFeatures(imgref, i, cell / 2, lbpstride / 2);
		 }
		 // Compute the features...
		 //		imgref=pyobj.GetImageRef(i);
		 GetXYBlocks(imgref, xbmax, ybmax);
		 lbpfeatures.SetFeature(i + nlpyramid, xbmax, ybmax, pyobj.GetScale(
		 i));
		 ComputeLBPFeatures(imgref, i + nlpyramid, cell, lbpstride);
		 }*/
	}
}

void EncodedPatchFeatures::ExtractCodeInfo(Image&image, vector<Coord>& ainfo)
// annotation information
{
	// Crop the image according to each annotation information by adding some boundary
	// and then proceed in the normal way to generate the codebook, this generates
	// the code book using local scale space pyramid...
	if (sspace == SingleRes) { // build a scale-space for positive images as well...
		UINT imwidth = image.columns(), imheight = image.rows();
		int tx1, ty1, tx2, ty2, xmax, xmin, ymax, ymin, minarea = width
				* height;
		//	char tiname[350], anntext[35];
		Image timage = image;
		REAL pad;

		for (UINT i = 0; i < ainfo.size(); ++i) {
			if ((ainfo[i].GetWidth() + 1) * (ainfo[i].GetHeight() + 1)
					>= minarea) {
				// annotation must have atleast detection window area
				ainfo[i].GetCoord(xmin, ymin, xmax, ymax);
				pad = ((ymax - ymin + 1) + (xmax - xmin + 1)) * 0.5;
				tx1 = MAX(0, round(xmin - pad));
				ty1 = MAX(0, round(ymin - pad));
				tx2 = MIN(imwidth, round(xmax + pad));
				ty2 = MIN(imheight, round(ymax + pad));
				CropImage2(timage, tx1, ty1, tx2, ty2);
				ExtractCodeInfo(timage);
				timage.flop();
				ExtractCodeInfo(timage);
				timage = image;
			}
		}
	} else {/*Simply Crop the Positive Windows with the offset and generate code book*/
		UINT rarea = width * height;
		LBPMap map[3], nmap[3];
		int x1, y1, x2, y2, offset = GetMaxFeatureOffset();
		REAL xpad, ypad;
		Image im, oim;
		for (vector<Coord>::iterator iter = ainfo.begin(); iter != ainfo.end(); ++iter) {
			im = image;
			iter->GetCoord(x1, y1, x2, y2);
			if ((y2 - y1) * (x2 - x1) >= (long) rarea) {
				ypad = (((double) (y2 - y1)) / (height) * offset);
				xpad = (((double) (x2 - x1)) / (width) * offset);
				x1 = (int) round(x1 - xpad);
				y1 = (int) round(y1 - ypad);
				CropImage(im, x1, y1, round(x2 + xpad) - x1,
						round(y2 + ypad) - y1);
				ResizeImage(im, width + 2 * offset, height + 2 * offset,
						Bilinear, true); // 2 cell borders
				ExtractCodeBookInfo(im, map, nmap, true); // Flip the features once again so that in case of odd size it has the rig
				im.flop();
				ExtractCodeBookInfo(im, map, nmap, true);
			}
		}
	}
}
void EncodedPatchFeatures::ExtractCodeInfo(Image &image) {

	UINT pycount = pyobj.GeneratePyramid(image);
	LBPMap map[3], nmap[3];
	for (UINT i = 0; i < pycount; i++) {
		Image& imgref = pyobj.GetImageRef(i);
		// Compute the Gradient Image
		//		ComputeEncodedMap(imgref, map);
		ExtractCodeBookInfo(imgref, map, nmap, true);
	}
}
void EncodedPatchFeatures::ExtractNegCodeInfo(Image &image) {
	if (nwinsampled == 0) { // use complete image
		ExtractCodeInfo(image);
	} else { // else sample the windows
		UINT pycount = pyobj.GeneratePyramid(image);
		vector<WindowInfoFeat> winfo;/*Detector window's information*/
		// Generate the Pyramid and built the codeInfo of each level
		UINT wCount = 0;
		UINT ybmax, xbmax;
		REAL response = 0, overlap = 0, scale;
		LBPMap(*pmap)[3] = new LBPMap[pycount][3];
		LBPMap(*nmap)[3] = new LBPMap[pycount][3];
		UINT offset = GetMaxFeatureOffset(), minoffset = GetFeaturesOffset(),
				doffset = offset - minoffset;
		UINT ystride = ycell, xstride = xcell;

		for (UINT i = 0; i < pycount; i++) {
			Image &imgref = pyobj.GetImageRef(i);
			ExtractCodeBookInfo(imgref, pmap[i], nmap[i], false, true);
			wCount = 0;
			GetXYBlocks(imgref, xbmax, ybmax);
			ybmax = ybmax * cell - 2 * doffset;
			xbmax = xbmax * cell - 2 * doffset;
			if (ybmax >= height && xbmax >= width) {
				int trwindows = (int) floor((ybmax - height) / ystride) + 1,
						tcwindows = (int) floor((xbmax - width) / xstride) + 1;
				scale = pyobj.GetScale(i);
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
		vector<int> index(winfo.size(), 0);
		for (UINT i = 0; i < index.size(); ++i)
			index[i] = i;
		RandomPermute(index);
		int xmin, ymin, xmax, ymax, sindex;
		// Randomly Permute the windows...
		UINT thres = nwinsampled < winfo.size() ? nwinsampled : winfo.size();
		for (UINT iter = 0; iter < thres; ++iter) {
			winfo[index[iter]].GetCoord(xmin, ymin, xmax, ymax);
			sindex = winfo[index[iter]].GetScaleIndex();
			GenerateCodeBookInfo(pmap[sindex], nmap[sindex], xmin + doffset,
					ymin + doffset); // SameCodeBook For all the Cells
		}
		delete[] pmap;
		delete[] nmap;
	}
}
void EncodedPatchFeatures::GenerateCodeBookInfo(LBPMap *pmap, LBPMap *nmap,
		UINT xmin, UINT ymin) { ///Copy the Codes from temporary locations to the codeinfo and ncodeinfo objects...
	// for the speed modularization is sacrificed...
	UINT xbmax = xmin + width, ybmax = ymin + height, endx, endy, winstride =
			cell;
	int x, y, rval;
	for (UINT rind = ymin; rind < ybmax; rind += winstride) {
		y = rind + foffset - maxradius; // patch top left coordinates for which rind+foffset is center...
		// because the map contains an offset maxradius for concentic circles..
		// means map(0,0) contain value of map(maxradius,maxradius) so value of hist(offset,offset) can be found at
		// map(offset-maxradius,offset-maxradius)
		for (UINT colind = xmin; colind < xbmax; colind += winstride) {
			x = colind + foffset - maxradius;

			endx = MIN(x + winstride, pmap[0].nxpoints);
			endy = MIN(y + winstride, pmap[0].nypoints);

			for (UINT yiter = y; yiter < endy; ++yiter) // replace ++yiter with +=patch stride...
				for (UINT xiter = x; xiter < endx; ++xiter) {
					rval = pmap[0].GetValue(xiter, yiter);
					if (rval >= 0) { // to omit the boundary values...
						codeinfo->SetValue(rval,
								pmap[1].GetValue(xiter, yiter),
								pmap[2].GetValue(xiter, yiter));
						if (ncodeinfo)
							ncodeinfo->SetValue(nmap[0].GetValue(xiter, yiter),
									nmap[1].GetValue(xiter, yiter),
									nmap[2].GetValue(xiter, yiter));
					}
				}

		}
	}
}
void EncodedPatchFeatures::InitalizeMaps(Image&image, PyramidType sspace_) {
	UINT xbmax, ybmax, pycount = 0;
	sspace = sspace_;
	// Read the image....
	if (sspace == SingleRes) {
		pycount = pyobj.GeneratePyramid(image);
		lbpfeatures.Initialize(pycount, 0, celldim);
		for (UINT i = 0; i < pycount; i++) {
			Image& imgref = pyobj.GetImageRef(i);
			// Compute the Gradient Image
			GetXYBlocks(imgref, xbmax, ybmax);
			lbpfeatures.SetFeature(xbmax, ybmax, pyobj.GetScale(i));
			ComputeLBPFeatures(imgref, i, cell, lbpstride);
		}
	} else if (sspace == DoubleRes) {
		cout << "ERROR: Only LTP Not Implemented In Detail; ";
		exit(EXIT_FAILURE);

		/*		pycount = pyobj.GeneratePyramid(image);
		 UINT nlpyramid = pyobj.GetInitIndex();
		 lbpfeatures.Initialize(pycount + nlpyramid, nlpyramid, celldim);
		 // Number of level in 2x...
		 for (UINT i = 0; i < pycount; i++) {
		 Image &imgref = pyobj.GetImageRef(i);

		 if (i < nlpyramid) {
		 Get2XYBlocks(imgref, xbmax, ybmax);
		 lbpfeatures.SetFeature(i, xbmax, ybmax, pyobj.GetScale(i) / 2);
		 ComputeLBPFeatures(imgref, i, cell / 2, lbpstride / 2);
		 }
		 // Compute the features...
		 //		imgref=pyobj.GetImageRef(i);
		 GetXYBlocks(imgref, xbmax, ybmax);
		 lbpfeatures.SetFeature(i + nlpyramid, xbmax, ybmax, pyobj.GetScale(
		 i));
		 ComputeLBPFeatures(imgref, i + nlpyramid, cell, lbpstride);
		 }*/
	} else {
		lbpfeatures.Initialize(1, 0, celldim);
		GetXYBlocks(image, xbmax, ybmax);
		lbpfeatures.SetFeature(xbmax, ybmax, 1.0);
		ComputeLBPFeatures(image, 0, cell, lbpstride);
	}
}
void EncodedPatchFeatures::ExtractHaarCodeBookInfo(Image &image, LBPMap *map,
		LBPMap *nmap, const bool &gencodebook, bool onlycomputecodes) {
	UINT xdim = image.columns(), ydim = image.rows(), rcodeval, gcodeval,
			bcodeval, nrcodeval, ngcodeval, nbcodeval, reg = 1;
	int nypoints = ydim - 2 * maxradius, nxpoints = xdim - 2 * maxradius, cenx,
			ceny;
	vector<REAL> rimpix(xdim * ydim, 0), gimpix(xdim * ydim, 0), bimpix(
			xdim * ydim, 0);
	pim.Process(image, rimpix, gimpix, bimpix);
	REAL *impix[3];
	impix[0] = &rimpix[0];
	impix[1] = &gimpix[0];
	impix[2] = &bimpix[0];
	if (!gencodebook) // extract the map info
		for (int i = 0; i < 3; ++i) // for 3 channels RGB
		{
			map[i].Init(nxpoints, nypoints, -1);
			nmap[i].Init(nxpoints, nypoints, -1);
		}

	double featval[15] = { 0 }, codeval[3] = { 0 };
	for (int i = 0; i < nypoints; ++i) { // patch top-left coordinates
		for (int j = 0; j < nxpoints; ++j) {

			cenx = j + maxradius; // radius[k], maxradius
			ceny = i + maxradius;
			// computing 15 filter values...
			//				fill(featval, featval + 45, 0);
			codeval[0] = codeval[1] = codeval[2] = 0;
			for (int chcount = 0; chcount < 3; ++chcount) {

				//					for (int l = -2; l < 2; ++l)// average of 4x4 neighbourhood
				//						for (int k = -2; k < 2; ++k)
				//							featval[0] += impix[chcount][(ceny + l)
				//									* xdim + (cenx + k)];
				//					featval[chcount][0] /= 16;
				double b1, b2, b3, b4;
				// horizontal
				featval[0] = impix[chcount][(ceny - 2) * xdim + (cenx - 2)]
						+ impix[chcount][(ceny - 2) * xdim + (cenx - 1)]
						- impix[chcount][(ceny - 1) * xdim + (cenx - 2)]
						- impix[chcount][(ceny - 1) * xdim + (cenx - 1)]; //+1 +1 ;-1 -1 : top left

				b1 = (impix[chcount][(ceny - 2) * xdim + (cenx - 2)]
						+ impix[chcount][(ceny - 2) * xdim + (cenx - 1)]
						+ impix[chcount][(ceny - 1) * xdim + (cenx - 2)]
						+ impix[chcount][(ceny - 1) * xdim + (cenx - 1)]) / 4;

				featval[1] = impix[chcount][ceny * xdim + (cenx - 2)]
						+ impix[chcount][ceny * xdim + (cenx - 1)]
						- impix[chcount][(ceny + 1) * xdim + (cenx - 2)]
						- impix[chcount][(ceny + 1) * xdim + (cenx - 1)]; //+1 +1 ;-1 -1 : bottom left

				b3 = (impix[chcount][ceny * xdim + (cenx - 2)]
						+ impix[chcount][ceny * xdim + (cenx - 1)]
						+ impix[chcount][(ceny + 1) * xdim + (cenx - 2)]
						+ impix[chcount][(ceny + 1) * xdim + (cenx - 1)]) / 4;

				featval[2] = impix[chcount][(ceny - 2) * xdim + cenx]
						+ impix[chcount][(ceny - 2) * xdim + (cenx + 1)]
						- impix[chcount][(ceny - 1) * xdim + cenx]
						- impix[chcount][(ceny - 1) * xdim + (cenx + 1)]; //top right

				b2 = (impix[chcount][(ceny - 2) * xdim + cenx]
						+ impix[chcount][(ceny - 2) * xdim + (cenx + 1)]
						+ impix[chcount][(ceny - 1) * xdim + cenx]
						+ impix[chcount][(ceny - 1) * xdim + (cenx + 1)]) / 4;

				featval[3] = impix[chcount][ceny * xdim + cenx]
						+ impix[chcount][ceny * xdim + (cenx + 1)]
						- impix[chcount][(ceny + 1) * xdim + cenx]
						- impix[chcount][(ceny + 1) * xdim + (cenx + 1)]; //+1 +1 ;-1 -1 : bottom right

				b4 = (impix[chcount][ceny * xdim + cenx] + impix[chcount][ceny
						* xdim + (cenx + 1)] + impix[chcount][(ceny + 1) * xdim
						+ cenx]
						+ impix[chcount][(ceny + 1) * xdim + (cenx + 1)]) / 4;

				// vertical
				featval[4] = impix[chcount][(ceny - 2) * xdim + (cenx - 2)]
						- impix[chcount][(ceny - 2) * xdim + (cenx - 1)]
						+ impix[chcount][(ceny - 1) * xdim + (cenx - 2)]
						- impix[chcount][(ceny - 1) * xdim + (cenx - 1)]; //+1 +1 ;-1 -1 : top left

				featval[5] = impix[chcount][ceny * xdim + (cenx - 2)]
						- impix[chcount][ceny * xdim + (cenx - 1)]
						+ impix[chcount][(ceny + 1) * xdim + (cenx - 2)]
						- impix[chcount][(ceny + 1) * xdim + (cenx - 1)]; //+1 +1 ;-1 -1 : bottom left

				featval[6] = impix[chcount][(ceny - 2) * xdim + cenx]
						- impix[chcount][(ceny - 2) * xdim + (cenx + 1)]
						+ impix[chcount][(ceny - 1) * xdim + cenx]
						- impix[chcount][(ceny - 1) * xdim + (cenx + 1)]; //top right

				featval[7] = impix[chcount][ceny * xdim + cenx]
						- impix[chcount][ceny * xdim + (cenx + 1)]
						+ impix[chcount][(ceny + 1) * xdim + cenx]
						- impix[chcount][(ceny + 1) * xdim + (cenx + 1)]; //+1 +1 ;-1 -1 : bottom right

				// diagonal

				featval[8] = impix[chcount][(ceny - 2) * xdim + (cenx - 2)]
						- impix[chcount][(ceny - 2) * xdim + (cenx - 1)]
						- impix[chcount][(ceny - 1) * xdim + (cenx - 2)]
						+ impix[chcount][(ceny - 1) * xdim + (cenx - 1)]; //+1 +1 ;-1 -1 : top left

				featval[9] = impix[chcount][ceny * xdim + (cenx - 2)]
						- impix[chcount][ceny * xdim + (cenx - 1)]
						- impix[chcount][(ceny + 1) * xdim + (cenx - 2)]
						+ impix[chcount][(ceny + 1) * xdim + (cenx - 1)]; //+1 +1 ;-1 -1 : bottom left

				featval[10] = impix[chcount][(ceny - 2) * xdim + cenx]
						- impix[chcount][(ceny - 2) * xdim + (cenx + 1)]
						- impix[chcount][(ceny - 1) * xdim + cenx]
						+ impix[chcount][(ceny - 1) * xdim + (cenx + 1)]; //top right

				featval[11] = impix[chcount][ceny * xdim + cenx]
						- impix[chcount][ceny * xdim + (cenx + 1)]
						- impix[chcount][(ceny + 1) * xdim + cenx]
						+ impix[chcount][(ceny + 1) * xdim + (cenx + 1)]; //+1 +1 ;-1 -1 : bottom right

				// 2nd layer codes ...
				featval[12] = b1 + b2 - b3 - b4;
				featval[13] = b1 - b2 + b3 - b4;
				featval[14] = b1 - b2 - b3 + b4;
				if (petype == PET_LTP) {
					// generate ternary code
					for (int count = 0; count < 15; ++count) {
						UINT treg = 0;
						// posvalue and negvalue contains the thresholded values for r,g & b
						// tcount is used to run over the values...
						for (int tl = 0; tl < ltplevels; ++tl) {
							if (featval[count] > tollevels[tl])
								treg = tl + 1;
							else if (featval[count] < -tollevels[tl])
								treg = tl + 1 + ltplevels; // codes will be +ve while vector contains the -ve
						}

						codeval[chcount] += basevalues[count] * treg;

					}
				}
			}
			if (onlycomputecodes) {
				map[0](j, i) = (codeval[0]); // for r,g and b channels..
				map[1](j, i) = (codeval[1]);
				map[2](j, i) = (codeval[2]);
			} else if (gencodebook) {
				codeinfo->SetValue(codeval[0], codeval[1], codeval[2]);

			} else {

				if (!softq) {
					map[0](j, i) = codeinfo->GetClusterCenter(codeval[0]); // for r,g and b channels..
					map[1](j, i) = codeinfo->GetClusterCenter(codeval[1]);
					map[2](j, i) = codeinfo->GetClusterCenter(codeval[2]);

				} else {
					map[0](j, i) = (codeval[0]); // for r,g and b channels..
					map[1](j, i) = (codeval[1]);
					map[2](j, i) = (codeval[2]);
				}
			}

		}
	}
}
void EncodedPatchFeatures::ExtractCircAngRadCodeBookInfo(Image &image,
		LBPMap *map, LBPMap *nmap, const bool &gencodebook,
		bool computeonlycodes) {

	long xdim = image.columns(), ydim = image.rows(), rcodeval, gcodeval,
			bcodeval, nrcodeval, ngcodeval, nbcodeval, reg = 1;
	int nypoints = ydim - 2 * maxradius, fx, fy, cx, cy, nxpoints = xdim - 2
			* maxradius, cenx, ceny, count;
	REAL tmpx, tmpy, fracx, fracy, pixval[3], cenval[3];
	REAL x, y, w1, w2, w3, w4;
	vector<REAL> rimpix(xdim * ydim, 0), gimpix(xdim * ydim, 0), bimpix(
			xdim * ydim, 0);
	pim.Process(image, rimpix, gimpix, bimpix);
	if (!gencodebook) // extract the map info
		for (int i = 0; i < 3; ++i) // for 3 channels RGB
		{
			map[i].Init(nxpoints, nypoints, -1);
			nmap[i].Init(nxpoints, nypoints, -1);
		}
	vector<vector<REAL> > spvalues(3, vector<REAL> (npoints, 0));
	vector<REAL> pixelvalues(npoints * 3, 0);
	REAL *rpixval = &pixelvalues[0], *gpixval = rpixval + npoints, *bpixval =
			gpixval + npoints;
	UINT ppcircle = npoints / 3; ///> points per circle
	long reg2 = reg << ppcircle, reg3 = reg2 << ppcircle;
	///> use online mapping of codes,
	vector<REAL> npixelvalues(npoints * 3, 0);
	REAL *nrpixval = &npixelvalues[0], *ngpixval = nrpixval + npoints,
			*nbpixval = ngpixval + npoints;
	for (int i = 0; i < nypoints; i += patchstride) { // patch top-left coordinates
		for (int j = 0; j < nxpoints; j += patchstride) {
			rcodeval = gcodeval = bcodeval = 0;
			nrcodeval = ngcodeval = nbcodeval = 0;

			cenx = j + maxradius; // radius[k], maxradius
			ceny = i + maxradius;

			cenval[0] = rimpix[ceny * xdim + cenx];
			cenval[1] = gimpix[ceny * xdim + cenx];
			cenval[2] = bimpix[ceny * xdim + cenx];
			count = 0;
			///> use online mapping of codes,
			fill(pixelvalues.begin(), pixelvalues.end(), 0);
			fill(npixelvalues.begin(), npixelvalues.end(), 0);

			REAL mean[3] = { 0, 0, 0 };
			for (vector<vector<Point<REAL> > >::iterator oiter =
					spoints.begin(); oiter != spoints.end(); ++oiter)
				for (vector<Point<REAL> >::iterator iter = oiter->begin(); iter
						!= oiter->end(); ++iter, ++count) {
					iter->GetCoord(x, y);
					tmpy = ceny + y;
					tmpx = cenx + x;

					/// do bilinear interpolation

					fx = floor(tmpx);
					fy = floor(tmpy);
					cx = ceil(tmpx);
					cy = ceil(tmpy);
					fracx = tmpx - fx;
					fracy = tmpy - fy;

					if (ABS(fracx) < 1e-6 && ABS(fracy) < 1e-6) {
						// no interpolation needed
						spvalues[0][count] = rimpix[fy * xdim + fx];
						spvalues[1][count] = gimpix[fy * xdim + fx];
						spvalues[2][count] = bimpix[fy * xdim + fx];
					} else {

						w1 = (1 - fracx) * (1 - fracy);
						w2 = fracx * (1 - fracy);
						w3 = (1 - fracx) * fracy;
						w4 = fracx * fracy;

						spvalues[0][count] = w1 * rimpix[fy * xdim + fx] + w2
								* rimpix[fy * xdim + cx] + w3 * rimpix[cy
								* xdim + fx] + w4 * rimpix[cy * xdim + cx];

						spvalues[1][count] = w1 * gimpix[fy * xdim + fx] + w2
								* gimpix[fy * xdim + cx] + w3 * gimpix[cy
								* xdim + fx] + w4 * gimpix[cy * xdim + cx];

						spvalues[2][count] = w1 * bimpix[fy * xdim + fx] + w2
								* bimpix[fy * xdim + cx] + w3 * bimpix[cy
								* xdim + fx] + w4 * bimpix[cy * xdim + cx];
					}
					// posvalue and negvalue contains the thresholded values for r,g & b
					// tcount is used to run over the values...
					if (oiter == spoints.begin()) { ///> mean values for r,g,b points lying on the circle.
						mean[0] += spvalues[0][count];
						mean[1] += spvalues[1][count];
						mean[2] += spvalues[2][count];
					}
				}

			mean[0] /= ppcircle;
			mean[1] /= ppcircle;
			mean[2] /= ppcircle;

			/// Compare with the central pixel values instead of mean...
			mean[0] = cenval[0];
			mean[1] = cenval[1];
			mean[2] = cenval[2];
			count = 0;
			/// compare with the central pixel values...
			for (; count < ppcircle; ++count) {
				int tind = (count + 1) % ppcircle, ioffset = 0;
				long tval = reg << count, tval2 = reg2 << count, tval3 = reg3
						<< count;
				if (petype == PET_LBP) {
					if (spvalues[0][count] >= mean[0])
						rcodeval += tval;

					if (spvalues[1][count] >= mean[1])
						gcodeval += tval;

					if (spvalues[2][count] >= mean[2])
						bcodeval += tval;
					///
					if (spvalues[0][count] >= spvalues[0][count + ppcircle])
						rcodeval |= tval2;

					if (spvalues[1][count] >= spvalues[1][count + ppcircle])
						gcodeval |= tval2;

					if (spvalues[2][count] >= spvalues[2][count + ppcircle])
						bcodeval |= tval2;

					if (spvalues[0][count] >= spvalues[0][tind])
						rcodeval |= tval3;

					if (spvalues[1][count] >= spvalues[1][tind])
						gcodeval |= tval3;

					if (spvalues[2][count] >= spvalues[2][tind])
						bcodeval |= tval3;

				} else if (petype == PET_SplitLTP) {
					ioffset = npoints - count - 1;
					if (spvalues[0][count] > mean[0] + tollevels[0]) {
						rpixval[ioffset] = 1;
						rcodeval += tval;
					} else if (spvalues[0][count] < mean[0] - tollevels[0]) {
						nrpixval[ioffset] = 1;
						nrcodeval += tval;
					}

					if (spvalues[1][count] > mean[1] + tollevels[0]) {
						gpixval[ioffset] = 1;
						gcodeval += tval;
					} else if (spvalues[1][count] < mean[1] - tollevels[0]) {
						ngpixval[ioffset] = 1;
						ngcodeval += tval;
					}

					if (spvalues[2][count] > mean[2] + tollevels[0]) {
						bpixval[ioffset] = 1;
						bcodeval += tval;
					} else if (spvalues[2][count] < mean[2] - tollevels[0]) {
						nbpixval[ioffset] = 1;
						nbcodeval += tval;
					}

					///> comparison with the innner circle values
					ioffset = npoints - count - ppcircle - 1;
					if (spvalues[0][count] > spvalues[0][count + ppcircle]
							+ tollevels[0]) {
						rpixval[ioffset] = 1;
						rcodeval |= tval2;
					} else if (spvalues[0][count] < spvalues[0][count
							+ ppcircle] - tollevels[0]) {
						nrpixval[ioffset] = 1;
						nrcodeval |= tval2;
					}

					if (spvalues[1][count] > spvalues[1][count + ppcircle]
							+ tollevels[0]) {
						gpixval[ioffset] = 1;
						gcodeval |= tval2;
					} else if (spvalues[1][count] < spvalues[1][count
							+ ppcircle] - tollevels[0]) {
						ngpixval[ioffset] = 1;
						ngcodeval |= tval2;
					}

					if (spvalues[2][count] > spvalues[2][count + ppcircle]
							+ tollevels[0]) {
						bpixval[ioffset] = 1;
						bcodeval |= tval2;
					} else if (spvalues[2][count] < spvalues[2][count
							+ ppcircle] - tollevels[0]) {
						nbpixval[ioffset] = 1;
						nbcodeval |= tval2;
					}
					ioffset = npoints - count - 2 * ppcircle - 1;
					///> angular comparison among the pixel values
					if (spvalues[0][count] > spvalues[0][tind] + tollevels[0]) {
						rpixval[ioffset] = 1;
						rcodeval |= tval3;
					} else if (spvalues[0][count] < spvalues[0][tind]
							- tollevels[0]) {
						nrpixval[ioffset] = 1;
						nrcodeval |= tval3;
					}

					if (spvalues[1][count] > spvalues[1][tind] + tollevels[0]) {
						gpixval[ioffset] = 1;
						gcodeval |= tval3;
					} else if (spvalues[1][count] < spvalues[1][tind]
							- tollevels[0]) {
						ngpixval[ioffset] = 1;
						ngcodeval |= tval3;
					}

					if (spvalues[2][count] > spvalues[2][tind] + tollevels[0]) {
						bpixval[ioffset] = 1;
						bcodeval |= tval3;
					} else if (spvalues[2][count] < spvalues[2][tind]
							- tollevels[0]) {
						nbpixval[ioffset] = 1;
						nbcodeval |= tval3;
					}
				}
			}

			if (computeonlycodes) {
				map[0](j, i) = (rcodeval); // for r,g and b channels..
				map[1](j, i) = (gcodeval);
				map[2](j, i) = (bcodeval);

				if (ncodeinfo) {
					nmap[0](j, i) = (nrcodeval); // for r,g and b channels..
					nmap[1](j, i) = (ngcodeval);
					nmap[2](j, i) = (nbcodeval);
				}

			} else if (gencodebook) {
#ifdef ONLINE_MAPPING
				codeinfo->SetValue(pixelvalues, rcodeval, gcodeval, bcodeval);
				if (ncodeinfo)
				ncodeinfo->SetValue(npixelvalues, nrcodeval, ngcodeval,
						nbcodeval);
#else
				codeinfo->SetValue(rcodeval, gcodeval, bcodeval);
				if (ncodeinfo)
					ncodeinfo->SetValue(nrcodeval, ngcodeval, nbcodeval);
#endif
			} else {
				if (!softq) { // map to codebook centers...
#ifdef ONLINE_MAPPING
					map[0](j, i)
					= codeinfo->GetClusterCenter(rpixval, rcodeval); // for r,g and b channels..
					map[1](j, i)
					= codeinfo->GetClusterCenter(gpixval, gcodeval);
					map[2](j, i)
					= codeinfo->GetClusterCenter(bpixval, bcodeval);

					if (ncodeinfo) {
						nmap[0](j, i) = ncodeinfo->GetClusterCenter(nrpixval,
								nrcodeval); // for r,g and b channels..
						nmap[1](j, i) = ncodeinfo->GetClusterCenter(ngpixval,
								ngcodeval);
						nmap[2](j, i) = ncodeinfo->GetClusterCenter(nbpixval,
								nbcodeval);
					}
#else
					map[0](j, i) = codeinfo->GetClusterCenter(rcodeval); // for r,g and b channels..
					map[1](j, i) = codeinfo->GetClusterCenter(gcodeval);
					map[2](j, i) = codeinfo->GetClusterCenter(bcodeval);

					if (ncodeinfo) {
						nmap[0](j, i) = ncodeinfo->GetClusterCenter(nrcodeval); // for r,g and b channels..
						nmap[1](j, i) = ncodeinfo->GetClusterCenter(ngcodeval);
						nmap[2](j, i) = ncodeinfo->GetClusterCenter(nbcodeval);
					}
#endif
				} else {
					map[0](j, i) = (rcodeval); // for r,g and b channels..
					map[1](j, i) = (gcodeval);
					map[2](j, i) = (bcodeval);

					if (ncodeinfo) {
						nmap[0](j, i) = (nrcodeval); // for r,g and b channels..
						nmap[1](j, i) = (ngcodeval);
						nmap[2](j, i) = (nbcodeval);
					}
				}
			}
		}
	}
}
void EncodedPatchFeatures::ExtractCodeBookInfo(Image &image, LBPMap *map,
		LBPMap *nmap, const bool &gencodebook, bool computeonlycodes) {

	if (ptype == PT_HAARPyramid) {
		ExtractHaarCodeBookInfo(image, map, nmap, gencodebook, computeonlycodes);
		return;
	}
	if (ptype == PT_CircAngRad) {
		ExtractCircAngRadCodeBookInfo(image, map, nmap, gencodebook,
				computeonlycodes);
		return;
	}
	ExtractCodeBookInfo(image, map, nmap, 0, 0, image.columns(), image.rows(),
			gencodebook, computeonlycodes);
}
void EncodedPatchFeatures::ExtractCodeBookInfo(Image &image, LBPMap *map,
		LBPMap *nmap, int sx, int sy, int ex, int ey, const bool &gencodebook,
		bool computeonlycodes) {

	long xdim = image.columns(), ydim = image.rows(), rcodeval, gcodeval,
			bcodeval, nrcodeval, ngcodeval, nbcodeval, reg = 1;
	int nypoints = ydim - 2 * maxradius, nxpoints = xdim - 2 * maxradius, cenx,
			ceny, count, fx, fy, cx, cy;
	REAL tmpx, tmpy, fracx, fracy, pixval[3], cenval[3];
	REAL x, y, w1, w2, w3, w4;
	vector<REAL> rimpix(xdim * ydim, 0), gimpix(xdim * ydim, 0), bimpix(
			xdim * ydim, 0);
	pim.Process(image, rimpix, gimpix, bimpix);
	if (!gencodebook) // extract the map info
		for (int i = 0; i < 3; ++i) // for 3 channels RGB
		{
			map[i].Init(nxpoints, nypoints, -1);
			nmap[i].Init(nxpoints, nypoints, -1);
		}
	vector<REAL> pixelvalues(npoints * 3, 0);
	REAL *rpixval = &pixelvalues[0], *gpixval = rpixval + npoints, *bpixval =
			gpixval + npoints;
	vector<REAL> npixelvalues(npoints * 3, 0);
	REAL *nrpixval = &npixelvalues[0], *ngpixval = nrpixval + npoints,
			*nbpixval = ngpixval + npoints;
	/*Coordinates configurations for a window codebook*/
	nypoints = MIN(ey-maxradius+1, ydim - 2 * maxradius);
	nxpoints = MIN(ex-maxradius+1, xdim - 2 * maxradius);
	sx = MAX(sx-(int)maxradius,0);
	sy = MAX(sy-(int)maxradius,0);
	for (int i = sy; i < nypoints; i += patchstride) { // patch top-left coordinates
		for (int j = sx; j < nxpoints; j += patchstride) {
			rcodeval = gcodeval = bcodeval = 0;
			nrcodeval = ngcodeval = nbcodeval = 0;

			cenx = j + maxradius; // radius[k], maxradius
			ceny = i + maxradius;

			cenval[0] = rimpix[ceny * xdim + cenx];
			cenval[1] = gimpix[ceny * xdim + cenx];
			cenval[2] = bimpix[ceny * xdim + cenx];
			count = 0;
			fill(pixelvalues.begin(), pixelvalues.end(), 0);
			fill(npixelvalues.begin(), npixelvalues.end(), 0);
			for (vector<vector<Point<REAL> > >::iterator oiter =
					spoints.begin(); oiter != spoints.end(); ++oiter)
				for (vector<Point<REAL> >::iterator iter = oiter->begin(); iter
						!= oiter->end(); ++iter, ++count) {
					iter->GetCoord(x, y);
					tmpy = ceny + y;
					tmpx = cenx + x;

					if ((ptype == PT_Circular && gridtype == Circular)
							|| (ptype == PT_HVCStrip || ptype == PT_DACStrip)
							|| ptype == PT_CircularShallow) { // do bilinear interpolation

						fx = floor(tmpx);
						fy = floor(tmpy);
						cx = ceil(tmpx);
						cy = ceil(tmpy);
						fracx = tmpx - fx;
						fracy = tmpy - fy;

						if (ABS(fracx) < 1e-6 && ABS(fracy) < 1e-6) {
							// no interpolation needed
							pixval[0] = rimpix[fy * xdim + fx];
							pixval[1] = gimpix[fy * xdim + fx];
							pixval[2] = bimpix[fy * xdim + fx];
						} else {

							w1 = (1 - fracx) * (1 - fracy);
							w2 = fracx * (1 - fracy);
							w3 = (1 - fracx) * fracy;
							w4 = fracx * fracy;

							pixval[0] = w1 * rimpix[fy * xdim + fx] + w2
									* rimpix[fy * xdim + cx] + w3 * rimpix[cy
									* xdim + fx] + w4 * rimpix[cy * xdim + cx];

							pixval[1] = w1 * gimpix[fy * xdim + fx] + w2
									* gimpix[fy * xdim + cx] + w3 * gimpix[cy
									* xdim + fx] + w4 * gimpix[cy * xdim + cx];

							pixval[2] = w1 * bimpix[fy * xdim + fx] + w2
									* bimpix[fy * xdim + cx] + w3 * bimpix[cy
									* xdim + fx] + w4 * bimpix[cy * xdim + cx];
						}
						// posvalue and negvalue contains the thresholded values for r,g & b
						// tcount is used to run over the values...
					}
					// replace with else statement.
					if ((ptype == PT_Circular && gridtype == Rectangular)
							|| (ptype != PT_Circular && ptype
									!= PT_CircularShallow)) {
						pixval[0] = rimpix[tmpy * xdim + tmpx];
						pixval[1] = gimpix[tmpy * xdim + tmpx];
						pixval[2] = bimpix[tmpy * xdim + tmpx];
					}
					if (petype == PET_LTP || petype == PET_CSLTP_LTP) {
						UINT treg[3] = { 0, 0, 0 };
						// posvalue and negvalue contains the thresholded values for r,g & b
						// tcount is used to run over the values...
						for (int tl = 0, tcount = 0; tl < ltplevels; ++tl)
							for (int iter = 0; iter < 3; ++iter, ++tcount) {
								if (pixval[iter] > cenval[iter] + tollevels[tl]) {
									treg[iter] = tl + 1;
								} else if (pixval[iter] < cenval[iter]
										- tollevels[tl]) {
									treg[iter] = tl + 1 + ltplevels; // codes will be +ve while vector contains the -ve
								}
							}
						rcodeval += basevalues[count] * treg[0];
						gcodeval += basevalues[count] * treg[1];
						bcodeval += basevalues[count] * treg[2];

						if (petype == PET_CSLTP_LTP) {
							rpixval[count] = pixval[0];
							gpixval[count] = pixval[1];
							bpixval[count] = pixval[2];
						}

					} else if (petype == PET_LBP) {
						UINT tval = reg << count;
						if (pixval[0] >= cenval[0])
							rcodeval += tval;

						if (pixval[1] >= cenval[1])
							gcodeval += tval;

						if (pixval[2] >= cenval[2])
							bcodeval += tval;
					} else if (petype == PET_SplitLTP) {
						long tval = reg << count;
						if (pixval[0] > cenval[0] + tollevels[0]) {
							rpixval[count] = 1;
							rcodeval += tval;
						} else if (pixval[0] < cenval[0] - tollevels[0]) {
							nrpixval[count] = 1;
							nrcodeval += tval;
						}

						if (pixval[1] > cenval[1] + tollevels[0]) {
							gpixval[count] = 1;
							gcodeval += tval;
						} else if (pixval[1] < cenval[1] - tollevels[0]) {
							ngpixval[count] = 1;
							ngcodeval += tval;
						}

						if (pixval[2] > cenval[2] + tollevels[0]) {
							bpixval[count] = 1;
							bcodeval += tval;
						} else if (pixval[2] < cenval[2] - tollevels[0]) {
							nbpixval[count] = 1;
							nbcodeval += tval;
						}

					} else {
						rpixval[count] = pixval[0];
						gpixval[count] = pixval[1];
						bpixval[count] = pixval[2];
					}
				}

			if (petype == PET_CSLTP_LTP || petype == PET_CSLTP) {
				UINT npatches = spoints.size(), nspoints = midpoint * 2,
						bvcount = petype == PET_CSLTP_LTP ? count : 0;
				for (UINT piter = 0; piter < npatches; ++piter) { // patch siter
					UINT siter = piter * nspoints, eiter = siter + nspoints - 1;
					UINT treg[3] = { 0, 0, 0 };
					for (; siter < eiter; ++siter, --eiter, ++bvcount) {
						for (int tl = 0; tl < ltplevels; ++tl) {
							if (rpixval[siter] > rpixval[eiter] + tollevels[tl])
								treg[0] = tl + 1;
							else if (rpixval[siter] < rpixval[eiter]
									- tollevels[tl])
								treg[0] = tl + ltplevels + 1;
							if (gpixval[siter] > gpixval[eiter] + tollevels[tl])
								treg[1] = tl + 1;
							else if (gpixval[siter] < gpixval[eiter]
									- tollevels[tl])
								treg[1] = tl + ltplevels + 1;
							if (bpixval[siter] > bpixval[eiter] + tollevels[tl])
								treg[2] = tl + 1;
							else if (bpixval[siter] < bpixval[eiter]
									- tollevels[tl])
								treg[2] = tl + ltplevels + 1;
						}
						rcodeval += basevalues[bvcount] * treg[0];
						gcodeval += basevalues[bvcount] * treg[1];
						bcodeval += basevalues[bvcount] * treg[2];
					}
				}
			} else if (petype == PET_CSLBP) {

				if (issigned) { /*by default use the unsigned for cslbp,ltp*/
					for (UINT iter = 0; iter < midpoint; ++iter) {
						UINT wt = reg << iter;
						rcodeval += ((rpixval[iter] - rpixval[iter + midpoint])
								> tollevels[0] ? wt : 0);
						gcodeval += ((gpixval[iter] - gpixval[iter + midpoint])
								> tollevels[0] ? wt : 0);
						bcodeval += ((bpixval[iter] - bpixval[iter + midpoint])
								> tollevels[0] ? wt : 0);
					}
				} else {
					for (UINT iter = 0; iter < midpoint; ++iter) {
						UINT wt = reg << iter;
						rcodeval
								+= (ABS(rpixval[iter] - rpixval[iter + midpoint])
										> tollevels[0] ? wt : 0);
						gcodeval
								+= (ABS(gpixval[iter] - gpixval[iter + midpoint])
										> tollevels[0] ? wt : 0);
						bcodeval
								+= (ABS(bpixval[iter] - bpixval[iter + midpoint])
										> tollevels[0] ? wt : 0);
					}
				}

			}
			if (computeonlycodes) {
				map[0](j, i) = (rcodeval); // for r,g and b channels..
				map[1](j, i) = (gcodeval);
				map[2](j, i) = (bcodeval);

				/*if (ncodeinfo)*/{
					nmap[0](j, i) = (nrcodeval); // for r,g and b channels..
					nmap[1](j, i) = (ngcodeval);
					nmap[2](j, i) = (nbcodeval);
				}

			} else if (gencodebook) {
#ifdef ONLINE_MAPPING
				codeinfo->SetValue(pixelvalues, rcodeval, gcodeval, bcodeval);
				if (ncodeinfo)
				ncodeinfo->SetValue(npixelvalues, nrcodeval, ngcodeval,
						nbcodeval);
#else
				codeinfo->SetValue(rcodeval, gcodeval, bcodeval);
				if (ncodeinfo)
					ncodeinfo->SetValue(nrcodeval, ngcodeval, nbcodeval);
#endif
			} else {
				if (!softq) { // map to codebook centers...
#if defined ONLINE_MAPPING || defined LIVE_QUANTIZATION
					map[0](j, i)
					= codeinfo->GetClusterCenter(rpixval, rcodeval); // for r,g and b channels..
					map[1](j, i)
					= codeinfo->GetClusterCenter(gpixval, gcodeval);
					map[2](j, i)
					= codeinfo->GetClusterCenter(bpixval, bcodeval);

					if (ncodeinfo) {
						nmap[0](j, i) = ncodeinfo->GetClusterCenter(nrpixval,
								nrcodeval); // for r,g and b channels..
						nmap[1](j, i) = ncodeinfo->GetClusterCenter(ngpixval,
								ngcodeval);
						nmap[2](j, i) = ncodeinfo->GetClusterCenter(nbpixval,
								nbcodeval);
					}
#else
					map[0](j, i) = codeinfo->GetClusterCenter(rcodeval); // for r,g and b channels..
					map[1](j, i) = codeinfo->GetClusterCenter(gcodeval);
					map[2](j, i) = codeinfo->GetClusterCenter(bcodeval);

					if (ncodeinfo) {
						nmap[0](j, i) = ncodeinfo->GetClusterCenter(nrcodeval); // for r,g and b channels..
						nmap[1](j, i) = ncodeinfo->GetClusterCenter(ngcodeval);
						nmap[2](j, i) = ncodeinfo->GetClusterCenter(nbcodeval);
					}
#endif
				} else {
					map[0](j, i) = (rcodeval); // for r,g and b channels..
					map[1](j, i) = (gcodeval);
					map[2](j, i) = (bcodeval);

					if (ncodeinfo) {
						nmap[0](j, i) = (nrcodeval); // for r,g and b channels..
						nmap[1](j, i) = (ngcodeval);
						nmap[2](j, i) = (nbcodeval);
					}
				}
			}
		}
	}
}
/*
 * As this is used with combinations of other features, so zeros are padded to have the
 * same cell dimensions.*/
#ifdef PATCH_CODE_VERSION2
void EncodedPatchFeatures::ComputeTextureFeatures(Image &imgref,
		UINT xcellsize, UINT ycellsize, vector<REAL>&features) {
	if (petype == PET_SplitLTP) {
		ComputeSplitTextureFeatures(imgref, xcellsize, ycellsize, 0, 0,
				imgref.columns() - 1, imgref.rows() - 1, features);
		return;
	}
	ComputeTextureFeatures(imgref, xcellsize, ycellsize, 0, 0,
			imgref.columns() - 1, imgref.rows() - 1, features);
}
void EncodedPatchFeatures::ComputeTextureFeatures(Image &imgref,
		UINT xcellsize, UINT ycellsize, int sx, int sy, int ex, int ey,
		vector<REAL>&features) {
	if (petype == PET_SplitLTP) {
		ComputeSplitTextureFeatures(imgref, xcellsize, ycellsize, sx, sy, ex,
				ey, features);
		return;
	}

	LBPMap tposmap[3], tnegmap[3];
	ExtractCodeBookInfo(imgref, tposmap, tnegmap, sx, sy, ex, ey, false); //extract the codebook information

	UINT width = ex - sx + 1, height = ey - sy + 1;
	int xoffset = -maxradius /*+MAX(ceil((REAL) (width % xcellsize) / 2),maxradius)*/;
	int yoffset = -maxradius/*+ MAX(ceil((REAL) (height % ycellsize) / 2),maxradius)*/;
	UINT xbmax = width / xcellsize, ybmax = height / ycellsize;

	UINT offset = 0, endx, endy;
	features.resize(xbmax * ybmax * celldim, 0);

	vector<vector<REAL> > pltpfeat(3, vector<REAL> (cellhistdim));
	REAL possum[3]; // normalization constants...

	for (UINT rind = 0; rind < ybmax * ycellsize; rind += ycellsize) {
		for (UINT colind = 0; colind < xbmax * xcellsize; colind += xcellsize) {

			for (int yiter = rind; yiter < rind + ycellsize; ++yiter) {
				int y = MIN(MAX(yiter+yoffset+sy,0),tposmap[0].nypoints-1);
				for (int xiter = colind; xiter < colind + xcellsize; ++xiter) {
					int x = MIN(MAX(xiter+xoffset+sx,0),tposmap[0].nxpoints-1);
					assert(tposmap[0].GetValue(x, y)!=-1);
					pltpfeat[0][tposmap[0].GetValue(x, y)]++;
					pltpfeat[1][tposmap[1].GetValue(x, y)]++;
					pltpfeat[2][tposmap[2].GetValue(x, y)]++;
				}
			}

			{ // normalization

				REAL *posptr = &features[offset];

				possum[0] = EPS;

				if (norm == L2) {
					for (UINT biter = 0; biter < cellhistdim; ++biter) { // add the three channels
						if (pim.GetChannel() != RGB)
							posptr[biter] = pltpfeat[0][biter];
						else
							posptr[biter] = pltpfeat[0][biter]
									+ pltpfeat[1][biter] + pltpfeat[2][biter];

						possum[0] += posptr[biter] * posptr[biter];
					}
					possum[0] = sqrt(possum[0]);
					for (UINT biter = 0; biter < cellhistdim; ++biter) {
						posptr[biter] = (posptr[biter] / possum[0]);
					}
				} else if (norm == NT_Average) {

					for (UINT biter = 0; biter < cellhistdim; ++biter) { // add the three channels
						posptr[biter] = (pltpfeat[0][biter]
								+ pltpfeat[1][biter] + pltpfeat[2][biter])
								/ 3.0;
					}

				} else {
					for (UINT biter = 0; biter < cellhistdim; ++biter) { // add the three channels
						posptr[biter] = pltpfeat[0][biter] + pltpfeat[1][biter]
								+ pltpfeat[2][biter];

						possum[0] += posptr[biter];
					}
					if (norm == L1) {
						for (UINT biter = 0; biter < cellhistdim; ++biter)
							posptr[biter] = (posptr[biter] / possum[0]);
					} else {
						for (UINT biter = 0; biter < cellhistdim; ++biter)
							posptr[biter] = sqrt(posptr[biter] / possum[0]);
					}

				}
			}
			for (UINT i = 0; i < 3; ++i) { // Zero out the temporary  histograms...
				fill(pltpfeat[i].begin(), pltpfeat[i].end(), 0);
			}
			offset += celldim;
		}
	}
}
void EncodedPatchFeatures::ComputeSplitTextureFeatures(Image &imgref,
		UINT xcellsize, UINT ycellsize, int sx, int sy, int ex, int ey,
		vector<REAL>&features) {

	LBPMap tposmap[3], tnegmap[3];
	ExtractCodeBookInfo(imgref, tposmap, tnegmap, sx, sy, ex, ey, false); //extract the codebook information

	UINT width = ex - sx + 1, height = ey - sy + 1;
	int xoffset = -maxradius/*+MAX(ceil((REAL) (width % xcellsize) / 2),maxradius)*/;
	int yoffset = -maxradius /*+MAX(ceil((REAL) (height % ycellsize) / 2),maxradius)*/;
	UINT xbmax = width / xcellsize, ybmax = height / ycellsize;
	UINT offset = 0;
	features.resize(xbmax * ybmax * celldim, 0);

	vector<vector<REAL> > pltpfeat(3, vector<REAL> (cellhistdim)), nltpfeat(3,
			vector<REAL> (cellhistdim));
	REAL possum[3], negsum[3]; // normalization constants...
	for (UINT rind = 0; rind < ybmax * ycellsize; rind += ycellsize) {
		for (UINT colind = 0; colind < xbmax * xcellsize; colind += xcellsize) {

			for (int yiter = rind; yiter < rind + ycellsize; ++yiter) {
				int y = MIN(MAX(yiter+yoffset+sy,0),tposmap[0].nypoints-1);
				for (int xiter = colind; xiter < colind + xcellsize; ++xiter) {
					int x = MIN(MAX(xiter+xoffset+sx,0),tposmap[0].nxpoints-1);
					assert(tposmap[0].GetValue(x, y)!=-1);
					pltpfeat[0][tposmap[0].GetValue(x, y)]++;
					pltpfeat[1][tposmap[1].GetValue(x, y)]++;
					pltpfeat[2][tposmap[2].GetValue(x, y)]++;

					nltpfeat[0][tnegmap[0].GetValue(x, y)]++;
					nltpfeat[1][tnegmap[1].GetValue(x, y)]++;
					nltpfeat[2][tnegmap[2].GetValue(x, y)]++;
				}
			}
			{ // normalization

				REAL *posptr = &features[offset], *negptr = &features[offset
						+ cellhistdim];

				possum[0] = EPS;
				negsum[0] = EPS;

				if (norm == L2) {
					for (UINT biter = 0; biter < cellhistdim; ++biter) { // add the three channels
						if (pim.GetChannel() != RGB) {
							posptr[biter] = pltpfeat[0][biter];
							negptr[biter] = nltpfeat[0][biter];
						} else {
							posptr[biter] = pltpfeat[0][biter]
									+ pltpfeat[1][biter] + pltpfeat[2][biter];
							negptr[biter] = nltpfeat[0][biter]
									+ nltpfeat[1][biter] + nltpfeat[2][biter];
						}
						possum[0] += posptr[biter] * posptr[biter];
						negsum[0] += negptr[biter] * negptr[biter];
					}
					possum[0] = sqrt(possum[0]);
					negsum[0] = sqrt(negsum[0]);
					for (UINT biter = 0; biter < cellhistdim; ++biter) {
						posptr[biter] = (posptr[biter] / possum[0]);
						negptr[biter] = (negptr[biter] / negsum[0]);
					}
				} else if (norm == NT_Average) {

					for (UINT biter = 0; biter < cellhistdim; ++biter) { // add the three channels
						posptr[biter] = (pltpfeat[0][biter]
								+ pltpfeat[1][biter] + pltpfeat[2][biter])
								/ 3.0;
						negptr[biter] = (nltpfeat[0][biter]
								+ nltpfeat[1][biter] + nltpfeat[2][biter])
								/ 3.0;
					}

				} else {
					for (UINT biter = 0; biter < cellhistdim; ++biter) { // add the three channels
						posptr[biter] = pltpfeat[0][biter] + pltpfeat[1][biter]
								+ pltpfeat[2][biter];
						negptr[biter] = nltpfeat[0][biter] + nltpfeat[1][biter]
								+ nltpfeat[2][biter];

						possum[0] += posptr[biter];
						negsum[0] += negptr[biter];
					}
					if (norm == L1) {
						for (UINT biter = 0; biter < cellhistdim; ++biter) {
							posptr[biter] = (posptr[biter] / possum[0]);
							negptr[biter] = (negptr[biter] / negsum[0]);
						}
					} else {
						for (UINT biter = 0; biter < cellhistdim; ++biter) {
							posptr[biter] = sqrt(posptr[biter] / possum[0]);
							negptr[biter] = sqrt(negptr[biter] / negsum[0]);
						}
					}

				}
			}
			for (UINT i = 0; i < 3; ++i) { // Zero out the temporary  histograms...
				fill(pltpfeat[i].begin(), pltpfeat[i].end(), 0);
				fill(nltpfeat[i].begin(), nltpfeat[i].end(), 0);
			}
			offset += celldim;
		}
	}
}
void EncodedPatchFeatures::ComputeLQPHistogram(Image &image, UINT sbin,
		LBPMap tlbpmap[], REAL *feat) {
	int dims[2];
	dims[0] = image.rows();
	dims[1] = image.columns();
	// As histograms are only computed for the centeral part, rest is all zeros
	// so there is no need
	// memory for caching histograms & their norms
	int blocks[2], offset = 0, roffset = 0;
	blocks[0] = (int) round((double) dims[0] / (double) sbin);
	blocks[1] = (int) round((double) dims[1] / (double) sbin);
	UINT bdim = blocks[0] * blocks[1];
	double *rhist = new double[(bdim * celldim)], *ghist = new double[(bdim
			* celldim)], *bhist = new double[(bdim * celldim)], *hist =
			new double[(bdim * celldim)];
	double *lbpnorm = new double[bdim];
	fill(lbpnorm, lbpnorm + bdim, 0);
	fill(rhist, rhist + bdim * celldim, 0);
	fill(ghist, ghist + bdim * celldim, 0);
	fill(bhist, bhist + bdim * celldim, 0);
	fill(hist, hist + bdim * celldim, 0); // for pooling of RGB channels....

	int out[2];
	out[0] = max(blocks[0] - 2, 0);
	out[1] = max(blocks[1] - 2, 0);

	int visible[2];
	visible[0] = blocks[0] * sbin;
	visible[1] = blocks[1] * sbin;
	UINT tx, ty;
	for (int x = maxradius; x < visible[1] - maxradius; x++) {
		for (int y = maxradius; y < visible[0] - maxradius; y++) {

			// first color channel
			/*because the value of map is offset by radius, value of x at x-radius */
			/*-1 because of < ; MIN(x, dims[1] - radius-1) for which lbp map is computed;
			 * final -radius the offset in lbpmap*/
			tx = MIN(x, dims[1] - maxradius-1) - maxradius;

			ty = MIN(y, dims[0] - maxradius-1) - maxradius;

			int ro = tlbpmap[0].GetValue(tx, ty), // red, GREEN & BLUE
					go = tlbpmap[1].GetValue(tx, ty), //
					bo = tlbpmap[2].GetValue(tx, ty); //
			// add to 4 histograms around pixel using linear interpolation
			double xp = ((double) x + 0.5) / (double) sbin - 0.5;
			double yp = ((double) y + 0.5) / (double) sbin - 0.5;
			int ixp = (int) floor(xp);
			int iyp = (int) floor(yp);
			double vx0 = xp - ixp;
			double vy0 = yp - iyp;
			double vx1 = 1.0 - vx0;
			double vy1 = 1.0 - vy0;
			double v = 1; // pooling of gradient magnitude can be an option

			if (ixp >= 0 && iyp >= 0) {
				*(rhist + ixp * blocks[0] + iyp + ro * bdim) += vx1 * vy1 * v;
				*(ghist + ixp * blocks[0] + iyp + go * bdim) += vx1 * vy1 * v;
				*(bhist + ixp * blocks[0] + iyp + bo * bdim) += vx1 * vy1 * v;
			}

			if (ixp + 1 < blocks[1] && iyp >= 0) {
				*(rhist + (ixp + 1) * blocks[0] + iyp + ro * bdim) += vx0 * vy1
						* v;
				*(ghist + (ixp + 1) * blocks[0] + iyp + go * bdim) += vx0 * vy1
						* v;
				*(bhist + (ixp + 1) * blocks[0] + iyp + bo * bdim) += vx0 * vy1
						* v;
			}

			if (ixp >= 0 && iyp + 1 < blocks[0]) {
				*(rhist + ixp * blocks[0] + (iyp + 1) + ro * bdim) += vx1 * vy0
						* v;
				*(ghist + ixp * blocks[0] + (iyp + 1) + go * bdim) += vx1 * vy0
						* v;
				*(bhist + ixp * blocks[0] + (iyp + 1) + bo * bdim) += vx1 * vy0
						* v;
			}

			if (ixp + 1 < blocks[1] && iyp + 1 < blocks[0]) {
				*(rhist + (ixp + 1) * blocks[0] + (iyp + 1) + ro * bdim) += vx0
						* vy0 * v;
				*(ghist + (ixp + 1) * blocks[0] + (iyp + 1) + go * bdim) += vx0
						* vy0 * v;
				*(bhist + (ixp + 1) * blocks[0] + (iyp + 1) + bo * bdim) += vx0
						* vy0 * v;
			}
		}
	}
	// compute energy in each block by summing over orientations
	for (int o = 0; o < celldim; o++) {
		double *srcr = rhist + o * bdim, *srcg = ghist + o * bdim, *srcb =
				bhist + o * bdim, *src = hist + o * bdim;
		double *dst = lbpnorm;
		double *end = lbpnorm + bdim;
		while (dst < end) {

			if (norm == L2) {
				if (pim.GetChannel() != RGB)
					*src = *srcr; // pooling of rgb channels...
				else
					*src = *srcr + *srcg + *srcb; // pooling of rgb channels...
				*(dst++) += *src * *src;
			} else {
				*src = *srcr + *srcg + *srcb; // pooling of rgb channels...
				*(dst++) += *src;
			}
			srcr++;
			srcg++;
			srcb++;
			src++;
		}
	}
	// compute features
	REAL sc1, sc2, sc3, sc4; /*Extract only the centeral pixels*/
	for (int x = 0; x < out[1]; x++) {
		for (int y = 0; y < out[0]; y++) {
			REAL *dst = feat + celldim * (x + y * out[1]); // write in column major order
			double *p, n1, n2, n3, n4;
			double *src = hist + (x + 1) * blocks[0] + (y + 1);
			if (useblocknorm) {/*normlaize by lbp energy of four neighbouring blocks*/
				p = lbpnorm + (x + 1) * blocks[0] + y + 1;
				n1 = 1.0 / (*p + *(p + 1) + *(p + blocks[0]) + *(p + blocks[0]
						+ 1) + EPS);
				p = lbpnorm + (x + 1) * blocks[0] + y;
				n2 = 1.0 / (*p + *(p + 1) + *(p + blocks[0]) + *(p + blocks[0]
						+ 1) + EPS);
				p = lbpnorm + x * blocks[0] + y + 1;
				n3 = 1.0 / (*p + *(p + 1) + *(p + blocks[0]) + *(p + blocks[0]
						+ 1) + EPS);
				p = lbpnorm + x * blocks[0] + y;
				n4 = 1.0 / (*p + *(p + 1) + *(p + blocks[0]) + *(p + blocks[0]
						+ 1) + EPS);

				for (int o = 0; o < celldim; o++) {
					double sum = *src;
					double h1 = sqrt(sum * n1);
					double h2 = sqrt(sum * n2);
					double h3 = sqrt(sum * n3);
					double h4 = sqrt(sum * n4);
					*(dst + o) = h4;
					*(dst + celldim + o) = h2;
					*(dst + celldim * 2 + o) = h3;
					*(dst + celldim * 3 + o) = h1;
					src += bdim;
				}
			} else {
				// single cell normalization...
				// single cell normalization...
				if (norm == L2) {
					REAL normval = 1 / sqrt(
							*(lbpnorm + (x + 1) * blocks[0] + y + 1) + EPS);
					for (int o = 0; o < celldim; o++) {
						*dst = *src * normval;
						src += bdim;
						++dst;
					}
				} else {
					REAL normval = 1.0 / (*(lbpnorm + (x + 1) * blocks[0] + y
							+ 1) + EPS);
					for (int o = 0; o < celldim; o++) {
						*dst = sqrt(*src * normval);
						src += bdim;
						++dst;
					}
				}
			}
		}
	}
	delete[] hist;
	delete[] rhist;
	delete[] ghist;
	delete[] bhist;
	delete[] lbpnorm;
}

void EncodedPatchFeatures::ComputeSplitLQPHistogram(Image &image, UINT sbin,
		LBPMap tposmap[], LBPMap tnegmap[], REAL *feat) {
	int dims[2];
	dims[0] = image.rows();
	dims[1] = image.columns();

	// memory for caching orientation histograms & their norms
	int blocks[2], histdim, bdim;
	blocks[0] = (int) round((double) dims[0] / (double) sbin);
	blocks[1] = (int) round((double) dims[1] / (double) sbin);
	histdim = blocks[0] * blocks[1] * cellhistdim;
	bdim = blocks[0] * blocks[1];
	double *prhist = new double[histdim], *pghist = new double[(histdim)],
			*pbhist = new double[histdim], *phist = new double[(histdim)];
	double *nrhist = new double[histdim], *nghist = new double[(histdim)],
			*nbhist = new double[histdim], *nhist = new double[(histdim)];
	double *pnorm = new double[(bdim)], *nnorm = new double[(bdim)];
	fill(pnorm, pnorm + bdim, 0);
	fill(nnorm, nnorm + bdim, 0);

	fill(prhist, prhist + histdim, 0);
	fill(pghist, pghist + histdim, 0);
	fill(pbhist, pbhist + histdim, 0);
	fill(phist, phist + histdim, 0); // for pooling of RGB channels....

	fill(nrhist, nrhist + histdim, 0);
	fill(nghist, nghist + histdim, 0);
	fill(nbhist, nbhist + histdim, 0);
	fill(nhist, nhist + histdim, 0); // for pooling of RGB channels....

	int out[2];
	out[0] = max(blocks[0] - 2, 0);
	out[1] = max(blocks[1] - 2, 0);

	int visible[2];
	visible[0] = blocks[0] * sbin;
	visible[1] = blocks[1] * sbin;
	UINT tx, ty;
	for (int x = maxradius; x < visible[1] - maxradius; x++) {
		for (int y = maxradius; y < visible[0] - maxradius; y++) {

			// first color channel
			/*because the value of map is offset by radius, value of x at x-radius */
			/*-1 because of < ; MIN(x, dims[1] - radius-1) for which lbp map is computed;
			 * final -radius the offset in lbpmap*/
			tx = MIN(x, dims[1] - maxradius-1) - maxradius;

			ty = MIN(y, dims[0] - maxradius-1) - maxradius;

			UINT pro = tposmap[0].GetValue(tx, ty), // red, GREEN & BLUE
					pgo = tposmap[1].GetValue(tx, ty), //
					pbo = tposmap[2].GetValue(tx, ty); //
			UINT nro = tnegmap[0].GetValue(tx, ty), // red, GREEN & BLUE
					ngo = tnegmap[1].GetValue(tx, ty), //
					nbo = tnegmap[2].GetValue(tx, ty); //

			// add to 4 histograms around pixel using linear interpolation
			double xp = ((double) x + 0.5) / (double) sbin - 0.5;
			double yp = ((double) y + 0.5) / (double) sbin - 0.5;
			int ixp = (int) floor(xp);
			int iyp = (int) floor(yp);
			double vx0 = xp - ixp;
			double vy0 = yp - iyp;
			double vx1 = 1.0 - vx0;
			double vy1 = 1.0 - vy0;
			double pvr = 1, pvg = 1, pvb = 1, // pooling of gradient magnitude can be an option
					nvr = 1, nvg = 1, nvb = 1; // pooling of gradient magnitude can be an option

			if (softq == QT_Discrim) { /// Get the discirminant code count ratio information extracted from the learned clusters...
				codeinfo->GetClusterCenterWithRatio(pro, pvr);
				codeinfo->GetClusterCenterWithRatio(pgo, pvg);
				codeinfo->GetClusterCenterWithRatio(pbo, pvb);

				ncodeinfo->GetClusterCenterWithRatio(nro, nvr);
				ncodeinfo->GetClusterCenterWithRatio(ngo, nvg);
				ncodeinfo->GetClusterCenterWithRatio(nbo, nvb);
			}

			if (ixp >= 0 && iyp >= 0) {
				*(prhist + ixp * blocks[0] + iyp + pro * bdim) += vx1 * vy1
						* pvr;
				*(pghist + ixp * blocks[0] + iyp + pgo * bdim) += vx1 * vy1
						* pvg;
				*(pbhist + ixp * blocks[0] + iyp + pbo * bdim) += vx1 * vy1
						* pvb;

				*(nrhist + ixp * blocks[0] + iyp + nro * bdim) += vx1 * vy1
						* nvr;
				*(nghist + ixp * blocks[0] + iyp + ngo * bdim) += vx1 * vy1
						* nvg;
				*(nbhist + ixp * blocks[0] + iyp + nbo * bdim) += vx1 * vy1
						* nvb;

			}

			if (ixp + 1 < blocks[1] && iyp >= 0) {

				*(prhist + (ixp + 1) * blocks[0] + iyp + pro * bdim) += vx0
						* vy1 * pvr;
				*(pghist + (ixp + 1) * blocks[0] + iyp + pgo * bdim) += vx0
						* vy1 * pvg;
				*(pbhist + (ixp + 1) * blocks[0] + iyp + pbo * bdim) += vx0
						* vy1 * pvb;

				*(nrhist + (ixp + 1) * blocks[0] + iyp + nro * bdim) += vx0
						* vy1 * nvr;
				*(nghist + (ixp + 1) * blocks[0] + iyp + ngo * bdim) += vx0
						* vy1 * nvg;
				*(nbhist + (ixp + 1) * blocks[0] + iyp + nbo * bdim) += vx0
						* vy1 * nvb;
			}

			if (ixp >= 0 && iyp + 1 < blocks[0]) {

				*(prhist + ixp * blocks[0] + (iyp + 1) + pro * bdim) += vx1
						* vy0 * pvr;
				*(pghist + ixp * blocks[0] + (iyp + 1) + pgo * bdim) += vx1
						* vy0 * pvg;
				*(pbhist + ixp * blocks[0] + (iyp + 1) + pbo * bdim) += vx1
						* vy0 * pvb;

				*(nrhist + ixp * blocks[0] + (iyp + 1) + nro * bdim) += vx1
						* vy0 * nvr;
				*(nghist + ixp * blocks[0] + (iyp + 1) + ngo * bdim) += vx1
						* vy0 * nvg;
				*(nbhist + ixp * blocks[0] + (iyp + 1) + nbo * bdim) += vx1
						* vy0 * nvb;

			}

			if (ixp + 1 < blocks[1] && iyp + 1 < blocks[0]) {

				*(prhist + (ixp + 1) * blocks[0] + (iyp + 1) + pro * bdim)
						+= vx0 * vy0 * pvr;
				*(pghist + (ixp + 1) * blocks[0] + (iyp + 1) + pgo * bdim)
						+= vx0 * vy0 * pvg;
				*(pbhist + (ixp + 1) * blocks[0] + (iyp + 1) + pbo * bdim)
						+= vx0 * vy0 * pvb;

				*(nrhist + (ixp + 1) * blocks[0] + (iyp + 1) + nro * bdim)
						+= vx0 * vy0 * nvr;
				*(nghist + (ixp + 1) * blocks[0] + (iyp + 1) + ngo * bdim)
						+= vx0 * vy0 * nvg;
				*(nbhist + (ixp + 1) * blocks[0] + (iyp + 1) + nbo * bdim)
						+= vx0 * vy0 * nvb;
			}
		}
	}
	for (int o = 0; o < cellhistdim; o++) {

		double *psrcr = prhist + o * bdim, *psrcg = pghist + o * bdim, *psrcb =
				pbhist + o * bdim, *psrc = phist + o * bdim;

		double *nsrcr = nrhist + o * bdim, *nsrcg = nghist + o * bdim, *nsrcb =
				nbhist + o * bdim, *nsrc = nhist + o * bdim;

		double *pdst = pnorm, *ndst = nnorm;
		double *end = pnorm + blocks[1] * blocks[0];
		while (pdst < end) {

			if (norm == L2) {
				if (pim.GetChannel() != RGB) {
					*psrc = *psrcr; // pooling of rgb channels...
					*nsrc = *nsrcr; // pooling of rgb channels...
				} else {
					*psrc = *psrcr + *psrcg + *psrcb; // pooling of rgb channels...
					*nsrc = *nsrcr + *nsrcg + *nsrcb; // pooling of rgb channels...
				}

				*(pdst++) += *psrc * *psrc;
				*(ndst++) += *nsrc * *nsrc;
			} else {
				*psrc = *psrcr + *psrcg + *psrcb; // pooling of rgb channels...
				*nsrc = *nsrcr + *nsrcg + *nsrcb; // pooling of rgb channels...

				*(pdst++) += *psrc;
				*(ndst++) += *nsrc;
			}

			psrcr++;
			psrcg++;
			psrcb++;
			psrc++;

			nsrcr++;
			nsrcg++;
			nsrcb++;
			nsrc++;
		}
	}
	// compute features
	REAL sc1, sc2, sc3, sc4;
	for (int x = 0; x < out[1]; x++) {
		for (int y = 0; y < out[0]; y++) {
			REAL *pdst = feat + celldim * (x + y * out[1]); // write in column major order
			//			REAL *dst = feat + celldim * dimmult * (x + y * out[1]); // write in column major order
			REAL *ndst = pdst + cellhistdim; //	cellhistdim = nbins * (useblocknorm == true ? 4 : 1);
			double *p, pn1, pn2, pn3, pn4, nn1, nn2, nn3, nn4;
			double *psrc = phist + (x + 1) * blocks[0] + (y + 1);
			double *nsrc = nhist + (x + 1) * blocks[0] + (y + 1);
			// single cell normalization...
			if (norm == L2) {
				REAL pnormval = 1 / sqrt(
						*(pnorm + (x + 1) * blocks[0] + y + 1) + EPS);
				REAL nnormval = 1 / sqrt(
						*(nnorm + (x + 1) * blocks[0] + y + 1) + EPS);
				for (int o = 0; o < cellhistdim; o++) {
					*pdst = *psrc * pnormval;
					*ndst = *nsrc * nnormval;

					psrc += bdim;
					nsrc += bdim;

					++pdst;
					++ndst;
				}
			} else if (norm == LOG2) {

				for (int o = 0; o < cellhistdim; o++) {
					*pdst = log2(*psrc + 1);
					*ndst = log2(*nsrc + 1);
					psrc += bdim;
					nsrc += bdim;
					++pdst;
					++ndst;
				}
			} else {
				REAL pnormval = *(pnorm + (x + 1) * blocks[0] + y + 1) + EPS;
				REAL nnormval = *(nnorm + (x + 1) * blocks[0] + y + 1) + EPS;
				for (int o = 0; o < cellhistdim; o++) {
					*pdst = sqrt(*psrc / pnormval);
					*ndst = sqrt(*nsrc / nnormval);
					psrc += bdim;
					nsrc += bdim;
					++pdst;
					++ndst;
				}
			}
		}
	}

	delete[] phist;
	delete[] prhist;
	delete[] pghist;
	delete[] pbhist;
	delete[] pnorm;

	delete[] nhist;
	delete[] nrhist;
	delete[] nghist;
	delete[] nbhist;
	delete[] nnorm;
}

//#else
//re-test the function
void EncodedPatchFeatures::ComputeLQPHistogramDiscrete(Image& imgref,
		UINT index, UINT sbin, LBPMap map[3]) {
	UINT x, y, cxbmax = lbpfeatures.GetXBlock(index), xbmax = cxbmax * sbin,
			cybmax = lbpfeatures.GetYBlock(index), ybmax = cybmax * sbin, endx,
			endy, offset = 0, roffset;
	int rval;
	vector<REAL> &features = lbpfeatures.GetFeatureRef(index);
	vector<vector<REAL> > feat(3, vector<REAL> (celldim));
	REAL possum[3];
	// in all files..boundary .
	int boundary = lfeatoffset - sbin, bcell = boundary / sbin;
	xbmax -= boundary;
	ybmax -= boundary;

	for (UINT rind = boundary, ystart = bcell; rind < ybmax; rind += sbin, ++ystart) {
		y = rind + sbin - maxradius; // patch top left coordinates for which rind+foffset is center...
		// because the map contains an offset maxradius for concentic circles..
		// means map(0,0) contain value of map(maxradius,maxradius) so value of hist(offset,offset) can be found at
		// map(offset-maxradius,offset-maxradius)
		offset = (ystart * cxbmax + bcell) * celldim;

		for (UINT colind = boundary; colind < xbmax; colind += sbin) {
			x = colind + sbin - maxradius;

			endx = MIN(x + sbin, map[0].nxpoints);
			endy = MIN(y + sbin, map[0].nypoints);

			for (UINT yiter = y; yiter < endy; ++yiter) // replace ++yiter with +=patch stride...
				for (UINT xiter = x; xiter < endx; ++xiter) {
					if (!softq) {
						rval = map[0].GetValue(xiter, yiter);
						if (rval >= 0) // because if feature is computed for red -channel
						//!it is implicit that it is computed for green & blue
						{
							feat[0][rval]++;
							feat[1][map[1].GetValue(xiter, yiter)]++;
							feat[2][map[2].GetValue(xiter, yiter)]++;
						}
					} // else pool all the features...
					else {
						const vector<REAL> & rvec =
								codeinfo->GetSoftClusterCenter(
										map[0].GetValue(xiter, yiter)), &gvec =
								codeinfo->GetSoftClusterCenter(
										map[1].GetValue(xiter, yiter)), &bvec =
								codeinfo->GetSoftClusterCenter(
										map[2].GetValue(xiter, yiter));

						for (UINT siter = 0; siter < celldim; ++siter) { // pool features...
							feat[0][siter] += rvec[siter];
							feat[1][siter] += gvec[siter];
							feat[2][siter] += bvec[siter];
						}
					}
				}

			REAL *posptr = &features[offset]; // add offsets in features code ...
			possum[0] = EPS;
			if (norm == L2) {
				for (UINT biter = 0; biter < celldim; ++biter) { // add the three channels
					posptr[biter] = feat[0][biter] + feat[1][biter]
							+ feat[2][biter];
					possum[0] += posptr[biter] * posptr[biter];
				}
				possum[0] = sqrt(possum[0]);

				for (UINT biter = 0; biter < celldim; ++biter)
					posptr[biter] = posptr[biter] / possum[0];

			} else {
				for (UINT biter = 0; biter < celldim; ++biter) { // add the three channels
					posptr[biter] = feat[0][biter] + feat[1][biter]
							+ feat[2][biter];
					possum[0] += posptr[biter];
				}
				for (UINT biter = 0; biter < celldim; ++biter)
					posptr[biter] = sqrt(posptr[biter] / possum[0]);
			}

			for (UINT i = 0; i < 3; ++i) // Zero out the temporary  histograms...
				fill(feat[i].begin(), feat[i].end(), 0);
			offset += celldim;
		}
	}
}
void EncodedPatchFeatures::ComputeSplitLQPHistogramDiscrete(Image& imgref,
		UINT index, UINT sbin, LBPMap map[3], LBPMap nmap[3]) {
	UINT x, y, cxbmax = lbpfeatures.GetXBlock(index), xbmax = cxbmax * sbin,
			cybmax = lbpfeatures.GetYBlock(index), ybmax = cybmax * sbin, endx,
			endy, offset = 0, roffset;
	int rval;
	vector<REAL> &features = lbpfeatures.GetFeatureRef(index);
	vector<vector<REAL> > feat(3, vector<REAL> (cellhistdim));
	vector<vector<REAL> > nfeat(3, vector<REAL> (cellhistdim));
	REAL possum[3], negsum[3];
	;
	// in all files..boundary .
	int boundary = lfeatoffset - sbin, bcell = boundary / sbin;
	xbmax -= boundary;
	ybmax -= boundary;

	for (UINT rind = boundary, ystart = bcell; rind < ybmax; rind += sbin, ++ystart) {
		y = rind + sbin - maxradius; // patch top left coordinates for which rind+foffset is center...
		// because the map contains an offset maxradius for concentic circles..
		// means map(0,0) contain value of map(maxradius,maxradius) so value of hist(offset,offset) can be found at
		// map(offset-maxradius,offset-maxradius)
		offset = (ystart * cxbmax + bcell) * celldim;

		for (UINT colind = boundary; colind < xbmax; colind += sbin) {
			x = colind + sbin - maxradius;

			endx = MIN(x + sbin, map[0].nxpoints);
			endy = MIN(y + sbin, map[0].nypoints);

			for (UINT yiter = y; yiter < endy; ++yiter) // replace ++yiter with +=patch stride...
				for (UINT xiter = x; xiter < endx; ++xiter) {
					rval = map[0].GetValue(xiter, yiter);
					if (rval >= 0) // because if feature is computed for red -channel
					//!it is implicit that it is computed for green & blue
					{
						feat[0][rval]++;
						feat[1][map[1].GetValue(xiter, yiter)]++;
						feat[2][map[2].GetValue(xiter, yiter)]++;

						nfeat[0][nmap[0].GetValue(xiter, yiter)]++;
						nfeat[1][nmap[1].GetValue(xiter, yiter)]++;
						nfeat[2][nmap[2].GetValue(xiter, yiter)]++;
					}
				}

			REAL *posptr = &features[offset], *negptr = posptr + cellhistdim; // add offsets in features code ...
			possum[0] = EPS;
			negsum[0] = EPS;
			if (norm == L2) {
				for (UINT biter = 0; biter < cellhistdim; ++biter) { // add the three channels
					posptr[biter] = feat[0][biter] + feat[1][biter]
							+ feat[2][biter];
					possum[0] += posptr[biter] * posptr[biter];

					negptr[biter] = nfeat[0][biter] + nfeat[1][biter]
							+ nfeat[2][biter];
					negsum[0] += negptr[biter] * negptr[biter];
				}
				possum[0] = sqrt(possum[0]);
				negsum[0] = sqrt(negsum[0]);

				for (UINT biter = 0; biter < cellhistdim; ++biter) {
					posptr[biter] = posptr[biter] / possum[0];
					negptr[biter] = negptr[biter] / negsum[0];
				}

			} else {
				for (UINT biter = 0; biter < cellhistdim; ++biter) { // add the three channels
					posptr[biter] = feat[0][biter] + feat[1][biter]
							+ feat[2][biter];
					possum[0] += posptr[biter];

					negptr[biter] = nfeat[0][biter] + nfeat[1][biter]
							+ nfeat[2][biter];
					negsum[0] += negptr[biter];
				}
				for (UINT biter = 0; biter < cellhistdim; ++biter) {
					posptr[biter] = sqrt(posptr[biter] / possum[0]);
					negptr[biter] = sqrt(negptr[biter] / negsum[0]);
				}
			}

			for (UINT i = 0; i < 3; ++i) // Zero out the temporary  histograms...
			{
				fill(feat[i].begin(), feat[i].end(), 0);
				fill(nfeat[i].begin(), nfeat[i].end(), 0);
			}
			offset += celldim;
		}
	}
}

#endif

void EncodedPatchFeatures::GetFeatures(UINT index, int x, int y, int width_,
		int height_, vector<REAL>&feat) { // There is offset of one cell so, input x=0 means x=8 by taking
	//	into account offset and the removed border
	int xcell = x / cell, ycell = y / cell;
	width_ /= cell;
	height_ /= cell;
	lbpfeatures.GetSkippedHOGCells(index, xcell, ycell, width_, height_, skip,
			&feat[0]);
}
void EncodedPatchFeatures::GetFeatures(UINT index, int x, int y, int width_,
		int height_, REAL*feat) { // There is offset of one cell so, input x=0 means x=8 by taking
	//	into account offset and the removed border
	int xcell = x / cell, ycell = y / cell;
	width_ /= cell;
	height_ /= cell;
	lbpfeatures.GetSkippedHOGCells(index, xcell, ycell, width_, height_, skip,
			feat);
}
/*
 * Contains only folded HOG's....*/
void EncodedPatchFeatures::GetFoldedFeatures(UINT index, int x, int y,
		int twidth, int theight, vector<REAL>&feat) { // There is offset of one cell so, input x=0 means x=8 by taking
	//	into account offset and the removed border
	int xcell = x / cell, ycell = y / cell;
	twidth /= cell;
	theight /= cell;
	lbpfeatures.GetSkippedHOGCells(index, xcell, ycell, twidth, theight, skip,
			&feat[0]);
}

/*
 * Contains only folded HOG's....*/
void EncodedPatchFeatures::GetFoldedFeatures(UINT index, int x, int y,
		int twidth, int theight, REAL*feat) { // There is offset of one cell so, input x=0 means x=8 by taking
	//	into account offset and the removed border
	int xcell = x / cell, ycell = y / cell;
	twidth /= cell;
	theight /= cell;
	lbpfeatures.GetSkippedHOGCells(index, xcell, ycell, twidth, theight, skip,
			feat);
}

// Only Flipping hog images....
void EncodedPatchFeatures::GetFlippedFeatures(UINT index, int x, int y,
		int twidth, int theight, vector<REAL>&feat) { // There is offset of one cell so, input x=0 means x=8 by taking
	//	into account offset and the removed border
	int xcell = x / cell, ycell = y / cell;
	twidth /= cell;
	theight /= cell;
	lbpfeatures.GetSkippedHOGCells(index, xcell, ycell, twidth, theight, skip,
			&feat[0]);
}

void EncodedPatchFeatures::GetFlippedFeatures(int twidth, int theight,
		REAL *ifeat, REAL *ofeat) { // There is offset of one cell so, input x=0 means x=8 by taking //	into account offset and the removed border
	//	twidth /= cell;
	//	theight /= cell;
	//	HInfo::FlipFeatures(twidth / skip, theight / skip, nbins * 2 * nthrlevels,
	//			ifeat, ofeat);
}

void EncodedPatchFeatures::UnFoldWeights(REAL *input, UINT w, UINT h,
		vector<REAL>& weights) {

	//	UINT tcwidth = w / cell, tcheight = h / cell;
	//	HInfo hinfo((UINT) ceil((double) tcwidth / (2.0 * skip)),
	//			(tcheight / skip), 1, nbins * 2 * nthrlevels, input);
	//	hinfo.GetUnFoldedHOGCells(tcwidth / skip, &weights[0]);// Flip the features
}

void EncodedPatchFeatures::PadFeatureMap(UINT index, UINT padx, UINT pady) { // does padding in the pixel space
	// first pad the hogFeatures map...
	padx /= cell;
	pady /= cell;
	lbpfeatures.PadFeatureMap(index, padx, pady);
}

void EncodedPatchFeatures::DotProduct(UINT index, UINT fwidth, UINT fheight,
		vector<REAL>&filter, Store&response) {
	fwidth /= cell;
	fheight /= cell;
	lbpfeatures.SkippedDotProduct(index, fwidth, fheight, skip, &filter[0],
			response);
}

void EncodedPatchFeatures::DotProduct(UINT index, UINT fwidth, UINT fheight,
		REAL *filter, Store&response) {
	fwidth /= cell;
	fheight /= cell;
	lbpfeatures.SkippedDotProduct(index, fwidth, fheight, skip, filter,
			response);
}
void SoftLQPFeatures::ComputeLBPFeatures(Image &image, UINT index, UINT sbin,
		UINT lbpstride) {
	if (hmethod == HM_Discrete) {
		ComputeDiscreteLQPHistogram(image, index, sbin, lbpstride);
		return;
	}
	LBPMap tposmap[3], tnegmap[3];
	ExtractCodeBookInfo(image, tposmap, tnegmap, false); //extract the codebook information
	vector<REAL> &tvec = lbpfeatures.GetFeatureRef(index);
	REAL *feat = &tvec[0];
	int dims[2];
	dims[0] = image.rows();
	dims[1] = image.columns();

	// memory for caching orientation histograms & their norms
	int blocks[2], histdim, bdim;
	blocks[0] = (int) round((double) dims[0] / (double) sbin);
	blocks[1] = (int) round((double) dims[1] / (double) sbin);
	histdim = blocks[0] * blocks[1] * cellhistdim;
	bdim = blocks[0] * blocks[1];
	double *prhist = new double[histdim], *pghist = new double[(histdim)],
			*pbhist = new double[histdim], *phist = new double[(histdim)];
	double *nrhist = new double[histdim], *nghist = new double[(histdim)],
			*nbhist = new double[histdim], *nhist = new double[(histdim)];
	double *pnorm = new double[(bdim)], *nnorm = new double[(bdim)];
	fill(pnorm, pnorm + bdim, 0);
	fill(nnorm, nnorm + bdim, 0);

	fill(prhist, prhist + histdim, 0);
	fill(pghist, pghist + histdim, 0);
	fill(pbhist, pbhist + histdim, 0);
	fill(phist, phist + histdim, 0); // for pooling of RGB channels....

	fill(nrhist, nrhist + histdim, 0);
	fill(nghist, nghist + histdim, 0);
	fill(nbhist, nbhist + histdim, 0);
	fill(nhist, nhist + histdim, 0); // for pooling of RGB channels....

	// memory for HOG features
	int out[2];
	out[0] = max(blocks[0] - 2, 0);
	out[1] = max(blocks[1] - 2, 0);

	int visible[2];
	visible[0] = blocks[0] * sbin;
	visible[1] = blocks[1] * sbin;
	UINT tx, ty;
	for (int x = maxradius; x < visible[1] - maxradius; x++) {
		for (int y = maxradius; y < visible[0] - maxradius; y++) {

			// first color channel
			/*because the value of map is offset by radius, value of x at x-radius */
			/*-1 because of < ; MIN(x, dims[1] - radius-1) for which lbp map is computed;
			 * final -radius the offset in lbpmap*/
			tx = MIN(x, dims[1] - maxradius-1) - maxradius;

			ty = MIN(y, dims[0] - maxradius-1) - maxradius;

			//			int pro = tposmap[0].GetValue(tx, ty), // red, GREEN & BLUE
			//					pgo = tposmap[1].GetValue(tx, ty), //
			//					pbo = tposmap[2].GetValue(tx, ty); //
			//			int nro = tnegmap[0].GetValue(tx, ty), // red, GREEN & BLUE
			//					ngo = tnegmap[1].GetValue(tx, ty), //
			//					nbo = tnegmap[2].GetValue(tx, ty); //

			const vector<REAL> & rvec = codeinfo->GetSoftClusterCenter(
					tposmap[0].GetValue(tx, ty)),
					&gvec = codeinfo->GetSoftClusterCenter(
							tposmap[1].GetValue(tx, ty)), &bvec =
							codeinfo->GetSoftClusterCenter(
									tposmap[2].GetValue(tx, ty));
			const vector<REAL> & nrvec = ncodeinfo->GetSoftClusterCenter(
					tnegmap[0].GetValue(tx, ty)),
					&ngvec = ncodeinfo->GetSoftClusterCenter(
							tnegmap[1].GetValue(tx, ty)), &nbvec =
							ncodeinfo->GetSoftClusterCenter(
									tnegmap[2].GetValue(tx, ty));
			//			cout << nrvec << endl;
			// add to 4 histograms around pixel using linear interpolation
			double xp = ((double) x + 0.5) / (double) sbin - 0.5;
			double yp = ((double) y + 0.5) / (double) sbin - 0.5;
			int ixp = (int) floor(xp);
			int iyp = (int) floor(yp);
			double vx0 = xp - ixp;
			double vy0 = yp - iyp;
			double vx1 = 1.0 - vx0;
			double vy1 = 1.0 - vy0;

			for (UINT orien = 0; orien < ncenters; ++orien) {
				if (ixp >= 0 && iyp >= 0) {
					Pool((prhist + ixp * blocks[0] + iyp + orien * bdim),
							vx1 * vy1 * rvec[orien]);
					Pool((pghist + ixp * blocks[0] + iyp + orien * bdim),
							vx1 * vy1 * gvec[orien]);
					Pool((pbhist + ixp * blocks[0] + iyp + orien * bdim),
							vx1 * vy1 * bvec[orien]);

					Pool((nrhist + ixp * blocks[0] + iyp + orien * bdim),
							vx1 * vy1 * nrvec[orien]);
					Pool((nghist + ixp * blocks[0] + iyp + orien * bdim),
							vx1 * vy1 * ngvec[orien]);
					Pool((nbhist + ixp * blocks[0] + iyp + orien * bdim),
							vx1 * vy1 * nbvec[orien]);

				}

				if (ixp + 1 < blocks[1] && iyp >= 0) {

					Pool((prhist + (ixp + 1) * blocks[0] + iyp + orien * bdim),
							vx0 * vy1 * rvec[orien]);
					Pool((pghist + (ixp + 1) * blocks[0] + iyp + orien * bdim),
							vx0 * vy1 * gvec[orien]);
					Pool((pbhist + (ixp + 1) * blocks[0] + iyp + orien * bdim),
							vx0 * vy1 * bvec[orien]);

					Pool((nrhist + (ixp + 1) * blocks[0] + iyp + orien * bdim),
							vx0 * vy1 * nrvec[orien]);
					Pool((nghist + (ixp + 1) * blocks[0] + iyp + orien * bdim),
							vx0 * vy1 * ngvec[orien]);
					Pool((nbhist + (ixp + 1) * blocks[0] + iyp + orien * bdim),
							vx0 * vy1 * nbvec[orien]);
				}

				if (ixp >= 0 && iyp + 1 < blocks[0]) {

					Pool((prhist + ixp * blocks[0] + (iyp + 1) + orien * bdim),
							vx1 * vy0 * rvec[orien]);
					Pool((pghist + ixp * blocks[0] + (iyp + 1) + orien * bdim),
							vx1 * vy0 * gvec[orien]);
					Pool((pbhist + ixp * blocks[0] + (iyp + 1) + orien * bdim),
							vx1 * vy0 * bvec[orien]);

					Pool((nrhist + ixp * blocks[0] + (iyp + 1) + orien * bdim),
							vx1 * vy0 * nrvec[orien]);
					Pool((nghist + ixp * blocks[0] + (iyp + 1) + orien * bdim),
							vx1 * vy0 * ngvec[orien]);
					Pool((nbhist + ixp * blocks[0] + (iyp + 1) + orien * bdim),
							vx1 * vy0 * nbvec[orien]);

				}

				if (ixp + 1 < blocks[1] && iyp + 1 < blocks[0]) {

					Pool(
							(prhist + (ixp + 1) * blocks[0] + (iyp + 1) + orien
									* bdim), vx0 * vy0 * rvec[orien]);
					Pool(
							(pghist + (ixp + 1) * blocks[0] + (iyp + 1) + orien
									* bdim), vx0 * vy0 * gvec[orien]);
					Pool(
							(pbhist + (ixp + 1) * blocks[0] + (iyp + 1) + orien
									* bdim), vx0 * vy0 * bvec[orien]);

					Pool(
							(nrhist + (ixp + 1) * blocks[0] + (iyp + 1) + orien
									* bdim), vx0 * vy0 * nrvec[orien]);
					Pool(
							(nghist + (ixp + 1) * blocks[0] + (iyp + 1) + orien
									* bdim), vx0 * vy0 * ngvec[orien]);
					Pool(
							(nbhist + (ixp + 1) * blocks[0] + (iyp + 1) + orien
									* bdim), vx0 * vy0 * nbvec[orien]);
				}
			}
		}
	}
	//	Print(prhist, prhist + histdim);
	//	cout << endl;
	//	Print(nrhist, nrhist + histdim);
	//	cout << endl;
	//	Print(pghist, pghist + histdim);
	//	cout << endl;
	//	Print(nghist, nghist + histdim);
	//	cout << endl;
	//	Print(pbhist, pbhist + histdim);
	//	cout << endl;
	//	Print(nbhist, nbhist + histdim);
	//	cout << endl;

	// compute energy in each block by summing over orientations
	for (int o = 0; o < cellhistdim; o++) {

		double *psrcr = prhist + o * bdim, *psrcg = pghist + o * bdim, *psrcb =
				pbhist + o * bdim, *psrc = phist + o * bdim;

		double *nsrcr = nrhist + o * bdim, *nsrcg = nghist + o * bdim, *nsrcb =
				nbhist + o * bdim, *nsrc = nhist + o * bdim;

		double *pdst = pnorm, *ndst = nnorm;
		double *end = pnorm + blocks[1] * blocks[0];
		while (pdst < end) {
			if (norm == L2) {
				if (pim.GetChannel() != RGB) {
					*psrc = *psrcr; // pooling of rgb channels...
					*nsrc = *nsrcr; // pooling of rgb channels...
				} else {
					*psrc = *psrcr + *psrcg + *psrcb; // pooling of rgb channels...
					*nsrc = *nsrcr + *nsrcg + *nsrcb; // pooling of rgb channels...
				}

				*(pdst++) += *psrc * *psrc;
				*(ndst++) += *nsrc * *nsrc;
			} else {
				*psrc = *psrcr + *psrcg + *psrcb; // pooling of rgb channels...
				*nsrc = *nsrcr + *nsrcg + *nsrcb; // pooling of rgb channels...

				*(pdst++) += *psrc;
				*(ndst++) += *nsrc;
			}
			psrcr++;
			psrcg++;
			psrcb++;
			psrc++;

			nsrcr++;
			nsrcg++;
			nsrcb++;
			nsrc++;
		}
	}
	// compute features
	REAL sc1, sc2, sc3, sc4;
	for (int x = 0; x < out[1]; x++) {
		for (int y = 0; y < out[0]; y++) {
			REAL *pdst = feat + celldim * (x + y * out[1]); // write in column major order
			//			REAL *dst = feat + celldim * dimmult * (x + y * out[1]); // write in column major order
			REAL *ndst = pdst + cellhistdim; //	cellhistdim = nbins * (useblocknorm == true ? 4 : 1);
			double *p, pn1, pn2, pn3, pn4, nn1, nn2, nn3, nn4;
			double *psrc = phist + (x + 1) * blocks[0] + (y + 1);
			double *nsrc = nhist + (x + 1) * blocks[0] + (y + 1);
			// single cell normalization...
			if (norm == L2) {
				REAL pnormval = 1 / sqrt(
						*(pnorm + (x + 1) * blocks[0] + y + 1) + EPS);
				REAL nnormval = 1 / sqrt(
						*(nnorm + (x + 1) * blocks[0] + y + 1) + EPS);
				for (int o = 0; o < cellhistdim; o++) {
					*pdst = *psrc * pnormval;
					*ndst = *nsrc * nnormval;

					psrc += bdim;
					nsrc += bdim;

					++pdst;
					++ndst;
				}
			} else if (norm == LOG2) {

				for (int o = 0; o < cellhistdim; o++) {
					*pdst = log2(*psrc + 1);
					*ndst = log2(*nsrc + 1);
					psrc += bdim;
					nsrc += bdim;
					++pdst;
					++ndst;
				}
			} else if (norm == L1) {
				REAL pnormval = *(pnorm + (x + 1) * blocks[0] + y + 1) + EPS;
				REAL nnormval = *(nnorm + (x + 1) * blocks[0] + y + 1) + EPS;
				for (int o = 0; o < cellhistdim; o++) {
					*pdst = (*psrc / pnormval);
					*ndst = (*nsrc / nnormval);
					psrc += bdim;
					nsrc += bdim;
					++pdst;
					++ndst;
				}
			} else {
				REAL pnormval = *(pnorm + (x + 1) * blocks[0] + y + 1) + EPS;
				REAL nnormval = *(nnorm + (x + 1) * blocks[0] + y + 1) + EPS;
				for (int o = 0; o < cellhistdim; o++) {
					*pdst = sqrt(*psrc / pnormval);
					*ndst = sqrt(*nsrc / nnormval);
					psrc += bdim;
					nsrc += bdim;
					++pdst;
					++ndst;
				}
			}
		}
	}
	//	cout << endl << tvec << endl << flush;
	delete[] phist;
	delete[] prhist;
	delete[] pghist;
	delete[] pbhist;
	delete[] pnorm;

	delete[] nhist;
	delete[] nrhist;
	delete[] nghist;
	delete[] nbhist;
	delete[] nnorm;
}

void SoftLQPFeatures::ComputeDiscreteLQPHistogram(Image &image, UINT index,
		UINT cellsize, UINT lbpstride) {

	LBPMap tposmap[3], tnegmap[3];
	ExtractCodeBookInfo(image, tposmap, tnegmap, false); //extract the codebook information
	vector<REAL> &features = lbpfeatures.GetFeatureRef(index);

	UINT x, y, cxbmax = lbpfeatures.GetXBlock(index),
			xbmax = cxbmax * cellsize, cybmax = lbpfeatures.GetYBlock(index),
			ybmax = cybmax * cellsize, endx, endy, offset = 0;

	vector<vector<REAL> > pltpfeat(3, vector<REAL> (cellhistdim)), nltpfeat(3,
			vector<REAL> (cellhistdim));
	REAL possum[3], negsum[3]; // normalization constants...
	int boundary = lfeatoffset - foffset, bcell = boundary / cell;
	xbmax -= boundary;
	ybmax -= boundary;

	for (UINT rind = boundary, ystart = bcell; rind < ybmax; rind += cellsize, ++ystart) {
		y = rind + foffset - maxradius; // patch top left coordinates for which rind+foffset is center...
		// because the map contains an offset maxradius for concentic circles..
		// means map(0,0) contain value of map(maxradius,maxradius) so value of hist(offset,offset) can be found at
		// map(offset-maxradius,offset-maxradius)
		offset = (ystart * cxbmax + bcell) * celldim;

		for (UINT colind = boundary; colind < xbmax; colind += cellsize) {
			x = colind + foffset - maxradius;

			endx = MIN(x + cellsize, tposmap[0].nxpoints);
			endy = MIN(y + cellsize, tposmap[0].nypoints);

			for (UINT yiter = y; yiter < endy; ++yiter) // replace ++yiter with +=patch stride...
				for (UINT xiter = x; xiter < endx; ++xiter) {
					/*	}
					 }
					 /*--------------------------------
					 for (UINT rind = 0; rind < ybmax; rind += cellsize) {
					 y = rind + yoffset /*foffset/- maxradius;
					 for (UINT colind = 0; colind < xbmax; colind += cellsize) {
					 x = colind + xoffset/*foffset/- maxradius;

					 endx = MIN(x+cellsize,tposmap[0].nxpoints);
					 endy = MIN(y+cellsize,tposmap[0].nypoints);

					 for (UINT yiter = y; yiter < endy; ++yiter) /// build the histogram using bilinear interpolation
					 for (UINT xiter = x; xiter < endx; ++xiter) {*/

					const vector<REAL> &rvec = codeinfo->GetSoftClusterCenter(
							tposmap[0].GetValue(xiter, yiter)), &gvec =
							codeinfo->GetSoftClusterCenter(
									tposmap[1].GetValue(xiter, yiter)), &bvec =
							codeinfo->GetSoftClusterCenter(
									tposmap[2].GetValue(xiter, yiter));
					const vector<REAL> &nrvec =
							ncodeinfo->GetSoftClusterCenter(
									tnegmap[0].GetValue(xiter, yiter)), &ngvec =
							ncodeinfo->GetSoftClusterCenter(
									tnegmap[1].GetValue(xiter, yiter)), &nbvec =
							ncodeinfo->GetSoftClusterCenter(
									tnegmap[2].GetValue(xiter, yiter));

					//					cout << endl << xiter << " ," << yiter << ", "
					//							<< tposmap[0].GetValue(xiter, yiter) << ", "
					//							<< tposmap[1].GetValue(xiter, yiter) << ", "
					//							<< tposmap[2].GetValue(xiter, yiter)
					//							<< tnegmap[0].GetValue(xiter, yiter) << ", "
					//							<< tnegmap[1].GetValue(xiter, yiter) << ", "
					//							<< tnegmap[2].GetValue(xiter, yiter) << endl;
					//					cout << rvec << endl << flush;
					for (UINT siter = 0; siter < cellhistdim; ++siter) { // pool features...
						if (pmethod == PM_Sum) {
							pltpfeat[0][siter] += rvec[siter];
							pltpfeat[1][siter] += gvec[siter];
							pltpfeat[2][siter] += bvec[siter];

							nltpfeat[0][siter] += nrvec[siter];
							nltpfeat[1][siter] += ngvec[siter];
							nltpfeat[2][siter] += nbvec[siter];
						} else /// Max pooling of features...
						{

							pltpfeat[0][siter]
									= MAX(pltpfeat[0][siter] ,rvec[siter]);
							pltpfeat[1][siter]
									= MAX(pltpfeat[1][siter] ,gvec[siter]);
							pltpfeat[1][siter]
									= MAX(pltpfeat[2][siter] ,bvec[siter]);

							nltpfeat[0][siter]
									= MAX(nltpfeat[0][siter] ,nrvec[siter]);
							nltpfeat[1][siter]
									= MAX(nltpfeat[1][siter] ,ngvec[siter]);
							nltpfeat[2][siter]
									= MAX(nltpfeat[2][siter] ,nbvec[siter]);
						}

					}
				}
			//			cout << pltpfeat[0] << endl;
			//			cout << pltpfeat[1] << endl;
			//			cout << pltpfeat[2] << endl;
			//
			//			cout << nltpfeat[0] << endl;
			//			cout << nltpfeat[1] << endl;
			//			cout << nltpfeat[2] << endl << flush;
			{ // normalization

				REAL *posptr = &features[offset], *negptr = &features[offset
						+ cellhistdim];

				possum[0] = EPS;
				negsum[0] = EPS;

				if (norm == L2) {
					for (UINT biter = 0; biter < cellhistdim; ++biter) { // add the three channels
						if (pim.GetChannel() != RGB) {
							posptr[biter] = pltpfeat[0][biter];
							negptr[biter] = nltpfeat[0][biter];
						} else {
							posptr[biter] = pltpfeat[0][biter]
									+ pltpfeat[1][biter] + pltpfeat[2][biter];
							negptr[biter] = nltpfeat[0][biter]
									+ nltpfeat[1][biter] + nltpfeat[2][biter];
						}
						possum[0] += posptr[biter] * posptr[biter];
						negsum[0] += negptr[biter] * negptr[biter];
					}
					possum[0] = sqrt(possum[0]);
					negsum[0] = sqrt(negsum[0]);
					for (UINT biter = 0; biter < cellhistdim; ++biter) {
						posptr[biter] = (posptr[biter] / possum[0]);
						negptr[biter] = (negptr[biter] / negsum[0]);
					}
				} else if (norm == NT_Average) {
					for (UINT biter = 0; biter < cellhistdim; ++biter) { // add the three channels
						posptr[biter] = (pltpfeat[0][biter]
								+ pltpfeat[1][biter] + pltpfeat[2][biter])
								/ 3.0;
						negptr[biter] = (nltpfeat[0][biter]
								+ nltpfeat[1][biter] + nltpfeat[2][biter])
								/ 3.0;
					}
				} else {
					for (UINT biter = 0; biter < cellhistdim; ++biter) { // add the three channels
						posptr[biter] = pltpfeat[0][biter] + pltpfeat[1][biter]
								+ pltpfeat[2][biter];
						negptr[biter] = nltpfeat[0][biter] + nltpfeat[1][biter]
								+ nltpfeat[2][biter];

						possum[0] += posptr[biter];
						negsum[0] += negptr[biter];
					}
					if (norm == L1) {
						for (UINT biter = 0; biter < cellhistdim; ++biter) {
							posptr[biter] = (posptr[biter] / possum[0]);
							negptr[biter] = (negptr[biter] / negsum[0]);
						}
					} else {
						for (UINT biter = 0; biter < cellhistdim; ++biter) {
							posptr[biter] = sqrt(posptr[biter] / possum[0]);
							negptr[biter] = sqrt(negptr[biter] / negsum[0]);
						}
					}
				}
				//				Print(posptr, posptr + cellhistdim);
				//				Print(negptr, negptr + cellhistdim);
			}
			for (UINT i = 0; i < 3; ++i) { // Zero out the temporary  histograms...
				fill(pltpfeat[i].begin(), pltpfeat[i].end(), 0);
				fill(nltpfeat[i].begin(), nltpfeat[i].end(), 0);
			}
			offset += celldim;
		}
	}
	//	cout << "final feature vecotrs" << features << endl << flush;
}

void GaborLQPFeatures::ExtractCodeBookInfo(vector<Image> &image, LBPMap *map,
		LBPMap *nmap, int sx, int sy, int ex, int ey, const bool &gencodebook,
		bool computeonlycodes) {
	UINT nimages = image.size(); // number of images
	long xdim = image[0].columns(), ydim = image[0].rows(), rcodeval, gcodeval,
			bcodeval, nrcodeval, ngcodeval, nbcodeval, reg = 1;
	int nypoints = ydim - 2 * maxradius, nxpoints = xdim - 2 * maxradius, cenx,
			ceny, count, fx, fy, cx, cy;
	REAL tmpx, tmpy, fracx, fracy, pixval[3], cenval[3];
	REAL x, y, w1, w2, w3, w4;
	vector<vector<REAL> > rimpix(nimages, vector<REAL> (xdim * ydim, 0)),
			gimpix(nimages, vector<REAL> (xdim * ydim, 0)), bimpix(nimages,
					vector<REAL> (xdim * ydim, 0));

	// Load the images
	for (UINT count = 0; count < nimages; ++count) {
		pim.Process(image[count], rimpix[count], gimpix[count], bimpix[count]);
	}

	if (!gencodebook) // extract the map info
		for (int i = 0; i < 3; ++i) // for 3 channels RGB
		{
			map[i].Init(nxpoints, nypoints, -1);
			nmap[i].Init(nxpoints, nypoints, -1);
		}
	vector<REAL> pixelvalues(npoints * 3, 0);
	REAL *rpixval = &pixelvalues[0], *gpixval = rpixval + npoints, *bpixval =
			gpixval + npoints;
	//#ifdef ONLINE_MAPPING ///> use online mapping of codes,
	vector<REAL> npixelvalues(npoints * 3, 0);
	REAL *nrpixval = &npixelvalues[0], *ngpixval = nrpixval + npoints,
			*nbpixval = ngpixval + npoints;
	//#endif
	/*Coordinates configurations for a window codebook*/
	nypoints = MIN(ey-maxradius+1, ydim - 2 * maxradius);
	nxpoints = MIN(ex-maxradius+1, xdim - 2 * maxradius);
	sx = MAX(sx-(int)maxradius,0);
	sy = MAX(sy-(int)maxradius,0);
	for (int i = sy; i < nypoints; i += patchstride) { // patch top-left coordinates
		for (int j = sx; j < nxpoints; j += patchstride) {
			rcodeval = gcodeval = bcodeval = 0;
			nrcodeval = ngcodeval = nbcodeval = 0;

			cenx = j + maxradius; // radius[k], maxradius
			ceny = i + maxradius;

			cenval[0] = rimpix[nglevels][ceny * xdim + cenx];
			cenval[1] = gimpix[nglevels][ceny * xdim + cenx];
			cenval[2] = bimpix[nglevels][ceny * xdim + cenx];
			count = 0;
			//#ifdef ONLINE_MAPPING ///> use online mapping of codes,
			fill(pixelvalues.begin(), pixelvalues.end(), 0);
			fill(npixelvalues.begin(), npixelvalues.end(), 0);
			//#endif
			for (int imgcount = 0; imgcount < nimages; ++imgcount) {

				if (imgcount != nglevels) { // for central pixels

					// replace with else statement.
					pixval[0] = rimpix[imgcount][ceny * xdim + cenx];
					pixval[1] = gimpix[imgcount][ceny * xdim + cenx];
					pixval[2] = bimpix[imgcount][ceny * xdim + cenx];
					if (petype == PET_LBP) {
						UINT tval = reg << count;
						if (pixval[0] >= cenval[0])
							rcodeval += tval;

						if (pixval[1] >= cenval[1])
							gcodeval += tval;

						if (pixval[2] >= cenval[2])
							bcodeval += tval;
					} else if (petype == PET_SplitLTP) {

						//#ifdef ONLINE_MAPPING
						long tval = reg << count;
						if (pixval[0] > cenval[0] + tollevels[0]) {
							rpixval[count] = 1;
							rcodeval += tval;
						} else if (pixval[0] < cenval[0] - tollevels[0]) {
							nrpixval[count] = 1;
							nrcodeval += tval;
						}

						if (pixval[1] > cenval[1] + tollevels[0]) {
							gpixval[count] = 1;
							gcodeval += tval;
						} else if (pixval[1] < cenval[1] - tollevels[0]) {
							ngpixval[count] = 1;
							ngcodeval += tval;
						}

						if (pixval[2] > cenval[2] + tollevels[0]) {
							bpixval[count] = 1;
							bcodeval += tval;
						} else if (pixval[2] < cenval[2] - tollevels[0]) {
							nbpixval[count] = 1;
							nbcodeval += tval;
						}

					} else {
						rpixval[count] = pixval[0];
						gpixval[count] = pixval[1];
						bpixval[count] = pixval[2];
					}
					count++;
				}
				for (vector<vector<Point<REAL> > >::iterator oiter =
						spoints.begin(); oiter != spoints.end(); ++oiter)
					for (vector<Point<REAL> >::iterator iter = oiter->begin(); iter
							!= oiter->end(); ++iter, ++count) {
						iter->GetCoord(x, y);
						tmpy = ceny + y;
						tmpx = cenx + x;

						if ((ptype == PT_Circular && gridtype == Circular)
								|| (ptype == PT_HVCStrip || ptype
										== PT_DACStrip) || ptype
								== PT_CircularShallow) { // do bilinear interpolation

							fx = floor(tmpx);
							fy = floor(tmpy);
							cx = ceil(tmpx);
							cy = ceil(tmpy);
							fracx = tmpx - fx;
							fracy = tmpy - fy;

							if (ABS(fracx) < 1e-6 && ABS(fracy) < 1e-6) {
								// no interpolation needed
								pixval[0] = rimpix[imgcount][fy * xdim + fx];
								pixval[1] = gimpix[imgcount][fy * xdim + fx];
								pixval[2] = bimpix[imgcount][fy * xdim + fx];
							} else {

								w1 = (1 - fracx) * (1 - fracy);
								w2 = fracx * (1 - fracy);
								w3 = (1 - fracx) * fracy;
								w4 = fracx * fracy;

								pixval[0] = w1 * rimpix[imgcount][fy * xdim
										+ fx] + w2 * rimpix[imgcount][fy * xdim
										+ cx] + w3 * rimpix[imgcount][cy * xdim
										+ fx] + w4 * rimpix[imgcount][cy * xdim
										+ cx];

								pixval[1] = w1 * gimpix[imgcount][fy * xdim
										+ fx] + w2 * gimpix[imgcount][fy * xdim
										+ cx] + w3 * gimpix[imgcount][cy * xdim
										+ fx] + w4 * gimpix[imgcount][cy * xdim
										+ cx];

								pixval[2] = w1 * bimpix[imgcount][fy * xdim
										+ fx] + w2 * bimpix[imgcount][fy * xdim
										+ cx] + w3 * bimpix[imgcount][cy * xdim
										+ fx] + w4 * bimpix[imgcount][cy * xdim
										+ cx];
							}
							// posvalue and negvalue contains the thresholded values for r,g & b
							// tcount is used to run over the values...
						}
						// replace with else statement.
						if ((ptype == PT_Circular && gridtype == Rectangular)
								|| (ptype != PT_Circular && ptype
										!= PT_CircularShallow)) {
							pixval[0] = rimpix[imgcount][tmpy * xdim + tmpx];
							pixval[1] = gimpix[imgcount][tmpy * xdim + tmpx];
							pixval[2] = bimpix[imgcount][tmpy * xdim + tmpx];
						}
						if (petype == PET_LBP) {
							UINT tval = reg << count;
							if (pixval[0] >= cenval[0])
								rcodeval += tval;

							if (pixval[1] >= cenval[1])
								gcodeval += tval;

							if (pixval[2] >= cenval[2])
								bcodeval += tval;
						} else if (petype == PET_SplitLTP) {

							//#ifdef ONLINE_MAPPING
							long tval = reg << count;
							if (pixval[0] > cenval[0] + tollevels[0]) {
								rpixval[count] = 1;
								rcodeval += tval;
							} else if (pixval[0] < cenval[0] - tollevels[0]) {
								nrpixval[count] = 1;
								nrcodeval += tval;
							}

							if (pixval[1] > cenval[1] + tollevels[0]) {
								gpixval[count] = 1;
								gcodeval += tval;
							} else if (pixval[1] < cenval[1] - tollevels[0]) {
								ngpixval[count] = 1;
								ngcodeval += tval;
							}

							if (pixval[2] > cenval[2] + tollevels[0]) {
								bpixval[count] = 1;
								bcodeval += tval;
							} else if (pixval[2] < cenval[2] - tollevels[0]) {
								nbpixval[count] = 1;
								nbcodeval += tval;
							}

						} else {
							rpixval[count] = pixval[0];
							gpixval[count] = pixval[1];
							bpixval[count] = pixval[2];
						}
					}
			}
			assert(count==npoints);
			if (computeonlycodes) {
				map[0](j, i) = (rcodeval); // for r,g and b channels..
				map[1](j, i) = (gcodeval);
				map[2](j, i) = (bcodeval);

				/*if (ncodeinfo)*/{
					nmap[0](j, i) = (nrcodeval); // for r,g and b channels..
					nmap[1](j, i) = (ngcodeval);
					nmap[2](j, i) = (nbcodeval);
				}

			} else if (gencodebook) {
#ifdef ONLINE_MAPPING
				codeinfo->SetValue(pixelvalues, rcodeval, gcodeval, bcodeval);
				if (ncodeinfo)
				ncodeinfo->SetValue(npixelvalues, nrcodeval, ngcodeval,
						nbcodeval);
#else
				codeinfo->SetValue(rcodeval, gcodeval, bcodeval);
				if (ncodeinfo)
					ncodeinfo->SetValue(nrcodeval, ngcodeval, nbcodeval);
#endif
			} else {
				if (!softq) { // map to codebook centers...
#if defined ONLINE_MAPPING || defined LIVE_QUANTIZATION
					map[0](j, i)
					= codeinfo->GetClusterCenter(rpixval, rcodeval); // for r,g and b channels..
					map[1](j, i)
					= codeinfo->GetClusterCenter(gpixval, gcodeval);
					map[2](j, i)
					= codeinfo->GetClusterCenter(bpixval, bcodeval);

					if (ncodeinfo) {
						nmap[0](j, i) = ncodeinfo->GetClusterCenter(nrpixval,
								nrcodeval); // for r,g and b channels..
						nmap[1](j, i) = ncodeinfo->GetClusterCenter(ngpixval,
								ngcodeval);
						nmap[2](j, i) = ncodeinfo->GetClusterCenter(nbpixval,
								nbcodeval);
					}
#else
					map[0](j, i) = codeinfo->GetClusterCenter(rcodeval); // for r,g and b channels..
					map[1](j, i) = codeinfo->GetClusterCenter(gcodeval);
					map[2](j, i) = codeinfo->GetClusterCenter(bcodeval);

					if (ncodeinfo) {
						nmap[0](j, i) = ncodeinfo->GetClusterCenter(nrcodeval); // for r,g and b channels..
						nmap[1](j, i) = ncodeinfo->GetClusterCenter(ngcodeval);
						nmap[2](j, i) = ncodeinfo->GetClusterCenter(nbcodeval);
					}
#endif
				} else {
					map[0](j, i) = (rcodeval); // for r,g and b channels..
					map[1](j, i) = (gcodeval);
					map[2](j, i) = (bcodeval);

					if (ncodeinfo) {
						nmap[0](j, i) = (nrcodeval); // for r,g and b channels..
						nmap[1](j, i) = (ngcodeval);
						nmap[2](j, i) = (nbcodeval);
					}
				}
			}
		}
	}
}
void SeparateGaborLQPFeatures::ExtractCodeBookInfo(vector<Image> &image,
		LBPMap *map, LBPMap *nmap, int sx, int sy, int ex, int ey,
		const bool &gencodebook, bool computeonlycodes) {
	UINT nimages = image.size(); // number of images
	long xdim = image[0].columns(), ydim = image[0].rows(), rcodeval, gcodeval,
			bcodeval, nrcodeval, ngcodeval, nbcodeval, reg = 1;
	int nypoints = ydim - 2 * maxradius, nxpoints = xdim - 2 * maxradius, cenx,
			ceny, count, fx, fy, cx, cy;
	REAL tmpx, tmpy, fracx, fracy, pixval[3], cenval[3];
	REAL x, y, w1, w2, w3, w4;
	vector<vector<REAL> > rimpix(nimages, vector<REAL> (xdim * ydim, 0)),
			gimpix(nimages, vector<REAL> (xdim * ydim, 0)), bimpix(nimages,
					vector<REAL> (xdim * ydim, 0));

	// Load the images
	for (UINT count = 0; count < nimages; ++count) {
		pim.Process(image[count], rimpix[count], gimpix[count], bimpix[count]);
	}

	if (!gencodebook) // extract the map info
		for (int i = 0; i < 3; ++i) // for 3 channels RGB
		{
			map[i].Init(nxpoints, nypoints, -1);
			nmap[i].Init(nxpoints, nypoints, -1);
		}
	vector<REAL> pixelvalues(npoints * 3, 0);
	REAL *rpixval = &pixelvalues[0], *gpixval = rpixval + npoints, *bpixval =
			gpixval + npoints;
	//#ifdef ONLINE_MAPPING ///> use online mapping of codes,
	vector<REAL> npixelvalues(npoints * 3, 0);
	REAL *nrpixval = &npixelvalues[0], *ngpixval = nrpixval + npoints,
			*nbpixval = ngpixval + npoints;
	//#endif
	/*Coordinates configurations for a window codebook*/
	nypoints = MIN(ey-maxradius+1, ydim - 2 * maxradius);
	nxpoints = MIN(ex-maxradius+1, xdim - 2 * maxradius);
	sx = MAX(sx-(int)maxradius,0);
	sy = MAX(sy-(int)maxradius,0);
	for (int i = sy; i < nypoints; i += patchstride) { // patch top-left coordinates
		for (int j = sx; j < nxpoints; j += patchstride) {
			rcodeval = gcodeval = bcodeval = 0;
			nrcodeval = ngcodeval = nbcodeval = 0;

			cenx = j + maxradius; // radius[k], maxradius
			ceny = i + maxradius;

			//#ifdef ONLINE_MAPPING ///> use online mapping of codes,
			fill(pixelvalues.begin(), pixelvalues.end(), 0);
			fill(npixelvalues.begin(), npixelvalues.end(), 0);
			//#endif

			count = 0;
			for (int imgcount = 0; imgcount < nimages; ++imgcount) {
				cenval[0] = rimpix[imgcount][ceny * xdim + cenx];
				cenval[1] = gimpix[imgcount][ceny * xdim + cenx];
				cenval[2] = bimpix[imgcount][ceny * xdim + cenx];

				for (vector<vector<Point<REAL> > >::iterator oiter =
						spoints.begin(); oiter != spoints.end(); ++oiter)
					for (vector<Point<REAL> >::iterator iter = oiter->begin(); iter
							!= oiter->end(); ++iter, ++count) {
						iter->GetCoord(x, y);
						tmpy = ceny + y;
						tmpx = cenx + x;

						if ((ptype == PT_Circular && gridtype == Circular)
								|| (ptype == PT_HVCStrip || ptype
										== PT_DACStrip) || ptype
								== PT_CircularShallow) { // do bilinear interpolation

							fx = floor(tmpx);
							fy = floor(tmpy);
							cx = ceil(tmpx);
							cy = ceil(tmpy);
							fracx = tmpx - fx;
							fracy = tmpy - fy;

							if (ABS(fracx) < 1e-6 && ABS(fracy) < 1e-6) {
								// no interpolation needed
								pixval[0] = rimpix[imgcount][fy * xdim + fx];
								pixval[1] = gimpix[imgcount][fy * xdim + fx];
								pixval[2] = bimpix[imgcount][fy * xdim + fx];
							} else {

								w1 = (1 - fracx) * (1 - fracy);
								w2 = fracx * (1 - fracy);
								w3 = (1 - fracx) * fracy;
								w4 = fracx * fracy;

								pixval[0] = w1 * rimpix[imgcount][fy * xdim
										+ fx] + w2 * rimpix[imgcount][fy * xdim
										+ cx] + w3 * rimpix[imgcount][cy * xdim
										+ fx] + w4 * rimpix[imgcount][cy * xdim
										+ cx];

								pixval[1] = w1 * gimpix[imgcount][fy * xdim
										+ fx] + w2 * gimpix[imgcount][fy * xdim
										+ cx] + w3 * gimpix[imgcount][cy * xdim
										+ fx] + w4 * gimpix[imgcount][cy * xdim
										+ cx];

								pixval[2] = w1 * bimpix[imgcount][fy * xdim
										+ fx] + w2 * bimpix[imgcount][fy * xdim
										+ cx] + w3 * bimpix[imgcount][cy * xdim
										+ fx] + w4 * bimpix[imgcount][cy * xdim
										+ cx];
							}
							// posvalue and negvalue contains the thresholded values for r,g & b
							// tcount is used to run over the values...
						}
						// replace with else statement.
						if ((ptype == PT_Circular && gridtype == Rectangular)
								|| (ptype != PT_Circular && ptype
										!= PT_CircularShallow)) {
							pixval[0] = rimpix[imgcount][tmpy * xdim + tmpx];
							pixval[1] = gimpix[imgcount][tmpy * xdim + tmpx];
							pixval[2] = bimpix[imgcount][tmpy * xdim + tmpx];
						}
						if (petype == PET_LBP) {
							UINT tval = reg << count;
							if (pixval[0] >= cenval[0])
								rcodeval += tval;

							if (pixval[1] >= cenval[1])
								gcodeval += tval;

							if (pixval[2] >= cenval[2])
								bcodeval += tval;
						} else if (petype == PET_SplitLTP) {

							//#ifdef ONLINE_MAPPING
							long tval = reg << count;
							if (pixval[0] > cenval[0] + tollevels[0]) {
								rpixval[count] = 1;
								rcodeval += tval;
							} else if (pixval[0] < cenval[0] - tollevels[0]) {
								nrpixval[count] = 1;
								nrcodeval += tval;
							}

							if (pixval[1] > cenval[1] + tollevels[0]) {
								gpixval[count] = 1;
								gcodeval += tval;
							} else if (pixval[1] < cenval[1] - tollevels[0]) {
								ngpixval[count] = 1;
								ngcodeval += tval;
							}

							if (pixval[2] > cenval[2] + tollevels[0]) {
								bpixval[count] = 1;
								bcodeval += tval;
							} else if (pixval[2] < cenval[2] - tollevels[0]) {
								nbpixval[count] = 1;
								nbcodeval += tval;
							}

						} else {
							rpixval[count] = pixval[0];
							gpixval[count] = pixval[1];
							bpixval[count] = pixval[2];
						}
					}
			}
			assert(count==npoints);
			if (computeonlycodes) {
				map[0](j, i) = (rcodeval); // for r,g and b channels..
				map[1](j, i) = (gcodeval);
				map[2](j, i) = (bcodeval);

				/*if (ncodeinfo)*/{
					nmap[0](j, i) = (nrcodeval); // for r,g and b channels..
					nmap[1](j, i) = (ngcodeval);
					nmap[2](j, i) = (nbcodeval);
				}

			} else if (gencodebook) {
#ifdef ONLINE_MAPPING
				codeinfo->SetValue(pixelvalues, rcodeval, gcodeval, bcodeval);
				if (ncodeinfo)
				ncodeinfo->SetValue(npixelvalues, nrcodeval, ngcodeval,
						nbcodeval);
#else
				codeinfo->SetValue(rcodeval, gcodeval, bcodeval);
				if (ncodeinfo)
					ncodeinfo->SetValue(nrcodeval, ngcodeval, nbcodeval);
#endif
			} else {
				if (!softq) { // map to codebook centers...
#if defined ONLINE_MAPPING || defined LIVE_QUANTIZATION
					map[0](j, i)
					= codeinfo->GetClusterCenter(rpixval, rcodeval); // for r,g and b channels..
					map[1](j, i)
					= codeinfo->GetClusterCenter(gpixval, gcodeval);
					map[2](j, i)
					= codeinfo->GetClusterCenter(bpixval, bcodeval);

					if (ncodeinfo) {
						nmap[0](j, i) = ncodeinfo->GetClusterCenter(nrpixval,
								nrcodeval); // for r,g and b channels..
						nmap[1](j, i) = ncodeinfo->GetClusterCenter(ngpixval,
								ngcodeval);
						nmap[2](j, i) = ncodeinfo->GetClusterCenter(nbpixval,
								nbcodeval);
					}
#else
					map[0](j, i) = codeinfo->GetClusterCenter(rcodeval); // for r,g and b channels..
					map[1](j, i) = codeinfo->GetClusterCenter(gcodeval);
					map[2](j, i) = codeinfo->GetClusterCenter(bcodeval);

					if (ncodeinfo) {
						nmap[0](j, i) = ncodeinfo->GetClusterCenter(nrcodeval); // for r,g and b channels..
						nmap[1](j, i) = ncodeinfo->GetClusterCenter(ngcodeval);
						nmap[2](j, i) = ncodeinfo->GetClusterCenter(nbcodeval);
					}
#endif
				} else {
					map[0](j, i) = (rcodeval); // for r,g and b channels..
					map[1](j, i) = (gcodeval);
					map[2](j, i) = (bcodeval);

					if (ncodeinfo) {
						nmap[0](j, i) = (nrcodeval); // for r,g and b channels..
						nmap[1](j, i) = (ngcodeval);
						nmap[2](j, i) = (nbcodeval);
					}
				}
			}
		}
	}
}
