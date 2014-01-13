/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/

#include "multipleCBQuantizedLQPFeatures.h"
void MultipleCBQuantizeLQPFeatures::InitalizeMaps(Image &image) {
	sspace = NoPyramid;
	UINT xbmax, ybmax;
	lbpfeatures.Initialize(1, 0, featdim);
	GetXYBlocks(image, xbmax, ybmax);
	lbpfeatures.SetFeature(xbmax, ybmax, 1.0);
	ComputeLBPFeatures(image, 0, cell, lbpstride);
}

void MultipleCBQuantizeLQPFeatures::InitalizeMaps(Pyramid &pyobj,
		PyramidType sspace_) {
	UINT pycount = 0;
	sspace = sspace_;
	// Read the image....
	if (sspace == SingleRes) {
		pycount = pyobj.GetNLevels();
		lbpfeatures.Initialize(pycount, 0, featdim);
#ifdef PATCH_THREADS
		boost::thread_group tgroup;
		for (UINT i = 0; i < pycount; i++) {

			Image& imgref = pyobj.GetImageRef(i);
			UINT xbmax, ybmax;
			// Compute the Gradient Image
			GetXYBlocks(imgref, xbmax, ybmax);
			lbpfeatures.SetFeature(xbmax, ybmax, pyobj.GetScale(i));
			vector<REAL>& features = lbpfeatures.GetFeatureRef(i);
			tgroup .add_thread(new boost::thread(boost::bind(
									&MultipleCBQuantizeLQPFeatures::ComputeLBPFeatures,
									boost::ref(*this), boost::ref(imgref),
									boost::ref(features), xbmax, ybmax, cell, lbpstride)));
			//			tgroup.add_thread(new boost::thread(boost::bind(
			//									&MultipleCBQuantizeLQPFeatures::ComputeLBPFeatures,
			//									boost::ref(*this), boost::ref(imgref), i, cell, lbpstride)));

		}
		tgroup.join_all();
#ifdef TMPDEBUG // for testing of features
		HOGInfo tmplbpfeatures;
		tmplbpfeatures.Initialize(pycount,0,featdim);
		for (UINT i = 0; i < pycount; i++) {
			Image& imgref = pyobj.GetImageRef(i);
			// Compute the Gradient Image
			UINT xbmax, ybmax;
			GetXYBlocks(imgref, xbmax, ybmax);
			tmplbpfeatures.SetFeature(xbmax, ybmax, pyobj.GetScale(i));
			vector<REAL>& features = tmplbpfeatures.GetFeatureRef(i);
			ComputeLBPFeatures(imgref,features,xbmax,ybmax, cell, lbpstride);
		}

		for (UINT i = 0; i < pycount; i++) {
			vector<REAL> &tfeat=lbpfeatures.GetFeatureRef(i),
			&feat=tmplbpfeatures.GetFeatureRef(i);

			for(UINT j=0; j < tfeat.size(); ++j)
			if(tfeat[j]!=feat[j])
			{
				cout<<" Discrepency caught";

			}
		}
#endif

#else

		for (UINT i = 0; i < pycount; i++) {
			Image& imgref = pyobj.GetImageRef(i);
			// Compute the Gradient Image
			UINT xbmax, ybmax;
			GetXYBlocks(imgref, xbmax, ybmax);
			lbpfeatures.SetFeature(xbmax, ybmax, pyobj.GetScale(i));
			ComputeLBPFeatures(imgref, i, cell, lbpstride);
		}
#endif

	} else /*if (sspace == DoubleRes)*/{
		cout << "ERROR: Only LTP Not Implemented In Detail; ";
		exit(EXIT_FAILURE);
	}
}
void MultipleCBQuantizeLQPFeatures::InitalizeMaps(Image&image,
		PyramidType sspace_) {
	UINT pycount = 0;
	sspace = sspace_;
	// Read the image....
	if (sspace == SingleRes) {
		pycount = pyobj.GeneratePyramid(image);
		lbpfeatures.Initialize(pycount, 0, featdim);
#ifdef PATCH_THREADS
		boost::thread_group tgroup;
		for (UINT i = 0; i < pycount; i++) {

			Image& imgref = pyobj.GetImageRef(i);
			UINT xbmax, ybmax;
			// Compute the Gradient Image
			GetXYBlocks(imgref, xbmax, ybmax);
			lbpfeatures.SetFeature(xbmax, ybmax, pyobj.GetScale(i));
			vector<REAL>& features = lbpfeatures.GetFeatureRef(i);
			tgroup .add_thread(new boost::thread(boost::bind(
									&MultipleCBQuantizeLQPFeatures::ComputeLBPFeatures,
									boost::ref(*this), boost::ref(imgref),
									boost::ref(features), xbmax, ybmax, cell, lbpstride)));
			//			tgroup.add_thread(new boost::thread(boost::bind(
			//									&MultipleCBQuantizeLQPFeatures::ComputeLBPFeatures,
			//									boost::ref(*this), boost::ref(imgref), i, cell, lbpstride)));

		}
		tgroup.join_all();
#ifdef TMPDEBUG // for testing of features
		HOGInfo tmplbpfeatures;
		tmplbpfeatures.Initialize(pycount,0,featdim);
		for (UINT i = 0; i < pycount; i++) {
			Image& imgref = pyobj.GetImageRef(i);
			// Compute the Gradient Image
			UINT xbmax, ybmax;
			GetXYBlocks(imgref, xbmax, ybmax);
			tmplbpfeatures.SetFeature(xbmax, ybmax, pyobj.GetScale(i));
			vector<REAL>& features = tmplbpfeatures.GetFeatureRef(i);
			ComputeLBPFeatures(imgref,features,xbmax,ybmax, cell, lbpstride);
		}

		for (UINT i = 0; i < pycount; i++) {
			vector<REAL> &tfeat=lbpfeatures.GetFeatureRef(i),
			&feat=tmplbpfeatures.GetFeatureRef(i);

			for(UINT j=0; j < tfeat.size(); ++j)
			if(tfeat[j]!=feat[j])
			{
				cout<<" Discrepency caught";
				assert(tfeat[j]!=feat[j]);
			}
		}
#endif

#else

		for (UINT i = 0; i < pycount; i++) {
			Image& imgref = pyobj.GetImageRef(i);
			// Compute the Gradient Image
			UINT xbmax, ybmax;
			GetXYBlocks(imgref, xbmax, ybmax);
			lbpfeatures.SetFeature(xbmax, ybmax, pyobj.GetScale(i));
			ComputeLBPFeatures(imgref, i, cell, lbpstride);
		}
#endif
	} else if (sspace == DoubleRes) {
		cout << "ERROR: Only LTP Not Implemented In Detail; ";
		exit(EXIT_FAILURE);
	} else {
		UINT xbmax, ybmax;
		lbpfeatures.Initialize(1, 0, featdim);
		GetXYBlocks(image, xbmax, ybmax);
		lbpfeatures.SetFeature(xbmax, ybmax, 1.0);
		ComputeLBPFeatures(image, 0, cell, lbpstride);
	}
}
void MultipleCBQuantizeLQPFeatures::ExtractPosCodeInfo(Image&image,
		vector<Coord>& ainfo)
// annotation information
{
	// Crop the image according to each annotation information by adding some boundary
	// and then proceed in the normal way to generate the codebook, this generates
	// the code book using local scale space pyramid...
	if (sspace == SingleRes) {
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
				//				ExtractCodeInfo(timage);
				ExtractCodeInfo(image, NULL, pcbooks[0]); // only 2 codebooks are used
				timage.flop();
				ExtractCodeInfo(image, NULL, pcbooks[0]);
				timage = image;
			}
		}
	} else {/*Simply Crop the Positive Windows with the offset and generate code book*/
		UINT rarea = width * height;
		LBPMap map[3];
		int x1, y1, x2, y2, offset = GetMaxFeatureOffset(), minoffset =
				GetFeaturesOffset(), doffset = offset - minoffset;
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
						Bilinear, true);// 2 cell borders
				//				ExtractCodeBookInfo(im, map, true); // Flip the features once again so that in case of odd size it has the rig
				ExtractCodeInfo(im, map, NULL);
				if (ncbs == 2)
					GenerateCodeBookInfo(map, doffset, doffset, pcbooks[0]);
				else
					GenerateCodeBookInfo(map, doffset, doffset, pcbooks);
				im.flop();

				ExtractCodeInfo(im, map, NULL);
				if (ncbs == 2)
					GenerateCodeBookInfo(map, doffset, doffset, pcbooks[0]);
				else
					GenerateCodeBookInfo(map, doffset, doffset, pcbooks);
			}
		}
	}
}

// to use either the sampling from the negative images or the complete image..
void MultipleCBQuantizeLQPFeatures::ExtractNegCodeInfo(Image &image) {

	vector<WindowInfoFeat> winfo;/*Detector window's information*/
	// Generate the Pyramid and built the codeInfo of each level
	UINT wCount = 0, pycount = pyobj.GeneratePyramid(image);

	UINT ybmax, xbmax;
	REAL response = 0, overlap = 0, scale;
	LBPMap (*map)[3] = new LBPMap[pycount][3];
	UINT offset = GetMaxFeatureOffset(), minoffset = GetFeaturesOffset(),
			doffset = offset - minoffset;
	UINT ystride = ycell, xstride = xcell;
	if (nwinsampled == 0) // Use Complete Images for the codeBook
	{
		for (UINT i = 0; i < pycount; i++) {
			Image &imgref = pyobj.GetImageRef(i);
			ExtractCodeInfo(imgref, NULL, pcbooks[ncbs - 1]);
		}
	} else { // Sample Images for the generation of codebook...
		for (UINT i = 0; i < pycount; i++) {
			Image &imgref = pyobj.GetImageRef(i);
			ExtractCodeInfo(imgref, map[i], NULL);
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
			GenerateCodeBookInfo(map[sindex], xmin + doffset, ymin + doffset,
					pcbooks[ncbs - 1]); // SameCodeBook For all the Cells
		}
	}
	delete[] map;
}
// To Use the Same CodeBook for all the cells...

void MultipleCBQuantizeLQPFeatures::GenerateCodeBookInfo(LBPMap *map,
		UINT xmin, UINT ymin, ComputeCode* cbook) {
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

			endx = MIN(x + winstride, map[0].nxpoints);
			endy = MIN(y + winstride, map[0].nypoints);

			for (UINT yiter = y; yiter < endy; ++yiter) // replace ++yiter with +=patch stride...
				for (UINT xiter = x; xiter < endx; ++xiter) {
					rval = map[0].GetValue(xiter, yiter);
					if (rval >= 0) {// to omit the boundary values...
						cbook->SetValue(rval, map[1].GetValue(xiter, yiter),
								map[2].GetValue(xiter, yiter));
					}
				}
		}
	}
}
void MultipleCBQuantizeLQPFeatures::GenerateCodeBookInfo(LBPMap *map,
		UINT xmin, UINT ymin, vector<ComputeCode*> &codebooks) {

	UINT xbmax = xmin + width, ybmax = ymin + height, endx, endy, winstride =
			cell, count = 0;
	ComputeCode *cbook;
	int x, y, rval;
	for (UINT rind = ymin; rind < ybmax; rind += winstride) {
		y = rind + foffset - maxradius; // patch top left coordinates for which rind+foffset is center...
		// because the map contains an offset maxradius for concentic circles..
		// means map(0,0) contain value of map(maxradius,maxradius) so value of hist(offset,offset) can be found at
		// map(offset-maxradius,offset-maxradius)
		for (UINT colind = xmin; colind < xbmax; colind += winstride) {
			x = colind + foffset - maxradius;

			endx = MIN(x + winstride, map[0].nxpoints);
			endy = MIN(y + winstride, map[0].nypoints);
			cbook = codebooks[count++]; /*Separate code book for each cell...*/

			for (UINT yiter = y; yiter < endy; ++yiter) // replace ++yiter with +=patch stride...
				for (UINT xiter = x; xiter < endx; ++xiter) {
					rval = map[0].GetValue(xiter, yiter);
					if (rval >= 0) // to omit the boundary values...
						// because if feature is computed for red -channel
						//!it is implicit that it is computed for green & blue
						cbook->SetValue(rval, map[1].GetValue(xiter, yiter),
								map[2].GetValue(xiter, yiter));
				}

		}
	}
	assert(count==ncbs);
}
void MultipleCBQuantizeLQPFeatures::GenerateCodeBookInfo(LBPMap *pmap,
		LBPMap *nmap, UINT xmin, UINT ymin) {

	UINT xbmax = xmin + width, ybmax = ymin + height, endx, endy, winstride =
			cell, count = 0;
	ComputeCode *pcbook, *ncbook;
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
			pcbook = pcbooks[count]; /*Separate code book for each cell...*/
			ncbook = ncbooks[count++]; /*Separate code book for each cell...*/

			for (UINT yiter = y; yiter < endy; ++yiter) // replace ++yiter with +=patch stride...
				for (UINT xiter = x; xiter < endx; ++xiter) {
					rval = pmap[0].GetValue(xiter, yiter);
					if (rval >= 0) // to omit the boundary values...
					// because if feature is computed for red -channel
					//!it is implicit that it is computed for green & blue
					{
						pcbook->SetValue(rval, pmap[1].GetValue(xiter, yiter),
								pmap[2].GetValue(xiter, yiter));
						ncbook->SetValue(nmap[0].GetValue(xiter, yiter),
								nmap[1].GetValue(xiter, yiter),
								nmap[2].GetValue(xiter, yiter));
					}
				}

		}
	}
	assert(count==ncbs);
}
// For Image level code book, specially if only single code book is to be learned
// Extract Image level codes and then use them for the generation of local cell based code books...
void MultipleCBQuantizeLQPFeatures::ExtractCodeInfo(Image &image, LBPMap *map,
		ComputeCode *cbook) {

	UINT xdim = image.columns(), ydim = image.rows(), rcodeval, gcodeval,
			bcodeval, reg = 1;
	int nypoints = ydim - 2 * maxradius, fx, fy, cx, cy, nxpoints = xdim - 2
			* maxradius, cenx, ceny, count;
	REAL tmpx, tmpy, fracx, fracy, pixval[3], cenval[3];
	REAL x, y, w1, w2, w3, w4;
	vector<REAL> rimpix(xdim * ydim, 0), gimpix(xdim * ydim, 0), bimpix(
			xdim * ydim, 0);
	pim.Process(image, rimpix, gimpix, bimpix);
	if (map) {
		for (int i = 0; i < 3; ++i)// for 3 channels RGB
			map[i].Init(nxpoints, nypoints); // foffset = 1 Cell of
	}
	for (int i = 0; i < nypoints; i += patchstride) { // patch top-left coordinates
		for (int j = 0; j < nxpoints; j += patchstride) {
			rcodeval = gcodeval = bcodeval = 0;

			cenx = j + maxradius;// radius[k], maxradius
			ceny = i + maxradius;

			cenval[0] = rimpix[ceny * xdim + cenx];
			cenval[1] = gimpix[ceny * xdim + cenx];
			cenval[2] = bimpix[ceny * xdim + cenx];
			count = 0;

			for (vector<vector<Point<REAL> > >::iterator oiter =
					spoints.begin(); oiter != spoints.end(); ++oiter)
				for (vector<Point<REAL> >::iterator iter = oiter->begin(); iter
						!= oiter->end(); ++iter, ++count) {
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
					} else {
						pixval[0] = rimpix[tmpy * xdim + tmpx];
						pixval[1] = gimpix[tmpy * xdim + tmpx];
						pixval[2] = bimpix[tmpy * xdim + tmpx];
					}
					if (petype == PET_LBP) {
						int tval = reg << count;
						if (pixval[0] >= cenval[0])
							rcodeval += tval;

						if (pixval[1] >= cenval[1])
							gcodeval += tval;

						if (pixval[2] >= cenval[2])
							bcodeval += tval;
					} else {
						UINT treg[3] = { 0, 0, 0 };
						// posvalue and negvalue contains the thresholded values for r,g & b
						// tcount is used to run over the values...
						for (int tl = 0, tcount = 0; tl < ltplevels; ++tl)
							for (int iter = 0; iter < 3; ++iter, ++tcount) {

								if (pixval[iter] >= cenval[iter]
										+ tollevels[tl]) {
									treg[iter] = tl + 1;
								} else if (pixval[iter] <= cenval[iter]
										- tollevels[tl]) {
									treg[iter] = tl + 1 + ltplevels; // codes will be +ve while vector contains the -ve
								}
							}
						rcodeval += basevalues[count] * treg[0];
						gcodeval += basevalues[count] * treg[1];
						bcodeval += basevalues[count] * treg[2];
					}
				}
			if (map) {
				map[0](j, i) = rcodeval; // for r,g and b channels..
				map[1](j, i) = gcodeval;
				map[2](j, i) = bcodeval;
			}
			if (cbook) {
				cbook->SetValue(rcodeval, gcodeval, bcodeval);
			}
		}
	}
}

/*
 * As this is used with combinations of other features, so zeros are padded to have the
 * same cell dimensions.*/
void MultipleCBQuantizeLQPFeatures::ComputeLBPFeatures(Image& imgref,
		UINT index, UINT winstride, UINT tlbpstride) {
	// for the speed modularization is sacrificed...
	LBPMap map[3];
	// because map is not computed for the

	ExtractCodeInfo(imgref, map, NULL);

	UINT cxbmax = lbpfeatures.GetXBlock(index), cybmax = lbpfeatures.GetYBlock(
			index);
	vector<REAL> &features = lbpfeatures.GetFeatureRef(index);
	ComputeLBPFeatures(imgref, features, cxbmax, cybmax, winstride, tlbpstride);
}

// For patch_threads
void MultipleCBQuantizeLQPFeatures::ComputeLBPFeatures(Image& imgref,
		vector<REAL> &features, UINT cxbmax, UINT cybmax, UINT winstride,
		UINT tlbpstride) {
	// for the speed modularization is sacrificed...
	if (petype == PET_SplitLTP) {
		ComputeSplitLBPFeatures(imgref, features, cxbmax, cybmax, winstride,
				tlbpstride);
		return;
	}

	LBPMap map[3];
	// because map is not computed for the
	ExtractCodeInfo(imgref, map, NULL);
	UINT x, y, xbmax = cxbmax * winstride, ybmax = cybmax * winstride, endx,
			endy, offset = 0;
	int rval;
	vector<vector<REAL> > feat(3, vector<REAL> (ncenters));
	REAL possum[3];
	//FIXME: Very Important for Parts based detector this /*rind + foffset - maxradius */ in featureComputation should be replace by WinStride...
	//TODO: Very Important for Parts based detector this /*rind + foffset - maxradius */  should be replace by WinStride...
	// in all files..boundary .
	int boundary = lfeatoffset - foffset, bcell = boundary / cell;
	xbmax -= boundary;
	ybmax -= boundary;
	ComputeCode *cb;
	for (UINT rind = boundary, ystart = bcell; rind < ybmax; rind += winstride, ++ystart) {
		y = rind + foffset - maxradius; // patch top left coordinates for which rind+foffset is center...
		// because the map contains an offset maxradius for concentic circles..
		// means map(0,0) contain value of map(maxradius,maxradius) so value of hist(offset,offset) can be found at
		// map(offset-maxradius,offset-maxradius)
		offset = (ystart * cxbmax + bcell) * featdim;

		for (UINT colind = boundary; colind < xbmax; colind += winstride) {
			x = colind + foffset - maxradius;
			endx = MIN(x + tlbpstride, map[0].nxpoints);
			endy = MIN(y + tlbpstride, map[0].nypoints);

			// Go over the codebooks one by one... loop over CodeBooks
			for (vector<ComputeCode*>::iterator iter = pcbooks.begin(); iter
					!= pcbooks.end(); ++iter) {
				cb = *iter;
				REAL *posptr = &features[offset];// add offsets in features code ...

				// for each cell
				for (UINT yiter = y; yiter < endy; ++yiter) // replace ++yiter with +=patch stride...
					for (UINT xiter = x; xiter < endx; ++xiter) {

						rval = cb->GetClusterCenter(
								map[0].GetValue(xiter, yiter));
						if (rval >= 0) // because if feature is computed for red -channel
						//!it is implicit that it is computed for green & blue
						{
							feat[0][rval]++;
							feat[1][cb->GetClusterCenter(
									map[1].GetValue(xiter, yiter))]++;
							feat[2][cb->GetClusterCenter(
									map[2].GetValue(xiter, yiter))]++;
						}
					}
				possum[0] = EPS;
				if (norm == L2) {
					for (UINT biter = 0; biter < ncenters; ++biter) { // add the three channels
						posptr[biter] = feat[0][biter] + feat[1][biter]
								+ feat[2][biter];
						possum[0] += posptr[biter] * posptr[biter];
					}
					possum[0] = sqrt(possum[0]);

					for (UINT biter = 0; biter < ncenters; ++biter)
						posptr[biter] = posptr[biter] / possum[0];

				} else {
					for (UINT biter = 0; biter < ncenters; ++biter) { // add the three channels
						posptr[biter] = feat[0][biter] + feat[1][biter]
								+ feat[2][biter];
						possum[0] += posptr[biter];
					}
					for (UINT biter = 0; biter < ncenters; ++biter)
						posptr[biter] = sqrt(posptr[biter] / possum[0]);
				}

				for (UINT i = 0; i < 3; ++i) // Zero out the temporary  histograms...
					fill(feat[i].begin(), feat[i].end(), 0);
				offset += ncenters;
			}
		}
	}
}
void MultipleCBQuantizeLQPFeatures::ComputeSplitLBPFeatures(Image& imgref,
		vector<REAL> &features, UINT cxbmax, UINT cybmax, UINT winstride,
		UINT tlbpstride) {
	// for the speed modularization is sacrificed...

	LBPMap pmap[3], nmap[3];
	// because map is not computed for the
	ExtractCodeBookInfo(imgref, pmap, nmap, false, true);
	UINT x, y, xbmax = cxbmax * winstride, ybmax = cybmax * winstride, endx,
			endy, offset = 0;
	int rval;
	vector<vector<REAL> > pfeat(3, vector<REAL> (ncenters)), nfeat(3,
			vector<REAL> (ncenters));
	REAL possum[3],negsum[3];
	//FIXME: Very Important for Parts based detector this /*rind + foffset - maxradius */ in featureComputation should be replace by WinStride...
	//TODO: Very Important for Parts based detector this /*rind + foffset - maxradius */  should be replace by WinStride...
	// in all files..boundary .
	int boundary = lfeatoffset - foffset, bcell = boundary / cell;
	xbmax -= boundary;
	ybmax -= boundary;
	ComputeCode *pcb, *ncb;
	for (UINT rind = boundary, ystart = bcell; rind < ybmax; rind += winstride, ++ystart) {
		y = rind + foffset - maxradius; // patch top left coordinates for which rind+foffset is center...
		// because the map contains an offset maxradius for concentic circles..
		// means map(0,0) contain value of map(maxradius,maxradius) so value of hist(offset,offset) can be found at
		// map(offset-maxradius,offset-maxradius)
		offset = (ystart * cxbmax + bcell) * featdim;

		for (UINT colind = boundary; colind < xbmax; colind += winstride) {
			x = colind + foffset - maxradius;
			endx = MIN(x + tlbpstride, pmap[0].nxpoints);
			endy = MIN(y + tlbpstride, pmap[0].nypoints);

			// Go over the codebooks one by one... loop over CodeBooks
			for (vector<ComputeCode*>::iterator piter = pcbooks.begin(), niter =
					ncbooks.begin(); piter != pcbooks.end(); ++piter, ++niter) {
				pcb = *piter;
				ncb = *niter;
				REAL *posptr = &features[offset], *negptr = posptr + ncenters;

				// for each cell
				for (UINT yiter = y; yiter < endy; ++yiter) // replace ++yiter with +=patch stride...
					for (UINT xiter = x; xiter < endx; ++xiter) {

						rval = pcb->GetClusterCenter(
								pmap[0].GetValue(xiter, yiter));
						if (rval >= 0) // because if feature is computed for red -channel
						//!it is implicit that it is computed for green & blue
						{
							pfeat[0][rval]++;
							pfeat[1][pcb->GetClusterCenter(
									pmap[1].GetValue(xiter, yiter))]++;
							pfeat[2][pcb->GetClusterCenter(
									pmap[2].GetValue(xiter, yiter))]++;
						}

						rval = ncb->GetClusterCenter(
								nmap[0].GetValue(xiter, yiter));
						if (rval >= 0) // because if feature is computed for red -channel
						//!it is implicit that it is computed for green & blue
						{
							nfeat[0][rval]++;
							nfeat[1][ncb->GetClusterCenter(
									nmap[1].GetValue(xiter, yiter))]++;
							nfeat[2][ncb->GetClusterCenter(
									nmap[2].GetValue(xiter, yiter))]++;
						}
					}
				possum[0] = EPS;
				negsum[0] = EPS;
				if (norm == L2) {
					for (UINT biter = 0; biter < ncenters; ++biter) { // add the three channels
						posptr[biter] = pfeat[0][biter] + pfeat[1][biter]
								+ pfeat[2][biter];
						possum[0] += posptr[biter] * posptr[biter];

						negptr[biter] = nfeat[0][biter] + nfeat[1][biter]
								+ nfeat[2][biter];
						negsum[0] += negptr[biter] * negptr[biter];
					}
					possum[0] = sqrt(possum[0]);
					negsum[0] = sqrt(negsum[0]);

					for (UINT biter = 0; biter < ncenters; ++biter) {
						posptr[biter] = posptr[biter] / possum[0];
						negptr[biter] = negptr[biter] / negsum[0];
					}

				} else {
					for (UINT biter = 0; biter < ncenters; ++biter) { // add the three channels
						posptr[biter] = pfeat[0][biter] + pfeat[1][biter]
								+ pfeat[2][biter];
						possum[0] += posptr[biter];

						negptr[biter] = nfeat[0][biter] + nfeat[1][biter]
								+ nfeat[2][biter];
						negsum[0] += negptr[biter];
					}
					for (UINT biter = 0; biter < ncenters; ++biter) {
						posptr[biter] = sqrt(posptr[biter] / possum[0]);
						negptr[biter] = sqrt(negptr[biter] / negsum[0]);
					}
				}

				for (UINT i = 0; i < 3; ++i) // Zero out the temporary  histograms...
				{
					fill(pfeat[i].begin(), pfeat[i].end(), 0);
					fill(nfeat[i].begin(), nfeat[i].end(), 0);
				}
				offset += 2 * ncenters;
			}
		}
	}
}

void MultipleCBQuantizeLQPFeatures::GetFeatures(UINT index, int x, int y,
		int width_, int height_, vector<REAL>&feat) {// There is offset of one cell so, input x=0 means x=8 by taking
//	into account offset and the removed border
	int xcell = x / cell, ycell = y / cell;
	width_ /= cell;
	height_ /= cell;
	lbpfeatures.GetHOGCells(index, xcell, ycell, width_, height_, &feat[0]);
}
void MultipleCBQuantizeLQPFeatures::GetFeatures(UINT index, int x, int y,
		int width_, int height_, REAL*feat) {// There is offset of one cell so, input x=0 means x=8 by taking
//	into account offset and the removed border
	int xcell = x / cell, ycell = y / cell;
	width_ /= cell;
	height_ /= cell;
	lbpfeatures.GetHOGCells(index, xcell, ycell, width_, height_, feat);
}
/*
 * Contains only folded HOG's....*/
void MultipleCBQuantizeLQPFeatures::GetFoldedFeatures(UINT index, int x, int y,
		int twidth, int theight, vector<REAL>&feat) {// There is offset of one cell so, input x=0 means x=8 by taking
//	into account offset and the removed border
	int xcell = x / cell, ycell = y / cell;
	twidth /= cell;
	theight /= cell;
	lbpfeatures.GetHOGCells(index, xcell, ycell, twidth, theight, &feat[0]);
}

/*
 * Contains only folded HOG's....*/
void MultipleCBQuantizeLQPFeatures::GetFoldedFeatures(UINT index, int x, int y,
		int twidth, int theight, REAL*feat) {// There is offset of one cell so, input x=0 means x=8 by taking
//	into account offset and the removed border
	int xcell = x / cell, ycell = y / cell;
	twidth /= cell;
	theight /= cell;
	lbpfeatures.GetHOGCells(index, xcell, ycell, twidth, theight, feat);
}

void MultipleCBQuantizeLQPFeatures::PadFeatureMap(UINT index, UINT padx,
		UINT pady) { // does padding in the pixel space
// first pad the hogFeatures map...
	padx /= cell;
	pady /= cell;
	lbpfeatures.PadFeatureMap(index, padx, pady);
}

void MultipleCBQuantizeLQPFeatures::DotProduct(UINT index, UINT fwidth,
		UINT fheight, vector<REAL>&filter, Store&response) {
	fwidth /= cell;
	fheight /= cell;
	lbpfeatures.DotProduct(index, fwidth, fheight, &filter[0], response);
}

void MultipleCBQuantizeLQPFeatures::DotProduct(UINT index, UINT fwidth,
		UINT fheight, REAL *filter, Store&response) {
	fwidth /= cell;
	fheight /= cell;
	lbpfeatures.DotProduct(index, fwidth, fheight, filter, response);
}
/*
 * void MultipleCBQuantizeLQPFeatures::ExtractCodeInfo(Image &image, LBPMap *map) {

 UINT xdim = image.columns(), ydim = image.rows(), rcodeval, gcodeval,
 bcodeval, reg = 1;
 int nypoints = ydim - 2 * maxradius, fx, fy, cx, cy, nxpoints = xdim - 2
 * maxradius, cenx, ceny, count;
 REAL tmpx, tmpy, fracx, fracy, pixval[3], cenval[3];
 REAL x, y, w1, w2, w3, w4;
 vector<REAL> rimpix(xdim * ydim, 0), gimpix(xdim * ydim, 0), bimpix(
 xdim * ydim, 0);
 pim.Process(image, rimpix, gimpix, bimpix);
 for (int i = 0; i < 3; ++i)// for 3 channels RGB
 map[i].Init(nxpoints, nypoints); // foffset = 1 Cell of

 for (int i = 0; i < nypoints; i += patchstride) { // patch top-left coordinates
 for (int j = 0; j < nxpoints; j += patchstride) {
 rcodeval = gcodeval = bcodeval = 0;

 cenx = j + maxradius;// radius[k], maxradius
 ceny = i + maxradius;

 cenval[0] = rimpix[ceny * xdim + cenx];
 cenval[1] = gimpix[ceny * xdim + cenx];
 cenval[2] = bimpix[ceny * xdim + cenx];
 count = 0;

 for (vector<Point<REAL> >::iterator iter = spoints[0].begin(); iter
 != spoints[0].end(); ++iter, ++count) {
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
 } else {
 pixval[0] = rimpix[tmpy * xdim + tmpx];
 pixval[1] = gimpix[tmpy * xdim + tmpx];
 pixval[2] = bimpix[tmpy * xdim + tmpx];
 }
 if (nlevels == 2) {
 int tval = reg << count;
 if (pixval[0] >= cenval[0])
 rcodeval += tval;

 if (pixval[1] >= cenval[1])
 gcodeval += tval;

 if (pixval[2] >= cenval[2])
 bcodeval += tval;
 } else {
 UINT treg[3] = { 0, 0, 0 };
 // posvalue and negvalue contains the thresholded values for r,g & b
 // tcount is used to run over the values...
 for (int tl = 0, tcount = 0; tl < ltplevels; ++tl)
 for (int iter = 0; iter < 3; ++iter, ++tcount) {

 if (pixval[iter] >= cenval[iter] + tollevels[tl]) {
 treg[iter] = tl + 1;
 } else if (pixval[iter] <= cenval[iter]
 - tollevels[tl]) {
 treg[iter] = tl + 1 + ltplevels; // codes will be +ve while vector contains the -ve
 }
 }
 rcodeval += basevalues[count] * treg[0];
 gcodeval += basevalues[count] * treg[1];
 bcodeval += basevalues[count] * treg[2];

 }
 }
 map[0](j, i) = rcodeval; // for r,g and b channels..
 map[1](j, i) = gcodeval;
 map[2](j, i) = bcodeval;
 }
 }
 }
 void MultipleCBQuantizeLQPFeatures::GenerateCodeBookInfo(Image &image,
 ComputeCode *cbook) {

 UINT xdim = image.columns(), ydim = image.rows(), rcodeval, gcodeval,
 bcodeval, reg = 1;
 int nypoints = ydim - 2 * maxradius, fx, fy, cx, cy, nxpoints = xdim - 2
 * maxradius, cenx, ceny, count;
 REAL tmpx, tmpy, fracx, fracy, pixval[3], cenval[3];
 REAL x, y, w1, w2, w3, w4;
 vector<REAL> rimpix(xdim * ydim, 0), gimpix(xdim * ydim, 0), bimpix(
 xdim * ydim, 0);
 pim.Process(image, rimpix, gimpix, bimpix);
 for (int i = 0; i < nypoints; i += patchstride) { // patch top-left coordinates
 for (int j = 0; j < nxpoints; j += patchstride) {
 rcodeval = gcodeval = bcodeval = 0;

 cenx = j + maxradius;// radius[k], maxradius
 ceny = i + maxradius;

 cenval[0] = rimpix[ceny * xdim + cenx];
 cenval[1] = gimpix[ceny * xdim + cenx];
 cenval[2] = bimpix[ceny * xdim + cenx];
 count = 0;

 for (vector<Point<REAL> >::iterator iter = spoints[0].begin(); iter
 != spoints[0].end(); ++iter, ++count) {
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
 } else {
 pixval[0] = rimpix[tmpy * xdim + tmpx];
 pixval[1] = gimpix[tmpy * xdim + tmpx];
 pixval[2] = bimpix[tmpy * xdim + tmpx];
 }
 if (petype == PET_LBP) {
 int tval = reg << count;
 if (pixval[0] >= cenval[0])
 rcodeval += tval;

 if (pixval[1] >= cenval[1])
 gcodeval += tval;

 if (pixval[2] >= cenval[2])
 bcodeval += tval;
 } else {
 UINT treg[3] = { 0, 0, 0 };
 // posvalue and negvalue contains the thresholded values for r,g & b
 // tcount is used to run over the values...
 for (int tl = 0, tcount = 0; tl < ltplevels; ++tl)
 for (int iter = 0; iter < 3; ++iter, ++tcount) {

 if (pixval[iter] >= cenval[iter] + tollevels[tl]) {
 treg[iter] = tl + 1;
 } else if (pixval[iter] <= cenval[iter]
 - tollevels[tl]) {
 treg[iter] = tl + 1 + ltplevels; // codes will be +ve while vector contains the -ve
 }
 }
 rcodeval += basevalues[count] * treg[0];
 gcodeval += basevalues[count] * treg[1];
 bcodeval += basevalues[count] * treg[2];
 }
 }
 cbook->SetValue(rcodeval, gcodeval, bcodeval);
 }
 }
 }
 */
