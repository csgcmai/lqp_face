/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/

#ifndef ENCODEDPATCHFEATURES_H_
#define ENCODEDPATCHFEATURES_H_
#include "ltpFeaturesRGB.h"
#include "encodeFeatures.h"
#include "../util/hogInfo.h"
#include "../util/classInfo.h"
#include "definitions.h"
#define EPS 0.0001
//#include "util.hpp"
#define PATCH_CODE_VERSION2 // with interpolation of histograms;
// Current Constraints include that
// if ptype == PT_Circular and nlevels > 2 then build code book for each level separately...
// just to be to handle the multitude of levels...
// if ptype != PT_Circular Number of patch stride is limited by number of levels of codes...
//#define LIVE_QUANTIZATION
class EncodedPatchFeatures: public Features {
public:
	EncodedPatchFeatures(const LBPParams &lbpparam, UINT width_, UINT height_,
			UINT nplevels, ProcessImage &pim_, const PatchParams& pparam) :
			Features(width_, height_, nplevels, false), pim(pim_) {
		// Histogram Parameters
		useblocknorm = false;

		issigned = false;
		cbparams = pparam.cbparams;
		nwinsampled = cbparams->nwinsampled;
		cellsize = lbpparam.cellsize;
		add = lbpparam.addfeat;
		nseparate = lbpparam.normfeatsep;
		gridtype = lbpparam.gridtype;
		lbpstride = lbpparam.stride;
		norm = lbpparam.featnorm;
		// PATCH Parameters
		ptype = pparam.type;
		patchsize = pparam.size;
		cpatchsize = pparam.cpatchsize;
		hmethod = lbpparam.hmethod;
		UINT ncpoints = 0, nspoints = 0; // used for the circular
		if (cpatchsize > 0) {
			cpatchsize = cpatchsize % 2 == 0 ? cpatchsize - 1 : cpatchsize;
			ncpoints = cpatchsize * cpatchsize - 1;
		}
		if ((ptype == PT_HVCStrip || ptype == PT_DACStrip) && cpatchsize <= 0) {
			cout << " Circular Patch Size = " << cpatchsize << endl;
			PrintErrorExit(" Circular patch size <=0 ");
		}
		petype = pparam.petype;
		cbtype = cbparams->cbtype;
		softq = cbparams->softq;
		meanp = cbparams->meanp;
		if (patchsize % 2 == 0 && ptype >= 0 && ptype <= 9) {
			patchsize -= 1;
			cout
					<< "\n !=!=!=!=!=!=!=!= § Even Patch Size § !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!= \n"
					<< " Converting it to Odd Size =" << patchsize
					<< "\n !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!= \n";
		}

		patchstride = pparam.stride;
		ncenters = cbparams->codebooksize;
		celldim = ncenters;
		cellhistdim = ncenters;
		if (pparam.pstep > 2) {
			pixelstep = 2;
			cout
					<< "\n !=!=!=!=!=!=!=!= § Attention § !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!= \n"
					<< " Using Pixel Step of 2 Instead of " << pparam.pstep
					<< "\n !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!= \n";
		} else
			pixelstep = pparam.pstep; // Maximum of 2 pixel is allowed in the current implementation
		nlevels = pparam.ncodelevels;
		nspoints = GetNPoints(lbpparam);
		UINT *npmult = pparam.npmult; //number of patch multiplier for multiple orientation see initialization in PatchParam constructor.
		npoints = nspoints * npmult[ptype] + ncpoints;
		icneighbours = 8;
		cell = lbpparam.cellsize;
		foffset = lbpparam.cellsize;
		UINT tvar = GetMaxOffset();
		lfeatoffset = (UINT) (ceil((REAL) tvar / cell) * cell);
		maxradius = tvar; // maximum of one cell is allowed..
		if (ptype == PT_HAARPyramid && maxradius != 2)
			maxradius = 2;

		pparam.PrintParamInfo(); // Print Patch Parameters
		cout << endl << " -------------- Histogram Info ----------------------"
				<< endl << " Width = " << width << " Height = " << height

				<< " CellStride= " << lbpstride << endl << " CellSize= " << cell
				<< endl << " OffsetSize= " << foffset << endl
				<< " Feature Normalization Type = " << norm << endl
				<< " Number of Points per Patch = " << npoints / npmult[ptype]
				<< " Number Code Books = " << 1 << endl << " Max Radius = "
				<< maxradius
				<< (add == true ? "  & Adding RGB Channels at Cell levels" :
						"  Concatenating RGB Channels at Cell levels ")
				<< " Local Feature Offset = " << lfeatoffset << endl << flush;
		cout << endl << " -----------------------------------------------------"
				<< endl;
		// Window Parameters

		cwidth = width / cell;
		cheight = height / cell;
		cout << endl << " Original Feature Dim = " << GetDim(width_, height_);

		skip = lbpstride / cell; // how many cells to skip
		// initialization of the variables
		ComputeFileNames();

		if (ptype == PT_HVCStrip || ptype == PT_DACStrip) {
			int nr = 3;
			spoints.resize(nr);
			spoints[0].resize(ncpoints);
			for (int tvar = 1; tvar < nr; ++tvar) {
				spoints[tvar].resize(nspoints);
			}

		} else
			spoints.resize(npmult[ptype], vector<Point<REAL> >(nspoints));
		ComputeSamplingPoints();
		ncodeinfo = NULL;
		cout << " Number of Points = " << npoints << " Sampling Points = "
				<< nspoints;
		if (ncpoints != 0)
			cout << " Circular Points =" << ncpoints << endl << flush;
		InitializeCodeInfo(pparam, nspoints, npmult);

		// whether their are multiple codebooks or not all will have be in same base...
		UINT base = codeinfo->GetBase();
		UINT tpoints =
				petype == PET_CSLBP || petype == PET_CSLBP ? npoints / 2 :
						npoints;
		basevalues.resize(tpoints, 0);
		cout << " \n Base Levels are = ";
		for (REAL i = 0; i < tpoints; ++i) {
			basevalues[i] = pow((REAL) base, i);
			cout << basevalues[i] << " , ";
		}
		cout << endl << "Finished Initialization " << endl << flush;

		//TODO : current implementation requires same cell & stride size, change this
		assert(lbpstride == cellsize);
		char ncenstr[256];
		sprintf(ncenstr, "%d", ncenters);
		cbfilename.append(ncenstr);
		midpoint = nspoints / 2;
	}
	void InitializeCodeInfo(const PatchParams &pparam, UINT nspoints,
			UINT *npmult) {
		int discinfo = cbtype == CBT_NegativePosDiscriminant ? 1 : 0;

		switch (petype) {
		case PET_LBP:
			assert(npoints <= 28);
			codeinfo = new ComputeCodeLBP(1, npoints, ncenters, softq, meanp,
					pparam.pcount, cbparams->knn, discinfo);
			break;
		case PET_CSLBP:
			assert(npoints <= 48);
			codeinfo = new ComputeCodeLBP(1, npoints / 2, ncenters, softq,
					meanp, pparam.pcount, cbparams->knn, discinfo);
			ltplevels = 1;
			tollevels = pparam.tolval;
			break;
		case PET_LTP:
			assert(npoints <= 18);
			ltplevels = pparam.tolval.size();
			codeinfo = new ComputeCodeLTP(ltplevels, npoints, ncenters, softq,
					meanp, pparam.pcount, cbparams->knn, discinfo);
			tollevels = pparam.tolval;
			cout << " Tolerance Levels are " << tollevels << endl;
			break;
		case PET_CSLTP:
			ltplevels = pparam.tolval.size();
			codeinfo = new ComputeCodeLTP(ltplevels, npoints / 2, ncenters,
					softq, meanp, pparam.pcount, cbparams->knn, discinfo);
			tollevels = pparam.tolval;
			cout << " Tolerance Levels are " << tollevels << endl;
			break;
		case PET_SplitLTP:
			ltplevels = pparam.tolval.size();
			assert(ltplevels == 1);
			tollevels = pparam.tolval;
			codeinfo = new ComputeCodeLBP(1, npoints, ncenters, softq, meanp,
					pparam.pcount, cbparams->knn, discinfo);
			ncodeinfo = new ComputeCodeLBP(1, npoints, ncenters, softq, meanp,
					pparam.pcount, cbparams->knn, discinfo);
			celldim = 2 * ncenters;
			cellhistdim = ncenters;
			cout << " Tolerance Levels are " << tollevels << endl;
			break;
		case PET_CSLTP_LTP:
			//			spoints.resize(npmult[ptype] * 2);// total number of patches...
			//			for (int tvar = 0; tvar < npmult[ptype] * 2; ++tvar) {
			//				if (tvar < npmult[ptype])
			//					spoints[tvar].resize(nspoints);
			//				else
			//					spoints[tvar].resize(nspoints / 2);
			//			}
			npoints = (nspoints + nspoints / 2) * npmult[ptype];
			assert(npoints<= 18);
			ltplevels = pparam.tolval.size();
			codeinfo = new ComputeCodeLTP(ltplevels, npoints, ncenters, softq,
					meanp, pparam.pcount, cbparams->knn, discinfo);
			tollevels = pparam.tolval;
			cout << " Tolerance Levels are " << tollevels << endl;
			break;
		default:
			cerr << "\n Should Use LQPPatch Features";
			cerr << " Nlevels =  " << nlevels << " NPoints= " << npoints
					<< endl;
			exit(EXIT_FAILURE);
		}
	}
	UINT GetMaxOffset() {

		if (ptype == PT_CircAngRad) {
			return (patchsize == 3 ? 2 : (patchsize / 2)) * pixelstep; // 2 to consider outer layer of pixels
		}
		return ((ptype >= 0 && ptype <= 9 || ptype == PT_CircularShallow) ? (patchsize
				/ 2) * pixelstep :
				1); // comparing neighbouring pixels
	}
	void ComputeFileNames() {
		char *ptypename[] = { "HStrip", "VStrip", "Circular", "Diagonal",
				"ADiagonal", "HVStrip", "DAStrip", "HVCStrip", "DACStrip",
				"HVDAStrip", "HVNStrip", "HaarPyramid", "CircularShallow",
				"CircularAngularRadial" };
		cbfilename = ptypename[ptype]
				+ PatchParams::GetPatchEncodingType(petype);
		//		cbfilename.append("_Nlevels_CodeBook .txt");
		cifilename = ptypename[ptype]
				+ PatchParams::GetPatchEncodingType(petype);

		char spsize[100];
		sprintf(spsize, "-patchsize%d-ncenters%d", patchsize, ncenters);
		cbfilename += spsize;
		cbfilename.append("-CodeBook-.txt");
		cifilename += spsize;
		cifilename.append("-CodeInfo-.txt");

		ncbfilename = ptypename[ptype]
				+ PatchParams::GetPatchEncodingType(petype) + "-neg-";
		//		cbfilename.append("_Nlevels_CodeBook .txt");
		ncifilename = ptypename[ptype]
				+ PatchParams::GetPatchEncodingType(petype) + "-neg-";

		sprintf(spsize, "-patchsize%d-ncenters%d", patchsize, ncenters);
		ncbfilename += spsize;
		ncbfilename.append("-CodeBook-.txt");
		ncifilename += spsize;
		ncifilename.append("-CodeInfo-.txt");
	}
	// GetDim returns the size of feature in a window without considering offset...
	UINT GetNPoints(const LBPParams & lbpparams) {
		if (ptype == PT_Circular || ptype == PT_HAARPyramid) {
			//			if (samptype == PST_Fixed)
			return patchsize * patchsize - 1;
			//			else {
			//				// Not doing dense points, limited by number of codes we can have in the memory...
			//				int tpsize = patchsize / 2;
			//				UINT sum = 0;
			//				for (int i = 0; i < tpsize; ++i)
			//					sum += (i + 1) * pixelstep * icneighbours;
			//				return sum;
			//			}
		} else if (ptype == PT_CircularShallow || ptype == PT_CircAngRad)
			return (patchsize / 2) * lbpparams.npoints;
		else
			return (patchsize - 1);
		//ptype	== PT_Circular ? GetNPoints(patchsize, pixelstep) /*patchsize * patchsize - 1 */
		//	: patchsize - 1
	}
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
	virtual ~EncodedPatchFeatures() {
		if (codeinfo)
			delete codeinfo;
		if (ncodeinfo)
			delete ncodeinfo;
	}

	virtual void InitalizeMaps(Image &, PyramidType);
	virtual void InitalizeMaps(Pyramid & pyobj_, PyramidType);
	virtual void InitalizeMaps(Image &image);

	virtual void GetFeatures(UINT, int, int, int, int, vector<REAL>&);
	virtual void GetFeatures(UINT, int, int, int, int, REAL*);
	virtual void GetFoldedFeatures(UINT, int, int, int, int, vector<REAL>&);
	virtual void GetFoldedFeatures(UINT, int, int, int, int, REAL*);
	virtual void
	GetFlippedFeatures(UINT, int, int, int, int, vector<REAL>&);
	virtual void UnFoldWeights(REAL*, UINT, UINT, vector<REAL>&);

	virtual UINT GetXSize(UINT index) {
		return lbpfeatures.GetXBlock(index) * cell;
	}
	virtual UINT GetYSize(UINT index) {
		return lbpfeatures.GetYBlock(index) * cell;
	}
	virtual REAL GetScale(UINT index) {
		return lbpfeatures.GetScale(index);
	}
	virtual UINT GetInitIndex() {
		return sspace == DoubleRes ? pyobj.GetInitIndex() : 0;
	}
	virtual void PadFeatureMap(UINT index, UINT padx, UINT pady);

	virtual void DotProduct(UINT index, UINT width, UINT height, vector<REAL>&,
			Store&response);
	virtual void DotProduct(UINT index, UINT width, UINT height, REAL *filter,
			Store&);
	virtual void DotProductWithQuadratic(UINT index, UINT width, UINT height,
			REAL, vector<REAL>&, vector<REAL>&, Store&) {
	}
	virtual void DotProductWithProjection(UINT index, UINT xsize, UINT ysize,
			const VectorXf &sigma, const MatrixXf &projMat, vector<REAL>&,
			Store&inhog) const {

	}
	void GetBaseCode(Image &image, vector<UINT>&); // Map the pixels to (2*number_of_levels+1)^number of neighbours
	void ComputeBinaryMap(Image &image, vector<vector<UINT> >& bmap);
	void ComputeLTPMap(Image &image, LBPMap *tpmap, LBPMap *tnmap);
	//	virtual void ComputeEncodedMap(Image &image, LBPMap *map);
	//	void ComputeLTPMap(Image &image, vector<vector<UINT> >&);
	//	virtual void ComputeLBPFeatures(Image& imgref, UINT index, UINT winstride,
	//			UINT tlbpstride);
	//	virtual void ComputeSplitLBPFeatures(Image &image, UINT index, UINT sbin,
	//			UINT lbpstride);
	virtual void ComputeLBPFeatures(Image &image, UINT index, UINT sbin,
			UINT lbpstride) {
		LBPMap plbpmap[3], nlbpmap[3];
		ExtractCodeBookInfo(image, plbpmap, nlbpmap, false); //extract the codebook information

		if (hmethod == HM_Discrete) {
			if (petype == PET_SplitLTP) {
				ComputeSplitLQPHistogramDiscrete(image, index, sbin, plbpmap,
						nlbpmap);
			} else
				ComputeLQPHistogramDiscrete(image, index, sbin, plbpmap);
		} else {
			vector<REAL> &tvec = lbpfeatures.GetFeatureRef(index);
			REAL *feat = &tvec[0];
			if (petype == PET_SplitLTP) {
				ComputeSplitLQPHistogram(image, sbin, plbpmap, nlbpmap, feat);
			} else
				ComputeLQPHistogram(image, sbin, plbpmap, feat);
		}
	}
	virtual void ComputeLQPHistogram(Image &image, UINT sbin, LBPMap plbpmap[],
			REAL *feat);
	virtual void ComputeSplitLQPHistogram(Image &image, UINT sbin,
			LBPMap plbpmap[], LBPMap nlbpmap[], REAL *feat);
	virtual void ComputeLQPHistogramDiscrete(Image& imgref, UINT index,
			UINT winstride, LBPMap map[3]);
	virtual void ComputeSplitLQPHistogramDiscrete(Image& imgref, UINT index,
			UINT winstride, LBPMap pmap[3], LBPMap nmap[3]);

	/*****Test Code ****/
	// All the inputs are in Pixels .....
	virtual void GetHOGFeatures(UINT index, int x, int y, int width_,
			int height_, REAL*feat) { // There is offset of one cell so, input x=0 means x=8 by taking
		//	into account offset and the removed border
		//		hoginfo.GetHOGCells(index, x / cell, y / cell, width_ / cell, height_
		//				/ cell, feat);
	}

	virtual void GetFoldedHOGFeatures(UINT index, int x, int y, int twidth,
			int theight, REAL*feat) { // There is offset of one cell so, input x=0 means x=8 by taking
		//	into account offset and the removed border
		//		hoginfo.GetFoldedHOGCells(index, x / cell, y / cell, twidth / cell,
		//				theight / cell, feat);
	}

	virtual void GetFlippedHOGFeatures(UINT index, int x, int y, int twidth,
			int theight, vector<REAL>&feat) { // There is offset of one cell so, input x=0 means x=8 by taking
		//	into account offset and the removed border
		//		hoginfo.GetFlippedFeatures(index, x / cell, y / cell, twidth / cell,
		//				theight / cell, &feat[0]);
	}
	virtual void GetFlippedFeatures(int twidth, int theight, REAL *ifeat,
			REAL *ofeat);

	virtual void HOGDotProduct(UINT index, UINT fwidth, UINT fheight,
			vector<REAL>&filter, Store&response) {
		//		hoginfo.DotProduct(index, fwidth / cell, fheight / cell, &filter[0],
		//				response);
	}
	virtual void LBPDotProduct(UINT index, UINT fwidth, UINT fheight,
			vector<REAL>&filter, Store&response) {
		lbpfeatures.SkippedDotProduct(index, fwidth / cell, fheight / cell,
				skip, &filter[0], response);
	}
	virtual void GetLBPFeatures(UINT index, int x, int y, int w, int h,
			REAL* feat) {
		lbpfeatures.GetSkippedHOGCells(index, x / cell, y / cell, w / cell,
				h / cell, skip, feat);
	}
	virtual void ComputeTextureFeatures(Image &imgref, UINT xcellsize,
			UINT ycellsize, vector<REAL>&features);
	virtual void ComputeTextureFeatures(Image &imgref, UINT xcellsize,
			UINT ycellsize, int sx, int sy, int ex, int ey,
			vector<REAL>&features);
	void ComputeSplitTextureFeatures(Image &imgref, UINT xcellsize,
			UINT ycellsize, int sx, int sy, int ex, int ey,
			vector<REAL>&features);
	void GenerateCodeBookInfo(LBPMap *pmap, LBPMap *nmap, UINT xmin, UINT ymin);
	void ComputeCodes() {
		UINT nclasses = 2;
		cout << "\n Computing Codes from Discriminative Multiclass clustering"
				<< endl;
		cout << "\n Number of Classes =" << nclasses << endl;
		UINT npexamples = 0, nnexamples = 0;
		// write file with the code info...
		ostringstream oss;
		oss << "codeinfo" << "-class-" << 2 << "-pos";

		if (!cbparams->cboffset.empty()) {
			cout << " Loading LUT from the positive LQP from file "
					<< oss.str() + cbparams->cboffset << endl;
			oss << "-" << cbparams->cboffset;
			if (!codeinfo->LoadLUT(oss.str())) {
				cerr << " Could't load the LUT for Positive LQP file" << endl;
				exit(EXIT_FAILURE);
			}
			if (ncodeinfo) {
				oss.clear();
				oss.str("");
				oss << "codeinfo" << "-class-" << 2 << "-neg-";
				cout << " Loading LUT for the negative LQP from file "
						<< oss.str() + cbparams->cboffset << endl;
				if (!ncodeinfo->LoadLUT(oss.str() + cbparams->cboffset)) {
					cerr << " Could't load the LUT from Positive LQP file"
							<< endl;
					exit(EXIT_FAILURE);
				}
			}
			cout << endl << " File(s) loaded Successfully" << endl << flush;
			return;
		}
		ofstream pofile(oss.str().c_str(), ios::out | ios::binary), nofile;
		codeinfo->WriteHeader(pofile);
		CodeVecType *pcodevec, *ncodevec; //	apply the delete operation as well...
		codeinfo->GenerateCodeVectors(pcodevec); // to initializate the codes ids
		delete[] pcodevec;
		pofile.write((char*) &nclasses, sizeof(UINT));

		if (ncodeinfo) {
			oss.clear();
			oss.str("");
			oss << "codeinfo" << "-class-" << 2 << "-neg";
			nofile.open(oss.str().c_str(), ios::out | ios::binary);
			ncodeinfo->WriteHeader(nofile);
			nofile.write((char*) &nclasses, sizeof(UINT));
			ncodeinfo->GenerateCodeVectors(ncodevec); // to initialize the codes ids
			delete[] ncodevec;
		}
		ClearData();
		cbparams->cbtype = CBT_PosCrop;
		ExtractPosCBInfo();
		SaveCodeInfo(pofile, nofile, 1, npexamples, nnexamples);
		cout << "\n Class " << 1 << " Has " << npexamples
				<< " Non-zero codes for Positive and " << nnexamples
				<< " for negative Ternary Images " << flush;

		ClearData();
		cbparams->cbtype = CBT_MultiClass;
		ExtractNegCBInfo();
		SaveCodeInfo(pofile, nofile, 0, npexamples, nnexamples);
		cout << "\n Class " << 0 << " Has " << npexamples
				<< " Non-zero codes for Positive and " << nnexamples
				<< " for negative Images " << flush;

		pofile.write((char*) &npexamples, sizeof(UINT));
		pofile.close();
		if (nofile.is_open()) {
			nofile.write((char*) &nnexamples, sizeof(UINT));
			nofile.close();
		}
	}
	virtual void ExtractCodeInfo(Image&image, vector<Coord>& ainfo);
	virtual void ExtractCodeInfo(Image&image);
	void ExtractNegCodeInfo(Image &image);
	virtual void ExtractPosCBInfo() {
		if (codeinfo->ExistCodeInfo(cifilename)) {
			if (ncodeinfo && ncodeinfo->ExistCodeInfo(ncifilename))
				return;
		}
		if (cbparams->cbtype == CBT_PosPyramid)
			sspace = SingleRes;
		else if (cbparams->cbtype == CBT_Complete
				|| cbparams->cbtype == CBT_PosCrop
				|| cbparams->cbtype == CBT_NegativeThenPosCrop)
			sspace = NoPyramid;
		ClassInfo *ocinfo = cbparams->ocinfo;
		if (ocinfo && cbparams->cbtype != CBT_Complete) { // use annotation information for the cropping of windows
			string imname;
			Image image;
			vector<Coord> vcoord; // to contain the flipped features coordinates
			ocinfo->ResetCounter();
			SetPyramidInterval(5);
			UINT count = 0;
#ifdef TMPDEBUG
			while (ocinfo->GetNextImageAnnotation(imname, vcoord) && count < 20) {
#else
			while (ocinfo->GetNextImageAnnotation(imname, vcoord)) {
#endif
				image.read(imname);
				cout << "Image " << ++count << ". Number of Annotations = "
						<< vcoord.size() << endl;
				ExtractCodeInfo(image, vcoord);
			}
		} else { // simply read the  images and generate codebook
			LBPMap map[3], nmap[3];
			cout
					<< " \n Reading list of images for codebook generation from validation file < "
					<< cbparams->vimgfname << endl;
			vector<string> imnames;
			ParseListFile(cbparams->vimgfname, imnames);
			Image image;
			for (UINT i = 0; i < imnames.size(); ++i) {
				cout << " Processing Cropped Image # = " << (i + 1) << " "
						<< imnames[i] << endl;
				image.read(imnames[i]);
				if (cbparams->sx < 0) {
					ExtractCodeBookInfo(image, map, nmap, true);
					image.flop();
					ExtractCodeBookInfo(image, map, nmap, true);
				} else {
					ExtractCodeBookInfo(image, map, nmap, cbparams->sx,
							cbparams->sy, cbparams->ex, cbparams->ey, true);
					image.flop();
					int sx = image.columns() - cbparams->ex, ex =
							image.columns() - cbparams->sx + 1;
					ExtractCodeBookInfo(image, map, nmap, sx, cbparams->sy, ex,
							cbparams->ey, true);
				}
			}
		}
	}
	void ExtractNegCBInfo() {
		if (codeinfo->ExistCodeInfo(cifilename)
				|| (ncodeinfo && ncodeinfo->ExistCodeInfo(ncifilename)))
			return;
		sspace = SingleRes;
		SetPyramidInterval(2);
		vector<string> imnames;
		if (cbparams->cbtype == CBT_Complete)
			ParseListFile(cbparams->vimgfname, imnames);
		else
			ParseListFile(cbparams->nimgfname, imnames);
		Image image;
		if (nwinsampled != 0) {
			cout << "\n Sampling Windows from Negative Images " << endl;
		}
#ifdef TMPDEBUG
		for (UINT i = 0; i < 20; ++i) {
#else
		for (UINT i = 0; i < imnames.size(); ++i) {
#endif
			cout << " Processing Image # = " << (i + 1) << " " << imnames[i]
					<< endl;
			image.read(imnames[i]);
			ExtractNegCodeInfo(image);
		}
	}
	void ExtractCBInfo(const string &pcb, const vector<string> &imnames) {
		cifilename = pcb + "-pos-codeinfo";
		ncifilename = pcb + "-neg-codeinfo";
		if (codeinfo->ExistCodeInfo(cifilename))
			return;
		//				vector<string> imnames;
		//				ParseListFile(pcb, imnames);
		Image image;
		cout << " Processing Images: ";
#ifdef TMPDEBUG
		for (UINT i = 0; i < 20; ++i) {
#else
		for (UINT i = 0; i < imnames.size(); ++i) {
#endif
			cout << ".+";
			image.read(imnames[i]);

			LBPMap map[3], nmap[3];
			ExtractCodeBookInfo(image, map, nmap, true);
		}
	}

	virtual void ExtractCodeBookInfo(Image &image, LBPMap*, LBPMap*, int, int,
			int, int, const bool&, bool = false); // last to compute only the codeinfo
	virtual void ExtractCodeBookInfo(Image &image, LBPMap*, LBPMap*,
			const bool&, bool = false); // last to compute only the codeinfo

	void ExtractCircAngRadCodeBookInfo(Image &image, LBPMap *map, LBPMap *nmap,
			const bool &gencodebook, bool computeonlycodes);
	void ExtractHaarCodeBookInfo(Image &image, LBPMap *map, LBPMap *nmap,
			const bool &gencodebook, bool = false);
	//*****/
	virtual void GenerateCodeBook() {
		if (ncodeinfo && codeinfo->ExistCodeBook(cbfilename)
				&& ncodeinfo->ExistCodeBook(ncbfilename))
			return;
		else if (codeinfo->ExistCodeBook(cbfilename) && !ncodeinfo)
			return;
		if ((cbparams->cbtype == CBT_Negative)) {
			ExtractNegCBInfo(); // for cbt_complete see function...
			GenerateCenters(NULL, NULL);
		} else if ((cbparams->cbtype == CBT_Complete
				|| cbparams->cbtype == CBT_PosCrop
				|| cbparams->cbtype == CBT_PosPyramid)) {
			ExtractPosCBInfo();
			GenerateCenters(NULL, NULL);
		} else if ((cbparams->cbtype == CBT_NegativeThenPosCrop)) {
			string cifn = cifilename;
			cifilename += "_Neg";
			ExtractNegCBInfo();
			GenerateCenters(NULL, NULL);
			REAL
					*pinitclusters = new REAL[ncenters
							* codeinfo->GetCodeLength()], *ninitclusters =
							new REAL[ncenters * ncodeinfo->GetCodeLength()];
			codeinfo->CopyClusters(pinitclusters);
			codeinfo->ClearData();
			if (ncodeinfo) {
				ncodeinfo->CopyClusters(ninitclusters);
				ncodeinfo->ClearData();
			}
			cifilename = cifn + "_Pos";
			ncifilename += "_Pos";
			ExtractPosCBInfo();
			GenerateCenters(pinitclusters, ninitclusters);
		} else if ((cbparams->cbtype == CBT_NegativePosDiscriminant)) {
			string cifn = cifilename;
			cifilename += "_Pos";
			ExtractPosCBInfo();
			vector<REAL> neg_pcodecount, neg_ncodecount, pos_pcodecount,
					pos_ncodecount; // negative and positive class code counts
			codeinfo->GetCodeCount(pos_pcodecount); /// count for codes for negative class
			codeinfo->ClearData();
			if (ncodeinfo) {
				ncodeinfo->GetCodeCount(pos_ncodecount);
				ncodeinfo->ClearData();
			}
			cifilename = cifn + "_Neg";
			ncifilename += "_Neg";
			ExtractNegCBInfo();
			codeinfo->GetCodeCount(neg_pcodecount); /// count for codes for negative class
			REAL * pratio, *nratio; // ratio information for both positive and negative
			FindRatio(pos_pcodecount, neg_pcodecount, pratio);
			codeinfo->AddCodeCount(pos_pcodecount);
			if (ncodeinfo) {
				ncodeinfo->GetCodeCount(neg_ncodecount);
				FindRatio(pos_ncodecount, neg_ncodecount, nratio);
				ncodeinfo->AddCodeCount(pos_ncodecount);
			}
			GenerateCenters(NULL, NULL, pratio, nratio);
		} else if (cbparams->cbtype == CBT_MultiClass) {
			ComputeCodes();
		}
	}
	virtual void GenerateTextureCodeBook(const string& pcbfname,
			const vector<string>&imgnames, const vector<string>& labels) {
		cbfilename = pcbfname + "-pos-codebook";
		ncbfilename = pcbfname + "-neg-codebook";
		if (ncodeinfo && codeinfo->ExistCodeBook(cbfilename)
				&& ncodeinfo->ExistCodeBook(ncbfilename))
			return;
		else if (codeinfo->ExistCodeBook(cbfilename) && !ncodeinfo)
			return;
		if ((cbparams->cbtype == CBT_Complete
				|| cbparams->cbtype == CBT_Negative)) {
			ClearData();
			ExtractCBInfo(pcbfname, imgnames); // for cbt_complete see function...
			GenerateCenters(NULL, NULL);
		} else if (cbparams->cbtype == CBT_MultiClass) {
			// simply generate different classes information
			// find set of labels..

			set<string> setlabel(labels.begin(), labels.end());
			UINT nclasses = setlabel.size();
			cout << "\n Computing Codes from Multiclass " << endl;
			cout << "\n Number of Classes =" << nclasses << endl;
			UINT clscount = 0, npexamples = 0, nnexamples = 0;
			vector<string> trnimages;
			// write file with the code info...
			ostringstream oss;
			oss << pcbfname << "codeinfo" << "-class-" << setlabel.size()
					<< "-pos";
			if (!cbparams->cboffset.empty()) {
				cout << " Loading LUT from the positive LQP from file "
						<< oss.str() + cbparams->cboffset << endl;
				oss << "-";
				if (!codeinfo->LoadLUT(oss.str() + cbparams->cboffset)) {
					cerr << " Could't load the LUT for Positive LQP file"
							<< endl;
					exit(EXIT_FAILURE);
				}
				if (ncodeinfo) {
					oss.clear();
					oss.str("");
					oss << pcbfname << "codeinfo" << "-class-"
							<< setlabel.size() << "-neg-";
					cout << " Loading LUT for the negative LQP from file "
							<< oss.str() + cbparams->cboffset << endl;
					if (!ncodeinfo->LoadLUT(oss.str() + cbparams->cboffset)) {
						cerr << " Could't load the LUT from Positive LQP file"
								<< endl;
						exit(EXIT_FAILURE);
					}
				}
				cout << endl << " File(s) loaded Successfully" << endl << flush;
				return;
			}

			ofstream pofile(oss.str().c_str(), ios::out | ios::binary), nofile;
			codeinfo->WriteHeader(pofile);
			CodeVecType *pcodevec, *ncodevec; //	apply the delete operation as well...
			//	vector<vector<CodeVecType> > codevec;
			//	codevec.clear();
			codeinfo->GenerateCodeVectors(pcodevec); // to initializate the codes ids
			delete[] pcodevec;
			pofile.write((char*) &nclasses, sizeof(UINT));

			if (ncodeinfo) {
				oss.clear();
				oss.str("");
				oss << pcbfname << "codeinfo" << "-class-" << setlabel.size()
						<< "-neg";
				nofile.open(oss.str().c_str(), ios::out | ios::binary);
				ncodeinfo->WriteHeader(nofile);
				nofile.write((char*) &nclasses, sizeof(UINT));
				ncodeinfo->GenerateCodeVectors(ncodevec); // to initialize the codes ids
				delete[] ncodevec;
			}
			for (set<string>::iterator iter = setlabel.begin();
					iter != setlabel.end(); ++iter, ++clscount) {
				// find all the images having same label
				FindImagesWithSameLabels(imgnames, labels, *iter, trnimages);
				cout << "\n Class " << clscount << " Has " << trnimages.size()
						<< "Images " << flush;
				cout << "\n Generating Codes";
				ClearData();
				ExtractCBInfo(pcbfname, trnimages); // for cbt_complete see function...

				SaveCodeInfo(pofile, nofile, clscount, npexamples, nnexamples);
				trnimages.clear();
			}
			cout << "\n Number of Positive Examples = " << npexamples;
			pofile.write((char*) &npexamples, sizeof(UINT));
			pofile.close();
			if (nofile.is_open()) {
				cout << "\n Number of Negative Examples = " << nnexamples;
				nofile.write((char*) &nnexamples, sizeof(UINT));
				nofile.close();

			}

		} else {
			cout
					<< " Error! Texture codebooks can be only generated from CBT_Complete type "
					<< endl;
		}
	}
	void GenerateCenters(REAL *initclusters, REAL *ninitclusters, REAL *pratio =
			NULL, REAL *nratio = NULL) {
		UINT nclusrounds = cbparams->nclusrounds;
		ClusteringDistanceMetric dmetric = cbparams->dmetric;
		long initime = time(0);
		if (ncenters == 59 && cbparams->inituniform && !initclusters) //  Initialize Cluster centers from the uniform codes
				{
			vector<REAL> indeces(58, 0);
			FindUniformIndeces(indeces);
			codeinfo->GenerateCodeBook(cifilename, cbfilename, dmetric,
					nclusrounds, &indeces[0]); //
			if (ncodeinfo)
				ncodeinfo->GenerateCodeBook(ncifilename, ncbfilename, dmetric,
						nclusrounds, &indeces[0]); //
		} else {
			codeinfo->GenerateCodeBook(cifilename, cbfilename, dmetric,
					nclusrounds, initclusters, pratio);
			cout << " \n Time Taken For "
					<< (ncodeinfo ? "Positive" : "Complete")
					<< " Code Book Generation =" << time(0) - initime << endl
					<< flush;
			if (ncodeinfo) {
				ncodeinfo->GenerateCodeBook(ncifilename, ncbfilename, dmetric,
						nclusrounds, ninitclusters, nratio);
				cout
						<< " \n Time Taken For Negative Code Book Number Generation ="
						<< time(0) - initime << endl << flush;
			}
		}
	}
	void SaveCodeInfo(ofstream &pfile, ofstream &nfile, UINT clslabel,
			UINT &npexamples, UINT &nnexamples) {
		/*codeinfo->WriteLQPCodes(fname + "-pos");
		 if (ncodeinfo)
		 ncodeinfo->WriteLQPCodes(fname + "-neg");*/
		pfile.write((char*) &clslabel, sizeof(UINT));
		npexamples += codeinfo->WriteLQPCodes(pfile);
		if (ncodeinfo) {
			nfile.write((char*) &clslabel, sizeof(UINT));
			nnexamples += ncodeinfo->WriteLQPCodes(nfile);
		}
	}
	void ClearData() {
		codeinfo->ClearData();
		if (ncodeinfo)
			ncodeinfo->ClearData();
	}
	template<class T>
	void FindUniformIndeces(vector<T> &indeces) {
		UINT tnpoints = 8;
		UINT niter = unsigned(1 << tnpoints);
		int index = 0;
		for (UINT i = 0; i < niter; i++) {
			if (Transitions(i, tnpoints) <= 2) //uniform
					{
				indeces[index] = i;
				++index;
			}
		}
	}

	void ComputeCircAngRadPoints() {
		REAL x, y;
		spoints.clear();
		spoints.resize(2); ///> first row contains samples from the reference circle, 2nd row contains the points for difference
		int hpsize = patchsize / 2;
		if (hpsize == 1) {
			for (int k = 0; k < 2; ++k) {
				REAL step = (M_PI * 2.0) / icneighbours;
				for (int i = 0; i < icneighbours; ++i) {
					x = (k + 1) * cos((double) (i * step)); // where (k+1) is radius
					y = -(k + 1) * sin((double) (i * step));
					x = ABS(x) <= 1e-6 ? 0 : x;
					y = ABS(y) <= 1e-6 ? 0 : y;
					spoints[k].push_back(Point<REAL>(x, y));
				}
			}
		} else {
			int nspoints = icneighbours * hpsize;
			for (int k = hpsize - 1; k <= hpsize; ++k) {
				REAL step = (M_PI * 2.0) / nspoints;
				for (int i = 0; i < nspoints; ++i) {
					x = k * cos((double) (i * step)); // where (k+1) is radius
					y = -k * sin((double) (i * step));
					x = ABS(x) <= 1e-6 ? 0 : x;
					y = ABS(y) <= 1e-6 ? 0 : y;
					spoints[hpsize - k].push_back(Point<REAL>(x, y));
				}
			}
		}
		cout << "\n Sample Points = ";
		for (int j = 0; j < spoints.size(); ++j) {
			cout << "Row " << j + 1 << " has -->";
			for (int i = 0; i < spoints[j].size(); ++i) {
				cout << i + 1 << ". (" << spoints[j][i].GetX() << " , "
						<< spoints[j][i].GetY() << "); ";
			}
			cout << endl;
		}
	}
	void ComputeSamplingPoints() {
		REAL x, y;
		int cpoints, count = 0;
		if (ptype == PT_HAARPyramid)
			return;
		if (ptype == PT_CircAngRad)
			return ComputeCircAngRadPoints();
		if (ptype == PT_CircularShallow) {

			int tradius = patchsize / 2;
			REAL step = (M_PI * 2.0) / npoints;

			for (int i = 0; i < npoints; ++i, ++count) {
				x = tradius * cos((double) (i * step)); // where (k+1) is radius
				y = -tradius * sin((double) (i * step));
				x = ABS(x) <= 1e-6 ? 0 : x;
				y = ABS(y) <= 1e-6 ? 0 : y;
				spoints[0][count].SetCoord(x, y);
			}

		} else if ((ptype == PT_Circular && gridtype == Circular)
				|| (ptype == PT_HVCStrip || ptype == PT_DACStrip) /*|| ptype
				 == PT_CircularShallow*/) {
			int hpsize = (
					(ptype == PT_HVCStrip || ptype == PT_DACStrip) ? cpatchsize
							/ 2 :
							patchsize / 2);

			int tradius = pixelstep;
			int k = 0;
			/*if (ptype == PT_CircularShallow) {
			 k = hpsize - 1;
			 tradius = patchsize / 2;
			 }*/

			for (; k < hpsize; ++k, tradius += pixelstep) {
				cpoints = ((k + 1) * icneighbours); // points on the ring to be sample will be multiple of ncneighbours
				REAL step = (M_PI * 2.0) / cpoints;

				for (int i = 0; i < cpoints; ++i, ++count) {
					x = tradius * cos((double) (i * step)); // where (k+1) is radius
					y = -tradius * sin((double) (i * step));
					x = ABS(x) <= 1e-6 ? 0 : x;
					y = ABS(y) <= 1e-6 ? 0 : y;
					spoints[0][count].SetCoord(x, y);
				}
			}
		} else if (ptype == PT_Circular && gridtype == Rectangular) {
			int hpsize = patchsize / 2;
			for (int x = -(int) (pixelstep * hpsize);
					x <= (int) (pixelstep * hpsize); x += pixelstep)
				for (int y = -(int) (pixelstep * hpsize);
						y <= (int) (pixelstep * hpsize); y += pixelstep)
					if (x != 0 || y != 0)
						spoints[0][count++].SetCoord(x, y);

		}
		if (ptype != PT_Circular && ptype != PT_CircularShallow) {
			count = 0;
			int hpsize = patchsize / 2;
			for (int x = -(int) (hpsize * pixelstep);
					x <= (int) (pixelstep * hpsize); x += pixelstep)
				if (x != 0) {
					switch (ptype) {
					case PT_HStrip:
						spoints[0][count].SetCoord(x, 0);
						break;
					case PT_VStrip:
						spoints[0][count].SetCoord(0, x);
						break;
					case PT_Diag:
						spoints[0][count].SetCoord(x, x);
						break;
					case PT_ADiag:
						spoints[0][count].SetCoord(-x, x);
						break;
					case PT_HVStrip:
					case PT_HVNStrip:
						spoints[0][count].SetCoord(x, 0);
						spoints[1][count].SetCoord(0, x);
						break;
					case PT_DAStrip:
						spoints[0][count].SetCoord(x, x);
						spoints[1][count].SetCoord(-x, x);
						break;
					case PT_HVDAStrip:
						spoints[0][count].SetCoord(x, 0);
						spoints[1][count].SetCoord(0, x);
						spoints[2][count].SetCoord(x, x);
						spoints[3][count].SetCoord(-x, x);
						break;
					case PT_HVCStrip:
						spoints[1][count].SetCoord(x, 0);
						spoints[2][count].SetCoord(0, x);
						break;
					case PT_DACStrip:
						spoints[1][count].SetCoord(x, x);
						spoints[2][count].SetCoord(-x, x);
						break;
					}
					++count;
				}
		}
		cout << "\n Number of Sampling Points = " << count
				<< "; They are as follows: \n ";
		int pcounter = 0;
		for (int j = 0; j < spoints.size(); ++j) {
			for (int i = 0; i < spoints[j].size(); ++i, ++pcounter) {
				cout << pcounter + 1 << " = (" << spoints[j][i].GetX() << ", "
						<< spoints[j][i].GetY() << "); ";
			}
			cout << endl;
		}
	}
	int Transitions(UINT c, int nbits) {
		int base = 1;
		int current = c & base, current2, changes = 0;
		for (int i = 1; i < nbits; i++) {
			base <<= 1;
			current2 = (c & base) >> i;
			if (current ^ current2)
				changes++;
			current = current2;
		}
		return changes; //(changes <= 2)? 1 : 0;
	}
	virtual UINT GetMaxFeatureOffset() {
		return lfeatoffset;
	}
	ComputeCode & GetCodeObjRef() {
		return *codeinfo;
	}
	void FindImagesWithSameLabels(const vector<string>&imgnames,
			const vector<string>&labels, string targetlabel,
			vector<string>& trnimages) {
		for (UINT i = 0; i < labels.size(); ++i)
			if (labels[i] == targetlabel)
				trnimages.push_back(imgnames[i]);
	}
	void FindRatio(vector<REAL>& pcount, vector<REAL>& ncount, REAL *&ratio) {
		long tpcount = 0, tncount = 0;
		assert(pcount.size()==ncount.size());
		// find the total sum
		for (vector<REAL>::iterator iter = pcount.begin(); iter != pcount.end();
				++iter)
			tpcount += *iter;

		for (vector<REAL>::iterator iter = ncount.begin(); iter != ncount.end();
				++iter)
			tncount += *iter;

		// normalize the pos and neg counts;
		for (vector<REAL>::iterator iter = pcount.begin(); iter != pcount.end();
				++iter)
			*iter /= tpcount;

		for (vector<REAL>::iterator iter = ncount.begin(); iter != ncount.end();
				++iter)
			*iter /= tncount;
		ratio = new REAL[pcount.size()];
		for (UINT i = 0; i < pcount.size(); ++i)
			ratio[i] = ncount[i] == 0 ? pcount[i] : pcount[i] / ncount[i];
	}
private:
	UINT GetFeatDim(UINT width_, UINT height_) const {
		return (width_ / cellsize) * (height_ / cellsize) * celldim;
	}
protected:
	static double ipercent; // store the percentage of images scanned from database, used for online-pruning of codes...
	bool useblocknorm;
	QuantizationType softq; // softq
	REAL meanp; // emean percentage
	/*cpatchsize is only used for the circular + otherstrip features...*/
	UINT npoints, patchsize, cpatchsize, nlevels, patchstride, icneighbours,
			pixelstep;
	PatchEncodingType petype;
	HistogramMethod hmethod;
	PatchType ptype;
	CBParams *cbparams;
	UINT midpoint;
	bool issigned;
	vector<vector<Point<REAL> > > spoints; // sampling points on the circle...
	vector<REAL> basevalues;
	NORM norm;
	UINT lbpstride;
	ProcessImage &pim;
	bool add, nseparate; // if add add all the color-channels else concatenates
	HOGInfo lbpfeatures;
	GridType gridtype;
	CBType cbtype;
	UINT cell, cwidth, cheight, skip;
	ComputeCode * codeinfo; // See constraint detail above
	UINT ncenters;
	UINT maxradius, ltplevels;
	string cbfilename, cifilename; //code information filename codebook file name....
	vector<REAL> tollevels; // tolerance levels used to compute the code
	//// local feat offset specially used because patch strip can have size of upto 3 cells, while foffset is fixed because all
	// other features are computed with an offset of 8 pixels, so in this case a boundary will be padded
	// whose difference is equal to featoffset-offset;
	UINT lfeatoffset; // local feat offset
	UINT celldim;
	// For split case
	string ncbfilename, ncifilename; //code information filename codebook file name....
	ComputeCode * ncodeinfo; // See constraint detail above
	UINT cellhistdim;
	UINT nwinsampled;
	void Compute(const Image&, UINT, REAL*);
	void GetXYBlocks(const Image&image, UINT&xbmax, UINT&ybmax) { // look in the code of Compute as the last two boundary blocks are not used
		xbmax = (UINT) (round((REAL) image.columns() / (REAL) cell) - 2);
		ybmax = (UINT(round((REAL) image.rows() / (REAL) cell) - 2));
	}

	void Get2XYBlocks(const Image&image, UINT&xbmax, UINT&ybmax) { // look in the code of Compute as the last two boundary blocks are not used

		xbmax = (UINT) (round((REAL) (image.columns() * 2) / (REAL) cell) - 2);
		ybmax = (UINT) (round((REAL) (image.rows() * 2) / (REAL) cell) - 2);
	}

	void
	GetLBPFeatures(UINT, UINT, UINT, UINT, UINT, REAL*); // scale index with starting location and size also size
};
/*
 class LTPSplitEncodedPatchFeatures: public EncodedPatchFeatures {
 public:

 protected:

 };*/
class SoftLQPFeatures: public EncodedPatchFeatures {
public:
	SoftLQPFeatures(const LBPParams &lbpparam, UINT width_, UINT height_,
			UINT nplevels, ProcessImage &pim_, const PatchParams& pparam) :
			EncodedPatchFeatures(lbpparam, width_, height_, nplevels, pim_,
					pparam) {
		pmethod = pparam.cbparams->pmethod;
	}
	virtual void ComputeLBPFeatures(Image &image, UINT index, UINT sbin,
			UINT lbpstride);
	void ComputeDiscreteLQPHistogram(Image &image, UINT index, UINT cellsize,
			UINT lbpstride);
	template<class T>
	void Pool(T *tar, T val) {
		if (pmethod == PM_Sum)
			*tar += val;
		else if (pmethod == PM_Max)
			*tar = MAX(*tar,val);
	}
private:
	PoolingMethod pmethod;
};
/**
 *  Gabor-LQP Features try to record the co-occurence statistics over different orientation and scales
 *  i.e LQP Features are computed over 3-D image planes, these planes include spatial, orientation and scale dimension..
 *  Input list files contain the list of neigbour (in scale and orientation) Gabor-images.
 */

class GaborLQPFeatures: public EncodedPatchFeatures {
public:
	GaborLQPFeatures(const LBPParams &lbpparam, UINT width_, UINT height_,
			UINT nplevels, ProcessImage &pim_, const PatchParams& pparam,
			UINT nglevels_) :
			EncodedPatchFeatures(lbpparam, width_, height_, nplevels, pim_,
					pparam) {
		nglevels = nglevels_;
		UINT nspoints = GetNPoints(lbpparam);
		UINT *npmult = pparam.npmult; //number of patch multiplier for multiple orientation see initialization in PatchParam constructor.
		npoints = (nglevels * 2 + 1) * nspoints * npmult[ptype] + nglevels * 2;
		cout << endl
				<< " -------------- GaborLBP Configuration Info ----------------------"
				<< " Number of Points per Patch = "
				<< (nglevels * 2 * npoints / npmult[ptype] + (nglevels - 1))
				<< endl
				<< " -----------------------------------------------------"
				<< endl;
		// Window Parameters

		if (ncodeinfo)
			delete ncodeinfo;
		ncodeinfo = NULL;
		delete codeinfo;
		cout << " Number of Points = " << npoints << " Sampling Points = "
				<< nspoints << endl << flush;
		InitializeCodeInfo(pparam, nspoints, npmult);
	}
	virtual void InitalizeMaps(vector<Image> &image) {
		sspace = NoPyramid;
		UINT xbmax, ybmax;
		lbpfeatures.Initialize(1, 0, celldim);
		GetXYBlocks(image[nglevels], xbmax, ybmax);
		lbpfeatures.SetFeature(xbmax, ybmax, 1.0);
		ComputeLBPFeatures(image, 0, cell, lbpstride);
	}
	virtual void ComputeLBPFeatures(vector<Image> &image, UINT index, UINT sbin,
			UINT lbpstride) {
		LBPMap tlbpmap[3], nlbpmap[3];
		ExtractCodeBookInfo(image, tlbpmap, nlbpmap, false); //extract the codebook information
		vector<REAL> &tvec = lbpfeatures.GetFeatureRef(index);
		REAL *feat = &tvec[0];
		if (petype == PET_SplitLTP) {
			ComputeSplitLQPHistogram(image[nglevels], sbin, tlbpmap, nlbpmap,
					feat);
		} else
			ComputeLQPHistogram(image[nglevels], sbin, tlbpmap, feat);
	}

	void ExtractCodeBookInfo(vector<Image> &image, LBPMap *map, LBPMap *nmap,
			const bool &gencodebook, bool computeonlycodes = false) {
		ExtractCodeBookInfo(image, map, nmap, 0, 0, image[0].columns(),
				image[0].rows(), gencodebook, computeonlycodes);
	}
	void ExtractCodeBookInfo(vector<Image> &image, LBPMap *map, LBPMap *nmap,
			int sx, int sy, int ex, int ey, const bool &gencodebook,
			bool computeonlycodes = false);

	virtual void ExtractPosCBInfo() {
		if (codeinfo->ExistCodeInfo(cifilename)) {
			if (ncodeinfo && ncodeinfo->ExistCodeInfo(ncifilename))
				return;
		}
		sspace = NoPyramid;
		/// Validation files contain Gabor images of neigbour(scale + orientation) in linear order with a step size of nglevels,
		LBPMap map[3], nmap[3];
		cout
				<< " \n Reading list of images for codebook generation from validation file";
		vector<string> imnames;
		ParseListFile(cbparams->vimgfname, imnames);
		vector<Image> image(nglevels * 2 + 1), fimage(nglevels * 2 + 1);
		for (UINT i = nglevels; i < imnames.size(); i += (nglevels * 2 + 1)) {
			cout << " Processing Cropped Image # = " << (i / nglevels + 1)
					<< " " << imnames[i] << endl;
			for (UINT j = i - nglevels, count = 0; j <= i + nglevels;
					++j, ++count) {
				cout << "Adding Image" << imnames[j] << endl << flush;
				image[count].read(imnames[j]);
				fimage[count] = image[count];
				fimage[count].flop();
			}
			if (cbparams->sx < 0) {
				ExtractCodeBookInfo(image, map, nmap, true);
				ExtractCodeBookInfo(fimage, map, nmap, true);
			} else {
				ExtractCodeBookInfo(image, map, nmap, cbparams->sx,
						cbparams->sy, cbparams->ex, cbparams->ey, true);
				int sx = image[0].columns() - cbparams->ex, ex =
						image[0].columns() - cbparams->sx + 1;
				ExtractCodeBookInfo(fimage, map, nmap, sx, cbparams->sy, ex,
						cbparams->ey, true);
			}
		}
	}
private:
	UINT nglevels; //+,- orientation and scale levels..
};

/**
 *  SeparateGaborLQP Features try to record the co-occurence statistics over different orientation and scales
 *  However in contrast to GaborLQPFeatures, SeparateGaborLQPFeatures just combines the different orientation and
 *  scales code at codebook level. In simple words it combines multiple lqp vector computed separately over neighbour orientations
 *  and scale at codebook level
 */

class SeparateGaborLQPFeatures: public EncodedPatchFeatures {
public:
	SeparateGaborLQPFeatures(const LBPParams &lbpparam, UINT width_,
			UINT height_, UINT nplevels, ProcessImage &pim_,
			const PatchParams& pparam, UINT nglevels_) :
			EncodedPatchFeatures(lbpparam, width_, height_, nplevels, pim_,
					pparam) {
		nglevels = nglevels_;
		UINT nspoints = GetNPoints(lbpparam);
		UINT *npmult = pparam.npmult; //number of patch multiplier for multiple orientation see initialization in PatchParam constructor.
		npoints = (nglevels * 2 + 1) * nspoints * npmult[ptype];
		cout << endl
				<< " -------------- Separate GaborLBP Configuration Info ----------------------\n"
				<< " Number of Points per Patch = " << npoints << endl
				<< " -----------------------------------------------------"
				<< endl;
		// Window Parameters

		if (ncodeinfo)
			delete ncodeinfo;
		ncodeinfo = NULL;
		delete codeinfo;
		cout << " Number of Points = " << npoints << " Sampling Points = "
				<< nspoints << endl << flush;
		InitializeCodeInfo(pparam, nspoints, npmult);
	}
	virtual void InitalizeMaps(vector<Image> &image) {
		sspace = NoPyramid;
		UINT xbmax, ybmax;
		lbpfeatures.Initialize(1, 0, celldim);
		GetXYBlocks(image[nglevels], xbmax, ybmax);
		lbpfeatures.SetFeature(xbmax, ybmax, 1.0);
		ComputeLBPFeatures(image, 0, cell, lbpstride);
	}
	virtual void ComputeLBPFeatures(vector<Image> &image, UINT index, UINT sbin,
			UINT lbpstride) {
		LBPMap tlbpmap[3], nlbpmap[3];
		ExtractCodeBookInfo(image, tlbpmap, nlbpmap, false); //extract the codebook information
		vector<REAL> &tvec = lbpfeatures.GetFeatureRef(index);
		REAL *feat = &tvec[0];
		if (petype == PET_SplitLTP) {
			ComputeSplitLQPHistogram(image[nglevels], sbin, tlbpmap, nlbpmap,
					feat);
		} else
			ComputeLQPHistogram(image[nglevels], sbin, tlbpmap, feat);
	}

	void ExtractCodeBookInfo(vector<Image> &image, LBPMap *map, LBPMap *nmap,
			const bool &gencodebook, bool computeonlycodes = false) {
		ExtractCodeBookInfo(image, map, nmap, 0, 0, image[0].columns(),
				image[0].rows(), gencodebook, computeonlycodes);
	}
	void ExtractCodeBookInfo(vector<Image> &image, LBPMap *map, LBPMap *nmap,
			int sx, int sy, int ex, int ey, const bool &gencodebook,
			bool computeonlycodes = false);

	virtual void ExtractPosCBInfo() {
		if (codeinfo->ExistCodeInfo(cifilename)) {
			if (ncodeinfo && ncodeinfo->ExistCodeInfo(ncifilename))
				return;
		}
		sspace = NoPyramid;
		/// Validation files contain Gabor images of neigbour(scale + orientation) in linear order with a step size of nglevels,
		LBPMap map[3], nmap[3];
		cout
				<< " \n Reading list of images for codebook generation from validation file";
		vector<string> imnames;
		ParseListFile(cbparams->vimgfname, imnames);
		vector<Image> image(nglevels * 2 + 1), fimage(nglevels * 2 + 1);
		for (UINT i = nglevels; i < imnames.size(); i += (nglevels * 2 + 1)) {
			cout << " Processing Cropped Image # = " << (i / nglevels + 1)
					<< " " << imnames[i] << endl;
			for (UINT j = i - nglevels, count = 0; j <= i + nglevels;
					++j, ++count) {
				cout << "Adding Image" << imnames[j] << endl << flush;
				image[count].read(imnames[j]);
				fimage[count] = image[count];
				fimage[count].flop();
			}
			if (cbparams->sx < 0) {
				ExtractCodeBookInfo(image, map, nmap, true);
				ExtractCodeBookInfo(fimage, map, nmap, true);
			} else {
				ExtractCodeBookInfo(image, map, nmap, cbparams->sx,
						cbparams->sy, cbparams->ex, cbparams->ey, true);
				int sx = image[0].columns() - cbparams->ex, ex =
						image[0].columns() - cbparams->sx + 1;
				ExtractCodeBookInfo(fimage, map, nmap, sx, cbparams->sy, ex,
						cbparams->ey, true);
			}
		}
	}
private:
	UINT nglevels; //+,- orientation and scale levels..
};

#endif
