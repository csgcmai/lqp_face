/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/

#ifndef ENCODEDMULTIPLECBPATCHFEATURES_H_
#define ENCODEDMULTIPLECBPATCHFEATURES_H_
#include "encodedPatchFeatures.h"
// Quantization aganist multiple codebooks
// currently, Only one codebook for the negative class is generated. Whereas cell based codebooks are generated
// for the +Ve classes...

class EncodedMultipleCBPatchFeatures: public EncodedPatchFeatures {
public:
	EncodedMultipleCBPatchFeatures(const LBPParams &lbpparam, UINT width_,
			UINT height_, UINT nplevels, ProcessImage &pim_,
			const PatchParams& pparam) :
		EncodedPatchFeatures(lbpparam, width_, height_, nplevels, pim_, pparam) {
		// Histogram Parameters
		delete codeinfo;
		codeinfo = NULL; // to be used for dellocation of memory
		char tvar[2525];
		string patchtype = PatchParams::GetPatchType(ptype);
		sprintf(tvar, "MC-CBS-CI-%s-%d-%d", patchtype.c_str(), patchsize,
				ncenters);
		cifilename = tvar;
		cifilename += "-%d.txt";
		sprintf(tvar, "MC-CBS-CB-%s-%d-%d", patchtype.c_str(), patchsize,
				ncenters);
		cbfilename = tvar;
		cbfilename += "-%d.txt";
		switch (cbtype) {
		case CBT_PosCropNegSampled:
		case CBT_PosPyramidNegSampled:
		case CBT_MultiClass:
			ncbs = 2;
			featdim = ncbs * ncenters;
			break;
		case CBT_ConstantPosCellNeg:
		case CBT_PosCellNeg:
		case CBT_PosCellNegSep:
		case CBT_ConstantPosCellNegSingleQuantization: // see the inherited class EncodedMultipleCBPatchFeaturesHistogram
			// featdim = ncenters * 2;/*dim of each cell because it is mapped against its own codebook &-ve cb*/
			// only doing this during returning the features to classifiers, while doing computation
			// each cell is quantized against all the codebooks., so featdim for each cell is
			ncbs = cwidth * cheight + 1;
			featdim = ncenters * ncbs;
		}
		//		ltplevels = (UINT) (nlevels / 2);
		UINT *npmult = pparam.npmult;
		cbooks.resize(ncbs, 0); /*Only LTP based CodeBooks are implemented...*/
		for (vector<ComputeCode*>::iterator iter = cbooks.begin(); iter
				!= cbooks.end(); ++iter)
			*iter = new ComputeCodeLTP(ltplevels, npoints * npmult[ptype],
					ncenters);
		nwinsampled = cbparams->nwinsampled;
		//
		//----------------- Print Parameter Info -----------
		cout << "\n ----------------- Mulitple CB Encoded Parameter Info "
			"----------- " << endl << "Number of CodeBooks = " << ncbs << endl
				<< " Feature Dim = " << featdim << endl << " CodeBookFileName"
				<< cbfilename << endl << " Number of Window Sampled ="
				<< nwinsampled
				<< "\n--------------------------------------------\n";

	}
	// GetDim returns the size of feature in a window without considering offset...
	virtual ~EncodedMultipleCBPatchFeatures() {
		for (vector<ComputeCode*>::iterator iter = cbooks.begin(); iter
				!= cbooks.end(); ++iter)
			delete *iter;
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
	virtual void InitalizeMaps(Image &, PyramidType);
	virtual void InitalizeMaps(Image &image);
	virtual void InitalizeMaps(Pyramid & pyobj_, PyramidType);

	virtual void GetFeatures(UINT, int, int, int, int, vector<REAL>&);
	virtual void GetFeatures(UINT, int, int, int, int, REAL*);
	virtual void GetFoldedFeatures(UINT, int, int, int, int, vector<REAL>&);
	virtual void GetFoldedFeatures(UINT, int, int, int, int, REAL*);

	virtual UINT GetInitIndex() {
		return sspace == DoubleRes ? pyobj.GetInitIndex() : 0;
	}
	virtual void PadFeatureMap(UINT index, UINT padx, UINT pady);

	virtual void DotProduct(UINT index, UINT width, UINT height, vector<REAL>&,
			Store&response);
	virtual void DotProduct(UINT index, UINT width, UINT height, REAL *filter,
			Store&);
	void ComputeLTPMap(Image &image, LBPMap *tpmap, LBPMap *tnmap);
	void ComputeLBPFeatures(Image& imgref, UINT index, UINT winstride,
			UINT tlbpstride);
	void ComputeLBPFeatures(Image& imgref, vector<REAL> &features, UINT cxbmax,
			UINT cybmax, UINT winstride, UINT tlbpstride);

	/*****Test Code ****/
	// All the inputs are in Pixels .....

	void ExtractPosCodeInfo(Image&image, vector<Coord>& ainfo);
	void ExtractNegCodeInfo(Image&image);
	void ExtractCodeInfo(Image &image, LBPMap *map, ComputeCode *cbook);

	void GenerateCodeBookInfo(LBPMap *map, UINT xmin, UINT ymin,
			ComputeCode* cbook); // for a single code book
	void GenerateCodeBookInfo(LBPMap *map, UINT xmin, UINT ymin,
			vector<ComputeCode*> & cbook); // for complete codebook
	//	void GenerateCodeBookInfo(Image & image, ComputeCode * cbook);


	bool ExistPosCodeInfo() {
		char tvar[2000];
		UINT count = 1;
		for (vector<ComputeCode*>::iterator iter = cbooks.begin(); iter
				!= cbooks.end() - 1; ++iter) {
			sprintf(tvar, cifilename.c_str(), count);
			if (!(*iter)->ExistCodeBook(tvar))
				return false;
		}
		return true;
	}
	bool ExistPosCodeBook() {
		char tvar[2000];
		UINT count = 1;
		for (vector<ComputeCode*>::iterator iter = cbooks.begin(); iter
				!= cbooks.end() - 1; ++iter, ++count) {
			sprintf(tvar, cbfilename.c_str(), count);
			if (!(*iter)->ExistCodeBook(tvar))
				return false;
		}
		return true;
	}
	virtual void GenerateCodeBook() {
		char ncbfname[2560];
		sprintf(ncbfname, cbfilename.c_str(), ncbs);

		if (ExistPosCodeBook() && cbooks[ncbs - 1]->ExistCodeBook(ncbfname))
			return;

		if (cbparams->cbtype == CBT_NegativeThenPosCrop) { // first generate the negative code book
			// and use it as initialization for the +ve cell based codebooks

			string cifname = cifilename + "Neg";
			ExtractNegCBInfo(cifname);
			ConstructNegCodeBook(NULL);
			REAL *ncluster = new REAL[ncenters
					* cbooks[ncbs - 1]->GetCodeLength()];
			cbooks[ncbs - 1]->CopyClusters(ncluster);
			//			ExtractPosCodeInfo();
			ExtractPosCBInfo();
			ConstructPosCodeBook(ncluster);
		} else {
			ExtractPosCBInfo();
			ConstructPosCodeBook(NULL);
			ExtractNegCBInfo(cifilename);
			ConstructNegCodeBook(NULL);
		}

	}

	void ExtractPosCBInfo() {
		if (ExistPosCodeInfo())
			return;
		ClassInfo *ocinfo = cbparams->ocinfo;
		if (cbtype == CBT_PosPyramid || cbtype == CBT_PosPyramidNegSampled) {
			sspace = SingleRes;
			SetPyramidInterval(5);
		} else
			sspace = NoPyramid; /*if to use local space pyramid*/
		cout << "\n Using  Scale Space = " << GetPyramidType(sspace) << endl;
		string imname;
		Image image;
		vector<Coord> vcoord;// to contain the flipped features coordinates
		ocinfo->ResetCounter();

		UINT count = 0;
		while (ocinfo->GetNextImageAnnotation(imname, vcoord)) {
			image.read(imname);
			cout << "Image " << ++count << ". Number of Annotations = "
					<< vcoord.size() << endl;
			ExtractPosCodeInfo(image, vcoord);
		}
	}
	void ExtractNegCBInfo(const string &cifn) {
		char cifname[2560];
		sprintf(cifname, cifn.c_str(), ncbs);
		if (cbooks[ncbs - 1]->ExistCodeInfo(cifname))
			return;
		sspace = SingleRes;
		SetPyramidInterval(2);
		vector<string> imnames;
		if (cbtype == CBT_Complete)/*Use complete validation set*/
			ParseListFile(cbparams->vimgfname, imnames);
		else
			ParseListFile(cbparams->nimgfname, imnames);
		Image image;
		for (UINT i = 0; i < imnames.size(); ++i) {
			cout << " Processing Image # = " << (i + 1) << " " << imnames[i]
					<< endl;
			image.read(imnames[i]);
			ExtractNegCodeInfo(image);
		}
	}

	void ConstructPosCodeBook(REAL *initclusters) {
		char ci[2000], cb[2000];
		UINT count = 1, nclusrounds = cbparams->nclusrounds;
		ClusteringDistanceMetric dmetric = cbparams->dmetric;

		if (ExistPosCodeBook())
			return;

		for (vector<ComputeCode*>::iterator iter = cbooks.begin(); iter
				!= (cbooks.end() - 1); ++iter, ++count) {
			sprintf(ci, cifilename.c_str(), count);
			sprintf(cb, cbfilename.c_str(), count);
			long initime = time(0);
			(*iter)->GenerateCodeBook(ci, cb, dmetric, nclusrounds,
					initclusters);
			cout << " \n Time Taken For Code Book Number " << count
					<< "Generation =" << time(0) - initime << endl << flush;
		}
		assert(count== ncbs);
	}
	void ConstructNegCodeBook(REAL *initclusters) {
		char ci[2000], cb[2000];
		UINT count = ncbs, nclusrounds = cbparams->nclusrounds;
		ClusteringDistanceMetric dmetric = cbparams->dmetric;
		sprintf(ci, cifilename.c_str(), count);
		sprintf(cb, cbfilename.c_str(), count);
		long initime = time(0);
		cbooks[count - 1]->GenerateCodeBook(ci, cb, dmetric, nclusrounds,
				initclusters);
		cout << " \n Time Taken For Code Book Number " << count
				<< "Generation =" << time(0) - initime << endl << flush;
	}

private:
	UINT GetFeatDim(UINT width_, UINT height_) const {
		return (width_ / cellsize) * (height_ / cellsize) * featdim;
	}
protected:
	UINT featdim;
	vector<ComputeCode*> cbooks;
	UINT ncbs;
	string imgfname;
};

#endif
