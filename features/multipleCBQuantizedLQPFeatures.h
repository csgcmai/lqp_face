/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/

#ifndef MULTIPLECBQUANTIZEDLQPFEATURES_H_
#define MULTIPLECBQUANTIZEDLQPFEATURES_H_
#include "encodedPatchFeatures.h"
#include "omp.h"
/// LQP Quantization, Features are quantized against each cell codebook

class MultipleCBQuantizeLQPFeatures: public EncodedPatchFeatures {
public:
	MultipleCBQuantizeLQPFeatures(const LBPParams &lbpparam, UINT width_,
			UINT height_, UINT nplevels, ProcessImage &pim_,
			const PatchParams& pparam) :
		EncodedPatchFeatures(lbpparam, width_, height_, nplevels, pim_, pparam) {
		// Histogram Parameters
		delete codeinfo;
		codeinfo = NULL; // to be used for dellocation of memory
		char tvar[2525];
		string patchtype = PatchParams::GetPatchType(ptype);
		ostringstream oss;
		oss << "MCBQ-CI-Pos" << patchtype << "-" << patchsize << "-"
				<< ncenters;
		cifilename = oss.str();
		oss.str("");
		oss << "MCBQ-CB-Pos" << patchtype << "-" << patchsize << "-"
				<< ncenters;
		cbfilename = oss.str();

		// for negative side lqp
		oss.str("");
		oss << "MCBQ-CI-Neg" << patchtype << "-" << patchsize << "-"
				<< ncenters;
		ncifilename = oss.str();

		oss.str("");
		oss << "MCBQ-CB-Neg" << patchtype << "-" << patchsize << "-"
				<< ncenters;
		ncbfilename = oss.str();

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
		case CBT_ConstantPosCellNegSingleQuantization: // see the inherited class MultipleCBQuantizeLQPFeaturesHistogram
			// featdim = ncenters * 2;/*dim of each cell because it is mapped against its own codebook &-ve cb*/
			// only doing this during returning the features to classifiers, while doing computation
			// each cell is quantized against all the codebooks., so featdim for each cell is
			ncbs = cwidth * cheight + 1;
			featdim = ncenters * ncbs;
		case CBT_CellQuantizeAll: /// cell quantization doing against all the cells codebooks
			ncbs = cwidth * cheight;
			featdim = ncenters * ncbs; // each cell is quantized against all the codebooks
		}
		UINT *npmult = pparam.npmult;
		pcbooks.resize(ncbs, 0); /*Only LTP based CodeBooks are implemented...*/

		if (petype == PET_SplitLTP) {
			ncbooks.resize(ncbs, 0); /*Only LTP based CodeBooks are implemented...*/
			featdim = featdim * 2;
			for (vector<ComputeCode*>::iterator iter = pcbooks.begin(), niter =
					ncbooks.begin(); iter != pcbooks.end(); ++iter, ++niter) {
				*iter = new ComputeCodeLBP(ltplevels, npoints * npmult[ptype],
						ncenters, softq, meanp, pparam.pcount);
				*niter = new ComputeCodeLBP(ltplevels, npoints * npmult[ptype],
						ncenters, softq, meanp, pparam.pcount);
			}

		} else {
			for (vector<ComputeCode*>::iterator iter = pcbooks.begin(); iter
					!= pcbooks.end(); ++iter)
				*iter = new ComputeCodeLTP(ltplevels, npoints * npmult[ptype],
						ncenters, softq, meanp, pparam.pcount);
		}

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
	virtual ~MultipleCBQuantizeLQPFeatures() {
		for (vector<ComputeCode*>::iterator iter = pcbooks.begin(), niter =
				ncbooks.begin(); iter != pcbooks.end(); ++iter, ++niter) {
			delete *iter;
			delete *niter;
		}
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
	void ComputeSplitLBPFeatures(Image& imgref, vector<REAL> &features,
			UINT cxbmax, UINT cybmax, UINT winstride, UINT tlbpstride);

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
	void GenerateCodeBookInfo(LBPMap *pmap, LBPMap *nmap, UINT xmin, UINT ymin);

	bool ExistPosCodeInfo() {
		UINT count = 1;
		for (vector<ComputeCode*>::iterator iter = pcbooks.begin(), niter =
				ncbooks.begin(); iter != pcbooks.end(); ++iter, ++niter) {
			ostringstream oss;
			oss << cifilename << "-" << count << ".txt";
			if (!(*iter)->ExistCodeBook(oss.str()))
				return false;

			oss.str("");
			oss << ncifilename << "-" << count << ".txt";
			if (!(*niter)->ExistCodeBook(oss.str()))
				return false;
		}

		return true;
	}
	bool ExistPosCodeBook() {
		UINT count = 1;
		for (vector<ComputeCode*>::iterator iter = pcbooks.begin(), niter =
				ncbooks.begin(); iter != pcbooks.end(); ++iter, ++count, ++niter) {
			ostringstream oss;
			oss << cbfilename << "-" << count << ".txt";
			if (!(*iter)->ExistCodeBook(oss.str()))
				return false;

			oss.str("");
			oss << ncbfilename << "-" << count << ".txt";
			if (!(*niter)->ExistCodeBook(oss.str()))
				return false;
		}
		return true;
	}
	virtual void GenerateCodeBook() {

		if (ExistPosCodeBook())
			return;

		if (cbparams->cbtype == CBT_NegativeThenPosCrop) { // first generate the negative code book
			// and use it as initialization for the +ve cell based codebooks
			string cifname = cifilename + "Neg";
			ExtractNegCBInfo(cifname);
			ConstructNegCodeBook(NULL);
			REAL *ncluster = new REAL[ncenters
					* pcbooks[ncbs - 1]->GetCodeLength()];
			pcbooks[ncbs - 1]->CopyClusters(ncluster);
			//			ExtractPosCodeInfo();
			ExtractPosCBInfo();
			ConstructPosCodeBook(ncluster);
		} else if (cbparams->cbtype == CBT_CellQuantizeAll) {

			if (!ExistPosCodeInfo()) {
				vector<string> imnames;
				ParseListFile(cbparams->vimgfname, imnames);
				Image image;
				for (UINT i = 0; i < imnames.size(); ++i) {
					cout << "\n Processing Image # = " << (i + 1) << " "
							<< imnames[i];
					image.read(imnames[i]);
					ExtractImageCodeInfo(image);
				}
			}
			ConstructPosCodeBook(NULL);
		} else {
			ExtractPosCBInfo();
			ConstructPosCodeBook(NULL);
			ExtractNegCBInfo(cifilename);
			ConstructNegCodeBook(NULL);
		}

	}
	void ExtractImageCodeInfo(Image &image) {
		LBPMap pmap[3], nmap[3];
		ExtractCodeBookInfo(image, pmap, nmap, false, true);
		UINT offset = GetMaxFeatureOffset(), minoffset = GetFeaturesOffset(),
				doffset = offset - minoffset;
		GenerateCodeBookInfo(pmap, nmap, doffset, doffset);

		image.flop(); // flipped image
		cout << "   ------ Processing Flipped Version  ";
		ExtractCodeBookInfo(image, pmap, nmap, false, true);
		GenerateCodeBookInfo(pmap, nmap, doffset, doffset);
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
		ostringstream oss;
		oss << cifn << "-" << ncbs << ".txt";
		if (pcbooks[ncbs - 1]->ExistCodeInfo(oss.str()))
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
		UINT nclusrounds = cbparams->nclusrounds;
		ClusteringDistanceMetric dmetric = cbparams->dmetric;

		if (ExistPosCodeBook())
			return;
#pragma omp parallel for
		for (UINT count = 0; count < pcbooks.size(); ++count) {
			ostringstream ciss, cbss;
			ciss << cifilename << "-" << count + 1 << ".txt";
			cbss << cbfilename << "-" << count + 1 << ".txt";
			long initime = time(0);
			pcbooks[count]->GenerateCodeBook(ciss.str(), cbss.str(), dmetric,
					nclusrounds, initclusters);
			cout << " \n Time Taken For Positive Code Book Number " << count
					<< "Generation =" << time(0) - initime << endl << flush;
		}

#pragma omp parallel for
		for (UINT count = 0; count < pcbooks.size(); ++count) {
			ostringstream ciss, cbss;
			ciss << ncifilename << "-" << count + 1 << ".txt";
			cbss << ncbfilename << "-" << count + 1 << ".txt";
			long initime = time(0);
			ncbooks[count]->GenerateCodeBook(ciss.str(), cbss.str(), dmetric,
					nclusrounds, initclusters);
			cout << " \n Time Taken For Negative Code Book Number " << count
					<< "Generation =" << time(0) - initime << endl << flush;
		}

	}
	void ConstructNegCodeBook(REAL *initclusters) {
		char ci[2000], cb[2000];
		UINT count = ncbs, nclusrounds = cbparams->nclusrounds;
		ClusteringDistanceMetric dmetric = cbparams->dmetric;
		ostringstream ciss, cbss;
		ciss << cifilename << "-" << count << ".txt";
		cbss << cbfilename << "-" << count << ".txt";
		long initime = time(0);
		pcbooks[count - 1]->GenerateCodeBook(ciss.str(), cbss.str(), dmetric,
				nclusrounds, initclusters);
		cout << " \n Time Taken For Code Book Number " << count
				<< "Generation =" << time(0) - initime << endl << flush;
	}

private:
	UINT GetFeatDim(UINT width_, UINT height_) const {
		return (width_ / cellsize) * (height_ / cellsize) * featdim;
	}
protected:
	UINT featdim; /// single cell feature dimension
	vector<ComputeCode*> pcbooks, ncbooks;
	UINT ncbs;
	string imgfname;
};

#endif
