/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/

#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_
#include "../util/util.hpp"
#include "../util/classInfo.h"
enum PyramidType {
	NoPyramid, SingleRes, DoubleRes
};
enum NORM {
	L1 = 0,
	L1sqrt,
	L2,
	LOG2/*features are translated by 1*/,
	LOG2Norm,
	L2Hystersis,
	NT_Average,
// for texture features...
};
enum GridType {
	Circular, Rectangular
};
enum HistogramMethod {
	HM_Discrete, ///> discrete binning of the codes
	HM_Bilinear
///> histogram binning of codes.
};
enum PatchType { // See the variable npmult in the structure
	PT_HStrip,
	PT_VStrip,
	PT_Circular,
	PT_Diag,
	PT_ADiag,
	PT_HVStrip,
	PT_DAStrip,
	PT_HVCStrip,
	PT_DACStrip,
	PT_HVDAStrip,
	PT_HVNStrip,
	PT_HAARPyramid, // multi-dimensional Haar organization...
	PT_CircularShallow, /// only the outer ring of pixels is sampled
	PT_CircAngRad,
///> combine the circular, angular and radial arrangements, circular same as LBP, Angular = diff between two pixels, radial= difference over ring of pixels
// compare neighbour pixels instead of comparison with the centeral pixels
//	PT_HVCNStrip
// compare center & neighbour pixels
};

enum PatchEncodingType {
	PET_LBP, PET_LTP, PET_CSLBP, PET_CSLTP, PET_SplitLTP, PET_CSLTP_LTP
};
enum CBType {
	/*	"0. Complete - Validation set to generate CB\n"
	 "1. PosCrop - CB from positive crop images\n"
	 "2. PosCropNegSampled - quantize against both poscrop and -ve CBs\n"
	 "3. PosPyramid\n - quantize against scale-space positive images CBs\n"
	 "4. PosPyramidNegSampled - quantize against both scale-space positive images and -ve CBs \n"
	 "5. ConstantPosCellNeg - quantize each cell against all +ve and -ve CBss\n"
	 "6. ConstantPosCellNegSingleQuantization - Window level histogram of cell CBss (obtained via option ConstantPosCellNeg.)\n"
	 "7. PosCellNeg - Each Cell is mapped only against its own CB and -ve class CB; negatives cells are accumulated at window level \n"
	 "8. PosCellNegSep - Same as PosCellNeg except negatives are not accumuated at windows levels\n"
	 "9. Negative - quantize only against negative class codebook\n"
	 "10.NegativeThenPosCrop - first initialize from negative class then purify with pos class\n"*/
	CBT_Complete, // complete training data is used with scale-space implementation
	CBT_PosCrop, // Use only Positive Cropped windows
	CBT_PosCropNegSampled, // Use two different codebooks  composed of  Positive Cropped windows & Sampled Negative Windows
	CBT_PosPyramid, // Use Positive Cropped windows With Scale Space
	CBT_PosPyramidNegSampled, // Use two different codebooks composed of  Positive Scale-space windows & Sampled Negative Windows
	CBT_ConstantPosCellNeg, // CodeBook for each cell based feature & one for the -ve examples.
	// Each cell is mapped against all the code books
	CBT_ConstantPosCellNegSingleQuantization, // all the codebooks from the above are accumulated at window level..
	//codebook of cell number 1 are accumulated into a window level histogram...
	CBT_PosCellNeg, // Each Cell is mapped against only its respective codeBook with the inclusion of negative codebook
	// Each cell is mapped against just its own cell & negative code book centers which are accumulated into final window
	// level histogram...
	CBT_PosCellNegSep, // Negative are not accumulated at window level.
	CBT_Negative,
	CBT_NegativeThenPosCrop,
	CBT_MultiClass, // To generate codes from each class examples(1 vs all) and then to use discriminative clustering...
	//	CBT_PosCropThenNegative, // Use complete -ve training set in scale space .... // Same as CBT_PosCellNeg except the -ve codebooks are not accumulated into window level histogram
	CBT_NegativePosDiscriminant,
	///> Find the ratio of positive to negative counts for each code and append this information during clustering
	// then either build the histogram using the ratio (bin) or simple count
	CBT_CellQuantizeAll,
///Quantize all cell against the codebooks of each cell, for each cell a separate codebook is learned.

};
enum WindowFeatureTypes // backward compatibility earlier versions///
{
	WFT_LBP = 1,
	WFT_HOG,
	WFT_HOG_LBP,
	WFT_HOG_COLOR,
	WFT_LBP_LTP,
	WFT_HOG_LBP_LTP,
	WFT_PROJ_HOG,
	WFT_LTPX,
	WFT_LBP256,
	WFT_CSLBP,
	WFT_CSLTP,
	WFT_LTP,
	WFT_PROJ_HOG_LBP_LTP_COLOR, // BUG IN COMPUTATION OF DET, NO PROBLEM FOR PR,
	WFT_PROJ_HOG_LBP_COLOR,
	WFT_LBP_COLOR, // these are the new code
	WFT_LTP_COLOR,
	WFT_LBP_LTP_COLOR,
	WFT_HOG_LBP_COLOR,
	WFT_HOG_LTP_COLOR,
	WFT_HOG_LBP_LTP_COLOR,
	WFT_CSLBP_COLOR,
	WFT_CSLTP_COLOR,
};
enum HistFeatures {
	HF_BIN, HF_BIN_TRAIN, HF_LHRATIO,
// liklihood ratio
};
enum FeatureTypes {
	FT_HOG,
	FT_LBP,
	FT_LTP,
	FT_LBPLTP,
	FT_HOGLBP,
	FT_HOGLTP,
	FT_HOGLBPLTP,
	FT_MultiLevelLBP,
	FT_MultiLevelLTP,
	FT_MultiLevelLBPLTP,
	FT_EncodedLTP,
	FT_EncodedLBP,
	FT_EncodedLBPPatch,
	FT_EncodedPatchPlusOther,
	FT_EncodedGenericPatch,
	FT_MultipleEncodedGenericPatchFeatures,
	FT_EncodedLTPPatch,
	FT_EncodedRFPatchFeatures,
	FT_LBPBlock,
	FT_CSLBP,
	FT_CSLTP,
	FT_LQPSoft,
	FT_GaborLQP,
	FT_SeparateGaborLQP,
// Change the function below
};
enum PoolingMethod {
	PM_Sum, PM_Max
};
enum EncodedHistogramNormType {
	EHNT_Local, // to do normalization of codebook quantization at cell level
	EHNT_Global, // to do normalization of codebook quantization at Window level
	EHNT_Both
// both at window and cell level (slowest == default)
};

static WindowFeatureTypes GetWindowFeatureTypes(const string &fname) {

	const int nftype = 22;
	string featname[] = { "LBP", "HOG", "HOG_LBP", "HOG_COLOR", "LBP_LTP",
			"HOG_LBP_LTP", "PROJ_HOG", "LTPX", "LBP256", "CSLBP", "CSLTP",
			"LTP", "PROJ_HOG_LBP_LTP_COLOR", "PROJ_HOG_LBP_COLOR", "LBP_COLOR",
			"LTP_COLOR", "LBP_LTP_COLOR", "HOG_LBP_COLOR", "HOG_LTP_COLOR",
			"HOG_LBP_LTP_COLOR", "CSLBP_COLOR", "CSLTP_COLOR" };
	int i = 0;
	for (; i < nftype && featname[i] != fname; ++i)
		;
	if (i == nftype) {
		cerr << " \n Error! UnRecognized Feature Type " << fname << endl;
		exit(EXIT_FAILURE);
	}
	cout << "\n ____Using Feature Type = " << featname[i] << "_______________"
			<< endl;
	return (WindowFeatureTypes) i;

}
static FeatureTypes GetFeatureType(const string &fname) {
	const int nftype = 21;
	string featname[] = { "HOG", "LBP", "LTP", "LBPLTP", "HOGLBP", "HOGLTP",
			"HOGLBPLTP", "MultiLevelLBP", "MultiLevelLTP", "MultiLevelLBPLTP",
			"EncodedLTP", "EncodedLBP", "EncodedLBPPatch",
			"EncodedPatchPlusOther", "EncodedGenericPatch",
			"MultipleEncodedGenericPatchFeatures", "EncodedLTPPatch",
			"EncodedRFPatchFeatures", "LBPBlock", "FT_CSLBP", "FT_CSLTP" };
	int i = 0;
	for (; i < nftype && featname[i] != fname; ++i)
		;
	if (i == nftype) {
		cerr << " \n Error! UnRecognized Feature Type " << fname << endl;
		exit(EXIT_FAILURE);
	}
	cout << "\n Using Feature Type = " << featname[i] << endl;
	return (FeatureTypes) i;
}
static string GetPyramidType(PyramidType ptype) {
	const int ptlen = 3;
	string pstring[] = { "NoPyramid", "SingleRes", "DoubleRes" };
	return pstring[MIN(ptype,ptlen-1)];
}

struct HOGParams {
	HOGParams() {
		cell = 8;
		norient = 36;
		chan = RGB;
	}
	UINT cell;
	UINT norient;
	Channel chan;
};
struct LBPParams { //
	LBPParams() // default parameters
	{
		cellsize = stride = 8; // 8 pixels cell stride size
		nbins = 59; // uniform bins
		radius = 1; // radius =1;
		npoints = 8; // 8 surrounding points in circle of radius r
		addfeat = true; // add feat channels R,G,B
		normfeatsep = false; // normalize each channel separately or combined
		featnorm = L1sqrt;
		nflevels = 1; // number of feature levels;
		folded = false; // obselete
		gridtype = Circular;
		tolerance = (5.0 / 255.0); // threshold range...
		lbptolerance = 0; // threshold range...
		nthrlevels = 1; // number of threshold levels
		ptype = PT_Circular;
		extent = 16;
		useblocknorm = false;
		issigned = false;
		usegradmag = false;
		usemeanval = false;
		hmethod = HM_Bilinear;
		xcellsize = ycellsize = -1;
		rinvariant = false;
	}
	static string GetHistogramMethod(const HistogramMethod &hm) {
		string shm[] = { "HM_Discrete-Binning", "HM_Bilinear-Interpolation" };
		return shm[hm];
	}
	void PrintParamInfo() const {
		cout << "-------------------LBPParams---------------------" << endl
				<< " cellsize = " << cellsize << ", stride =" << stride << endl
				<< " nbins = " << nbins << ", radius = " << radius << endl
				<< " npoints =" << npoints << ", nflevels=" << nflevels << endl
				<< " folded =" << folded << ", norm =" << featnorm << endl
				<< " addfeat = " << addfeat << ", normfeatsep=" << normfeatsep
				<< endl << "LTP tolerance =" << tolerance << "LBTP tolerance ="
				<< lbptolerance << ", extent =" << extent << endl
				<< " Patchtype = " << ptype << ", nthrlevels=" << nthrlevels
				<< endl << "useblocknorm = " << useblocknorm << " usesigned= "
				<< issigned << " usegradimage= " << usegradmag << endl
				<< " usemeanval= " << usemeanval << endl << " Using "
				<< GetHistogramMethod(hmethod)
				<< " method for generating histograms"
				<< (rinvariant == true ? "Using Rotation Invariant Features" :
						"") << "----------------------------------------"
				<< endl;
	}
	UINT cellsize, stride, nbins, radius, npoints;
	int xcellsize, ycellsize; /// redundant (mostly used only in the case of features computation for texture classification)...
	bool addfeat, normfeatsep;
	NORM featnorm;
	UINT nflevels;
	bool folded;
	GridType gridtype;
	REAL tolerance; // used for the ltp features
	REAL lbptolerance; // used for some types of lbp features
	UINT nthrlevels;
	PatchType ptype;
	UINT extent;
	bool useblocknorm;
	bool usegradmag;
	bool issigned;
	bool usemeanval; // use mean of the
	bool rinvariant; /// use rotation invarinat features
	HistogramMethod hmethod;
};
/**
 * QuantizationType represent type of method used for encoding the input codes
 * QT_VQ: Hard KMeans coding to a single cluster
 * QT_Median:
 * QT_LLC: Locality-Constrained linear coding
 */
enum QuantizationType {
	QT_VQ, QT_Median, QT_LLC, QT_Discrim
};
enum ClusteringDistanceMetric {
	CDM_L2, CDM_L1
};
struct CBParams {
	CBParams() {
		ocinfo = NULL;
		cbtype = CBT_PosCrop;
		codebooksize = 100;
		vimgfname = "trainval.txt";
		nimgfname = vimgfname;
		nwinsampled = 10;
		nclusrounds = 10;
		dmetric = CDM_L2;
		inituniform = false;
		softq = QT_VQ;
		meanp = 1;
		cboffset = "";
		knn = 5;
		pmethod = PM_Sum;
		sx = sy = ex = ey = -1;
	}
	static string GetCodeBookType(const CBType cbtype) {
		const int ncbtypes = 14;
		string scbtype[] = { "CBT_Complete", "CBT_PosCrop",
				"CBT_PosCropNegSampled", "CBT_PosPyramid",
				"CBT_PosPyramidNegSampled", "CBT_ConstantPosCellNeg",
				"CBT_ConstantPosCellNegSingleQuantization", "CBT_PosCellNeg",
				"CBT_PosCellNegSep", "CBT_Negative", "CBT_NegativeThenPosCrop",
				"CBT_MultiClass", "CBT_NegativePosDiscriminant",
				"CBT_CellQuantizeAll" };
		return scbtype[MIN(ncbtypes - 1, cbtype)];
	}
	static string GetQuantizationType(QuantizationType ctype) {
		string cnames[] = { "QT_VQ", "QT_Median", "QT_LLC", "QT_Discriminant" };
		return cnames[MIN(ctype, 3)];
	}
	static string GetPoolingMethodType(PoolingMethod ctype) {
		string cnames[] = { "PM_Sum", "PM_Max" };
		return cnames[MIN(ctype, 1)];
	}
	void PrintParamInfo() const {
		cout << " With Code Book Size = " << codebooksize << endl
				<< " Code Book Type = " << GetCodeBookType(cbtype) << endl;
		if (!nimgfname.empty())
			cout << " Negative Images are read from the file " << nimgfname
					<< endl;
		if (!vimgfname.empty())
			cout << " Validation Images are from the file " << vimgfname
					<< endl;

		cout
				<< (nwinsampled != 0 ? " Number of Negative Windows Sampled " :
						"Pixels from the complete Image are used for CodeBook Construction")
				<< nwinsampled << " Number of Clustering Rounds ="
				<< nclusrounds << endl << " Soft Quantization Type = "
				<< GetQuantizationType(softq) << endl;
		if ((softq) != QT_VQ)
			cout << "  Threshold (which percentage of mean distance) " << meanp
					<< "  KNN used for LCC Mapping  " << knn << "  and Using "
					<< GetPoolingMethodType(pmethod) << "   " << endl;

		cout
				<< (cboffset.empty() ? " CodeBook will be learned and lookup table will be built" :
						" LookUp Table will be loaded from the file with suffix")
				<< cboffset << endl;
		if (sx >= 0 && sy >= 0)
			cout << " Patch(Part) Geometry = " << sx << " , " << sy << " , "
					<< ex << " , " << ey << endl;

	}
	UINT codebooksize;
	CBType cbtype;
	ClassInfo *ocinfo; // for positive class ...
	string nimgfname; // either negative file name
	string vimgfname;
	UINT nwinsampled; //
	UINT nclusrounds;
	ClusteringDistanceMetric dmetric;
	bool inituniform; // initalize from uniform codes...
	QuantizationType softq; // do soft quantization of codes
	REAL meanp; // which percentage of an example mean distance from all the clusters to use as threshold
//\ie - code_i=max(0-(mean_dis_from_all_clusters - dist_i));
	UINT knn;
	string cboffset; // flag used to read the data-type., if mentioned a file containing the cluster information is read from hard-disk
	PoolingMethod pmethod;
	int sx, sy, ex, ey; /// geometry used for specify the patch (part) coordinates where the codes will be extracted for codebook learning
};
struct PatchParams {
	PatchParams() {
		size = 5;
		cpatchsize = 0;
		stride = 1;
		type = PT_Circular;
		ncodelevels = 2; // Binary LBP, else LTP...
		hnormtype = EHNT_Both;
		cbparams = NULL;
		tolval.resize(1, 5.0 / 255.0);
		petype = PET_LTP;
		const int nptypes = 14;
		UINT tvar[] = { 1, 1, 1, 1, 1, 2, 2, 2, 2, 4, 4, 1, 1, 3 };
		npmult = new UINT[nptypes];
		copy(tvar, tvar + nptypes, npmult);
		pcount = 0;

	}
	PatchParams(const PatchParams& pparam) {
		*this = pparam; /*call by default operator*/
	}
	~PatchParams() {
		//		delete[] npmult;
	}

	void PrintParamInfo() const {
		string ehnt[] = { "EHNT_Local", "EHNT_Global", // to do normalization of codebook quantization at Window level
				"EHNT_Both" };
		string pet[] = { "PET_LBP", "PET_LTP", "PET_CSLBP", "PET_CSLTP",
				"PET_SplitLTP", "PET_CSLTP_LTP" };
		cout
				<< "\n------------------------------------ LQP(Patch) Parameters Info"
				<< "------------------------------------\n" << " Patch Type ="
				<< GetPatchName(type) << ", Patch Size = " << size
				<< ", Circular Patch Size = " << cpatchsize
				<< ", Patch Stride = " << stride << " (Considering Every "
				<< pstep << " Pixel" << ")\n Number of Code Levels = "
				<< ncodelevels << "\n Tolerance values = " << tolval
				<< " PatchEncodingType = " << pet[petype] << " Pruning Count = "
				<< pcount << endl;

		if (cbparams)
			cbparams->PrintParamInfo();
		cout << " Cell Based CodeBook Histogarm Normalization Type ="
				<< ehnt[hnormtype] << endl
				<< "------------------------------------------------------------------"
				<< endl;
	}
	static string GetPatchEncodingType(const PatchEncodingType &petype) {
		string pet[] = { "PET_LBP", "PET_LTP", "PET_CSLBP", "PET_CSLTP",
				"PET_SplitLTP", "PET_CSLTP_LTP" };
		return pet[petype];
	}
	static string GetPatchName(const PatchType& ptype) {
		string ptypename[] = { "HorizontalStrip", "VerticalStrip",
				"CircularPatch", "DiagonalStrip", "ADiagonalStrip", "HVStrip",
				"DAStrip", "HVCStrip", "DACStrip", "HVDAStrip", "HVNStrip",
				"HaarPyramid", "CircularShallow", "Circular+Angular+Radial" };
		//		PT_HStrip,
		//		PT_VStrip,
		//		PT_Circular,
		//		PT_Diag,
		//		PT_ADiag,
		//		PT_HVStrip,
		//		PT_DAStrip,
		//		PT_HVDAStrip,
		//		PT_HVNStrip, // compare neighbour pixels instead of comparison with the centeral pixels
		//		PT_HVCNStrip
		return ptypename[ptype];
	}
	static string GetPatchType(const PatchType &ptype) { // input the patch number output string...
		const int nptype = 5;
		string patchname[] = { "HorizontalStrip", "VerticalStrip", "Circular",
				"DiagonalStrip", "ADiagonalStrip" };

		return patchname[MIN(nptype-1,ptype)];
	}
	UINT size, cpatchsize;
	UINT stride;
	PatchType type;
	UINT pstep; //pixel step to be used in horizontal,vertical &n square
	UINT ncodelevels; // 2->binary = lbp, 3->ternary = ltp, 5->quinary =2 different ltp's with different threshold levels...
	// EncodedMultipleCBPatchFeaturesHistogram
	EncodedHistogramNormType hnormtype;
	vector<REAL> tolval;
	CBParams *cbparams;
	PatchEncodingType petype;
	UINT *npmult;
	UINT pcount; /// count used for pruning the codes with small counts
};
enum FeatVectorNormalType {
	FVNT_range, FVNT_standard, FVNT_midrange, FVNT_none, FVNT_range2
};

struct NLFParams {
	NLFParams() {
		lplsdim = qplsdim = 14;
		ctype = QuadEmbedSVM; //
		fvntype = FVNT_range;
		polydegree = 0; // if not zero than add that power polynomials to Full Quadratic features...
		qcost = 1;
	}
	static string GetFeatureVectorNormalType(
			const FeatVectorNormalType& fvntype) {
		string val[] = { "FVNT_range", "FVNT_standard", "FVNT_midrange",
				"FVNT_none", "FVNT_range2" };
		return val[MIN(fvntype, 4)];
	}
	void PrintParams() {
		cout << " --------------- NLF-Params----------------" << endl
				<< " PLS Linear Dim" << lplsdim << " PLS Quadratic Dim = "
				<< qplsdim << " Classifier Type=" << GetClassifierType(ctype)
				<< " Feature vector Normalization type = "
				<< GetFeatureVectorNormalType(fvntype) << " PolyNomial Degree ="
				<< polydegree << " QCost =" << qcost << endl
				<< " --------------- NLF-Params----------------" << endl;
	}
	void WriteFile(ofstream &ofile) {
		ofile.write((char*) &ctype, sizeof(ClassifierType));
		ofile.write((char*) &lplsdim, sizeof(UINT));
		ofile.write((char*) &qplsdim, sizeof(UINT));
		ofile.write((char*) &polydegree, sizeof(UINT));
		ofile.write((char*) &fvntype, sizeof(FeatVectorNormalType));
		ofile.write((char*) &qcost, sizeof(REAL));
	}
	void ReadFile(ifstream &ifile) {
		ifile.read((char*) &ctype, sizeof(ClassifierType));
		ifile.read((char*) &lplsdim, sizeof(UINT));
		ifile.read((char*) &qplsdim, sizeof(UINT));
		ifile.read((char*) &polydegree, sizeof(UINT));
		ifile.read((char*) &fvntype, sizeof(FeatVectorNormalType));
		ifile.read((char*) &qcost, sizeof(REAL));
	}
	UINT lplsdim, qplsdim; // plsdim2 to be used in cascade cases...
	ClassifierType ctype; // if true use linear+ non-linear filter ...
	FeatVectorNormalType fvntype; // feature vector normalization type...
	UINT polydegree; // polynomial degree for coordinate-quad..
	REAL qcost; // Quadratic coefficent scaling factor...
};
#endif /* DEFINITIONS_H_ */
