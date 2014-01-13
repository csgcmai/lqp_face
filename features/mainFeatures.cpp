/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/
#include "../features/includeHeaders.h"
#include "../util/processImage.h"
#include "../features/customValidators.h"
#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <iterator>

UINT Features::xcell = 8, Features::ycell = 8, Features::norient = 36,
		Features::foffset = 8, Features::cellsize = 8;
//TODO: Remove it
//#include "../util/boosting.h"
//UINT DecisionStumps::dimfeature = 0;

void WriteFeatures(const string & fname, const UINT & xbmax, const UINT & ybmax,
		vector<REAL> &feat) {
	ofstream ofile(fname.c_str(), ios::out | ios::binary);
	if (!ofile) {
		cout << " Couldn't open the file" << fname << endl;
		exit(EXIT_FAILURE);
	}
	UINT version = 001, rsize = feat.size() / ybmax;
	ofile.write((char*) &version, sizeof(UINT));
	ofile.write((char*) &ybmax, sizeof(UINT));
	ofile.write((char*) &rsize, sizeof(UINT));
	REAL *featptr = &feat[0];
	ofile.write((char*) featptr, sizeof(REAL) * feat.size());
	ofile.close();
}

class Process {
public:
	Process(po::variables_map &obj, LBPParams &lbpparam, PatchParams &pparam,
			CBParams &cbparams) {
		fobj = NULL;
		pim = NULL;
		ProcessInput(obj, lbpparam, pparam, cbparams);
	}
	void
	ProcessInput(po::variables_map&, LBPParams &, PatchParams &, CBParams&);
	void DumpFeatures(const string &lstfile, const string &offile, UINT width,
			UINT height);
	void DumpGaborLQPFeatures(const string &lstfile, UINT nglevels, UINT width,
			UINT height);
	void DumpPartFeatures(const string &lstfile, const string &feattype,
			int xcellsize, int ycellsize, int sx, int sy, int ex, int ey);

	void InitializeImageProcessor(po::variables_map&);
	void ParseLocalPatternParams(po::variables_map&, LBPParams&);
	void ParsePatchParams(po::variables_map&, LBPParams&, CBParams&);
	string InitializeLocalPatternFeatures(const FeatureTypes &ftype,
			LBPParams &lbpparam, UINT width, UINT height, UINT nplevels);
	~Process() {
		delete fobj;
		delete pim;
	}
private:
	Features *fobj;
	ProcessImage *pim;
	string dir;
};
void Process::InitializeImageProcessor(po::variables_map&obj)
// initialize different Image Processing Parameters...
		{
	//§§§§§§§§§§§§§§§§§§§§§§§    Process Image Parameter       ////////////////////:::
	Channel ch = obj["Color-Channel"].as<Channel>(); // which color channel to use..
	REAL gnorm = obj["Gamma-Normalization"].as<REAL>();
	REAL alpha = obj["Cont-Alpha"].as<REAL>();
	REAL tou = (REAL) obj["Cont-Tou"].as<REAL>() / 255;
	REAL k1 = obj["DoG-K1"].as<REAL>();
	REAL k2 = obj["DoG-K2"].as<REAL>();
	bool dog = obj["DoG"].as<bool>();
	bool coneq = obj["Cont-Equal"].as<bool>();
	pim = new ProcessImage(gnorm, dog, k1, k2, coneq, tou, alpha, ch);
}

void Process::DumpFeatures(const string &lstfile, const string &feattype,
		UINT width, UINT height) {
	vector<REAL> feat(fobj->GetDim(width, height));
	UINT xbmax = (UINT) (round((REAL) width / (REAL) fobj->cellsize)), ybmax =
			(UINT) (round((REAL) height / (REAL) fobj->cellsize)), wcount = 0;
	cout << " \n Feature Size for Window " << width << " X " << height << " ="
			<< feat.size() << endl << " Xbmax=" << xbmax << " , Ybmax= "
			<< ybmax << " , Cell Size " << fobj->cellsize << endl;

	vector<string> lstimages;
	ParseListFile(lstfile, lstimages);
	Image im;
	stringstream ss;
	for (UINT i = 0; i < lstimages.size(); ++i) {
		im.read(lstimages[i]);
		cout << "\n Processing Image " << lstimages[i] << flush << endl;
		fobj->InitalizeMaps(im, NoPyramid);
		fobj->GetFeatures(0, 0, 0, width, height, feat);
		ss.str("");
		ss << SplitFilenameFromDir(lstimages[i]) << "-" << feattype;
		WriteFeatures(ss.str(), xbmax, ybmax, feat);
		cout << " Number = " << wcount << " / " << lstimages.size() * 2;
		++wcount;
		// Flip the image
		im.flop();
		fobj->InitalizeMaps(im, NoPyramid);
		fobj->GetFeatures(0, 0, 0, width, height, feat);
		ss.str("");
		ss << SplitFilenameFromDir(lstimages[i]) << "-" << feattype << "-f";
		WriteFeatures(ss.str(), xbmax, ybmax, feat);
		cout << " Number = " << wcount << " / " << lstimages.size() * 2;
		++wcount;
	}
}
void Process::DumpGaborLQPFeatures(const string &lstfile, UINT nglevels,
		UINT width, UINT height) {
	vector<REAL> feat(fobj->GetDim(width, height));
	UINT xbmax = (UINT) (round((REAL) width / (REAL) fobj->cellsize)), ybmax =
			(UINT) (round((REAL) height / (REAL) fobj->cellsize)), wcount = 0;
	cout << " \n Feature Size for Window " << width << " X " << height << " ="
			<< feat.size() << endl << " Xbmax=" << xbmax << " , Ybmax= "
			<< ybmax << " , Cell Size " << fobj->cellsize << endl;

	vector<string> lstimages;
	ParseListFile(lstfile, lstimages);
	stringstream ss;

	if (lstimages.size() % (nglevels * 2 + 1) != 0) {
		cout
				<< "\n Error! Number of images in the list file is not factor of (Number of Gabor-LQP Levels)"
				<< endl;
		exit(EXIT_FAILURE);
	}
	vector<Image> image(nglevels * 2 + 1), fimage(nglevels * 2 + 1);
	string feattype = "glqp";
	for (UINT i = nglevels; i < lstimages.size(); i += (nglevels * 2 + 1)) {
		cout << " Processing Image # = " << (i / nglevels + 1) << " "
				<< lstimages[i] << endl;
		for (UINT j = i - nglevels, count = 0; j <= i + nglevels;
				++j, ++count) {
			image[count].read(lstimages[j]);
			fimage[count] = image[count];
			fimage[count].flop();
		}

		fobj->InitalizeMaps(image);
		fobj->GetFeatures(0, 0, 0, width, height, feat);
		ss.str("");
		ss << SplitFilenameFromDir(lstimages[i]) << "-" << feattype;
		WriteFeatures(ss.str(), xbmax, ybmax, feat);
		cout << " Number = " << wcount << " / " << lstimages.size() * 2;
		++wcount;

		fobj->InitalizeMaps(fimage);
		fobj->GetFeatures(0, 0, 0, width, height, feat);
		ss.str("");
		ss << SplitFilenameFromDir(lstimages[i]) << "-" << feattype << "-f";
		WriteFeatures(ss.str(), xbmax, ybmax, feat);
		cout << " Number = " << wcount << " / " << lstimages.size() * 2;
		++wcount;
	}
}
void Process::DumpPartFeatures(const string &lstfile, const string &feattype,
		int xcellsize, int ycellsize, int sx, int sy, int ex, int ey) {
	int width = ex - sx + 1, height = ey - sy + 1;
	vector<REAL> feat(fobj->GetDim(width, height));
	UINT xbmax = (UINT) (round((REAL) width / (REAL) xcellsize)), ybmax =
			(UINT) (round((REAL) height / (REAL) ycellsize)), wcount = 0;
	cout << " \n Feature Size for Window " << width << " X " << height << " ="
			<< feat.size() << endl << " Xbmax=" << xbmax << " , Ybmax= "
			<< ybmax << " , xCell Size " << xcellsize << " , yCell Size "
			<< ycellsize << endl;

	vector<string> lstimages;
	ParseListFile(lstfile, lstimages);
	Image im;
	stringstream ss;
	for (UINT i = 0; i < lstimages.size(); ++i) {
		im.read(lstimages[i]);
		cout << "\n Processing Image " << lstimages[i] << flush << endl;
		fobj->ComputeTextureFeatures(im, xcellsize, ycellsize, sx, sy, ex, ey,
				feat);
		ss.str("");
		ss << SplitFilenameFromDir(lstimages[i]) << "-" << feattype;
		WriteFeatures(ss.str(), xbmax, ybmax, feat);
		cout << " Number = " << wcount << " / " << lstimages.size() * 2;
		++wcount;
		// Flip the image
		im.flop();
		int tsx = im.columns() - ex, tex = im.columns() - sx + 1;
		fobj->ComputeTextureFeatures(im, xcellsize, ycellsize, tsx, sy, tex, ey,
				feat);
		ss.str("");
		ss << SplitFilenameFromDir(lstimages[i]) << "-" << feattype << "-f";
		WriteFeatures(ss.str(), xbmax, ybmax, feat);
		cout << " Number = " << wcount << " / " << lstimages.size() * 2;
		++wcount;
	}
}

string Process::InitializeLocalPatternFeatures(const FeatureTypes &ftype,
		LBPParams &lbpparam, UINT width, UINT height, UINT nplevels) {

	string feattype;
#ifndef WINDOW_CLASSIFIER
	/*For different Colour Channels, and Gray Scale with L1, L1Sqrt, L2, Normalization...*/
	/*There are two different types of implementations (i.e With Cell & Without Cell Suffix)
	 1. Cell Suffix: Features are computed for the complete image at Once
	 (image level) and features can be extracted with step size=single cell.
	 Note that number of pixels = single cell size are ignored around the four sides to avoid the boundary artifacts.
	 2. Without Cell Suffix: Features are computed at a given window (window level) thus can be computed with any step size.
	 */
	switch (ftype) {
		case FT_LBP:
		feattype = "LBP";
		fobj = new LBPFeaturesRGBBlock(lbpparam, width, height, nplevels,
				*pim);
		break;
		case FT_LTP:
		feattype = "LTP";
		fobj = new OnlyLTPFeaturesRGBCell(lbpparam, width, height,
				nplevels, *pim);
		break;
		case FT_LBPLTP:
		feattype = "LBPLTP";
		fobj = new LTPFeaturesRGBCell(lbpparam, width, height, nplevels,
				*pim);
		break;
	}
#else
	/*For Backward Compatibility, Works at window levels*/
	switch (ftype) {
	case FT_LBP:
		feattype = "LBP";
		fobj = new LBPFeatures(lbpparam, width, height, nplevels, *pim);
		break;
	case FT_LTP:
		feattype = "LTP";
		fobj = new OnlyLTPFeatures(lbpparam, width, height, nplevels, *pim);
		break;
	case FT_LBPLTP:
		feattype = "LBPLTP";
		fobj = new LTPFeatures(lbpparam, width, height, nplevels, *pim);
		break;
	}
#endif
	return feattype;
}
void Process::ProcessInput(po::variables_map &obj, LBPParams &lbpparam,
		PatchParams &pparam, CBParams &cbparam) {
	string feattype;
	UINT width, height;
	InitializeImageProcessor(obj); // initialize Image Processor...

	pparam.cbparams = &cbparam;
	pparam.ncodelevels = pparam.petype;

	width = obj["Win-Width"].as<int>(); // Reference Width & Height of images or image windows
	height = obj["Win-Height"].as<int>();

	string tfile;
	if (obj.count("TrainingFile"))

		tfile = obj["TrainingFile"].as<string>();
	else
		tfile = cbparam.vimgfname;
	//	UINT xstride = lbpparam.cellsize, ystride = xstride;
	UINT nplevels = 1; // obj["Number-Levels"].as<int> (); // number of pyramid levels

	// Parse Miscelleneous parameters
	FeatureTypes ftype = obj["FeatureType"].as<FeatureTypes>();
	int nglevels = obj["GaborLQPLevels"].as<int>();

	/* Region window coordinates used to learn Region level codebook (these coordinates can be used to compute part-level features
	 * or features for a given region)
	 */
	int sx, sy, ex, ey, xcellsize, ycellsize;
	sx = obj["WinX1"].as<int>(); /* */
	sy = obj["WinY1"].as<int>();
	ex = obj["WinX2"].as<int>();
	ey = obj["WinY2"].as<int>();
	xcellsize = lbpparam.cellsize;
	ycellsize = lbpparam.cellsize;

	cbparam.sx = sx;
	cbparam.sy = sy;
	cbparam.ex = ex;
	cbparam.ey = ey;

	/* Normalize the tolerance Values to be in the range [0,1] */
	lbpparam.tolerance = lbpparam.tolerance / 255.0;
	lbpparam.lbptolerance = lbpparam.lbptolerance / 255.0;

	if (obj.count("ToleranceLevels")) {
		pparam.tolval = obj["ToleranceLevels"].as<vector<REAL> >();
		for (vector<REAL>::iterator iter = pparam.tolval.begin();
				iter != pparam.tolval.end(); ++iter)
			*iter = *iter / 255;

	} else {
		pparam.tolval.resize(1, 0);
		pparam.tolval[0] = lbpparam.tolerance;
	}

	/* Now Compute the features */
	switch (ftype) {
	case FT_LBP:
	case FT_LTP:
	case FT_LBPLTP:
		feattype = InitializeLocalPatternFeatures(ftype, lbpparam, width,
				height, nplevels);
		break;

	case FT_GaborLQP:
	case FT_LQPSoft:
	case FT_SeparateGaborLQP:
	case FT_EncodedLTP: // Generate Dictionary of the features...
	case FT_EncodedLBP:
	case FT_EncodedLBPPatch:
	case FT_EncodedLTPPatch:
	case FT_EncodedGenericPatch:
	case FT_EncodedPatchPlusOther:
	case FT_EncodedRFPatchFeatures:
	case FT_MultipleEncodedGenericPatchFeatures:
		// always uses  non-symmetric, current implementation does not learn symmetric vocabularies
		feattype = "LQP";
		EncodedPatchFeatures *tobj;
		switch (ftype) {
		case FT_EncodedGenericPatch:
		case FT_LQPSoft:
		case FT_SeparateGaborLQP:
		case FT_GaborLQP:
			switch (cbparam.cbtype) {
			case CBT_CellQuantizeAll:
				fobj = new MultipleCBQuantizeLQPFeatures(lbpparam, width,
						height, nplevels, *pim, pparam);
				break;
			case CBT_PosCrop:
			case CBT_PosPyramid:
			case CBT_Complete:
			case CBT_NegativePosDiscriminant:
			default:
				if (ftype == FT_LQPSoft)
					fobj = new SoftLQPFeatures(lbpparam, width, height,
							nplevels, *pim, pparam);
				else if (ftype == FT_GaborLQP)
					fobj = new GaborLQPFeatures(lbpparam, width, height,
							nplevels, *pim, pparam, nglevels);
				else if (ftype == FT_SeparateGaborLQP) {
					cout << " \n Featurey TYpe = FT_SeparateGaborLQP ";
					fobj = new SeparateGaborLQPFeatures(lbpparam, width, height,
							nplevels, *pim, pparam, nglevels);
				} else
					fobj = new EncodedPatchFeatures(lbpparam, width, height,
							nplevels, *pim, pparam);
				break;
			}
			break;
		}
		fobj->GenerateCodeBook(); //
		break;
	default:
		cout
				<< " \n UnRecognized Feature Type: Please Select a Valid Feature Type"
				<< endl;
		exit(EXIT_FAILURE);
	}
	if (obj.count("OnlyGenerateCodeBook")) {
		cout << "\n Program was called only for the generation of codebooks"
				<< ", It's finished so Quittttttttttting..... " << endl;
		exit(EXIT_SUCCESS);
	}
	if (ftype == FT_GaborLQP || ftype == FT_SeparateGaborLQP)
		DumpGaborLQPFeatures(tfile, nglevels, width, height);
	else if (cbparam.sx >= 0 && cbparam.sy >= 0)
		DumpPartFeatures(tfile, feattype, xcellsize, ycellsize, sx, sy, ex, ey);
	else
		DumpFeatures(tfile, feattype, width, height);
}
void AddLBPOptions(po::options_description &lbp_options, LBPParams &lbpparam) {

	//	lbpparam.gridtype = (GridType) obj["LBPGridType"].as<int> ();
	//	lbpparam.addfeat = (obj["Add-RGBLBP"].as<int> () == 1 ? true : false);
	//	lbpparam.cellsize = obj["Cell-Size"].as<int> ();
	//	lbpparam.radius = obj["Radius"].as<int> ();
	//	lbpparam.npoints = obj["Number-Points"].as<int> ();
	//	lbpparam.nflevels = obj["Feature-Levels"].as<int> ();
	//	lbpparam.stride = obj["LBP-Stride"].as<int> ();
	//	lbpparam.nbins = obj["LBP-NBins"].as<int> ();
	//
	//	lbpparam.tolerance = ((REAL) obj["LTP-Tolerance"].as<int> ()) / 255.0;
	//	lbpparam.lbptolerance = ((REAL) obj["LBP-Tolerance"].as<int> ()) / 255.0;
	//
	//	lbpparam.normfeatsep = (obj["Norm-Sep"].as<int> () == 1 ? true : false);
	//	lbpparam.featnorm = (NORM) obj["Normalization"].as<int> ();
	//	//	lbpparam.folded = obj.count("Fold-Features"); // obselete ...
	//	lbpparam.nthrlevels = obj["LTPLevels"].as<int> ();
	//	lbpparam.extent = obj["Extent"].as<int> ();
	//	lbpparam.hmethod = (HistogramMethod) obj["HistogramMethod"].as<int> ();
	//	lbpparam.useblocknorm = obj.count("BlockNorm");
	//	lbpparam.usemeanval = obj.count("UseMeanCenVal");
	//	lbpparam.usegradmag = obj.count("UseGradMag");
	lbp_options.add_options()("LBPGridType",
			po::value<GridType>(&lbpparam.gridtype)->default_value(Circular),
			"0. Circular\n 2. Rectangular")("Add-RGBLBP",
			po::value<bool>(&lbpparam.addfeat)->default_value(1),
			"Merge RGB Channels \n0. Concat RGB-Channels Hist\n1. Add RGB-Channels Hist  ")(
			"Norm-Sep",
			po::value<bool>(&lbpparam.normfeatsep)->default_value(0),
			"For Above Option 1: \n0. Add RGB-channels & then Normalize Once \n1. Normalize Each RGB channel & then add")/*
	 // Not used
	 ("Skip",
	 po::value<int>()->default_value(0),
	 "Skip 'Skip' Number of Pixels on Borders")*/("LBP-NBins",
			po::value<UINT>(&lbpparam.nbins)->default_value(59),
			"Number of LBPBins can be either 256 or 59")("Cell-Size",
			po::value<UINT>(&lbpparam.cellsize)->default_value(10),
			"LBP Feature: Cell size")("Radius",
			po::value<UINT>(&lbpparam.radius)->default_value(1),
			"Radius of LBP Features")("Number-Points",
			po::value<UINT>(&lbpparam.npoints)->default_value(8),
			"Number of Points around each point")("Feature-Levels",
			po::value<UINT>(&lbpparam.nflevels)->default_value(1),
			"Number of Feature Levels(default = 1)")("LTP-Tolerance",
			po::value<float>(&lbpparam.tolerance)->default_value(5.0),
			"+,- Threshold for the LTP Features (in the range [0,255])")(
			"LBP-Tolerance",
			po::value<float>(&lbpparam.lbptolerance)->default_value(0.0),
			"Threshold used during some types of LBP Features...")(
			"Normalization",
			po::value<NORM>(&lbpparam.featnorm)->default_value(L1sqrt),
			"Histogram Normalization :\n 0. L1 Norm  \n 1. L1Sqrt Norm"
					" \n 2. L2 Norm \n 3. LOG2(some cases) \n 6. Average(of RGB Channels for Texture Case)")(
			"LBP-Stride", po::value<UINT>(&lbpparam.stride)->default_value(16),
			"LBP Cell Stride-Size in Pixels")("BlockNorm",
			po::bool_switch(&lbpparam.useblocknorm),
			"Use Neigbhour block normalization for LBP Same as HOG(by default false)")(
			"UseMeanCenVal", po::bool_switch(&lbpparam.usemeanval),
			"Use Mean of Square box radius X radius instead of center pixel value(by default false)")(
			"UseGradMag", po::bool_switch(&lbpparam.usegradmag),
			"Use Gradient as voting for the LBP/LTP histogram  bins(by default false)")(
			"LTPLevels",
			po::value<UINT>(&lbpparam.nthrlevels)->default_value(1),
			"Number of LTP Level's(Threshold is used as base)")(
			"HistogramMethod",
			po::value<HistogramMethod>(&lbpparam.hmethod)->default_value(
					HM_Bilinear), "Binning Method for Histograms "
					"\n 0. Discrete Binning "
					"\n 1. Bilinear Interpolation")("Extent",
			po::value<UINT>(&lbpparam.extent)->default_value(16),
			"Cropped Image Border");

}
void AddMiscOptions(po::options_description & misc_options) {
	misc_options.add_options()("FeatureType",
			po::value<FeatureTypes>()->default_value(FT_EncodedGenericPatch),
			"Features Type :\n"
					"1. LBP\n"
					"2. LTP\n"
					"3. LBPLTP\n"
					"14.LQP (EncodedPatchFeatures -- covers all the encoded features)\n"
					"22.GaborLQP -- Multiple Orientation-Scale Co-occurence LQP\n"
					"23.SeparateGaborLQP -- Combines LQP Computed over multiple orientations & scales at codebook level")(
			"Win-Width", po::value<int>()->default_value(64), "Window Width")(
			"Win-Height", po::value<int>()->default_value(128), "Window Height")/*(
	 "X-Stride", po::value<int>()->default_value(8),
	 "X Stride in Pixels (Must be a multiple of cell size)")("Y-Stride",
	 po::value<int>()->default_value(8), "Y Stride in Pixels")*/(
			"OnlyGenerateCodeBook",
			"LQP: Only Generate CodeBook and then quit...")("GaborLQPLevels",
			po::value<int>()->default_value(1),
			"Number of Scale and Orientation to include in Gabor-LQP")
	///Window Specification for learning codebooks & extracting features from a given region of the training images...
	("WinX1", po::value<int>()->default_value(-1),
			"Xcoord of window top left coordinate\n (for learning codebooks & extracting features from only the given region)")(
			"WinY1", po::value<int>()->default_value(-1),
			"Ycoord of window top left coordinate")("WinX2",
			po::value<int>()->default_value(-1),
			"Ycoord of window bottom right coordinate")("WinY2",
			po::value<int>()->default_value(-1),
			"Ycoord of window bottom right coordinate")("TrainingFile",
			po::value<string>(),
			"File contains list of training images for which features need to be computed")/*(
			 "Number-Levels", po::value<int>()->default_value(10),
			 "Number of Levels in each Pyramid Octave")("Pos-Train-File",
			 po::value<string>(), "File containing list of positive images")(
			 "Neg-Train-File", po::value<string>(),
			 "File containing list of negative images")

			 ///##################### PATCH PLUS FEATURE TYPE
			 ("PatchPlusFeatureType", po::value<int>()->default_value(0),
			 "Patch Plus Feature Type, One of FeatureType")("MEPF-ConfigFile",
			 po::value<string>()->default_value(""),
			 "MultipleEncodedPatchFeatures Configuration File")("NumberCircles",
			 po::value<int>()->default_value(1),
			 "Number of Concentric Circles for Computation of LBP/LTP")
			 //++++++++++++++++++++++++ Combined Features++++++++++++++++++++++++
			 ("Combined-Feature-Type,F", po::value<vector<int> >(),
			 "Specify Multiple features to have combined features from the Feature-Type")(
			 "Feature-ConfigFile", po::value<string>(),
			 "Feature Configuration File To Specify Various Combinations")*/;

}
void AddCBOptions(po::options_description & cb_options, CBParams &cbparam) {
	//	CBParams cbparam;

	//	cbparam.cbtype = (CBType) obj["CodeBook-Type"].as<int> ();
	//	cbparam.codebooksize = obj["CodeBookSize"].as<int> ();
	//	cbparam.nwinsampled = obj["CodeBook-NegSampled"].as<int> ();
	//	cbparam.nclusrounds = obj["ClusteringRounds"].as<int> ();
	//	cbparam.inituniform = false; // argument not used..
	//	cbparam.pmethod = (PoolingMethod) obj["CodeBook-Pooling"].as<int> ();
	//	cbparam.softq = (QuantizationType) obj["SoftQuantize"].as<int> ();
	//	cbparam.meanp = obj["MeanDivisor"].as<REAL> ();
	//	cbparam.cboffset = obj["ReadCodeBook"].as<string> ();
	//	cbparam.sx = sx;
	//	cbparam.ex = ex;
	//	cbparam.sy = sy;
	//	cbparam.ey = ey;

	// CodeBook Parameters
	cb_options.add_options()("CodeBookSize",
			po::value<UINT>(&cbparam.codebooksize)->default_value(100),
			"Code Book Size (Number of visual words) for Mapping the features")(
			"CodeBook-Type",
			po::value<CBType>(&cbparam.cbtype)->default_value(CBT_Complete),
			"Different CodeBook Types "
					"0. Complete - use complete training set to generate CB\n"
					"1. PosCrop - CB from positive crop images\n"
					"2. PosCropNegSampled - quantize against both poscrop and -ve CBs\n"
					"3. PosPyramid\n - quantize against scale-space positive images CBs\n"
					"4. PosPyramidNegSampled - quantize against both scale-space positive images and -ve CBs \n"
					"5. ConstantPosCellNeg - quantize each cell against all +ve and -ve CBss\n"
					"6. ConstantPosCellNegSingleQuantization - Window level histogram of cell CBss (obtained via option ConstantPosCellNeg.)\n"
					"7. PosCellNeg - Each Cell is mapped only against its own CB and -ve class CB; negatives cells are accumulated at window level \n"
					"8. PosCellNegSep - Same as PosCellNeg except negatives are not accumulated at windows levels\n"
					"9. Negative - quantize only against negative class codebook\n"
					"10.NegativeThenPosCrop - first initialize from negative class then purify with pos class\n"
					"11.MultiClass \n"
					"12.Discriminative- CodeBook is constructed from both positive and negative classes with an extra bin for codecount_pos/codecount_neg ratio\n"
					"13. A Separate CodeBook for each cell and each cell is quantized against all cells")(
			"CodeBook-DMetric", po::value<int>()->default_value(0),
			"Distance to Compute Code Book \n 0.L2 1.Hamming or L1")(
			"ClusteringRounds",
			po::value<UINT>(&cbparam.nclusrounds)->default_value(3),
			"Number of KMeans Clustering Rounds for generation of CodeBook")(
			"SoftQuantize",
			po::value<QuantizationType>(&cbparam.softq)->default_value(QT_VQ),
			"Different Vector Quantization Methods"
					"\n 0. VQ -- Normal KMeans"
					"\n 1. Mean "
					"\n 2. LLC Mapping"
					"\n 3. Discriminant")("CodeBook-Pooling",
			po::value<PoolingMethod>(&cbparam.pmethod)->default_value(PM_Sum),
			"Pooling Methods for codes"
					"\n 0. Sum Pooling"
					"\n 1. Max Pooling")("CodeBook-NegSampled",
			po::value<UINT>(&cbparam.nwinsampled)->default_value(40),
			"Number of -ve windows to Sampled For Negative CodeBook(0=Complete Negative Set)")(
			"MeanDivisor", po::value<REAL>(&cbparam.meanp)->default_value(1),
			"Divisor(md) of mean distance to use as distance limiting threshold ie max(0,mu/md-d_i)")(
			"ReadCodeBook",
			po::value<string>(&cbparam.cboffset)->default_value(""),
			"post-fix flag for file containing the "
					"\n LQP lookup table (if specified codebook is read from the given file otherwise no)")(
			"Validation-File",
			po::value<string>(&cbparam.vimgfname)->default_value(""),
			"Validation File Used for CodeBook Generation");
	;
}
void AddLQPOptions(po::options_description &lqp_options, PatchParams & pparam) {
	//	pparam.size = obj["Patch-PatchSize"].as<int> ();
	//	pparam.cpatchsize = obj["Patch-CircPatchSize"].as<int> ();
	//	pparam.stride = obj["Patch-PatchStride"].as<int> ();
	//	pparam.type = (PatchType) obj["PatchType"].as<int> ();
	//	pparam.pstep = obj["Patch-PixelStride"].as<int> ();
	//	pparam.ncodelevels = obj["Patch-CodingType"].as<int> ();
	//	pparam.pcount = obj["Patch-PruneCount"].as<int> ();
	//	pparam.petype = (PatchEncodingType) obj["Patch-CodingType"].as<int> ();
	//	if (obj.count("ToleranceLevels")) {
	//		pparam.tolval = obj["ToleranceLevels"].as<vector<REAL> > ();
	//		for (UINT i = 0; i < pparam.tolval.size(); ++i)
	//			pparam.tolval[i] /= 255;
	//	} else {
	//		pparam.tolval.resize(1, 0);
	//		pparam.tolval[0] = ((REAL) obj["LTP-Tolerance"].as<int> ()) / 255.0;
	//	}
	//	pparam.cbparams = &cbparam;

	// Image Processor
	lqp_options.add_options() // Image Processor
	// Patch Parameters
	("PatchType",
			po::value<PatchType>(&pparam.type)->default_value(PT_Circular),
			"Patch Type \n "
			//PT_HStrip, PT_VStrip, PT_Circular, PT_Diag, PT_ADiag
					"0. Horizontal Strip\n"
					"1. Vertical Strip\n"
					"2. Circular/Rect Patch\n"
					"3. Diagonal Strip\n"
					"4. ADiagonal Strip\n"
					"5. Hor+Vert Strip\n"
					"6. Diag+ADiag Strip\n"
					"7. Hor+Vert+Circ Strip\n"
					"8. Diag+ADiag+Circ Strip\n"
					"9. Hor+Vert+Diag+ADiag Strip\n"
					"10. Hor+Vert+Diag+ADiag Neighbour/comparsion Strip\n"
					"11. Haar Pyramid With 4x4 Neighbourhood\n"
					"12. Circular Shallow Only Outer Ring is considered\n"
					"13. Circular+Angular+ Radial difference among the pixels\n")(
			"Patch-PatchSize", po::value<UINT>(&pparam.size)->default_value(5),
			"Patch size (diameter) for Patch Features")("Patch-CodingType",
			po::value<PatchEncodingType>(&pparam.petype)->default_value(
					PET_LTP), "Patch Coding Type\n"
					"0. PET_LBP\n"
					"1. PET_LTP\n"
					"2. PET_CSLBP\n"
					"3. PET_CSLTP\n"
					"4. PET_SplitLTP\n")("Patch-CircPatchSize",
			po::value<UINT>(&pparam.cpatchsize)->default_value(0),
			"Circular Patch size works only with patch type 7&8")(
			"Patch-PatchStride",
			po::value<UINT>(&pparam.stride)->default_value(1),
			"Distance Between LQP Patches(2 means skip every 2nd Pixel)")(
			"Patch-PixelStride",
			po::value<UINT>(&pparam.pstep)->default_value(1),
			"Pixel Stride inside the patch")("ToleranceLevels",
			po::value<vector<REAL> >(&pparam.tolval)/*->default_value(
			 vector<REAL>(1, 5.0), "default")*/,
			" Tolerance"
					"Values for MultiLevelLQPs i.e Quinary (Input=Vector) in the range [0,255]")(
			"Patch-PruneCount",
			po::value<UINT>(&pparam.pcount)->default_value(0),
			"Prune the codes having count < PruneCount")

	("CodeBook-HistogramNormalization",
			po::value<EncodedHistogramNormType>(&pparam.hnormtype)->default_value(
					EHNT_Global),
			"CodeBook Histogram Normalization For Cell Based CodeBooks \n "
					"0. EHNT_Local\n 1. EHNT_Global\n2. EHNT_Both");

}
void AddImageProcessOptions(po::options_description &im_options
/*,ProcessImage *pim*/) {

	// Image Processor
	im_options.add_options()("Color-Channel",
			po::value<Channel>()->default_value(RGB),
			"Image Color Channel to use "
					"0. Red,\n"
					"1. Green,\n"
					"2. Blue,\n"
					"3. Gray,\n"
					"4. RGB,\n"
					"5. O1,\n"
					"6. O2,\n"
					"7. O3,\n"
					"8. HSL,\n"
					"9. LAB,\n"
					"10. WandellOCS,\n"
					"11. SandeOCS,\n"
					"12. ABSGradientMag\n")("Gamma-Normalization",
			po::value<REAL>()->default_value(1),
			"Perform Gamma Normalization (for LBP & LTP)")("DoG",
			po::bool_switch(), "Apply DoG Filter")("DoG-K1",
			po::value<REAL>()->default_value(1),
			"DoG filter Inner kernel StDev")("DoG-K2",
			po::value<REAL>()->default_value(3),
			"DoG filter outer kernel StDev")("Cont-Equal", po::bool_switch(),
			"Perform Contrast Equalization")("Cont-Alpha",
			po::value<REAL>()->default_value(0.1),
			"Contrast Equalization:Alpha value")("Cont-Tou",
			po::value<REAL>()->default_value(10),
			"Contrast Equalization:Threshold value in the range [0,255]");
}
void AddImageProcessOptions(po::options_description &im_options,
		ProcessImage *pim) {

	//	Channel ch = (Channel) obj["Color-Channel"].as<int> (); // which color channel to use..
	//	REAL gnorm = obj["Gamma-Normalization"].as<REAL> ();
	//	REAL alpha = obj["Cont-Alpha"].as<REAL> ();
	//	REAL tou = (REAL) obj["Cont-Tou"].as<int> () / 255;
	//	REAL k1 = obj["DoG-K1"].as<REAL> ();
	//	REAL k2 = obj["DoG-K2"].as<REAL> ();
	//	bool dog = obj.count("DoG");
	//	bool coneq = obj.count("Cont-Equal");

	Channel ch;
	bool dog, coneq;
	REAL gnorm, alpha, tou, k1, k2;
	// Image Processor
	im_options.add_options()("Color-Channel",
			po::value<Channel>(&ch)->default_value(RGB),
			"Image Color Channel to use "
					"0. Red,\n"
					"1. Green,\n"
					"2. Blue,\n"
					"3. Gray,\n"
					"4. RGB,\n")("Gamma-Normalization",
			po::value<REAL>(&gnorm)->default_value(1),
			"Perform Gamma Normalization (for LBP & LTP)")("DoG",
			po::bool_switch(&dog), "Apply DoG Filter")("DoG-K1",
			po::value<REAL>(&k1)->default_value(1),
			"DoG filter Inner kernel StDev")("DoG-K2",
			po::value<REAL>(&k2)->default_value(3),
			"DoG filter outer kernel StDev")("Cont-Equal",
			po::bool_switch(&coneq), "Perform Contrast Equalization")(
			"Cont-Alpha", po::value<REAL>(&alpha)->default_value(0.1),
			"Contrast Equalization:Alpha value")("Cont-Tou",
			po::value<REAL>(&tou)->default_value(10),
			"Contrast Equalization:Threshold value");

	pim = new ProcessImage(gnorm, dog, k1, k2, coneq, tou, alpha, ch);

}

int main(int ac, char* av[]) {

	// Set Program options
	po::options_description generic("Generic Options");
	generic.add_options()("help", "Produce help message");

	po::options_description input_options("Input/Output/Miscellaneous Options");
	po::options_description im_options("Image Processor (IP) Options");
	po::options_description lbp_options("Local Patterns (LBP/LTP) Options");
	po::options_description lqp_options("Local Patterns (LQP) Options");
	po::options_description cb_options("LQP CodeBook Options");

	LBPParams lbpparam;
	PatchParams lqpparam;
	CBParams cbparam;

	AddMiscOptions(input_options);
	AddImageProcessOptions(im_options);
	AddLBPOptions(lbp_options, lbpparam);
	AddLQPOptions(lqp_options, lqpparam);
	AddCBOptions(cb_options, cbparam);
	po::options_description all_options;
	all_options.add(generic).add(input_options).add(im_options).add(lbp_options).add(
			lqp_options).add(cb_options);
	//	po::variables_map vm;
	//	po::store(po::command_line_parser(argc, argv).options(all_options).run(),
	//			vm);
	//	po::notify(vm);

	try {
		po::variables_map vm;
		po::store(
				po::parse_command_line(ac, av, all_options,
						(po::command_line_style::default_style
								^ po::command_line_style::allow_guessing)), vm);

		po::notify(vm);

		if (vm.count("help")) {
			cout << all_options << "\n";
			return 1;
		} else
			Process pobj(vm, lbpparam, lqpparam, cbparam);
	} catch (exception& e) {
		cerr << "error: " << e.what() << "\n";
		return 1;
	} catch (...) {
		cerr << "Exception of unknown type!\n";
	}
	return 0;
}
