/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/
#include "includeHeaders.h"

#include "processImage.h"
#include <boost/program_options.hpp>
namespace po = boost::program_options;
UINT Features::xcell = 8, Features::ycell = 8, Features::norient = 36,
		Features::foffset = 8, Features::cellsize = 8;
//Pyramid Features::pyobj(64,128,10,Bilinear);
//PyramidType Features::sspace = NoPyramid;
#include <iterator>
class Process {
public:
	Process(po::variables_map &obj, HistParams &histparams) {
		ProcessInput(obj, histparams);
	}
	void ProcessInput(po::variables_map &obj, HistParams &histparams);
	~Process() {
		delete fobj;
		delete hogroot;
	}
private:
	string cname;
	Features *fobj;
	HOGRoot *hogroot;
	ProcessImage *pim;
};


void Process::ProcessInput(po::variables_map &obj, HistParams &histparams) {
	cname = obj["Class"].as<string> ();
	UINT height = obj["Win-Height"].as<int> (), width =
			obj["Win-Width"].as<int> (),
			cellsize = obj["Cell-Size"].as<int> (), radius = obj["Radius"].as<
					int> (), npoints = obj["Number-Points"].as<int> (),
			xstride = obj["X-Stride"].as<int> (), ystride = obj["Y-Stride"].as<
					int> (), nplevels = obj["Number-Levels"].as<int> (),
			nflevels = obj["Feature-Levels"].as<int> (), lbpstride =
					obj["LBP-Stride"].as<int> (), ftype =
					obj["FeatureType"].as<int> (), extent = obj["Extent"].as<
					int> (), norient = obj["NORIENT"].as<int> (), ltol =
					obj["LTP-TLevels"].as<int> (), plsdim = obj["PLS-Dim"].as<
					int> (), nbins = obj["LBP-NBins"].as<int> (), skip =
					obj["Skip"].as<int> (), detvalue =
					obj["Detector"].as<int> (), csize =
					obj["CacheSize"].as<int> ();
	plsdim = (plsdim > 0 ? plsdim : 0);
	csize = pow(2.0, 20.0) * csize;
	NORM norm = (NORM) obj["Normalization"].as<int> ();
	Channel ch = (Channel) obj["Color-Channel"].as<int> ();

	REAL svmj = obj["SVM-J"].as<REAL> (), svmc = obj["SVM-C"].as<REAL> (),
			gnorm = obj["Gamma-Normalization"].as<REAL> (), tolerance =
					((REAL) obj["LTP-Tolerance"].as<int> ()) / 255.0, alpha =
					obj["Cont-Alpha"].as<REAL> (), tou =
					(REAL) obj["Cont-Tou"].as<int> () / 255, k1 =
					obj["DoG-K1"].as<REAL> (), k2 = obj["DoG-K2"].as<REAL> ();

	bool dog = obj.count("DoG"), coneq = obj.count("Cont-Equal"), foldfeat =
			obj.count("Fold-Features"), embedSVM = (obj["Embed-SVM"].as<int> ()
			== 1 ? true : false), addrgb =
			(obj["Add-RGBLBP"].as<int> () == 1 ? true : false), nseparately =
			(obj["Norm-Sep"].as<int> () == 1 ? true : false), qda = (detvalue
			== 2 ? true : false);

	pim = new ProcessImage(gnorm, dog, k1, k2, coneq, tou, alpha, ch);
	cout << "\n Folding Parameter = " << foldfeat << endl;
	LBPParams lbpparam;
	lbpparam.cellsize = cellsize;
	lbpparam.stride = xstride;
	lbpparam.nbins = nbins;
	lbpparam.radius = radius;
	lbpparam.npoints = npoints;
	lbpparam.addfeat = addrgb;
	lbpparam.normfeatsep = nseparately;
	lbpparam.featnorm = norm;
	lbpparam.nflevels = nflevels;
	lbpparam.folded = foldfeat;
	lbpparam.gridtype = (GridType) obj["LBPGridType"].as<int> ();
	lbpparam.tolerance = tolerance;// used for the ltp features
	lbpparam.nthrlevels = 1;
	lbpparam.extent = extent;
	lbpparam.useblocknorm = obj.count("BlockNorm");
	lbpparam.issigned = obj.count("Signed");
	lbpparam.usegradmag = obj.count("UseGradMag");

	HOGParams hogparams;
	hogparams.cell = 8;
	hogparams.norient = norient;
	hogparams.chan = ch;

	if (obj.count("Neg-Train-File"))
		histparams.nfname = obj["Neg-Train-File"].as<string> ();
	else
		histparams.nfname = "";

	if (obj.count("Pos-Train-File"))
		histparams.pfname = obj["Pos-Train-File"].as<string> ();
	else
		histparams.pfname = "";

	int sparse = 0;
	if (obj.count("HF-Type")) {
		histparams.ftype = (WindowFeatureTypes) ftype;
		HistFeatures type = (HistFeatures) obj["HF-Type"].as<int> ();
		HistogramFeatures *histfeat;
		switch (type) {
		case HF_BIN:
			histfeat = new HistogramFeatures(lbpparam, hogparams, width,
					height, nplevels, *pim, histparams);
			histfeat->GenerateHistogramInformation();
			sparse = histparams.nbins; // use sparse svm learning mechanisms...
			break;
		case HF_BIN_TRAIN:
			histfeat = new HistogramFeaturesTrain(lbpparam, hogparams, width,
					height, nplevels, *pim, histparams);
			histfeat->GenerateHistogramInformation();
			break;
		case HF_LHRATIO:// liklihood ratio
		default:
			histfeat = new HistogramFeaturesLLRatio(lbpparam, hogparams, width,
					height, nplevels, *pim, histparams);
			histfeat->GenerateHistogramInformation();

		}
		fobj = histfeat;
	} else if (obj.count("Combined-Feature-Type") || obj.count(
			"Feature-ConfigFile")) {
		fobj = new MultipleWindowFeatures(lbpparam, width, height, nplevels,
				*pim, obj);

	} else {
		switch (ftype) {
		case WFT_LBP:
			fobj = new LBPFeatures(lbpparam, width, height, nplevels, *pim);
			break;
		case WFT_HOG:
			fobj = new HOGFeatures(width, height, nplevels);
			break;
		case WFT_HOG_LBP:
			fobj = new HOGLBPFeatures(cellsize, npoints, radius, width, height,
					nplevels, norm, nflevels, lbpstride, *pim);
			//		fobj =new ProjectedLBPFeatures(cellsize,npoints,radius,
			//				width,height,nplevels,norm,nflevels,lbpstride,*pim,foldfeat,
			//				norient);
			break;
		case WFT_HOG_COLOR:
			fobj = new HOGColorFeatures(width, height, nplevels);
			break;
		case WFT_LBP_LTP:
			fobj = new LTPFeatures(cellsize, npoints, radius, width, height,
					nplevels, norm, nflevels, lbpstride, *pim, foldfeat,
					tolerance);
			break;
		case WFT_HOG_LBP_LTP:
			fobj
					= new HOGLTPFeatures(cellsize, npoints, radius, width,
							height, nplevels, norm, nflevels, lbpstride, *pim,
							false, tolerance);
			break;
		case WFT_PROJ_HOG:
			fobj = new ProjectedFeatures(width, height, nplevels, norient, ch);
			break;
		case WFT_LTPX:
			fobj = new OnlyLTPFeaturesX(cellsize, npoints, radius, width,
					height, nplevels, norm, nflevels, lbpstride, *pim,
					foldfeat, tolerance, ltol);
			break;
		case WFT_LBP256:
			fobj = new LBPFeatures256(cellsize, npoints, radius, width, height,
					nplevels, norm, nflevels, lbpstride, *pim);
			break;

		case WFT_CSLTP:
			//			fobj = new CSLTPFeatures(cellsize, npoints, radius, width, height,
			//					nplevels, norm, nflevels, lbpstride, *pim, foldfeat,
			//					tolerance);
			cout << " Not Defined THis option CS_LTP";
			exit(EXIT_FAILURE);
			break;
		case WFT_LTP:
			fobj = new OnlyLTPFeatures(cellsize, npoints, radius, width,
					height, nplevels, norm, nflevels, lbpstride, *pim,
					foldfeat, tolerance);
			break;
		case WFT_PROJ_HOG_LBP_LTP_COLOR:
			fobj = new ProjectedLTPFeaturesRGB(cellsize, npoints, radius,
					width, height, nplevels, norm, nflevels, lbpstride, *pim,
					foldfeat, addrgb, nbins, nseparately, tolerance, norient,
					Circular);
			break;
		case WFT_PROJ_HOG_LBP_COLOR:
			fobj = new ProjectedLBPFeaturesRGB(cellsize, npoints, radius,
					width, height, nplevels, norm, nflevels, lbpstride, *pim,
					foldfeat, addrgb, nbins, nseparately, norient);
			break;
		case WFT_LBP_COLOR:
			if (obj.count("UseMeanValue"))
				fobj = new MeanLBPFeaturesRGB(lbpparam, width, height,
						nplevels, *pim);
			else
				fobj = new LBPFeaturesRGB(lbpparam, width, height, nplevels,
						*pim);
			break;
		case WFT_HOG_LBP_COLOR:
			fobj = new HOGLBPFeaturesRGB(lbpparam, width, height, nplevels,
					*pim, hogparams);
			break;
		case WFT_LTP_COLOR:
			fobj = new OnlyLTPFeaturesRGB(lbpparam, width, height, nplevels,
					*pim);
			break;
		case WFT_HOG_LTP_COLOR:
			fobj = new HOGOnlyLTPFeaturesRGB(lbpparam, width, height, nplevels,
					*pim, hogparams);
			break;
		case WFT_LBP_LTP_COLOR:
			fobj = new LTPFeaturesRGB(lbpparam, width, height, nplevels, *pim);
			break;

		case WFT_HOG_LBP_LTP_COLOR:
			fobj = new HOGLTPFeaturesRGB(lbpparam, width, height, nplevels,
					*pim, hogparams);
		case WFT_CSLBP:
			fobj = new CSLBPFeatures(lbpparam, width, height, nplevels, *pim);
			break;
		case WFT_CSLBP_COLOR:
			fobj
					= new CSLBPFeaturesRGB(lbpparam, width, height, nplevels,
							*pim);
			break;
		case WFT_CSLTP_COLOR:
			fobj
					= new CSLTPFeaturesRGB(lbpparam, width, height, nplevels,
							*pim);
			break;
		default:
			fobj = new LBPFeatures(lbpparam, width, height, nplevels, *pim);
		}
	}
#ifdef TMPDEBUG
	//	fobj->TestRoutine();
#endif
	hogroot = new HOGRoot(cname, width, height, xstride, ystride, extent,
			*fobj, plsdim, embedSVM, skip, qda, sparse, svmc);
	if (obj.count("Train")) {
		if (!obj.count("Pos-Train-File") || !obj.count("Neg-Train-File")) {
			cout
					<< "TRAIN : Positive or Negative Training List file is missing "
					<< endl << " So Aborting .... " << endl;
			return;
		}
		string plfile = obj["Pos-Train-File"].as<string> (), nlfile =
				obj["Neg-Train-File"].as<string> ();
		hogroot->Run(cname, plfile, nlfile, svmj, csize);

	} else if (obj.count("Plot")) {
		PLOTYPE ptype = (PLOTYPE) obj["Plot"].as<int> ();
		string lpimages, lnimages;
		if (ptype == 1) {
			if (!obj.count("Pos-Test-File") || !obj.count("Neg-Test-File")) {
				cout
						<< "DET Plot: Positive or Negative Test List fileM is missing "
						<< endl << " So Aborting .... " << endl;
				return;
			}
			lpimages = obj["Pos-Test-File"].as<string> ();
			lnimages = obj["Neg-Test-File"].as<string> ();
		}
		if (ptype == 2) {
			if (!obj.count("Test-File")) {
				cout << "PR Plot: Test Images List file is missing " << endl
						<< " So Aborting .... " << endl;
				return;
			}
			lpimages = obj["Test-File"].as<string> ();
		}
		string lmfile;
		if (!obj.count("Model-File"))
			lmfile = cname + ".mod";
		REAL rth = obj["ReduceThreshold"].as<REAL> ();
		hogroot->GeneratePlotData(cname, lpimages, lnimages, lmfile, ptype, rth);
	} else if (obj.count("Test")) {

	} else {
		cout << " Neither of Training, Plotting or Testing Option Selected  "
				<< endl << " So Aborting .... " << endl;
		return;
	}
}
int main(int ac, char* av[]) {
	try {

		po::options_description mainopt("Detector Options");
		mainopt.add_options()("help", "produce help message")("Class",
				po::value<string>()->default_value("person"), "Class Name")(
				"Train", "Train Detector")("Test", po::value<string>(),
				"Test Detector: Name of File")("Plot",
				po::value<int>()->default_value(1),
				"Generate Plot Data \n 1.DET Curve \n 2.PR Curve")("Detector",
				po::value<int>()->default_value(1),
				"Detector to be Used \n 1.Linear SVM \n 2.QDA Classifier")(
				"Color-Channel", po::value<int>()->default_value(4),
				"Image Color Channel to use "
					"\n 0. Red Channel \n"
					" 1. Green Channel \n"
					" 2. Blue Channel \n"
					" 3. Gray Scale \n"
					" 4. RGB (Default) \n"
					" 5. O1: luminance component...\n"
					" 6. O2 :Red-Green channel \n"
					" 7. O3 :Blue-Yellow channel "
					" 8. HSL ,\n"
					" 9. LAB ,\n"
					" 10. WandellOCS ,\n"
					" 11. SandeOCS ,\n"
					" 12. ABSGradientMag")(
				"FeatureType",
				po::value<int>()->default_value(1),
				"Features Type \n1. WFT_LBP,"
					"\n2. WFT_HOG,"
					"\n3. WFT_HOG_LBP,"
					"\n4. WFT_HOG_COLOR(HOG+Color Histogram),"
					"\n5. WFT_LBP_LTP,"
					"\n6. WFT_HOG_LBP_LTP,"
					"\n7. WFT_PROJ_HOG,"
					"\n8. WFT_LTPX,"
					"\n9. WFT_LBP256,"
					"\n10. WFT_CSLBP,"
					"\n11. WFT_CSLTP,"
					"\n12. WFT_LTP,"
					"\n13. WFT_PROJ_HOG_LBP_LTP_COLOR,  // BUG IN COMPUTATION OF DET, NO PROBLEM FOR PR,"
					"\n14. WFT_PROJ_HOG_LBP_COLOR	,"
					"\n15. WFT_LBP_COLOR,"
					"\n16. WFT_LTP_COLOR,"
					"\n17. WFT_LBP_LTP_COLOR,"
					"\n18. WFT_HOG_LBP_COLOR,"
					"\n19. WFT_HOG_LTP_COLOR,"
					"\n20. WFT_HOG_LBP_LTP_COLOR"
					"\n21. WFT_CSLBP_COLOR"
					"\n22. WFT_CSLTP_COLOR")("Add-RGBLBP",
				po::value<int>()->default_value(1),
				"1. Add RGB-Channels Hist\n 2. Concat RGB-Channels Hist")(
				"BlockNorm",
				"Use Neigbhour block normalization for LBP(by default false)")(
				"Norm-Sep", po::value<int>()->default_value(1),
				"1. Normalize Each RGB channel & then add"
					" \n 0.Add RGB channels and then Normalize")("Embed-SVM",
				po::value<int>()->default_value(1), "Include SVM training ...")(
				"Skip", po::value<int>()->default_value(0),
				"Skip Number of Pixels on Borders")("LBP-NBins",
				po::value<int>()->default_value(59),
				"Number of LBPBins can be either 256 or 59")("Fold-Features",
				"Valid for Feature-Type 1,5,7 & 8")("PLS-Dim",
				po::value<int>()->default_value(0), "PLS Projection Dimension")// features will be projected to that dimension
		("Cell-Size", po::value<int>()->default_value(16),
				"LBP Feature: Cell size")("Radius",
				po::value<int>()->default_value(1), "Radius of LBP Features")(
				"Number-Points", po::value<int>()->default_value(8),
				"Number of Points around each point")("Feature-Levels",
				po::value<int>()->default_value(1),
				"Number of Feature Levels(default = 1)")("LTP-Tolerance",
				po::value<int>()->default_value(5),
				"+,- Threshold for the LTP Features")("LTP-TLevels", po::value<
				int>()->default_value(1), "Number of Tolerance Levels for LTPs")(
				"NORIENT", po::value<int>()->default_value(36),
				"Number of Orientations for Projected HOG")(
				"Gamma-Normalization", po::value<REAL>()->default_value(1),
				"Gamma Normalization for LBP & LTP Features")("DoG",
				"Apply DoG Filter")("DoG-K1", po::value<REAL>()->default_value(
				1), "DoG filter Inner kernel Stdev")("DoG-K2",
				po::value<REAL>()->default_value(3),
				"DoG filter outer kernel Stdev")("Cont-Equal",
				"Perform Contrast Equalization")("Cont-Alpha",
				po::value<REAL>()->default_value(0.1),
				"Contrast Equalization:Alpha value")("Cont-Tou",
				po::value<int>()->default_value(10),
				"Contrast Equalization:Threshold value")("Normalization",
				po::value<int>()->default_value(L1sqrt),
				"Block Normalization:\n 0. L1 Norm  \n 1. L1Sqrt Norm"
					" \n 2. L2 Norm \n 3. LOG2 \n 4. LOG2 followed by L1-Norm")(
				"Number-Levels", po::value<int>()->default_value(10),
				"Number of Levels in each Pyramid Octave")("X-Stride",
				po::value<int>()->default_value(8), "X Stride in Pixels")(
				"Y-Stride", po::value<int>()->default_value(8),
				"Y Stride in Pixels")("LBP-Stride",
				po::value<int>()->default_value(16),
				"LBP Block Stride in Pixels")("Win-Width",
				po::value<int>()->default_value(64), "Detection Window Width")(
				"Win-Height", po::value<int>()->default_value(128),
				"Detection Window Height")("Extent",
				po::value<int>()->default_value(16), "Pos Cropped Image Border")(
				"SVM-J", po::value<REAL>()->default_value(3),
				"SVM J: Relative Weighting of +ve to -ve Examples")(
				"Pos-Train-File", po::value<string>(), "Pos Train:  List File")(
				"Neg-Train-File", po::value<string>(), "Neg Train: List File")(
				"Pos-Test-File", po::value<string>(),
				"Pos Test:  List File for DET Computation")("Neg-Test-File",
				po::value<string>(), "Neg Test:  List File for DET Computation")(
				"Test-File", po::value<string>(),
				"PR Curve: List File of Images")("Model-File",
				po::value<string>(), "Name of Learned Model file")("CacheSize",
				po::value<int>()->default_value(1024),
				"Cache Size for the 2nd Stage Training of the detector")(
				"SVM-C", po::value<REAL>()->default_value(1),
				"SVM-C Paramters(1 means select by Heuristic)")("LBPGridType",
				po::value<int>()->default_value(0),
				"0. Circular\n 1. Rectangular")("UseGradMag",
				"If Specified Use Gradient Magnitude for quantization in LBP/CSLBP/LTP")(
				"Signed", "Use Signed CSLBP/LTP")("UseMeanValue",
				"Compare with mean of patch instead of center pixel value")(
				"Combined-Feature-Type,F", po::value<vector<int> >(),
				"Specify Multiple features to have combined features from the Feature-Type")(
				"Feature-ConfigFile", po::value<vector<string> >(),
				"Feature Configuration File")("ReduceThreshold",
				po::value<REAL>()->default_value(0),
				"Reduce the Threshold by subtracting this value");

		HistParams histparams;

		po::options_description hfdesc("Histogram of Features:");
		hfdesc.add_options()("HF-Type", po::value<int>(),
				"Histogram Feature Options \n 1. BIN,"
					"\n 2.BIN_TRAIN,"
					"\n 3.LHRATIO,// liklihood ratio")("HF-Nbins", po::value<
				UINT>(&histparams.nbins)->default_value(30),
				"Histogram Features Number of Bins")(
				"HF-Sigma",
				po::value<REAL>(&histparams.sigma)->default_value(5),
				"Gaussian Sigma for Smoothing of Histograms(Only in LHRatio)(0 Means No Smoothing)")(
				"HF-Epsilon",
				po::value<REAL>(&histparams.lhratio_const)->default_value(1e-5),
				"Epsilon For likelihood ratio constant")("HF-Norm", po::value<
				int>()->default_value(0),
				"Same as Original Feature Normalization");
		po::options_description desc;

		desc.add(mainopt).add(hfdesc);
		po::variables_map vm;
		po::store(po::parse_command_line(ac, av, desc), vm);
		po::notify(vm);

		if (vm.count("help")) {
			cout << desc << "\n";
			return 1;
		} else {
			histparams.norm = (NORM) vm["HF-Norm"].as<int> ();
			Process pobj(vm, histparams);
		}
	} catch (exception& e) {
		cerr << "error: " << e.what() << "\n";
		return 1;
	} catch (...) {
		cerr << "Exception of unknown type!\n";
	}

	return 0;
}
