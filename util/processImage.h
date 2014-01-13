/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/

#ifndef PROCESSIMAGE_H_
#define PROCESSIMAGE_H_
#include "util.hpp"
// Process the Image with
// Perform a Set of Image Processing Operations on the input image...
// those include gamma correction followed by DoG filtering and in the
// end contrast equalization....
// with standard deviation of 1....

//@Warning: In Magick++ version 6.3.7 the default border adding method
// that is virtualPixelMethod type is EdgeType while in the recent
// version it is virtualPixelMethod by default is TileType

// contrast Equalization test ltp threshold deemed to set higher around can be set to 17.465
class ProcessImage {
public:
	ProcessImage(REAL gamma_, bool dog_ = false, REAL k1_ = 1, REAL k2_ = 3,
			bool coneq_ = false, REAL tou_ = 10/*255*/, REAL alpha_ = 0.1,
			Channel tche = Gray, bool constr_ = false) :
		gamma(gamma_), k1(k1_), k2(k2_), dog(dog_), coneq(coneq_), tou(tou_),
				alpha(alpha_) {
		ksize = 3;
		cout << endl << " ------------------------------------" << endl
				<< " Image Processing Parameters \n Gamma Normalization:  "
				<< gamma << endl << " DoG Filtering : "
				<< (dog == false ? "False" : "True") << " K1 = " << k1
				<< " K2 = " << k2 << " Kernel 1 Spatial Width " << ceil(
				2.0 * k1) * 2.0 + 1 << " Kernel 2 Spatial Width " << ceil(
				2.00 * k2) * 2 + 1 << endl << " Contrast Equalization : "
				<< (coneq == false ? "False" : "True") << " Tou = " << tou
				<< ": Alpha = " << alpha;
		channel = tche;
		uchstat = false;
		cout << ": Image Channel =  (" << channel << ")  ";
		minr = 0;
		maxr = 1;

		if (!dog && gamma == 1)
			constr = constr_;
		else
			constr = !coneq;
		if (constr)
			cout
					<< "\n-----------------------\n Using Contrast Stretching Procedure\n";
		cout << endl << " Using " << GetChannelName(channel) << " Channel \n";
		cout << endl << " ------------------------------------" << endl;
	}
	void Process(Image &image, vector<REAL>&);
	void ProcessGray(Image &image, vector<REAL>&);
	void Process(Image &image, vector<REAL>&, vector<REAL>&, vector<REAL>&);// Case RGB,OCS,HSL
	void ProcessRGB(Image &image, vector<REAL>&, vector<REAL>&, vector<REAL>&);// Case RGB
	void ProcessOCS(Image &image, vector<REAL>&, vector<REAL>&, vector<REAL>&);// Case RGB
	void ProcessSaveRes(Image &image, vector<REAL>&, vector<REAL>&,
			vector<REAL>&, const string&);// Case RGB
	void ProcessHSL(Image &image, vector<REAL>&, vector<REAL>&, vector<REAL>&);// Case HSL
	void ReadImage(Image& image, vector<REAL>&pixval);
	void ReadGrayImage(const Image& image, vector<REAL>&pixval);
	void ReadHSLImage(const Image& image, vector<REAL>&, vector<REAL>&,
			vector<REAL>&);
	void ReadImage(const Image& image, vector<REAL>&, vector<REAL>&,
			vector<REAL>&);
	void SubtractImage(const Image& image1, const Image& image2,
			vector<REAL>&pixval);
	void SubtractImage(const Image& image1, const Image& image2, vector<REAL>&,
			vector<REAL>&, vector<REAL>&);
	void SubtractImageHSL(const Image& image1, const Image& image2,
			vector<REAL>&, vector<REAL>&, vector<REAL>&);
	void DoContrastStretching(vector<REAL>&r, vector<REAL>&g, vector<REAL>&b,
			REAL imgstat[]);
	void DoContrastStretching(vector<REAL>&gray, REAL imgstat[]);
	void DoContrastEqualization(vector<REAL>&r, vector<REAL>&g, vector<REAL>&b,
			REAL imgstat[]);
	void DoContrastEqualization(vector<REAL>&gray, REAL imgstat[]);
	void ExtractChannelStats(vector<REAL>&r, vector<REAL>&g, vector<REAL>&b,
			REAL imgstat[]);
	void ExtractChannelStats(vector<REAL>&gray, REAL imgstat[]);
	REAL GetGamma() {
		return gamma;
	}
	void Gaussian(vector<REAL>& im, vector<REAL>&ker, UINT imwidth);
	void SetChannel(Channel che) {
		channel = che;
	}
	Channel GetChannel() {
		return channel;
	}
	static void Normalize(vector<REAL>& ch, REAL factor);
	static void Normalize(vector<REAL>& rch, vector<REAL>& gch,
			vector<REAL>& bch, REAL factor);

private:
	REAL gamma, k1, k2;
	UINT ksize;// width of kernels used for the DoG filtering...
	bool dog, coneq, constr;// perform the dog_filtering & contrast equalization or not....,if not contrast equalization
	bool uchstat;
	// then contrast stretch is used...
	REAL tou, alpha; // for contrast Equalization
	Channel channel;
	REAL minr, maxr;// min of range & max of output range...
};
#endif /* PROCESSIMAGE_H_ */
