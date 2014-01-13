/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/
//#include <boost/bind.hpp>
#include<vector>
#include<Magick++.h>
#include<dirent.h>
#include<iostream>
#include<string>
#include<cmath>
#include<cstdio>
#include<fstream>
#include <cstdlib>
#include<algorithm>
#include<numeric>
#include <sstream>
#include<set>
#include"Sorter.h"
#include <boost/thread.hpp>
#include <boost/program_options.hpp>
#ifdef WITH_BLITZ
#include <blitz/array.h>
#include <blitz/tinyvec.h>
#include <blitz/tinymat.h>
#include <blitz/tinyvec-et.h>
#endif
//#include "../Eigen/QR"
//using namespace Eigen;

//#include "../Eigen/Core"
//#include "../Eigen/Array"
//#define EIGEN2_SUPPORT
#define EIGEN_DEFAULT_TO_ROW_MAJOR // to use row major order for eigen class
#include "../Eigen/Core"
#include "../Eigen/Dense"
using namespace std;
using namespace Eigen;

#include "bips.h"

#ifndef UTIL_H_
#define UTIL_H_
using namespace Magick;
typedef unsigned int UINT;
typedef float REAL;
typedef unsigned char UCHAR;
#ifdef WITH_BLITZ
typedef blitz::Array<REAL, 2> Array2DReal;
#endif
enum FilterType {
	Gaussian, Bilinear
};
enum DistanceFunction {
	DF_NormalizedIntersection, DF_CHISquared, DF_L2,
};
#ifndef MIN
#define MIN(a,b) ( ( (a) < (b) ) ? (a) : (b) )
#endif

#ifndef MAX
#define MAX(a,b) ( ( (a) > (b) ) ? (a) : (b) )
#endif
// last three signifies the Opponent ColorSpace....
enum Channel {
	Red = 0,
	Green,
	Blue,
	Gray,
	RGB,
	O1,
	O2,
	O3,
	HSL,
	LAB,
	WandellOCS,
	SandeOCS,
	ABSGradientMag
};
// Where O1== luminance component...
// O2 is the red-green channel
// O3 is the blue-yellow channel....

enum ClassifierType {
	SVMLIGHT,
	L1LR,
	L2SVM_L1Reg,
	L1SVM_L2Reg,
	LatentSVM,
	PLSQDA,
	GaussianSVM,
	QuadEmbedSVM,
	LinearEmbedSVM,
	CubicEmbedSVM,
	IKSVM,
	LibLinearPhi2,
	L2LR,
	LinQuadEmbed,
	NormalIKSVM,
	PolyEmbedSVM,
	CT_RandomForests,
	CT_SMIDAS,
	LIBOCAS
};

#ifndef ABS
#define ABS(x) ( ( (x) < 0 ) ?  -(x)   : (x) )
#endif

static string GetClassifierType(const ClassifierType& ctype) {
	string val[] = { "SVMLIGHT", "L1LR", "L2SVM_L1Reg", "L1SVM_L2Reg",
			"LatentSVM", "PLSQDA", "GaussianSVM", "QuadEmbedSVM",
			"LinearEmbedSVM", "CubicEmbedSVM", "IKSVM", "LibLinearPhi2",
			"L2LR", "LinQuadEmbed", "NormalIKSVM", "PolyEmbedSVM",
			"CT_RandomForests" };
	return val[MIN(ctype,16)];
}

static string GetChannelName(Channel ctype) {
	const int nchn = 13;
	string chnames[] = { "Red", "Green", "Blue", "Gray", "RGB", "O1", "O2",
			"O3", "HSL", "LAB", "WandellOCS", "SandeOCS", "ABSGradientMag" };
	return chnames[MIN(nchn,ctype)];
}
void ReadDir(const string &dname, const string &, const string &ofname);
void ReadDir(const string &dname, const string &ofname);
void ReadDir(const string dname, vector<string> &lname);
void ReadImage(unsigned char *mat, const Image & image, unsigned int offset);
void ReadImage(const Image &image, vector<REAL> &pixval, REAL = 1);
void ReadColorImage(const Image&image, vector<REAL>&pixval, int = 1);
void Scale(Image &simage, UINT rows, UINT cols, bool aspect);
void Scale(Image &simage, REAL scale);
string int2str(int input);
void CopyBoundaries(Image& input, Image& output, char dir);
void Padding(Image& input, Image& output, UINT size, char dir);

void ResizeImage(Image & image, UINT width, UINT height, FilterType = Bilinear,
		bool = false);

void Reflect(Image& input, Image& output, UINT size, char dir);
void ReflectImage(Image& output, UINT size, char dir);
void BilinearInterpolation(vector<float>& input, int cols, int rows,
		float scale, vector<float>& output);
void BilinearInterpolation(vector<REAL>& input, UINT cols, UINT rows, UINT dim,
		REAL scale, vector<REAL>& output);
void Compare(const string& nfile1, const string& nfile2, vector<string>&diff);
void ReadFile(const string &nfile, vector<string>&vstr);

// Modified version of BinarySearch gives the closest element's
// index to the element being searched...
//
//static int BinarySearch(vector<BoostInfo> &data, int size, REAL searchElement) {
//	//	int size = data.size();
//	int low = 0; // low end of the search area
//	int high = size - 1; // high end of the search area
//	int middle = (low + high + 1) / 2; // middle element
//	int location = -1; // return value; -1 if not found
//
//	if (searchElement >= data[high].weight)
//		location = high;
//	else if (searchElement <= data[low].weight)
//		location = low;
//	else {
//		do // loop to search for element
//		{
//			// print remaining elements of vector to be searched
//			//				displaySubElements(data, low, high);
//
//			// output spaces for alignment
//			for (int i = 0; i < middle; i++)
//				cout << "   ";
//
//			//				cout << " * " << endl; // indicate current middle
//
//			// if the element is found at the middle
//			if (searchElement == data[middle].weight)
//				location = middle; // location is the current middle
//			else if (searchElement < data[middle].weight) // middle is too high
//				high = middle - 1; // eliminate the higher half
//			else
//				// middle element is too low
//				low = middle + 1; // eliminate the lower half
//
//			middle = (low + high + 1) / 2; // recalculate the middle
//		} while ((low <= high) && (location == -1));
//
//		if (location == -1) // no elements in the list, so return the next element in the ordered list.
//			location = MIN(low, size - 1);
//	}
//
//	return location; // return location of search key
//}

// Comparison Functors...
struct cmpInts {
	bool operator ()(int i1, int i2) {
		return i1 < i2;
	}
};
template<class T>
struct Cmp {
	bool operator ()(const T &i1, const T &i2) {
		return i1 < i2;
	}
};

template<class T>
struct CmpAbsDescend {
	bool operator ()(const T &i1, const T &i2) {
		return ABS(i1) > ABS(i2);
	}
};
template<class T> // utility function to display
ostream& operator<<(ostream& os, const vector<T>& v) {
	copy(v.begin(), v.end(), ostream_iterator<T> (cout, " "));
	return os;
}
template<class T>
void Print(T *start, T *end) {
	cout << endl;
	for (T *p = start; p < end; ++p)
		cout << *p << " ";
	cout << endl;
}
template<class T> // utility function to display
ofstream& operator<<(ofstream& os, const vector<T>& v) {
	copy(v.begin(), v.end(), ostream_iterator<T> (os, " "));
	return os;
}
template<class T>
void Print(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &mMatrix) {

	for (UINT i = 0; i < mMatrix.rows(); ++i) {
		for (UINT j = 0; j < mMatrix.cols(); ++j)
			cout << mMatrix(i, j) << " ";
		cout << endl;
	}

}
template<class T>
void WriteMatrix2File(const string &fname,
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &mMatrix) {

	ofstream ofile(fname.c_str(), ios::out);
	if (!ofile) {
		cout << " Couldn't Create the File " << fname << endl;
		return;
	}
	for (UINT i = 0; i < mMatrix.rows(); ++i) {
		for (UINT j = 0; j < mMatrix.cols(); ++j)
			ofile << mMatrix(i, j) << " ";
		ofile << endl;
	}
	ofile.close();
}
template<class T>
void WriteMatrix2File(const string &fname, Eigen::Map<T> &mMatrix) {

	ofstream ofile(fname.c_str(), ios::out);
	if (!ofile) {
		cout << " Couldn't Create the File " << fname << endl;
		return;
	}
	for (UINT i = 0; i < mMatrix.rows(); ++i) {
		for (UINT j = 0; j < mMatrix.cols(); ++j)
			ofile << mMatrix(i, j) << " ";
		ofile << endl;
	}
	ofile.close();
}
template<class T>
void Print(Eigen::Map<T> mMatrix) {

	for (UINT i = 0; i < mMatrix.rows(); ++i) {
		for (UINT j = 0; j < mMatrix.cols(); ++j)
			cout << mMatrix(i, j) << " ";
		cout << endl;
	}

}
template<class T>
void RandomPermute(vector<T> &array) {
	int n = array.size(); // The number of items left to shuffle (loop invariant).
	while (n > 1) {
		int k = rand() % n; // 0 <= k < n.
		n--; // n is now the last pertinent index;
		T temp = array[n]; // swap array[n] with array[k] (does nothing if k == n).
		array[n] = array[k];
		array[k] = temp;
	}
}
template<class T>
void RandomPermute(vector<T> &array, UINT limit)
// Samples the indeces from a given limit...
{
	int n;
	if (limit >= array.size()) {
		n = array.size();
	} else
		n = limit; // The number of items left to shuffle (loop invariant).
	while (n > 1) {
		int k = rand() % n; // 0 <= k < n.
		n--; // n is now the last pertinent index;
		T temp = array[n]; // swap array[n] with array[k] (does nothing if k == n).
		array[n] = array[k];
		array[k] = temp;
	}
}
void DrawRectangle(Image &image, int startx, int starty, int width, int height,
		const string &col);

void CropImage(Image &image, int xmin, int ymin, int width, int height);
void CropImage2(Image &image, int xmin, int ymin, int xmax, int ymax);
void PrintErrorExit(string msg);
void ParseListFile(const string&, vector<string>&);
void Dec2Bin(UINT, string&);
string SplitFilenameFromDir(const string& str);

//template <class T>
//UINT Bin2Dec(T *ptr,const UINT &size)
//{
//	UINT val = 0, count = 0;
//	for (UINT k = size;k >= 0; --k, ++count)
//		val += (UINT) pow(2.0, (double) count) * ptr[k];
//	return val;
//
//}
UINT Bin2Dec(const string&);
//template <class T>
//T FindMaximum(vector<T>& ivec,UINT &index);
//{
//	index=0;
//	T maxval=ivec[0];
//	UINT k=0;
////	for(vector<T>::iterator iter=ivec.begin();;)
////	for(vector<T>::iterator iter=ivec.begin(); iter!=ivec.end();++iter,++k)
//	for(int k=0; k < ivec.size();++k)
//		if( ivec[k] > maxval)
//		{
//			index = k;
//			maxval = ivec[k];
//		}
//	return maxval;
//}
// Read Matlab Matrix
#ifdef WITH_BLITZ
Array2DReal ReadData(const string &fname);
void Write2File(const string & ofname, const Array2DReal & flowu);
#endif
string SplitFilename(const string& str);
//{
//	std::ofstream ofile(ofname.c_str());
//	blitz::TinyVector<int, 2> extent = flowu.extent();
//	for (int i = 0; i < extent[0]; ++i) {
//		for (int j = 0; j < extent[1]; ++j) {
//			ofile << flowu(i, j) << " ";
//		}
//		ofile << endl;
//	}
//	ofile.close();
//}
unsigned long getActiveMemory();
UINT GetNumberLines(const string & fname);
void fwrite_u8_pnm(const char* file, int nx, int ny, int nc,
		const u8 im[/*0:nx-1@dx,0:ny-1@dy,0:nc-1@dc*/], int dx, int dy, int dc);
void CopyConstBoundaries(Image& input, UINT size);
template<class T>
void WriteVector2File(const string& fname, T * feat, int length) {
	ofstream ofile(fname.c_str(), ios::out);
	if (!ofile) {
		cout << " Couldn't Create the File " << fname << endl;
		return;
	}
	for (int i = 0; i < length; ++i)
		ofile << feat[i] << " ";
	ofile << endl;
	ofile.close();
}
template<class T>
void WriteVector2File(ofstream &ofile, T * feat, int length) {
	for (int i = 0; i < length; ++i)
		ofile << feat[i] << " ";
	ofile << endl;
	//	ofile.close();
}
template<class T>
void WriteVector2File(const string &fname, vector<T> &feat) {
	int length = feat.size();
	ofstream ofile(fname.c_str(), ios::out | ios::binary);
	if (!ofile) {
		cout << " Couldn't Create the File " << fname << endl;
		return;
	}
	ofile.write((char*) &length, sizeof(UINT));
	ofile.write((char*) (&feat[0]), sizeof(T) * length);
	ofile.close();
}
void WriteStrings2File(const string &fname, vector<string> &feat);
template<class T>
void WriteMatrix2File(const string& fname, T **feat, const UINT &rows,
		const UINT &cols) {
	ofstream ofile(fname.c_str(), ios::out);
	if (!ofile) {
		cout << " Couldn't Create the File " << fname << endl;
		return;
	}
	for (UINT i = 0; i < rows; ++i) {
		for (UINT j = 0; j < cols; ++j)
			ofile << feat[i][j] << " ";
		ofile << endl;
	}

	ofile.close();
}
template<class T>
void WriteImageFile(const string& fname, vector<T>& r, vector<T>&g,
		vector<T>&b, UINT rows, UINT cols) {
	ofstream ofile(fname.c_str(), ios::out);
	if (!ofile) {
		cout << " Couldn't Create the File " << fname << endl;
		return;
	}
	ofile.write((char*) &rows, sizeof(UINT));
	ofile.write((char*) &cols, sizeof(UINT));
	ofile.write((char*) &r[0], sizeof(T) * (rows * cols));
	ofile.write((char*) &g[0], sizeof(T) * (rows * cols));
	ofile.write((char*) &b[0], sizeof(T) * (rows * cols));
	ofile.close();
}
template<class T>
void WriteMatrix2File(const string& fname, vector<vector<T> >&feat) {
	ofstream ofile(fname.c_str(), ios::out);
	if (!ofile) {
		cout << " Couldn't Create the File " << fname << endl;
		return;
	}
	for (UINT i = 0; i < feat.size(); ++i) {
		for (UINT j = 0; j < feat[i].size(); ++j)
			ofile << feat[i][j] << " ";
		ofile << endl;
	}

	ofile.close();
}
//void WriteRealMatrix2BinFile(const string& fname, vector<vector<REAL> >&feat);
//void ReadRealMatrixFromBinFile(const string& fname, vector<vector<REAL> >&feat);
//void WriteRealMatrix2BinFile(const string& fname, vector<vector<REAL> >&feat,
//		const vector<string>&);
//void ReadRealMatrixFromBinFile(const string& fname, vector<vector<REAL> >&feat,
//		vector<string>&);

template<class T>
void WriteRealMatrix2BinFile(const string& fname, vector<vector<T> >&feat) {
	ofstream ofile(fname.c_str(), ios::out | ios::binary);
	if (!ofile) {
		cout << " Couldn't Create the File " << fname << endl;
		return;
	}
	UINT nfeat = feat.size(), featdim = feat[0].size();
	ofile.write((char*) &nfeat, sizeof(UINT));
	ofile.write((char*) &featdim, sizeof(UINT));
	for (UINT i = 0; i < feat.size(); ++i) {
		ofile.write((char*) (&feat[i][0]), sizeof(T) * featdim);
	}
	ofile.close();
}
template<class T>
void WriteRealMatrix2BinFile(const string& fname, vector<vector<T> >&feat,
		const vector<string>&labels) {
	assert(feat.size() == labels.size());
	ofstream ofile(fname.c_str(), ios::out | ios::binary);
	if (!ofile) {
		cout << " Couldn't Create the File " << fname << endl;
		return;
	}
	UINT nfeat = feat.size(), featdim = feat[0].size();
	ofile.write((char*) &nfeat, sizeof(UINT));
	ofile.write((char*) &featdim, sizeof(UINT));
	for (UINT i = 0; i < feat.size(); ++i) {
		UINT tvar = labels[i].length();
		ofile.write((char*) &tvar, sizeof(UINT));
		ofile.write((char*) labels[i].c_str(), sizeof(char) * tvar);
		ofile.write((char*) (&feat[i][0]), sizeof(T) * featdim);
	}
	ofile.close();
}
template<class T>
void ReadRealMatrixFromBinFile(const string& fname, vector<vector<T> >&feat) {
	ifstream ifile(fname.c_str(), ios::in | ios::binary);
	if (!ifile) {
		cout << " Couldn't Create the File " << fname << endl;
		return;
	}
	UINT nfeat, featdim;
	ifile.read((char*) &nfeat, sizeof(UINT));
	ifile.read((char*) &featdim, sizeof(UINT));
	feat.resize(nfeat, vector<T> (featdim, 0));
	for (UINT i = 0; i < feat.size(); ++i) {
		ifile.read((char*) (&feat[i][0]), sizeof(T) * featdim);
	}

	ifile.close();
}
template<class T>
void ReadRealMatrixFromBinFile(const string& fname, vector<vector<T> >&feat,
		vector<string>&labels) {
	ifstream ifile(fname.c_str(), ios::in | ios::binary);
	if (!ifile) {
		cout << " Couldn't Create the File " << fname << endl;
		return;
	}
	UINT nfeat, featdim;
	ifile.read((char*) &nfeat, sizeof(UINT));
	ifile.read((char*) &featdim, sizeof(UINT));
	feat.resize(nfeat, vector<T> (featdim, 0));
	char lab[1024];
	for (UINT i = 0; i < feat.size(); ++i) {
		UINT tvar;
		ifile.read((char*) &tvar, sizeof(UINT));
		ifile.read((char*) lab, sizeof(char) * tvar);
		labels.push_back(lab);
		ifile.read((char*) (&feat[i][0]), sizeof(T) * featdim);
	}
	ifile.close();
}

template<class T>
void WriteMatrix2File(const string& fname, T *feat, const UINT &rows,
		const UINT &cols) {
	ofstream ofile(fname.c_str(), ios::out);
	if (!ofile) {
		cout << " Couldn't Create the File " << fname << endl;
		return;
	}
	T *ptr = feat;
	for (UINT i = 0; i < rows; ++i) {
		for (UINT j = 0; j < cols; ++j)
			ofile << ptr[j] << " ";
		ofile << endl;
		ptr += cols;
	}

	ofile.close();
}
template<class T>
T Dot(vector<T>& feat) {
	T sum = 0;
	for (typename vector<T>::iterator iter = feat.begin(); iter != feat.end(); ++iter)
		sum += *iter * *iter;
	return sum;
}
template<class T> /*Can Use inner_product instead*/
T Dot(vector<T>& feat, vector<T>&filter) {
	T sum = 0;
	for (typename vector<T>::iterator iter = feat.begin(), fiter =
			filter.begin(); iter != feat.end(); ++iter, ++fiter)
		sum += *iter * *fiter;
	return sum;
}
template<class T> /*Can Use inner_product instead*/
T Dot(vector<T>& feat, T *filter) {
	T sum = 0;
	T *ptr = filter;
	for (typename vector<T>::iterator iter = feat.begin(); iter != feat.end(); ++iter, ++ptr)
		sum += *iter * *ptr;
	return sum;
}
template<class T>
void ReadTextFile(const string &fname, vector<vector<T> >&vec, UINT &nexamples,
		UINT &dimfeature) {
	ifstream ifile(fname.c_str(), ios::in);
	if (!ifile) {
		cout << endl << " Error Couldn't Load the file " << fname << endl;
		exit(EXIT_FAILURE);
	}
	ifile >> nexamples >> dimfeature;

	vec.resize(nexamples, vector<T> (dimfeature, 0));

	for (UINT i = 0; i < nexamples; ++i) {

		for (UINT j = 0; j < dimfeature && ifile; ++j)
			ifile >> vec[i][j];

	}
}
template<class T>
void ReadTextFile(const string &fname, T **&vec, UINT &nexamples,
		UINT &dimfeature) {
	ifstream ifile(fname.c_str(), ios::in);
	if (!ifile) {
		cout << endl << " Error Couldn't Load the file " << fname << endl;
		exit(EXIT_FAILURE);
	}
	ifile >> nexamples >> dimfeature;

	vec = new T*[nexamples];

	for (UINT i = 0; i < nexamples; ++i) {
		vec[i] = new T[dimfeature];

		for (UINT j = 0; j < dimfeature && ifile; ++j)
			ifile >> vec[i][j];
	}

}
template<class T>
void ReadTextFile(const string &fname, T *&vec, UINT &nexamples,
		UINT &dimfeature) {
	ifstream ifile(fname.c_str(), ios::in);
	if (!ifile) {
		cout << endl << " Error Couldn't Load the file " << fname << endl;
		exit(EXIT_FAILURE);
	}
	ifile >> nexamples >> dimfeature;

	vec = new T[nexamples * dimfeature];
	REAL *tvec = vec;
	for (UINT i = 0; i < nexamples * dimfeature; ++i, ++tvec) {
		ifile >> *tvec;
	}
	ifile.close();
}
template<class T>
void ReadFile(const string &fname, T **&vec, UINT &nexamples, UINT &dimfeature) {
	ifstream ifile(fname.c_str(), ios::in | ios::binary);
	if (!ifile) {
		cout << endl << " Error Couldn't Load the file " << fname << endl;
		exit(EXIT_FAILURE);
	}
	ifile.read((char*) &nexamples, sizeof(UINT));
	ifile.read((char*) &dimfeature, sizeof(UINT));
	cout << "\n File  " << fname << " has " << nexamples
			<< " Examples  with Dimension =" << dimfeature << endl;

	vec = new T*[nexamples];

	for (UINT i = 0; i < nexamples; ++i) {
		vec[i] = new T[dimfeature];

		ifile.read((char*) vec[i], sizeof(REAL) * dimfeature);
	}
	ifile.close();
}
template<class T>
void ReadFile(const string &fname, vector<vector<T> > &vec) {
	ifstream ifile(fname.c_str(), ios::in | ios::binary);
	if (!ifile) {
		cout << endl << " Error Couldn't Load the file " << fname << endl;
		exit(EXIT_FAILURE);
	}
	UINT nexamples, dimfeature;
	ifile.read((char*) &nexamples, sizeof(UINT));
	ifile.read((char*) &dimfeature, sizeof(UINT));
	cout << "\n File  " << fname << " has " << nexamples
			<< " Examples  with Dimension =" << dimfeature << endl;

	vec.resize(nexamples, vector<T> (dimfeature, 0));

	for (UINT i = 0; i < nexamples; ++i) {
		T *vecptr = &vec[i][0];
		ifile.read((char*) vecptr, sizeof(REAL) * dimfeature);
	}
	ifile.close();
}
void FindCrossValidationIndeces(UINT kfold, UINT tnexamples,
		vector<vector<UINT> >& trainidx, vector<vector<UINT> >& testidx);

template<class T>
UINT Bin2Dec(T *ptr, const UINT &size) {
	UINT val = 0, count = size - 1;
	for (UINT k = 0; k < size; ++k, --count)
		val += (UINT) pow(2.0, (double) count) * ptr[k];
	return val;
}
template<class Template>
void WriteStorage(ofstream &ofile, Template &tem) {
	long m = tem.rows(), n = tem.cols();
	ofile.write((char*) &m, sizeof(long));
	ofile.write((char*) &n, sizeof(long));
	double val;
	for (long i = 0; i < m; ++i)
		for (long j = 0; j < n; ++j) {
			val = tem(i, j);
			ofile.write((char*) &val, sizeof(double));
		}

}
template<class Template>
void ReadStorage(ifstream &ifile, Template &tem) {
	long m, n;
	ifile.read((char*) &m, sizeof(long));
	ifile.read((char*) &n, sizeof(long));
	if (m > 0 && n > 0) {
		tem = Template::Zero(m, n);
		double val;
		for (long i = 0; i < m; ++i)
			for (long j = 0; j < n; ++j) {
				ifile.read((char*) &val, sizeof(double));
				tem(i, j) = val;
			}
	}
}
template<class T1, class T2>
map<T1, T2> AddMaps(map<T1, T2> & m1, map<T1, T2> &m2) {
	map<T1, T2> tobj;
	for (typename map<T1, T2>::iterator p = m1.begin(); (p != m1.end()); ++p)
		tobj[p->first] = p->second;

	for (typename map<T1, T2>::iterator p = m2.begin(); (p != m2.end()); ++p)
		tobj[p->first] = tobj[p->first] + p->second;

	return tobj;
}
template<class T1, class T2>
void DisplayMaps(const map<T1, T2> & m1) {

	cout << "\n-------------------------------------------\n";
	for (typename map<T1, T2>::const_iterator p = m1.begin(); (p != m1.end()); ++p)
		cout << "Key = " << p->first << ", Value =" << p->second << ";" << endl;

}
bool CanReadFile(const string &fname);
//////////////////////////////////////////////////////////////////////////////
//
// process_mem_usage(double &, double &) - takes two doubles by reference,
// attempts to read the system-dependent data for a process' virtual memory
// size and resident set size, and return the results in KB.
//
// On failure, returns 0.0, 0.0

void ProcessMemoryUsage(double& vm_usage, double& resident_set);
void Split(int *start, int *end, int *&ls, int *&le, int *&rs, int *&re,
		int threshold);
void DisplayMemoryUsage();
template<class T, class C>
REAL L2Dist(vector<T> & vec1, vector<C> & vec2) {

	REAL dist = 0;
	typename vector<T>::iterator v1 = vec1.begin();
	typename vector<C>::iterator v2 = vec2.begin();
	for (; v1 != vec1.end(); ++v1, ++v2) {
		REAL d = (*v1 - *v2);
		dist += d * d;
	}
	return sqrt(dist);
}
template<class T, class C>
REAL L2Dist(T * vec1, C*vec2, UINT size) {

	REAL dist = 0;
	T *v1 = vec1;
	C *v2 = vec2;
	for (; v1 != vec1 + size; ++v1, ++v2) {
		REAL d = (*v1 - *v2);
		dist += d * d;
	}
	return sqrt(dist);
}
void MyGaussianFilter(REAL sigma, vector<REAL> &filter);
void PadImage(Image& out, UINT dim[]);
void GetCrossValidationIndeces(UINT kfold, UINT npos, UINT nneg,
		vector<vector<UINT> >&trainidx, vector<vector<UINT> >&testidx);
void PrintIdx(vector<UINT>&trainidx, vector<UINT>&testidx,
		UINT ptrainoffset = 0, UINT ptestoffset = 0);
template<class T>
void Convolve(vector<T> & data, vector<T>&filter) {
	UINT hfsize = filter.size() / 2;
	if (data.size() < filter.size())
		return;
	UINT dsize = data.size();
	vector<REAL> response(dsize + 2 * hfsize, 0);
	copy(data.begin(), data.end(), response.begin() + hfsize);
	// central part
	fill(data.begin(), data.end(), 0);
	for (UINT i = hfsize, count = 0; i < response.size() - hfsize; ++i, count++) {
		UINT j = count; // i-hfsize
		for (typename vector<T>::iterator iter = filter.begin(); iter
				!= filter.end(); ++iter, ++j) {
			data[count] += *iter * response[j];
		}
	}

}

template<class T>
void DisplayVector(vector<T>& data) {
	cout << endl;
	UINT count = 0;
	for (typename vector<T>::iterator iter = data.begin(); iter != data.end(); ++iter) {
		cout << *iter << " , ";
		++count;
		if (count % 10 == 0)
			cout << endl;
	}
	cout << endl;
}
template<typename T>
std::string ToString(const T& x) {
	std::ostringstream ss;
	ss << x;
	return ss.str();
}
string strtokenize(const string& str, char ch, bool linstance = false); /// can tokenize based on the last instance of the character, otherwise first instance is used.
void AppendPath(const string &pathdir, vector<string>&fnames);
template<class T>
T FindDistance(vector<T>& feat1, vector<T>& feat2,
		DistanceFunction dftype = DF_L2) {
	double distance = 0;
	double eps = 2.2204e-16;
	assert(feat1.size() == feat2.size());
	typename vector<T>::iterator v1 = feat1.begin();
	typename vector<T>::iterator v2 = feat2.begin();
	switch (dftype) {
	case DF_CHISquared:
		//		for (int i = 0; i < feat1.size(); ++i)
		//			distance += ((feat1[i] - feat2[i]) * (feat1[i] - feat2[i]))
		//					/ (feat1[i] + feat2[i] + 2.2204e-16);

		for (; v1 != feat1.end(); ++v1, ++v2) {
			distance += ((*v1 - *v2) * (*v1 - *v2)) / ((*v1 + *v2) + eps);
			//			distance += d * d;
		}
		distance /= 2;
		break;
	case DF_L2:
		for (; v1 != feat1.end(); ++v1, ++v2) {
			distance += (*v1 - *v2) * (*v1 - *v2);
		}
		// having sqrt doesn't have any impact
		break;
	case DF_NormalizedIntersection:
	default:
		for (; v1 != feat1.end(); ++v1, ++v2)
			distance += MIN( *v1 , *v2);
		distance = -distance;// because it measures the proximity
		break;
	}
	//		cout << endl << feat1 << endl;
	//		cout << feat2 << endl;
	//		cout << distance << endl;
	return distance;
}
template<class T>
T FindDistance(vector<T>& feat1, vector<T>& feat2, DistanceFunction dftype,
		UINT cellsize, vector<REAL>&dist) {
	double distance = 0;
	double eps = 2.2204e-16;
	assert(feat1.size() == feat2.size());
	typename vector<T>::iterator v1 = feat1.begin();
	typename vector<T>::iterator v2 = feat2.begin();
	switch (dftype) {
	case DF_CHISquared:
		for (; v1 != feat1.end(); ++v1, ++v2) {
			distance += ((*v1 - *v2) * (*v1 - *v2)) / ((*v1 + *v2) + eps);
		}
		distance /= 2;
		break;
	case DF_L2:
		for (; v1 != feat1.end(); ++v1, ++v2) {
			distance += (*v1 - *v2) * (*v1 - *v2);
		}
		// its is L2 Square distance
		break;
	case DF_NormalizedIntersection:
	default:
		for (; v1 != feat1.end(); ++v1, ++v2)
			distance += MIN( *v1 , *v2);
		distance = -distance;// because it measures the proximity
		break;
	}
	//		cout << endl << feat1 << endl;
	//		cout << feat2 << endl;
	//		cout << distance << endl;
	return distance;
}
// program_option validator...
void validate(boost::any& v, const std::vector<std::string>& values,
		Channel* target_type, int);
#endif /*UTIL_H_*/
