/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/
#ifndef FILEIO_H_
#define FILEIO_H_
// The header file used to read & write the different class informations
// Two different File Types 1). f: feature file type 2). l: SVM Light binary format file.
enum FTYPE {
	FeatureFile, SVMLight
};
#include "util.hpp"
#include<boost/filesystem/operations.hpp> // to use in merge function
namespace bf = boost::filesystem;
class FileIO {
public:

	FileIO(const string &, char, UINT, long = 0, FTYPE = FeatureFile);
	FileIO(const string &, char, FTYPE = FeatureFile);
	FileIO() {
		datapos = 0;
		dimfeature = 0;
		nexamples = 0;
	}
	~FileIO();
	void Initialize(const string&, char, FTYPE = FeatureFile);
	void OpenRead(const string&, FTYPE = FeatureFile);
	void OpenWrite(const string&, FTYPE = FeatureFile);
	void OpenWrite(const string&, UINT, long = 0);
	void WriteHeader();
	void WriteSVMLightHeader();
	void ReadSVMLightHeader();
	void ReadHeader();
	void WriteNExamples(long);
	virtual bool WriteFeature(float*);
	virtual bool WriteFeature(vector<REAL>&);
	//	template<class T>
	//	bool WriteFeature(Matrix<T, Dynamic, 1> &feat);
	//	template<class T>
	//	bool WriteFeature(Matrix<T, Dynamic, 1>&feat, int target);
	template<class T>
	bool WriteFeature(Matrix<T, Dynamic, 1> &feat) {

		for (int i = 0; i < dimfeature; ++i) {
			file.write((char*) &feat[i], sizeof(float));
			//		cout<<feat[i]<<" ";
		}
		//	file.write((char*)&feat[0],sizeof(float)*dimfeature);
		//	cout<<endl;
		if (file)
			return true;

		return false;
	}
	template<class T>
	bool WriteFeature(Matrix<T, Dynamic, 1>&feat, int target) {
		file.write((char*) &target, sizeof(int));
		for (int i = 0; i < dimfeature; ++i)
			file.write((char*) &feat[i], sizeof(float));

		if (file)
			return true;

		return false;
	}
	//	bool WriteFeature(VectorXd&);
	bool WriteFeature(float*, int);
	bool WriteFeature(vector<REAL>&, int);
	virtual bool ReadFeature(float*);
	bool ReadSVMFeature(float *feat, int& target);
	bool ReadFeature(VectorXf &);
	bool ReadSVMFeature(VectorXf& feat, int& target);
	static UINT MergeFeatureFiles(const string &, const string&);
	long GetNExamples() {
		return nexamples;
	}
	UINT GetDimFeature() {
		return dimfeature;
	}
	static void Convert2SVMLight(const string&, const string&, const string&);
	static void Convert2SVMLight(const string &pfname, const string &nfname,
			const string&, const string &nofile);
	static void Convert2SVMLight(const string &pfname,
			const vector<string>& nfnames, UINT, const string &nofile);
	static void WriteFile(const string&fname, UINT dfeature, REAL **pfeat,
			UINT npos, REAL **nfeat, UINT nneg);
	static void PermuteFeatures(const string&, const string&, FTYPE);
	static void Convert2SparseFormat(const string &pfname);
	static void ReadFile(const string&nffile, vector<vector<REAL> >& featmat,
			vector<int> &labels, FTYPE ftype = SVMLight);

	void GetClassCount(UINT &poscount, UINT &negcount);
	void CloseFile() {
		if (file.is_open())
			file.close();
	}
	bool IsOpen() {
		return file.is_open();
	}
	void ResetToDataPart() {
		if (IsOpen())
			file.seekg(datapos, ios::beg);
	}
	static void WriteFeatureInSVMFormat(const string&fname,
			vector<vector<REAL> >&features, vector<int>&target);
//	static void WriteFile(vector<vector<REAL> >& featmat, vector<int> &labels,
//			const string& nofile);
protected:
	UINT dimfeature;
	long nexamples;
	fstream file;
	bool append;
	FTYPE ftype;
	streampos datapos;
};
class SparseFileIO: public FileIO {
public:
	//	FileIO(const string &, char, UINT, long = 0, FTYPE = FeatureFile);
	//	FileIO(const string &, char, FTYPE = FeatureFile)
	SparseFileIO(const string &fname, char fun, UINT dfeature, long nexam = 0,
			FTYPE fftype = FeatureFile) :
		FileIO(fname, fun, dfeature, nexam, fftype) {
	}
	SparseFileIO(const string &fname, char fun, FTYPE fftype = FeatureFile) :
		FileIO(fname, fun, fftype) {
	}
	SparseFileIO() {
		datapos = 0;
		dimfeature = 0;
		nexamples = 0;
	}
	~SparseFileIO() {

	}
	virtual bool WriteFeature(float*);
	virtual bool WriteFeature(vector<REAL>&);
	virtual bool ReadFeature(float*);
	static UINT MergeFeatureFiles(const string &, const string&);
	static void Convert2SVMLight(const string&, const string&, const string&);
private:
	static bool WriteFeature(ofstream &file, UINT dimfeature, float *feat,
			int fclass);

};

#endif /*FILEIO_H_*/
