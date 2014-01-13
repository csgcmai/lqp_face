/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/
#include "fileIO.h"

FileIO::FileIO(const string &fname, char fun, UINT dfeature, long nexam,
		FTYPE fftype) {
	dimfeature = dfeature;
	nexamples = nexam;
	ftype = fftype;
	datapos = 0;
	if (fun == 'r')
		OpenRead(fname, fftype);
	else if (fun == 'w' || fun == 'a') {
		append = fun == 'a' ? true : false;
		OpenWrite(fname, fftype);
	} else {
		cout << "FileIO: Error in the Specification of File type";
		exit(EXIT_FAILURE);
	}

}

FileIO::FileIO(const string &fname, char fun, FTYPE fftype) {
	datapos = 0;
	Initialize(fname, fun, fftype);
}
void FileIO::Initialize(const string &fname, char fun, FTYPE fftype) {
	CloseFile();
	dimfeature = 0;
	nexamples = 0;
	ftype = fftype;
	if (fun == 'r')
		OpenRead(fname, fftype);
	else {
		cout << " Error in the Specification of Function type";
		exit(EXIT_FAILURE);
	}

}
FileIO::~FileIO() {
	CloseFile();
}
void FileIO::OpenRead(const string & fname, FTYPE fftype) {
	file.open(fname.c_str(), ios_base::in | ios_base::binary);
	if (!file) {
		cout << " Error in Opening for Reading File " << fname << endl;
		exit(EXIT_FAILURE);
	}
	if (fftype == SVMLight)
		ReadSVMLightHeader();
	else if (fftype == FeatureFile)
		ReadHeader();
	else {
		cout
				<< " Error Invalid File Type, Should be Either Feature file or SVM LIght file"
				<< endl;
		exit(EXIT_FAILURE);
	}
	datapos = file.tellg(); // to restart at the data position
}
void FileIO::OpenWrite(const string & fname, UINT dfeature, long nexam) {

	file.open(fname.c_str(), ios_base::out | ios_base::binary);
	if (!file) {
		cout << " Error in Opening for Writing File " << fname << endl;
		exit(EXIT_FAILURE);
	}
	dimfeature = dfeature;
	nexamples = nexam;
	WriteHeader();
}

void FileIO::GetClassCount(UINT &poscount, UINT &negcount) {
	assert(ftype == SVMLight);
	streampos gcpos = file.tellg();
	poscount = negcount = 0;
	vector<REAL> feat(dimfeature, 0);
	int target;
	for (UINT i = 0; i < nexamples; ++i) {
		ReadSVMFeature(&feat[0], target);
		if (target == 1)
			++poscount;
		else
			++negcount;
	}
	file.seekg(gcpos, ios::beg);
	assert(nexamples == (poscount + negcount));

}
void FileIO::OpenWrite(const string & fname, FTYPE fftype) {
	if (!append)
		file.open(fname.c_str(), ios_base::out | ios_base::binary);
	else
		file.open(fname.c_str(), ios_base::out | ios_base::app
				| ios_base::binary);
	if (!file) {
		cout << " Error in Opening for Writing File " << fname << endl;
		exit(EXIT_FAILURE);
	}
	if (fftype == SVMLight)
		WriteSVMLightHeader();
	else if (fftype == FeatureFile)
		WriteHeader();
	else {
		cout
				<< " Error Invalid File Type, Should be Either Feature file or SVM LIght file"
				<< endl;
		exit(EXIT_FAILURE);
	}

}
void FileIO::WriteHeader() {
	file.write((char*) &nexamples, sizeof(long));
	file.write((char*) &dimfeature, sizeof(UINT));
}
void FileIO::ReadHeader() {
	file.read((char*) &nexamples, sizeof(long));
	file.read((char*) &dimfeature, sizeof(UINT));
}

void FileIO::WriteNExamples(long nexam) {
	nexamples = nexam;
	file.seekp(0, ios::beg);
	file.write((char*) &nexamples, sizeof(long));
}
bool FileIO::ReadFeature(float *feat) {
	if (!file.eof() && file) {
		file.read((char*) feat, sizeof(float) * dimfeature);
		return true;
	}
	return false;
}
//bool FileIO::ReadFeature(VectorXf & feat) {
//	if (!file.eof() && file) {
//		for (int i = 0; i < dimfeature; ++i)
//			file.read((char*) &feat[i], sizeof(float));
//		//		file.read((char*) feat, sizeof(float) * dimfeature);
//
//		return true;
//	}
//	return false;
//}
bool FileIO::ReadSVMFeature(float *feat, int& target) {
	if (!file.eof() && file) {
		file.read((char*) &target, sizeof(int));
		file.read((char*) feat, sizeof(float) * dimfeature);
		return true;
	}
	return false;
}
//bool FileIO::ReadSVMFeature(VectorXf &feat, int& target) {
//	if (!file.eof() && file) {
//		file.read((char*) &target, sizeof(int));
//		for (int i = 0; i < dimfeature; ++i)
//			file.read((char*) &feat[i], sizeof(float));
//		return true;
//	}
//	return false;
//}
bool FileIO::WriteFeature(float *feat) {
#ifdef CHECK_NAN
	for(int i=0; i < dimfeature;++i)
	if(isnan(feat[i]) || isinf(feat[i]))
	{
		cout<< " Nan Value found : " << i <<
		" in feature vector"<<endl;
		exit(EXIT_FAILURE);
	}
#endif
	file.write((char*) feat, sizeof(float) * dimfeature);
	if (file)
		return true;

	return false;
}
bool FileIO::WriteFeature(vector<REAL>&feat) {

	//	for(int i=0; i < dimfeature;++i)
	//			{
	//	//			file.write((char*)&feat[i],sizeof(float));
	//				cout<<feat[i]<<" ";
	//			}
	file.write((char*) &feat[0], sizeof(float) * dimfeature);

	//	cout<<endl;
	if (file)
		return true;

	return false;
}

void FileIO::WriteSVMLightHeader() {
	int version = 6, data_typeid = 4, target_typeid = 3,// 4 means float, 3 means int,5 means double
			tfeature = nexamples, dfeature = dimfeature;
	{// Header for the output SVMLight format...

		file.write((char*) &version, sizeof(int));
		file.write((char*) &data_typeid, sizeof(int));
		file.write((char*) &target_typeid, sizeof(int));
		file.write((char*) &tfeature, sizeof(int));
		file.write((char*) &dfeature, sizeof(int));
	}

}
void FileIO::ReadSVMLightHeader() {
	int version = 6, data_typeid = 4, target_typeid = 3, tfeature;// 4 means float, 3 means int,5 means double

	{// Header for the output SVMLight format...

		file.read((char*) &version, sizeof(int));
		file.read((char*) &data_typeid, sizeof(int));
		file.read((char*) &target_typeid, sizeof(int));
		file.read((char*) &tfeature, sizeof(int));
		file.read((char*) &dimfeature, sizeof(int));
	}
	nexamples = tfeature;
}

bool FileIO::WriteFeature(float *feat, int fclass) {
	file.write((char*) &fclass, sizeof(int));
	file.write((char*) feat, sizeof(float) * dimfeature);
	if (file)
		return true;

	return false;
}
bool FileIO::WriteFeature(vector<REAL>&feat, int fclass) {
	file.write((char*) &fclass, sizeof(int));
	file.write((char*) &feat[0], sizeof(float) * feat.size());
	if (file)
		return true;

	return false;
}
void FileIO::Convert2SVMLight(const string &pfname, const string &nfname,
		const string &nofile) {
	FileIO pfile(pfname, 'r'), nfile(nfname, 'r');
	UINT dfeature = pfile.GetDimFeature(), nexamples = pfile.GetNExamples()
			+ nfile.GetNExamples();
	FileIO ofile(nofile, 'w', dfeature, nexamples, SVMLight);
	float *feat = new float[dfeature];
	cout << " Total Number of Examples Written = " << nexamples
			<< " Positive = " << pfile.GetNExamples() << " Negative = "
			<< nfile.GetNExamples() << " Dimension of features = " << dfeature
			<< endl;

	for (int i = 0; i < pfile.GetNExamples(); ++i) {
		pfile.ReadFeature(feat);
		ofile.WriteFeature(feat, 1);
	}

	// Writing the -Ve features
	for (int i = 0; i < nfile.GetNExamples(); ++i) {
		nfile.ReadFeature(feat);
		ofile.WriteFeature(feat, -1);
	}

}
void FileIO::Convert2SparseFormat(const string &pfname) {
	FileIO pfile(pfname, 'r', SVMLight);

	UINT dfeature = pfile.GetDimFeature(), nexamples = pfile.GetNExamples();
	string ofname = pfname + "_text";
	float *feat = new float[dfeature];
	int target;
	ofstream ofile(ofname.c_str(), ios::out);
	for (int i = 0; i < pfile.GetNExamples(); ++i) {
		pfile.ReadSVMFeature(feat, target);
		ofile << target << " ";
		for (UINT j = 0; j < dfeature; ++j) {
			if (feat[j] != 0) {
				ofile << j + 1 << ":" << feat[j] << " ";
			}
		}
		ofile << endl;
	}
	ofile.close();
	pfile.CloseFile();
	string oldfname = pfname + "_old";
	rename(pfname.c_str(), oldfname.c_str());
	rename(ofname.c_str(), pfname.c_str());

}
void FileIO::WriteFile(const string&fname, UINT dfeature, REAL **pfeat,
		UINT npos, REAL **nfeat, UINT nneg) {
	//FTYPE fftype,
	UINT nexamples = npos + nneg;
	FileIO ofile(fname, 'w', dfeature, nexamples, SVMLight);
	cout << " Total Number of Examples Written = " << nexamples
			<< " Positive = " << npos << " Negative = " << nneg
			<< " Dimension of features = " << dfeature << endl;

	for (UINT i = 0; i < npos; ++i) {
		ofile.WriteFeature(pfeat[i], 1);
	}
	for (UINT i = 0; i < nneg; ++i) {
		ofile.WriteFeature(nfeat[i], -1);
	}

	ofile.CloseFile();
}
void FileIO::Convert2SVMLight(const string &pfname,
		const vector<string> &nfname, UINT maxnn, const string &nofile) {
	FileIO pfile(pfname, 'r');
	FileIO *nfiles = new FileIO[(nfname.size())];
	UINT dfeature = pfile.GetDimFeature(), nexamples = pfile.GetNExamples();
	for (int i = 0; i < nfname.size(); ++i) {
		nfiles[i].Initialize(nfname[i], 'r');
		nexamples += nfiles[i].GetNExamples();
	}
	nexamples = nexamples < maxnn ? nexamples : maxnn;
	FileIO ofile(nofile, 'w', dfeature, nexamples, SVMLight);
	float *feat = new float[dfeature];
	cout << " Total Number of Examples Written = " << nexamples
			<< " Positive = " << pfile.GetNExamples() << " Negative = "
			<< (nexamples - pfile.GetNExamples())
			<< " Dimension of features = " << dfeature << endl;

	//	ofstream tofile("features.txt.tmp");
	//	tofile << " Number of features = " << dfeature<<endl;
	//	tofile << " Number of Total Examples = " << nexamples<<
	//	" Positive = "<< pfile.GetNExamples()<< " Negative = "<< nfile1.GetNExamples();
	// Writing the +Ve features
	for (int i = 0; i < pfile.GetNExamples(); ++i) {
		pfile.ReadFeature(feat);
		ofile.WriteFeature(feat, 1);
	}

	// Writing the -Ve features
	UINT excount = 0;
	for (int j = 0; j < nfname.size(); ++j) {
		for (int i = 0; i < nfiles[j].GetNExamples() && excount < maxnn; ++i) {
			nfiles[j].ReadFeature(feat);
			ofile.WriteFeature(feat, -1);
		}
	}
	delete[] feat;
	delete[] nfiles;
}
void FileIO::PermuteFeatures(const string& iffile, const string&offile,
		FTYPE filetype) {
	FileIO ifile(iffile, 'r', filetype);
	UINT dimfeature = ifile.GetDimFeature(), nexamples = ifile.GetNExamples();
	vector<vector<REAL> > features(nexamples, vector<REAL> (dimfeature + 1, 0));
	vector<UINT> index(nexamples, 0);
	REAL *tptr;
	int target;
	for (UINT i = 0; i < nexamples; ++i) {
		if (filetype == SVMLight) {
			tptr = &(features[i])[1];
			ifile.ReadSVMFeature(tptr, target);
			features[i][0] = target;
		} else
			ifile.ReadFeature(&(features[i])[0]);
		index[i] = i;
	}

	RandomPermute(index); // permute randomly example indices
	ifile.CloseFile();

	FileIO ofile(offile, 'w', dimfeature, nexamples, filetype);
	UINT tindex;

	for (UINT i = 0; i < nexamples; ++i) {
		tindex = index[i];
		if (filetype == SVMLight) {
			target = features[tindex][0];
			tptr = &(features[tindex])[1];
			ofile.WriteFeature(tptr, target);

		} else
			ofile.WriteFeature(&(features[tindex])[0]);
	}
	ofile.CloseFile();
}

UINT FileIO::MergeFeatureFiles(const string& nfile1, const string &nfile2) {
	bf::path p1(nfile1), p2(nfile2);
	if (!bf::exists(p1)) {
		FileIO file2(nfile2, 'r');
		long nexamples = file2.GetNExamples();
		file2.CloseFile();
		//		remove(nfile2.c_str());
		rename(nfile2.c_str(), nfile1.c_str());
		return nexamples;
	} else if (!bf::exists(p2)) {
		FileIO file1(nfile1, 'r');
		long nexamples = file1.GetNExamples();
		file1.CloseFile();
		return nexamples;
	}
	{

		FileIO file1(nfile1, 'r'), file2(nfile2, 'r');
		UINT dfeature = file1.GetDimFeature();
		long nexamples = file1.GetNExamples() + file2.GetNExamples();
		if (dfeature != file2.GetDimFeature()) {
			cout << " Error: The dimension of features doesn't match";
			exit(EXIT_FAILURE);
		}
		string ofname = "_tmp_.nfeat";
		FileIO ofile(ofname.c_str(), 'w', dfeature, nexamples);

		float *feat = new float[dfeature];
		for (int i = 0; i < file1.GetNExamples(); ++i) {
			file1.ReadFeature(feat);
			ofile.WriteFeature(feat);
		}
		for (int i = 0; i < file2.GetNExamples(); ++i) {
			file2.ReadFeature(feat);
			ofile.WriteFeature(feat);
		}
		file1.CloseFile();
		file2.CloseFile();
		ofile.CloseFile();

		remove(nfile1.c_str());
		rename(ofname.c_str(), nfile1.c_str());
		delete[] feat;

		return nexamples;
	}
}
void FileIO::WriteFeatureInSVMFormat(const string&fname,
		vector<vector<REAL> >&features, vector<int>&target) {
	UINT dfeature = features[0].size(), nexam = features.size();
	assert(target.size() == features.size());
	FileIO ofile(fname, 'w', dfeature, nexam, SVMLight);

	for (UINT i = 0; i < nexam; ++i) {
		REAL *feat = &(features[i][0]);
		ofile.WriteFeature(feat, target[i]);
	}
	ofile.CloseFile();
}

bool SparseFileIO::WriteFeature(float *feat) {
	UINT nzcount = 0;
	for (UINT i = 0; i < dimfeature; ++i)
		if (feat[i] != 0)
			nzcount++;

	file.write((char*) &nzcount, sizeof(UINT));
	for (UINT i = 0; i < dimfeature; ++i)
		if (feat[i] != 0) {
			file.write((char*) &i, sizeof(UINT));
			file.write((char*) (feat + i), sizeof(float));
		}
	if (file)
		return true;

	return false;
}

bool SparseFileIO::ReadFeature(float *feat) {
	if (!file.eof() && file) {
		UINT nzcount;
		UINT index;
		file.read((char*) &nzcount, sizeof(UINT));
		for (UINT j = 0; j < nzcount; ++j) {
			file.read((char*) &index, sizeof(UINT));
			file.read((char*) (feat + index), sizeof(float));
		}
		return true;
	}
	return false;
}
bool SparseFileIO::WriteFeature(vector<REAL>&feat) {
	WriteFeature(&feat[0]);
}
bool SparseFileIO::WriteFeature(ofstream &file, UINT dimfeature, float *feat,
		int fclass) {
	file << fclass << " ";

	for (UINT i = 0; i < dimfeature; ++i)
		if (feat[i] != 0) {
			file << i + 1 << ":" << feat[i] << " ";
		}
	file << endl;
	if (file)
		return true;
	return false;
}
void SparseFileIO::Convert2SVMLight(const string &pfname, const string &nfname,
		const string &nofile) {
	SparseFileIO pfile(pfname, 'r'), nfile(nfname, 'r');
	UINT dfeature = pfile.GetDimFeature(), nexamples = pfile.GetNExamples()
			+ nfile.GetNExamples();
	ofstream ofile(nofile.c_str(), ios::out);
	float *feat = new float[dfeature];
	cout << " Total Number of Examples Written = " << nexamples
			<< " Positive = " << pfile.GetNExamples() << " Negative = "
			<< nfile.GetNExamples() << " Dimension of features = " << dfeature
			<< endl;

	for (int i = 0; i < pfile.GetNExamples(); ++i) {
		pfile.ReadFeature(feat);
		WriteFeature(ofile, dfeature, feat, 1);
	}

	// Writing the -Ve features
	for (int i = 0; i < nfile.GetNExamples(); ++i) {
		nfile.ReadFeature(feat);
		WriteFeature(ofile, dfeature, feat, -1);
	}
	ofile.close();
}
UINT SparseFileIO::MergeFeatureFiles(const string& nfile1, const string &nfile2) {
	bf::path p1(nfile1), p2(nfile2);
	if (!bf::exists(p1)) {
		SparseFileIO file2(nfile2, 'r');
		long nexamples = file2.GetNExamples();
		file2.CloseFile();
		//		remove(nfile2.c_str());
		rename(nfile2.c_str(), nfile1.c_str());
		return nexamples;
	} else if (!bf::exists(p2)) {
		SparseFileIO file1(nfile1, 'r');
		long nexamples = file1.GetNExamples();
		file1.CloseFile();
		return nexamples;
	}
	{

		SparseFileIO file1(nfile1, 'r'), file2(nfile2, 'r');
		UINT dfeature = file1.GetDimFeature();
		long nexamples = file1.GetNExamples() + file2.GetNExamples();
		if (dfeature != file2.GetDimFeature()) {
			cout << " Error: The dimension of features doesn't match";
			exit(EXIT_FAILURE);
		}
		string ofname = "_tmp_.nfeat";
		SparseFileIO ofile(ofname.c_str(), 'w', dfeature, nexamples);

		float *feat = new float[dfeature];
		for (int i = 0; i < file1.GetNExamples(); ++i) {
			file1.ReadFeature(feat);
			ofile.WriteFeature(feat);
		}
		for (int i = 0; i < file2.GetNExamples(); ++i) {
			file2.ReadFeature(feat);
			ofile.WriteFeature(feat);
		}
		file1.CloseFile();
		file2.CloseFile();
		ofile.CloseFile();

		remove(nfile1.c_str());
		rename(ofname.c_str(), nfile1.c_str());
		delete[] feat;

		return nexamples;
	}
}

void FileIO::ReadFile(const string&nffile, vector<vector<REAL> >& featmat,
		vector<int> &labels, FTYPE ftype) {
	cout << " Reading Feature Files" << nffile << endl;
	FileIO ffile(nffile, 'r', ftype);
	UINT dfeature = ffile.GetDimFeature(), nexamples = ffile.GetNExamples();
	featmat.resize(nexamples, vector<REAL> (dfeature, 0));// read data in the form dimfeature x nexamples as column major order is read...
	labels.resize(nexamples, 0);
	UINT excount = 0;
	int target;
	// Read the features into memory...
	for (UINT i = 0; i < nexamples; ++i, ++excount) // Get the sum of Positive examples...
		if (ftype == SVMLight)
			ffile.ReadSVMFeature(&(featmat[i][0]), labels[i]);
		else
			ffile.ReadFeature(&(featmat[i][0]));
	assert(excount == nexamples);
	cout << "      Finished Reading Files " << endl;
}

//void FileIO::WriteFile(vector<vector<REAL> >& featmat, vector<int> &labels,
//		const string& nofile) {
//
//	UINT nexamples = featmat.size();
//	UINT dimfeature = featmat[0].size();
//	FileIO outfile(nofile, 'w', featmat[0].size(), featmat.size(), SVMLight);
//	for (int i = 0; i < featmat.size(); ++i) {
//		outfile.WriteFeature(&(featmat[i][0]), labels[i]);
//	}
//	outfile.CloseFile();
//}
