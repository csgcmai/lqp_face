/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/

#ifndef ENCODEFEATURES_H_
#define ENCODEFEATURES_H_
#include <unordered_map>
#include "../util/util.hpp"
//#include "dataPoints.hpp"
//#include "cluster.hpp"
#include "../MPI_KMeans/mpi_kmeans.h"
#include "../util/classInfo.h"
#include "definitions.h"
/*ONLINE_MAPPING used for large number of codes such as 2^48*/
//#define ONLINE_MAPPING // Uses Maps (and unordered maps) for the large number of codes, build a lookup table on the fly...
 /* For machines with small memory, codes are not stored they are only constructed at the time of codebook learning*/

/*Online-pruning can be also used deleting codes with small counts and thus saving memory*/
//#define PRUNE_ONLINE // Prunes the lively built lookup table to keep the memory footprint less than PRUNE_LIMIT
/*Note that STL Map has almost a factor of 2 overhead */
#define PRUNE_LIMIT 0.5 // Maximum size (in Gigabytes)to be used for each of the look table up (used during online pruning)...
#define MPI_KMEANS // Use MPI_KMEANS algorithm
#define WITH_PTR // removes the memory overhead of vectors (and pointers) allocations.
//#define ONLINE_MAPPING ///> online mapping of codes, relatively slow but allows much larger neighbourhood with more pixels to be encoded..
typedef char /*original int*/CodeVecType;
typedef unsigned int CodeCountType; //2^32 = 4e9 patches can
typedef short ClusterIdType;



#define MAX_NUMBER_CODES 12e6
/*
 * Compute Code for LQP's
 * */
//template<class T>
//void Print(vector<T>::iterator &iter, int length) {
//	cout << endl;
//	for (int i = 0; i < length; ++i)
//		cout << " " << *(iter + i);
//	cout << endl << flush;
//}
#ifdef ONLINE_MAPPING /*TODO: Rewrite using polymorphism*/
struct CodeCache {
	CodeCache(UINT count, ClusterIdType cid = -1) {
		codecount = count;
		clusid = cid;

		/*
		 * Printing for debugging
		 */
//		cout << endl;
//		for (int i = 0; i < length; ++i)
//			cout << " " << *(iter + i);
//		cout << endl << flush;
	}
	CodeCache(const CodeCache &cobj) {
		*this = cobj;
	}
	void operator=(const CodeCache &cobj) {
		codecount = cobj.codecount;
		clusid = cobj.clusid;
	}
	void WriteToFile(ofstream &ofile) {
		ofile.write((char*) &codecount, sizeof(UINT));
		ofile.write((char*) &clusid, sizeof(ClusterIdType));
	}
	void ReadFromFile(ifstream &ifile) {
		ifile.read((char*) &codecount, sizeof(UINT));
		ifile.read((char*) &clusid, sizeof(ClusterIdType));
	}
	static long GetSize() {
		return sizeof(UINT) + sizeof(ClusterIdType);
	}
	~CodeCache() {
	}
	UINT codecount;
	ClusterIdType clusid;
};
/*Can be a ordered_map with slow access time O(log(n)) or unordered_map with fast access time i.e.O(1)*/
typedef unordered_map<long, CodeCache*> CodeTable;
#else
struct CodeCache {
	CodeCache() {
		codecount = 0;
		code = new CodeVecType[length];
		clusid = -1;
	}
	template<class T>
	CodeCache(T *vec, UINT count, ClusterIdType cid) {
		codecount = count;
		clusid = cid;
		code = new CodeVecType[length];
		copy(vec, vec + length, code);
	}
	template<class T>
	CodeCache(const T &iter, UINT count) {
		codecount = count;
		clusid = -1;

		code = new CodeVecType[length];
		copy(iter, iter + length, code);
		/*
		 * Printing for debugging
		 */
//		cout << endl;
//		for (int i = 0; i < length; ++i)
//			cout << " " << *(iter + i);
//		cout << endl << flush;
	}
	CodeCache(const CodeCache &cobj) {
		*this = cobj;
	}
	void operator=(const CodeCache &cobj) {
		codecount = cobj.codecount;
		clusid = cobj.clusid;
		copy(code, code + length, cobj.code);
	}
	void WriteToFile(ofstream &ofile) {
		ofile.write((char*) &codecount, sizeof(UINT));
		ofile.write((char*) &clusid, sizeof(ClusterIdType));
		ofile.write((char*) code, sizeof(CodeVecType) * length);
	}
	void ReadFromFile(ifstream &ifile) {
		ifile.read((char*) &codecount, sizeof(UINT));
		ifile.read((char*) &clusid, sizeof(ClusterIdType));
		ifile.read((char*) code, sizeof(CodeVecType) * length);
	}
	static long GetSize() {
		return sizeof(CodeVecType) * length + sizeof(UINT)
		+ sizeof(ClusterIdType) + sizeof(CodeVecType*);
	}
	~CodeCache() {
		delete[] code;
	}
	CodeVecType *code;
	UINT codecount;
	ClusterIdType clusid;
	static int length;
};
#endif
/*tocheck-result between unordered_map and map...*/


class ComputeCode {
public:
	ComputeCode(UINT nlevels_ = 1, UINT neighbours_ = 8, UINT nclusters = 100,
			QuantizationType soft_ = QT_VQ, REAL mratio_ = 1, UINT pcode_ = 0,
			UINT knn_ = 5, int discinfo_ = 0) :
			nlevels(nlevels_), neighbours(neighbours_), ncenters(nclusters), soft(
					soft_), mratio(mratio_), pcount(pcode_), knn(knn_) {
		//		clusters = new REAL[ncenters];
		discinfo = discinfo_;
		codelen = neighbours + discinfo;
#ifndef ONLINE_MAPPING
		CodeCache::length = codelen; ///> number of codes = number of neighbours considered
		cout<<"\n ------------> Codes are Computed Offline and Stored in An Array<--------------- \n";
#else
		cout<<"\n ------------> Codes are Computed only at codebook learning stage and stored in Map (or Unordered)<--------------- \n";
#endif
		clusters = NULL;
		AllocateClustersMemory();
		cout << "\n Mapping Codes to Center via "
				<< CBParams::GetQuantizationType(soft) << endl;

	}
	virtual ~ComputeCode() {
		if (clusters)
			delete[] clusters;
	}
	void WriteChannel(const string &fname);
	REAL ComputeSSE(vector<vector<CodeVecType> > &codevec,
			vector<vector<REAL> > &clusters);
	REAL ComputeSSE(vector<vector<CodeVecType> > &codevec, REAL *clusters);
	REAL ComputeSSE(CodeVecType *codevec, REAL *clusters);
	void GenerateCodeVectors(CodeVecType *& codevec);
	bool ReadChannel(const string &fname);
	void WriteChannelBinary(const string &fname);
	bool ReadChannelBinary(const string &fname);
	void WriteClusters(const string &fname);
	UINT FindIndexNonZeroCodes(vector<UINT>& indexes, REAL *ratio = NULL);
	void GenerateCodeBook(const string & cifname, const string & cbfname,
			const ClusteringDistanceMetric &dmetric, UINT nclusrounds,
			REAL *initclusters = NULL, REAL *ratio = NULL);
	void WriteLQPCodes(const string &fname);
	UINT WriteLQPCodes(ofstream &ofile);
	void WriteHeader(ofstream &ofile) {
		int version = 8; // version 5 write in blocks
		//only in version 008
		ofile.write((char*) &discinfo, sizeof(UINT)); // length vector
		//version 007
		ofile.write((char*) &version, sizeof(int));
		ofile.write((char*) &ncodes, sizeof(UINT)); // ncodes
		ofile.write((char*) &lenvec, sizeof(UINT)); // base values
		ofile.write((char*) &midp, sizeof(UINT)); // midpoint
		ofile.write((char*) &neighbours, sizeof(UINT)); // length vector
		ofile.write((char*) &ncenters, sizeof(UINT)); // quantization
		ofile.write((char*) &clusind[0], sizeof(ClusterIdType) * ncodes); // lookup table...
	}
	void WriteCode(const string &fname) {
		WriteChannelBinary(fname);
	}
	bool ReadCode(const string &fname) {
		return ReadChannelBinary(fname);
	}
	bool ReadCodeText(const string&fname) {
		return ReadChannel(fname);
	}
	void WriteCodeText(const string &fname) {
		WriteChannel(fname);
	}
	void PruneCodes(); ///> Prune the codes having small number of counts i.e count < pcode
#ifdef ONLINE_MAPPING ///> use online mapping of codes,
	void SaveCodeStatistics(const string & fname);
	void SetValue(vector<REAL>& codevec, const long &rvalue, const long &gvalue,
			const long &bvalue
			/*,const float ipercen/* what percentage of images have been processed*/) {
		///> function used with online-mapping flag, initializes codeCache entries.

#ifdef PRUNE_ONLINE
		const double limit = pow(2, 30) * PRUNE_LIMIT / 2, usage =
		codecache.size() * CodeCache::GetSize();
//		if (codecache.size())
		if (pcount > 0 && usage >= limit) {
			cout << endl << "Usage = " << usage << endl << " Limit = " << limit
			<< endl
			<< " Doing Online Pruning of Codes that have counts less than Prune_LIMIT"
			<< endl;
			PruneCodes();
			cout << endl << " Finished Pruning of online codes " << endl;
		}

#endif
		if (codecache.count(rvalue))
			codecache[rvalue]->codecount++;
		else {
//For Testing purposes
			vector<REAL> tvar;
			GenerateCodeVector(rvalue, tvar);
//			for (vector<REAL>::iterator iter = tvar.begin(); iter != tvar.end();
//					++iter)
//				cout << " " << (*iter) << flush;
//			cout << endl << flush;
			long val = GetCode(vector<CodeVecType>(tvar.begin(), tvar.end()));
//			vector<REAL>::iterator rviter = codevec.begin();
//			for (vector<CodeVecType>::iterator iter = tvar.begin();
//					iter != tvar.end(); iter++, rviter++)
//				assert(*iter==*rviter);
			assert(val==rvalue);
			codecache[rvalue] = new CodeCache(1);
		}

		if (codecache.count(gvalue))
			codecache[gvalue]->codecount++;
		else
			codecache[gvalue] = new CodeCache(1);

		if (codecache.count(bvalue))
			codecache[bvalue]->codecount++;
		else
			codecache[bvalue] = new CodeCache(1);
		/*assert(codecache.size() < 30e6);
		 /*if (tmpcodecache.count(rvalue))
		 tmpcodecache[rvalue]++;
		 else
		 tmpcodecache[rvalue] = 1;

		 if (tmpcodecache.count(gvalue))
		 tmpcodecache[gvalue]++;
		 else
		 tmpcodecache[gvalue] = 1;

		 if (tmpcodecache.count(bvalue))
		 tmpcodecache[bvalue]++;
		 else
		 tmpcodecache[bvalue] = 1;*/
	}
	ClusterIdType GetClusterCenter(REAL *codevec, long code) { ///> function used with online-mapping
		if (codecache.count(code)) {
			//			cout << " " << codecache[code]->clusid << " " << flush;
			return codecache[code]->clusid;
			/*	if (tcid < 0) {
			 cout << " \n Found -ve cid= " << tcid << " " << code<< " ";
			 for (int i = 0; i < neighbours; ++i)
			 cout << codevec[i] << " ";
			 cout << endl;
			 for (CodeTable::iterator iter = codecache.begin();
			 iter != codecache.end(); ++iter) {
			 cout << "\n Code=" << iter->first << " "
			 << codecache.count(iter->first) << " " << flush;
			 }
			 cout << endl;
			 }
			 return tcid;*/
		}
		ClusterIdType cid = assign_point_to_cluster_ordinary(codevec, clusters,
				neighbours, ncenters);
		codecache[code] = new CodeCache(0, cid);
		//		cout << " " << cid << " " << flush;
		if (cid < 0 || cid >= ncenters) {
			vector<REAL> vec(codevec, codevec + neighbours);
			cout << " Glitch found" << cid << flush << endl;
			cout << codevec << flush;
		}
		return cid;
	}
#endif
	template<class T>
	ClusterIdType GetClusterCenter(REAL *codevec, T code) { ///> function used with online-mapping
		ClusterIdType cid = assign_point_to_cluster_ordinary(codevec, clusters,
				neighbours, ncenters);
		//		if (cid < 0 || cid >= ncenters) {
		//			vector<REAL> vec(codevec, codevec + neighbours);
		//			cout << " Glitch found" << cid << flush << endl;
		//			cout << codevec << flush;
		//		}
		assert(cid>=0 && cid < ncenters);
		return cid;
	}
	void SetValue(const UINT &rvalue, const UINT &gvalue, const UINT &bvalue) {
		// map the three value codes to respective
		codecount[rvalue]++;
		codecount[gvalue]++;
		codecount[bvalue]++;
	}
	void SetValue(const UINT &rvalue) {
		// map the three value codes to respective
		codecount[rvalue]++;
	}
	UINT GetBase() {
		return lenvec;
	}
	UINT GetNCodes() {
		return ncodes;
	}
	UINT GetCodeLength() {
		return codelen;
	}
	template<class T>
	void GetCodeCount(vector<T>&vec) {
		vec.resize(codecount.size(), 0);
		copy(codecount.begin(), codecount.end(), vec.begin());
	}
	bool ExistCodeBook(const string&fname) {
		if (ReadCode(fname)) {
			ClearCodeCount();
			return true;
		}
		return false;
	}
	bool ExistCodeInfo(const string&fname) {
		return ReadCode(fname);
	}
	void AllocateClustersMemory() {
		if (clusters)
			delete[] clusters;
		clusters = new REAL[ncenters * (codelen)];
		fill(clusters, clusters + ncenters * (codelen), 0);
	}
	UINT GetNCenters() {
		return ncenters;
	}
	template<class T>
	ClusterIdType GetClusterCenter(T code) {
		return clusind[code];
	}
	void GetClusterCenterWithRatio(UINT &code, double &value) {
		code = clusind[code];
		value = *(clusters + code * codelen + neighbours); //TODO: replace with a vector
	}
	const vector<REAL> & GetSoftClusterCenter(UINT code) {
		return softcodes[code];
	}
	void CopyClusters(REAL *outclus) {
		copy(clusters, clusters + ncenters * codelen, outclus);
	}
	void ClearData() {
#ifdef ONLINE_MAPPING ///> use online mapping of codes,
		for (CodeTable::iterator iter = codecache.begin();
				iter != codecache.end(); ++iter)
			delete iter->second;
		codecache.clear();
#else
		if (codecount.size() > 0) {
			codecount.clear();
			clusind.clear();
		}
		ComputeInfo();
#endif

	}
	virtual void ComputeInfo()=0;
	bool LoadLUT(const string & fname) {
		ifstream ifile(fname.c_str(), ios::in | ios::binary);
		if (!ifile) {
			cerr << endl << " Could't load " << fname << " file" << endl;
			exit(EXIT_FAILURE);
		}
		ifile.seekg(0, ios::end);
		UINT length = ifile.tellg();
		ifile.seekg(0, ios::beg);
		assert(length / sizeof(ClusterIdType) == ncodes);
		clusind.resize(ncodes, 0);
		fill(clusind.begin(), clusind.end(), 0);
		ifile.read((char*) (&clusind[0]), sizeof(ClusterIdType) * ncodes); //sizeof(ClusterIdType)*length / sizeof(ClusterIdType)
		cout << endl << clusind[0] << " -------------- last ="
				<< clusind[ncodes - 1] << flush;
		return true;
	}
	template<class T>
	void AddCodeCount(vector<T> &ncodecount) {
		typename vector<T>::iterator niter = ncodecount.begin();
		for (vector<UINT>::iterator iter = codecount.begin();
				iter != codecount.end(); ++iter, ++niter)
			*iter += *niter;
	}
protected:
#ifdef ONLINE_MAPPING ///> use online mapping of codes,
//	unordered_map<long, UINT> tmpcodecache; // to know total number of nonzero codes...
	CodeTable codecache;
	map<long, UINT> tmpcodecache; // to know total number of nonzero codes...

#endif
	int discinfo; /// whether discriminant class information will be included or not, variable contains number of bins including the discriminant information

	UINT pcount; // prune the codes with counts < pcount;
	UINT midp;
	UINT ncodes, lenvec, codelen; // offset required for the computation of multi-level lbp's...
	UINT nlevels, neighbours, ncenters; // ncenters is the size of code book
	QuantizationType soft; // use soft-quantization..
	UINT knn; // number of nearest neighbours for LLC approximation
	REAL mratio; // what percentage of mean to use as thresholding the distance... mean all distances falling below
// this threshold are considered and all which are above truncated...
// Code Information Values...
	vector<UINT> codeid; // code is equivalent to index... so duplication..
	vector<CodeCountType> codecount;
	vector<vector<REAL> > softcodes; /// softcodes used for mapping..

	vector<ClusterIdType> clusind; // cluster index for this code...
	REAL *clusters;
	vector<int> mapping; // used to calculate the initial codes.
	void ClearCodeCount() {
#ifdef ONLINE_MAPPING
		for (CodeTable::iterator iter = codecache.begin();
				iter != codecache.end(); ++iter) {
			iter->second->codecount = 0;
		}
#else
		codecount.clear();
#endif
	}
	virtual long GetCode(const vector<CodeVecType> & codevec)=0;
	virtual long GetCode(const CodeVecType *codevec)=0;
	virtual void InitializeMapping()=0;
	//	template<class T>
	//	void ExtractClusterInfo(const Clustering::PointsSpace<T> &ps);
	//	template<class T>
	//	void CopyData(vector<vector<CodeVecType> > &codevec,
	//			Clustering::PointsSpace<T> &ps);
	void CopyData(vector<vector<CodeVecType> > &codevec, vector<UINT> &nzidx,
			REAL *data, uInt*ptcount);
	void CopyData(CodeVecType *codevec, vector<UINT> &nzidx, REAL *data,
			uInt *ptcount, REAL *ratio = NULL);
	void InitializeClusters(REAL *data, UINT nnzcodes, REAL *clusters);
	void InitializeClusters(vector<vector<CodeVecType> > &codevec,
			vector<UINT> &indeces, REAL *clusters);
	void AllocateMemory();
	void GenerateCodeVectors(vector<vector<CodeVecType> > &codevec);
	void GenerateCodeVector(const long &value, vector<REAL>&codevec);/*Generate a single code vector from given code*/
	void GenerateCodeVector(const long &value, REAL *codevec);
	void MapCodes2Centers(REAL* x = NULL);
	/**
	 * This Method maps the codes to k nearest cluster centers using locality-constrained
	 * linear coding method..
	 * @param cbook: learned codebook
	 * @param x: new feature x which will be mapped
	 * @param knn: number of k-nearest-neighbours
	 * @param lambda: regularization constant
	 * @param mapping: the llc mapping for the given feature or code
	 */
	void LLCMapping(Map<MatrixXf> & cbook, REAL* tx, UINT knn, REAL lambda,
			REAL *mapping);
};
class ComputeCodeLBP: public ComputeCode {
public:
	ComputeCodeLBP(UINT nlevels_ = 1, UINT neighbours_ = 8,
			UINT nclusters = 100, QuantizationType softq = QT_VQ,
			REAL meanp = 1, UINT pcode_ = 0, UINT knn_ = 5, int discinfo_ = 0) :
			ComputeCode(nlevels_, neighbours_, nclusters, softq, meanp, pcode_,
					knn_, discinfo_) {
		ComputeInfo();

	}
	virtual ~ComputeCodeLBP() {

	}
	virtual void InitializeMapping() {
		mapping.resize(lenvec, 0);
		for (UINT i = 1; i <= lenvec; ++i) {
			mapping[i] = i;
		}
	}
	virtual void ComputeInfo();
	virtual long GetCode(const vector<CodeVecType> & codevec);
	virtual long GetCode(const CodeVecType *codevec);
};

class ComputeCodeLTP: public ComputeCode {
public:
	ComputeCodeLTP(UINT nlevels_ = 1, UINT neighbours_ = 8,
			UINT nclusters = 100, QuantizationType softq = QT_VQ,
			REAL meanp = 1, UINT pcode_ = 0, UINT knn_ = 5, int discinfo_ = 0) :
			ComputeCode(nlevels_, neighbours_, nclusters, softq, meanp, pcode_,
					knn_, discinfo_) {

		ComputeInfo();

	}
	virtual void InitializeMapping() {
		mapping.resize(lenvec, 0); // lenvec refers to base of vector..
		mapping[0] = 0;
		for (int i = 1; i <= nlevels; ++i) {
			mapping[i] = i;
			mapping[i + nlevels] = -i;
		}
	}
	virtual ~ComputeCodeLTP() {

	}
	virtual void ComputeInfo();
	virtual long GetCode(const vector<CodeVecType> & codevec);
	virtual long GetCode(const CodeVecType *codevec);
};

#endif /* ENCODEFEATURES_H_ */
