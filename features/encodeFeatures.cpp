/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/
#include "encodeFeatures.h"
#include "../Eigen/Dense"
/*
 * encodeFeatures.h
 *
 *  Created on: Sep 9, 2010
 *      Author: shussain
 */
#ifndef ONLINE_MAPPING
int CodeCache::length = 0;
#endif
#ifndef STATIC_CODES
void ComputeCode::WriteChannel(const string &fname) {
	ofstream ofile(fname.c_str(), ios::out);
	if (!ofile) {
		cout << " Error: Opening File " << fname;
		exit(EXIT_FAILURE);
	}
#ifdef TMPDEBUG
	vector < vector<CodeVecType> > codevec;
	GenerateCodeVectors(codevec);
#endif
	for (UINT i = 0; i < ncodes; ++i) {
		ofile << codeid[i] << " " << codecount[i] << " " << clusind[i] << " ";
#ifdef TMPDEBUG
		for (UINT j = 0; j < codelen; ++j)
		ofile << (int) codevec[i][j] << " ";
#endif
		ofile << endl;
	}
	ofile.close();
}
REAL ComputeCode::ComputeSSE(vector<vector<CodeVecType> > &codevec,
		vector<vector<REAL> > &clusters) {

	REAL sse = 0;
	for (UINT i = 0; i < ncodes; ++i) {
		if (codecount[i] > 0) {
			REAL d = L2Dist(codevec[i], clusters[clusind[i]]);
			sse += d * d * codecount[i];
		}
	}
	return sse;
}
REAL ComputeCode::ComputeSSE(vector<vector<CodeVecType> > &codevec,
		REAL *clusters) {

	REAL sse = 0;
	for (UINT i = 0; i < ncodes; ++i) {
		if (codecount[i] > 0) {
			REAL d = L2Dist(&(codevec[i][0]), clusters + codelen * clusind[i],
					codelen);
			sse += d * d * codecount[i];
		}
	}
	return sse;
}

REAL ComputeCode::ComputeSSE(CodeVecType *codevec, REAL *clusters) {

	REAL sse = 0;
	char *pcodevec = codevec;
	for (UINT i = 0; i < ncodes; ++i) {
		if (codecount[i] > pcount * 3) {
			REAL d = L2Dist(pcodevec, clusters + codelen * clusind[i], codelen);
			sse += d * d * codecount[i];
		}
		pcodevec += codelen;
	}
	return sse;
}
bool ComputeCode::ReadChannel(const string &fname) {
	ifstream ifile(fname.c_str(), ios::in);
	if (!ifile) {
		cout << " Error: Opening File " << fname;
		//			exit(EXIT_FAILURE);
		return false;
	}
	for (UINT i = 0; i < ncodes; ++i) {
		ifile >> codeid[i] >> codecount[i] >> clusind[i];
		//		for (UINT j = 0; j < neighbours; ++j)
		//			ifile >> codevec[i][j];

	}
	ifile.close();
	return true;
}
void ComputeCode::WriteChannelBinary(const string &fname) {
	ofstream ofile(fname.c_str(), ios::out | ios::binary);
	if (!ofile) {
		cout << " Error: Opening File " << fname;
		exit(EXIT_FAILURE);
	}
#ifdef ONLINE_MAPPING
	int version = 005; // version write in blocks
	ofile.write((char*) &version, sizeof(int));
	// only in version 4
	ofile.write((char*) &soft, sizeof(bool));
	ofile.write((char*) &mratio, sizeof(REAL));
	// version 3
	ncodes = codecache.size();
	ofile.write((char*) &ncodes, sizeof(UINT));
	ofile.write((char*) &lenvec, sizeof(UINT));
	ofile.write((char*) &nlevels, sizeof(UINT));
	ofile.write((char*) &neighbours, sizeof(UINT));
	ofile.write((char*) &ncenters, sizeof(UINT));
	ofile.write((char*) clusters, sizeof(REAL) * ncenters * neighbours);
	for (CodeTable::iterator iter = codecache.begin(); iter != codecache.end();
			++iter) {
		ofile.write((char*) &(iter->first), sizeof(long));
		iter->second->WriteToFile(ofile);
	}
#else
	int version = 007;
	ofile.write((char*) &version, sizeof(int));
	// only in version 7
	ofile.write((char*) &discinfo, sizeof(int));
	// only in version 5
	ofile.write((char*) &soft, sizeof(int));
	ofile.write((char*) &knn, sizeof(UINT));
	ofile.write((char*) &mratio, sizeof(REAL));

	// version 3
	ofile.write((char*) &ncodes, sizeof(UINT));
	ofile.write((char*) &lenvec, sizeof(UINT));
	ofile.write((char*) &nlevels, sizeof(UINT));
	ofile.write((char*) &neighbours, sizeof(UINT));
	ofile.write((char*) &ncenters, sizeof(UINT));
	ofile.write((char*) clusters, sizeof(REAL) * ncenters * codelen);
	for (UINT i = 0; i < ncodes; ++i)
	//			tchinfo[i].WriteTextBinary(ofile, false);
	{
		ofile.write((char*) &codeid[i], sizeof(UINT));
		ofile.write((char*) &codecount[i], sizeof(CodeCountType));
		ofile.write((char*) &clusind[i], sizeof(ClusterIdType));
	}
#endif
	ofile.close();
}
bool ComputeCode::ReadChannelBinary(const string &fname) {
	ifstream ifile(fname.c_str(), ios::in | ios::binary);
	if (!ifile) {
		cout << " Error: Opening File " << fname;
		//			exit(EXIT_FAILURE);
		return false;
	}
#ifdef ONLINE_MAPPING
	int version; // version write in blocks
	ifile.read((char*) &version, sizeof(int));
	ifile.read((char*) &soft, sizeof(bool));
	ifile.read((char*) &mratio, sizeof(REAL));
	ifile.read((char*) &ncodes, sizeof(UINT));
	ifile.read((char*) &lenvec, sizeof(UINT));
	ifile.read((char*) &nlevels, sizeof(UINT));
	ifile.read((char*) &neighbours, sizeof(UINT));
	ifile.read((char*) &ncenters, sizeof(UINT));
	ifile.read((char*) clusters, sizeof(REAL) * ncenters * neighbours);
//	CodeCache::length = neighbours;
	UINT tcodeid;
	for (UINT i = 0; i < ncodes; ++i) {
		ifile.read((char*) &tcodeid, sizeof(long));
		codecache[tcodeid] = new CodeCache(0);
		codecache[tcodeid]->ReadFromFile(ifile);
	}
#else
	int version;
	ifile.read((char*) &version, sizeof(int));
	if (version == 7) { // 5 used for online mapping
		ifile.read((char*) &discinfo, sizeof(int));
		ifile.read((char*) &soft, sizeof(int));
		ifile.read((char*) &knn, sizeof(UINT));
		ifile.read((char*) &mratio, sizeof(REAL));
	}
	if (version == 6) { // 5 used for online mapping
		ifile.read((char*) &soft, sizeof(int));
		ifile.read((char*) &knn, sizeof(UINT));
		ifile.read((char*) &mratio, sizeof(REAL));
	}
	if (version == 4) {
		ifile.read((char*) &soft, sizeof(bool));
		ifile.read((char*) &mratio, sizeof(REAL));
	}
	if (version >= 3) {
		ifile.read((char*) &ncodes, sizeof(UINT));
		ifile.read((char*) &lenvec, sizeof(UINT));
		ifile.read((char*) &nlevels, sizeof(UINT));
		ifile.read((char*) &neighbours, sizeof(UINT));
		codelen = neighbours + discinfo;
		AllocateMemory();
		ifile.read((char*) &ncenters, sizeof(UINT));
		AllocateClustersMemory();
		ifile.read((char*) clusters, sizeof(REAL) * ncenters * codelen);
		//			codevec.resize(ncodes, vector<CodeVecType>(neighbours));
		for (UINT i = 0; i < ncodes; ++i) {
			ifile.read((char*) &codeid[i], sizeof(UINT));
			ifile.read((char*) &codecount[i], sizeof(CodeCountType));
			ifile.read((char*) &clusind[i], sizeof(ClusterIdType));
		}
		if (version >= 4 && soft) {
			cout << "\n Mapping Codes to Center via "
			<< CBParams::GetQuantizationType(soft) << endl;
			//			if (fname
			//					== "CircularShallowPET_SplitLTP-patchsize5-ncenters150-CodeBook-.txt150") {
			//				ReadRealMatrixFromBinFile("Mapping-1328285267", softcodes);
			//			} else
			//				ReadRealMatrixFromBinFile("Mapping-1328285342", softcodes);

			MapCodes2Centers();
		}
	} else {
		vector<vector<CodeVecType> > codevec;
		codevec.resize(ncodes, vector<CodeVecType> (neighbours));
		if (version == 2) {
			ifile.read((char*) &ncodes, sizeof(UINT));
			ifile.read((char*) &neighbours, sizeof(UINT));
			AllocateMemory();
		}

		for (UINT i = 0; i < ncodes; ++i) {
			ifile.read((char*) &codeid[i], sizeof(UINT));
			ifile.read((char*) &codecount[i], sizeof(CodeCountType));
			ifile.read((char*) &clusind[i], sizeof(ClusterIdType));
			const CodeVecType* tcodevec = &(codevec[i][0]);
			ifile.read((char*) tcodevec, sizeof(CodeVecType) * neighbours);
		}
	}

#endif
	ifile.close();
	return true;
}
void ComputeCode::WriteLQPCodes(const string &fname) {
	ofstream ofile(fname.c_str(), ios::out | ios::binary);
	if (!ofile) {
		cout << " Error: Opening File " << fname;
		exit(EXIT_FAILURE);
	}
	WriteHeader(ofile);
	WriteLQPCodes(ofile);
	ofile.close();
}
UINT ComputeCode::WriteLQPCodes(ofstream &ofile) {
	vector<UINT> nzidx;
	UINT nnz = FindIndexNonZeroCodes(nzidx);
	ofile.write((char*) &nnz, sizeof(UINT));
	ofile.write((char*) (&nzidx[0]), sizeof(UINT) * nnz);
	for (UINT i = 0; i < nnz; ++i)
		nzidx[i] = codecount[nzidx[i]];
	// writing code-count...
	ofile.write((char*) (&nzidx[0]), sizeof(UINT) * nnz);
	return nnz;
}
void ComputeCode::WriteClusters(const string &fname) {
	ofstream ofile(fname.c_str(), ios::out | ios::binary);
	if (!ofile) {
		cout << " Error: Opening File " << fname;
		exit(EXIT_FAILURE);
	}
	int version = 10; // version 10 includes the cluster-index(lqp mapping)
	ofile.write((char*) &version, sizeof(int));
	// only in version 4
	ofile.write((char*) &soft, sizeof(bool));
	ofile.write((char*) &mratio, sizeof(REAL));
	// version 3
	ofile.write((char*) &ncodes, sizeof(UINT));
	ofile.write((char*) &lenvec, sizeof(UINT)); // lenvec base vector...
	ofile.write((char*) &nlevels, sizeof(UINT)); // 1, for binary et ternary, 2 for quinary approx=base/2;
	ofile.write((char*) &codelen, sizeof(UINT));
	ofile.write((char*) &ncenters, sizeof(UINT));
	ofile.write((char*) clusters, sizeof(REAL) * ncenters * codelen);
	ofile.write((char*) (&codecount[0]), sizeof(CodeCountType) * ncodes);
	// included only in version 10
	ofile.write((char*) (&clusind[0]), sizeof(ClusterIdType) * ncodes);
	ofile.close();
}
void ComputeCode::LLCMapping(Map<MatrixXf> & cbook, REAL* tx, UINT knn,
		REAL lambda, REAL* mapping) {
	/// first find the k-nn
	Map<MatrixXf> x(&tx[0], 1, cbook.cols());
	/// finding the distance to all the cluster centers
	MatrixXf t = (cbook - x.replicate(cbook.rows(), 1));
	t = t.cwiseProduct(t);
	MatrixXf v = t.rowwise().sum();
	typedef REAL * T3;
	vector<SorterElement<T3> > res3 = sort_elements(v.data(),
			v.data() + cbook.rows(), Cmp<REAL>());

	int j = 0;
	MatrixXf B(knn, cbook.cols());

	for (vector<SorterElement<T3> >::iterator i = res3.begin();
			i != res3.end() && j < knn; ++i, ++j) {
		B.row(j) = (cbook.row(i->_ind) - x);
	}
	/// Solve the linear problem of lcc
	B = B * B.transpose();
	B += MatrixXf::Identity(knn, knn) * lambda * B.trace();
	VectorXf w = B.jacobiSvd(ComputeThinU | ComputeThinV).solve(
			VectorXf::Ones(knn));
	w /= w.sum();
	fill(mapping, mapping + cbook.rows(), 0);
	j = 0;
	for (vector<SorterElement<T3> >::iterator i = res3.begin();
			i != res3.end() && j < knn; ++i, ++j) {
		mapping[i->_ind] = w(j);
	}
}
void ComputeCode::MapCodes2Centers(REAL *ratio) {
	cout << endl << "Mapping Codes to Clusters ------>" << endl << flush;
#ifndef WITH_PTR
	vector<vector<CodeVecType> > codevec;
	GenerateCodeVectors(codevec); // Regenerate Code Vectors...
	//#ifdef TMPDEBUG
	//		ReadCode(cifname); // to find the cluster assignments  if CodeVector is Cleared above then used it...
	for (UINT i = 0; i < ncodes; ++i)
	clusind[i] = assign_point_to_cluster_ordinary(&(codevec[i][0]),
			clusters, codelen, ncenters);
	if (soft) {
		softcodes.resize(ncodes, vector<REAL> (ncenters, 0));
		for (UINT i = 0; i < ncodes; ++i)
		assign_point_to_cluster_soft(&(codevec[i][0]), clusters,
				codelen, ncenters, &(softcodes[i][0]), mratio);
	}
	cout << " \n New SSE = " << ComputeSSE(codevec, clusters);
	codevec.clear();
#else
#ifdef ONLINE_MAPPING
	cout << endl;
	vector<REAL> code(neighbours, 0);
	for (CodeTable::iterator iter = codecache.begin(); iter != codecache.end();
			++iter) {
		GenerateCodeVector(iter->first, code);
		iter->second->clusid = assign_point_to_cluster_ordinary(&code[0],
				clusters, codelen, ncenters);
		//		cout << iter->second->clusid << " ";
		////		Print((int*)iter->second->code, (int*)(iter->second->code + codelen));
		//		for (int i = 0; i < codelen; ++i)
		//			cout << (int) iter->second->code[i] << " ";
		//		cout << endl;
	}
	//	cout << endl;
	//	for (CodeTable::iterator iter = codecache.begin();
	//			iter != codecache.end(); ++iter) {
	//		cout << "\n Code=" << iter->first << " " << codecache.count(iter->first)
	//				<< " " << flush;
	//	}
#else
	CodeVecType *codevec;
	GenerateCodeVectors(codevec);
	CodeVecType *pcodevec = codevec;
	if (soft == QT_VQ) {
		if (ratio == NULL) {
			for (UINT i = 0; i < ncodes; ++i) {
				clusind[i] = assign_point_to_cluster_ordinary(pcodevec,
						clusters, codelen, ncenters);
				pcodevec += codelen;
			}
		} else {
			for (UINT i = 0; i < ncodes; ++i) {
				pcodevec[neighbours] = ratio[i];
				clusind[i] = assign_point_to_cluster_ordinary(pcodevec,
						clusters, codelen, ncenters);
				pcodevec += codelen;
			}
		}
		cout << " \n New SSE = " << ComputeSSE(codevec, clusters);
	} else if (soft == QT_Median) {

		softcodes.resize(ncodes, vector<REAL> (ncenters, 0));
		for (UINT i = 0; i < ncodes; ++i) {
			assign_point_to_cluster_soft(pcodevec, clusters, codelen, ncenters,
					&(softcodes[i][0]), mratio);
			pcodevec += codelen;
		}
	} else {
		long initime = time(0);
		REAL *tvec = new REAL[codelen];
		Map<MatrixXf> cbook(clusters, ncenters, codelen);
		softcodes.resize(ncodes, vector<REAL> (ncenters, 0));
		//		WriteMatrix2File("codebook",cbook);
		for (UINT i = 0; i < ncodes; ++i) {
			copy(pcodevec, pcodevec + codelen, tvec);
			LLCMapping(cbook, tvec, knn, 1e-4, &(softcodes[i][0]));
			pcodevec += codelen;
		}
		cout << " \n Time Taken by Mapping = " << time(0) - initime << endl;
		delete[] tvec;
	}
	delete[] codevec;
#endif
#endif
}
UINT ComputeCode::FindIndexNonZeroCodes(vector<UINT>& indexes, REAL *ratio) {
	vector<UINT>::iterator cid = codeid.begin();
	UINT nnzcodes = 0;
	if (ratio == NULL) { /// No Discriminant information is available...
		for (vector<CodeCountType>::iterator iter = codecount.begin();
				iter != codecount.end(); ++iter, ++cid) {
			if (*iter > 0)
				nnzcodes++;
			if (*iter > pcount * 3)
				indexes.push_back(*cid);
		}
		cout << "\n Number of Non Zero Codes = " << nnzcodes << flush;
		cout << "\n Number of Codes Count > " << pcount << " = "
				<< indexes.size() << flush;
		return indexes.size();
	} else {
		UINT i = 0;
		for (vector<CodeCountType>::iterator iter = codecount.begin();
				iter != codecount.end(); ++iter, ++cid, ++i) {
			if (*iter > 0 || (*iter == 0 && ratio[i] != 0))
				nnzcodes++;
			if (*iter > pcount * 3 || (*iter == 0 && ratio[i] != 0)) /// to include the case when their exist a positive code but no negative one.
				indexes.push_back(*cid);
		}
		cout << "\n Number of Non Zero Codes = " << nnzcodes << flush;
		cout << "\n Number of Codes Count > " << pcount << " = "
				<< indexes.size() << flush;
		return indexes.size();
	}
}

#ifdef ONLINE_MAPPING
void ComputeCode::SaveCodeStatistics(const string & fname) {
	ofstream ofile(fname.c_str(), ios::out);
	if (!ofile) {
		cout << " Couldn't Create the File >>> " << fname << endl;
		exit(EXIT_FAILURE);
	}

	for (map<long, UINT>::iterator iter = tmpcodecache.begin();
			iter != tmpcodecache.end(); ++iter) {
		ofile << iter->first << " " << iter->second << endl;
	}
	ofile.close();
}
void ComputeCode::PruneCodes() {
	vector<long> relements;
	for (CodeTable::iterator iter = codecache.begin(); iter != codecache.end();
			++iter) {
		if (iter->second->codecount < pcount * 3)
			relements.push_back(iter->first);
	}
	cout << "\n Pruning Codes ........." << "\n Total Number of Codes = "
			<< codecache.size() << "\n Total Number of Codes having count < "
			<< pcount << "= " << relements.size();

	for (vector<long>::iterator iter = relements.begin();
			iter != relements.end(); ++iter) {
		codecache.erase(*iter);
	}
	cout << " \n Remaining Number of Codes = " << codecache.size();
}
void ComputeCode::GenerateCodeBook(const string & cifname,
		const string & cbfname, const ClusteringDistanceMetric &dmetric,
		UINT nclusrounds, REAL *initclusters, REAL *ratio) {
	//	using namespace Clustering;
	WriteCode(cifname);
	vector<UINT> nzindx;
	cout << " \n Before Clustering Size " << flush;
	DisplayMemoryUsage();
	//	SaveCodeStatistics("code-statistics");
	cout << "\n Number of Total Codes= " << codecache.size() << endl << flush;
	PruneCodes();
	UINT nnzcodes = codecache.size();
	assert(nnzcodes < MAX_NUMBER_CODES);
	REAL *data = new REAL[nnzcodes * codelen], *pdata = data;
	AllocateClustersMemory();
	uInt * ptcount = new uInt[nnzcodes], *assignment = new uInt[nnzcodes],
			lcount = 0;
	fill(assignment, assignment + nnzcodes, 0);
	vector<REAL> code(codelen, 0);
	for (CodeTable::iterator iter = codecache.begin(); iter != codecache.end();
			++iter) {
		GenerateCodeVector(iter->first, pdata);
		ptcount[lcount++] = iter->second->codecount;
		pdata += codelen;
	}

	if (initclusters == NULL) // Randomly Permute & Initialize
		InitializeClusters(data, nnzcodes, clusters);
	else
		// Initialize from the given centers...
		copy(initclusters, initclusters + ncenters * codelen, clusters);

	std::cout << "\nStarting Kmeans ..." << std::endl;
	std::cout << " ... with " << nnzcodes << " training points  of dimensions "
			<< codelen << std::endl << flush;
	std::cout << " ... for " << ncenters << " clusters " << std::endl;
	std::cout << " ... Doing " << nclusrounds << " Rounds" << std::endl;
	cout << " \n Memory Usage =";
	DisplayMemoryUsage();
	cout << " .... Running Kmeans \n" << flush;
	long initime = time(0);
	REAL sse = kmeans(clusters, data, assignment, codelen, nnzcodes, ncenters,
			10E3, nclusrounds, ptcount);

	cout << " \n Time Take By KMeans =" << time(0) - initime << endl << flush;
	std::cout << "Done!" << std::endl;
	std::cout << "Sum of Squared Error : " << sse << std::endl;
	initime = time(0);
	std::cout << "Mapping Codes to Clusters ... " << std::endl << flush;
	MapCodes2Centers();
	std::cout << "Time Taken By Mapping =... " << time(0) - initime << endl
			<< flush;
	//#endif
	std::cout << "Writing to File... " << std::endl << flush;
	initime = time(0);
	WriteMatrix2File("Clusters", clusters, ncenters, codelen);
	WriteCode(cbfname);
	std::cout << "Time Taken By File Writing =... " << time(0) - initime << endl
			<< flush;
	cout << " \n After Clearing Size = ";
	DisplayMemoryUsage();
	delete[] data;
	delete[] assignment;
	delete[] ptcount;
	ClearCodeCount();
}
#else
void ComputeCode::GenerateCodeBook(const string & cifname,
		const string & cbfname, const ClusteringDistanceMetric &dmetric,
		UINT nclusrounds, REAL *initclusters, REAL *ratio) {
//	using namespace Clustering;
	cout << " With PTR Generating CodeVectors" << endl << flush;
	char *codevec;
	GenerateCodeVectors(codevec);
	WriteCode(cifname);
	vector<UINT> nzindx;
	cout << " \n Before Clustering Size " << flush;
	//	DisplayMemoryUsage();
	UINT nnzcodes = FindIndexNonZeroCodes(nzindx, ratio);
	REAL * data = new REAL[nnzcodes * codelen];
	AllocateClustersMemory();
	WriteClusters(cifname + "-v2");
	uInt * ptcount = new uInt[nnzcodes], *assignment = new uInt[nnzcodes];
	//	fill(clusters, clusters + ncenters * neighbours, 0);
	fill(assignment, assignment + nnzcodes, 0);
	CopyData(codevec, nzindx, data, ptcount, ratio);
	delete[] codevec;

	if (initclusters == NULL)// Randomly Permute & Initialize
	InitializeClusters(data, nnzcodes, clusters);
	else
	// Initialize from the given centers...
	copy(initclusters, initclusters + ncenters * codelen, clusters);

	std::cout << "\nStarting Kmeans ..." << std::endl;
	std::cout << " ... with " << nnzcodes << " training points  of dimensions "
	<< codelen << std::endl << flush;
	std::cout << " ... for " << ncenters << " clusters " << std::endl;
	std::cout << " ... Doing " << nclusrounds << " Rounds" << std::endl;
	cout << " \n Memory Usage =";
	//	DisplayMemoryUsage();
	cout << " .... Running Kmeans \n" << flush;
	long initime = time(0);
	REAL sse = kmeans(clusters, data, assignment, codelen, nnzcodes, ncenters,
			10E3, nclusrounds, ptcount);

	cout << " \n Time Take By KMeans =" << time(0) - initime << endl << flush;
	std::cout << "Done!" << std::endl;
	std::cout << "Sum of Squared Error : " << sse << std::endl;
	initime = time(0);
	std::cout << "Mapping Codes to Clusters ... " << std::endl << flush;
	MapCodes2Centers(ratio);
	std::cout << "Time Taken By Mapping =... " << time(0) - initime << endl
	<< flush;
	//#endif
	std::cout << "Writing to File... " << std::endl << flush;
	initime = time(0);
	WriteMatrix2File("Clusters", clusters, ncenters, codelen);
	WriteVector2File("lookup.txt", &clusind[0], ncodes);
	WriteCode(cbfname);
	WriteClusters(cbfname + "-Clusters");
	std::cout << "Time Taken By File Writing =... " << time(0) - initime
	<< endl << flush;

	cout << " \n After Clearing Size = ";
	//	DisplayMemoryUsage();
	delete[] data;
	delete[] assignment;
	delete[] ptcount;

	ClearCodeCount();
}
#endif
//template<class T>
//void ComputeCode::ExtractClusterInfo(const Clustering::PointsSpace<T> &ps) {
//	// copy back the cluster mapping
//	using namespace Clustering;
//	const vector<DataPoints<T> > &chdata = ps.GetDataRef();
//	for (typename vector<DataPoints<T> >::const_iterator iter = chdata.begin(); iter
//			!= chdata.end(); ++iter)
//		clusind[iter->GetDataIndex()] = iter->GetCIndex();
//}
//template<class T>
//void ComputeCode::CopyData(vector<vector<CodeVecType> > &codevec,
//		Clustering::PointsSpace<T> &ps) {
//	using namespace Clustering;
//	ps.ResetCounter();
//	DataPoints<T> dpoint(codelen);
//	for (UINT i = 0; i < ncodes; ++i) {
//		if (codecount[i] > 0) {
//			copy(codevec[i].begin(), codevec[i].end(), dpoint.data.begin());
//			dpoint.ninstances = codecount[i];
//			dpoint.dindex = codeid[i];
//			dpoint.cindex = 0;
//			ps.AddData(dpoint);
//		}
//	}
//	cout << " Number of Non-Zeros Codes " << ps.GetNPoints() << flush;
//}
void ComputeCode::CopyData(vector<vector<CodeVecType> > &codevec,
		vector<UINT> &nzidx, REAL *data, uInt*ptcount) {
	REAL *ptr = data;
	UINT j = 0;
	for (UINT i = 0; i < ncodes; ++i) {
		if (codecount[i] > pcount * 3) {
			assert(codeid[i] == nzidx[j]);
			ptcount[j] = codecount[i];
			copy(codevec[i].begin(), codevec[i].end(), ptr);
			ptr += codelen;
			++j;
		}
	}
	cout << "\n Total Number of Codes having Count > " << pcount << " = " << j;
}
//void ComputeCode::CopyData(CodeVecType *codevec, vector<UINT> &nzidx,
//		REAL *data, uInt*ptcount, REAL *ratio) {
//	REAL *ptr = data;
//	CodeVecType *sptr = codevec;
//	UINT j = 0;
//	for (UINT i = 0; i < ncodes; ++i) {
//		if (codecount[i] > pcount * 3) {
//			assert(codeid[i] == nzidx[j]);
//			ptcount[j] = codecount[i];
//			copy(sptr, sptr + codelen, ptr);
//			ptr += codelen;
//			++j;
//		}
//		sptr += codelen;
//	}
//	cout << "\n Total Number of Codes having Count > " << pcount << " = " << j;
//}
void ComputeCode::CopyData(CodeVecType *codevec, vector<UINT> &nzidx,
		REAL *data, uInt*ptcount, REAL *ratio) {
	REAL *ptr = data;
	CodeVecType *sptr = codevec;
	UINT j = 0;
	for (UINT i = 0; i < ncodes; ++i) {
		if (codeid[i] == nzidx[j]) {
			ptcount[j] = codecount[i];
			copy(sptr, sptr + codelen, ptr);
			if (ratio != NULL)
				ptr[neighbours] = ratio[i];
			ptr += codelen;
			++j;
		}
		sptr += codelen;
	}
	cout << "\n Total Number of Codes having Count > " << pcount << " = " << j;
}
void ComputeCode::InitializeClusters(REAL *data, UINT nnzcodes,
		REAL *clusters) {
	vector<UINT> index(nnzcodes, 0);
	for (UINT i = 0; i < nnzcodes; ++i)
		index[i] = i;
	RandomPermute(index);
	REAL *init = clusters;
	for (UINT i = 0; i < ncenters; ++i) {
		copy(data + index[i] * codelen, data + (index[i] + 1) * codelen, init);
		init += codelen;
	}
}
void ComputeCode::InitializeClusters(vector<vector<CodeVecType> > &codevec,
		vector<UINT> &indeces, REAL *clusters) {
	REAL *ptr = clusters;
	for (vector<UINT>::iterator iter = indeces.begin(); iter != indeces.end();
			++iter) {
		//if (codecount[*iter] > 0) {
		vector<CodeVecType> & tcodevec = codevec[*iter];
		copy(tcodevec.begin(), tcodevec.end(), ptr);
		ptr += codelen;
		//				}
	}
}
void ComputeCode::AllocateMemory() {
	codecount.resize(ncodes, 0);
	clusind.resize(ncodes, 0);
	codeid.resize(ncodes, -1);
}
void ComputeCode::GenerateCodeVector(const long &value, vector<REAL>&codevec) {
	long var = value;
	codevec.resize(neighbours);
	fill(codevec.begin(), codevec.end(), 0);
#ifdef ONLINE_MAPPING
	for (int j = 0; j < neighbours; ++j) {
#else
		for (int j = neighbours - 1; j >= 0; --j) {

#endif
		codevec[j] = (var % lenvec - midp);
		var /= lenvec;
	}
}
void ComputeCode::GenerateCodeVector(const long &value, REAL *codevec) {
	long var = value;
	fill(codevec, codevec + neighbours, 0);
#ifdef ONLINE_MAPPING
	for (int j = 0; j < neighbours; ++j) {
#else
		for (int j = neighbours - 1; j >= 0; --j) {

#endif
		codevec[j] = (var % lenvec - midp);
		var /= lenvec;
	}
}
void ComputeCode::GenerateCodeVectors(vector<vector<CodeVecType> > &codevec) {
	// lenvec=base;
	codevec.resize(ncodes, vector<CodeVecType>(neighbours));

	for (UINT i = 0; i < ncodes; ++i) {
		UINT var = i;
		for (int j = neighbours - 1; j >= 0; --j) {
			codevec[i][j] = var % lenvec - midp;
			var /= lenvec;
		}
	}
	/*	InitializeMapping();
	 for (UINT i = 0; i < neighbours; ++i) {
	 long count = 0;
	 UINT ovar = pow((double) lenvec, (double) i);
	 for (long j = 0; j < ovar; ++j)
	 for (UINT k = 0; k < lenvec; ++k) {
	 UINT tvar = (UINT) (pow((double) lenvec,
	 (double) (neighbours - i - 1)));
	 for (UINT l = 0; l < tvar; ++l, ++count)
	 codevec[count][i] = mapping[k];
	 }
	 }*/
	UINT i = 0;
	for (vector<vector<CodeVecType> >::iterator iter = codevec.begin();
			iter != codevec.end(); ++iter, ++i)
		codeid[i] = GetCode(*iter);
}
void ComputeCode::GenerateCodeVectors(CodeVecType *& codevec) {
	long arraysize = (long) ncodes * (long) codelen;
	cout << "\n Array Size = " << arraysize << endl;
	codevec = new CodeVecType[arraysize];
	InitializeMapping();
	char *pcodevec = codevec;
	for (UINT i = 0; i < ncodes; ++i) {
		UINT var = i;
		for (int j = neighbours - 1; j >= 0; --j) {
			pcodevec[j] = var % lenvec - midp;
			var /= lenvec;
		}
		pcodevec += codelen;
	}
	//	for (UINT i = 0; i < neighbours; ++i) {
	//		long count = 0  ;
	//		pcodevec = codevec;
	//		cout << "i " << i << flush << endl;
	//		for (long j = 0; j < pow((double) lenvec, (double) i); ++j)
	//			for (UINT k = 0; k < lenvec; ++k) {
	//				UINT tvar = (UINT) (pow((double) lenvec,
	//						(double) (neighbours - i - 1)));
	//				for (UINT l = 0; l < tvar; ++l, ++count) {
	//					pcodevec[i] = mapping[k];
	//					pcodevec += neighbours;
	//				}
	//			}
	//	}
	UINT i = 0;
	cout << "\n Base to Decimal conversion = " << arraysize << endl;
	pcodevec = codevec;
	for (; i < ncodes; ++i) {
		codeid[i] = GetCode(pcodevec);
		pcodevec += codelen;
	}
}
void ComputeCodeLBP::ComputeInfo() {
	lenvec = (UINT) (pow((double) 2, (double) nlevels));
	ncodes = (UINT) (pow((double) lenvec, (double) neighbours));
	midp = 0;
	//		ComputeInfo();
	cout << " \n Total Number of Codes = " << ncodes << " Length of Vector = "
			<< neighbours << " Total Data Size ="
			<< (sizeof(clusind) + sizeof(codecount) + sizeof(codeid)) * ncodes
			<< endl << flush;
#ifndef ONLINE_MAPPING ///> use online mapping of codes,
	AllocateMemory();
#endif

}
long ComputeCodeLBP::GetCode(const vector<CodeVecType> & codevec) { // To Generate the code from the code vec
//		[0 1 2 -1 -2]
//	UINT codeval = 0;
//	for (int i = neighbours - 1, j = 0; i >= 0; --i, ++j)
//		codeval += (UINT) (pow((double) lenvec, (double) i) * codevec[j]);
//	return codeval;
	long codeval = 0;

#ifdef ONLINE_MAPPING
	for (int i = neighbours - 1; i >= 0; --i)

#else
		for (int i = 0; i < neighbours; ++i)

#endif
		codeval = codeval * lenvec + codevec[i]; // scaling of code
	return codeval;
}
long ComputeCodeLBP::GetCode(const CodeVecType *codevec) { // To Generate the code from the code vec
//		[0 1 2 -1 -2]
	long codeval = 0;
#ifdef ONLINE_MAPPING
	for (int i = neighbours - 1; i >= 0; --i)
#else
		for (int i = 0; i < neighbours; ++i)
#endif
		codeval = codeval * lenvec + codevec[i]; // scaling of code
	return codeval;
}

void ComputeCodeLTP::ComputeInfo() {
	lenvec = (nlevels * 2 + 1);
	midp = nlevels; // translation coefficient...
	ncodes = (UINT) (pow((double) lenvec, (double) neighbours));
	cout << " \n Total Number of Codes = " << ncodes << " Lenght of Vector = "
			<< neighbours << " Total Data Size ="
			<< (sizeof(clusind) + sizeof(codecount) + sizeof(codeid)) * ncodes
			<< endl << flush;

#ifndef ONLINE_MAPPING ///> use online mapping of codes,
	AllocateMemory();
#endif
}
long ComputeCodeLTP::GetCode(const vector<CodeVecType> & codevec) { // To Generate the code from the code vec
	long codeval = 0;
	/*for (int i = neighbours - 1, j = 0; i >= 0; --i, ++j)
	 codeval += (UINT) (pow((double) lenvec, (double) i)
	 * (codevec[j] >= 0 ? codevec[j] : nlevels + ABS(
	 codevec[j])));*/
#ifdef ONLINE_MAPPING
	for (int i = neighbours - 1; i >= 0; --i)
#else
		for (int i = 0; i < neighbours; ++i)
#endif
		codeval = codeval * lenvec + (codevec[i] + midp);

	return codeval;
}
long ComputeCodeLTP::GetCode(const CodeVecType *codevec) { // To Generate the code from the code vec
	long codeval = 0;
	/*for (int i = neighbours - 1, j = 0; i >= 0; --i, ++j)
	 codeval += (UINT) (pow((double) lenvec, (double) i)
	 * (codevec[j] >= 0 ? codevec[j] : nlevels + ABS(
	 codevec[j])));*/
#ifdef ONLINE_MAPPING
	for (int i = neighbours - 1; i >= 0; --i)
#else
		for (int i = 0; i < neighbours; ++i)
#endif
		codeval = codeval * lenvec + (codevec[i] + midp); // scaling of code
	return codeval;
}
#endif
