/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/
#include "classInfo.h"
bool ClassInfo::GetNextImageAnnotation(string &iname, vector<Coord>& vcoord) {
	if (counter < iinfo.size()) {
		iname = iinfo[counter].iname;
		vcoord = iinfo[counter++].coord;
		return true;
	}
	return false;
}

bool ClassInfo::GetImageAnnotation(string &iname, vector<Coord>& vcoord) {
	for (vector<ImageInfo>::iterator iter = iinfo.begin(); iter != iinfo.end(); ++iter)
		if (iname == iter->iname) {
			iname = iter->iname;
			vcoord = iter->coord;
			return true;
		}
	return false;
}

void ClassInfo::ReadInfo(const string & fname) {
	ifstream ifile(fname.c_str(), ios::in | ios::binary);
	if (!ifile) {
		cout << " Error: Couldn't read the Annotation Information : " << fname;
		cout << endl << " Quittttingggggggggggg. ... " << endl;
		exit(EXIT_FAILURE);
	}
	UINT niobjects = 0, snstr, version, ncomp, icounter = 0, x1, y1, x2, y2,
			nobjects, cid;// number of image objects
	ifile.read((char*) &version, sizeof(UINT));
	/*	ifile.read((char*)&width,sizeof(UINT));
	 ifile.read((char*)&height,sizeof(UINT));
	 ifile.read((char*)&nobjects,sizeof(UINT));

	 ifile.read((char*)&snstr,sizeof(UINT));
	 char *nstr = new char[snstr+1];
	 ocounter=0;
	 while(ifile && !ifile.eof() && ocounter < nobjects)
	 {

	 ifile.read((char*)nstr,sizeof(char)*snstr);
	 nstr[snstr]='\0';
	 iinfo.push_back(ImageInfo(nstr));
	 ifile.read((char*)&niobjects,sizeof(UINT));
	 iinfo[icounter].coord.resize(niobjects);
	 ocounter+=niobjects;
	 for(int i=0; i < niobjects;++i)
	 {
	 ifile.read((char*)&x1,sizeof(UINT));
	 ifile.read((char*)&y1,sizeof(UINT));
	 ifile.read((char*)&x2,sizeof(UINT));
	 ifile.read((char*)&y2,sizeof(UINT));
	 iinfo[icounter].coord[i].SetCoord(x1,y1,x2,y2);
	 }
	 icounter++;
	 }
	 cout<< "\n Number of Pos Images read = "<< icounter << " Number of Objects ="
	 << ocounter<<endl;
	 */
	ifile.read((char*) &ncomp, sizeof(UINT));
	ifile.read((char*) &cid, sizeof(UINT));// Component Id
	ifile.read((char*) &nobjects, sizeof(UINT));
	ifile.read((char*) &width, sizeof(UINT));
	ifile.read((char*) &height, sizeof(UINT));

	char nstr[255];
	ocounter = 0;
	while (ifile && !ifile.eof() && ocounter < nobjects) {
		ifile.read((char*) &snstr, sizeof(UINT));
		ifile.read((char*) nstr, sizeof(char) * snstr);
		nstr[snstr] = '\0';
		iinfo.push_back(ImageInfo(nstr));
		ifile.read((char*) &niobjects, sizeof(UINT));
		iinfo[icounter].coord.resize(niobjects);
		ocounter += niobjects;
		for (int i = 0; i < niobjects; ++i) {
			ifile.read((char*) &x1, sizeof(UINT));
			ifile.read((char*) &y1, sizeof(UINT));
			ifile.read((char*) &x2, sizeof(UINT));
			ifile.read((char*) &y2, sizeof(UINT));
			iinfo[icounter].coord[i].SetCoord(x1, y1, x2, y2);
		}
		icounter++;
	}
	cout << "\n Component # = " << cid + 1 << " Number of Pos Images read = "
			<< icounter << " Number of Objects =" << ocounter << endl;
}

/*int main()
 {
 ClassInfo cobj("/local_scratch2/shussain/datasets/VOC2006/VOC2006/person_pos.lst");
 }*/
