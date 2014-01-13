/*
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
*/
#include "util.hpp"

void ReadDir(const string &dname, const string &offset, const string &ofname)
// offset to add before the name of each file in the directory...
{

	ofstream ofile(ofname.c_str(), ios::out);
	DIR *fs = opendir(dname.c_str()); // Give Path directory to read
	if (fs == NULL) {
		cout << "Could not open dir \n" << endl;
		return;
	}

	struct dirent *dirstr;

	while ((dirstr = readdir(fs)) != NULL) {

		if (strncmp(dirstr->d_name, ".", 1) != 0 && strncmp(dirstr->d_name,
				"..", 2) != 0) { // read all the files except which start with . or ..
			ofile << offset + dirstr->d_name << endl;
			//			cout<< offset+dirstr->d_name<<endl;
		}

	}
}

void ReadDir(const string &dname, const string &ofname) {

	ofstream ofile(ofname.c_str(), ios::out);
	DIR *fs = opendir(dname.c_str()); // Give Path directory to read
	if (fs == NULL) {
		cout << "Could not open dir \n" << endl;
		return;
	}

	struct dirent *dirstr;

	while ((dirstr = readdir(fs)) != NULL) {

		if (strncmp(dirstr->d_name, ".", 1) != 0 && strncmp(dirstr->d_name,
				"..", 2) != 0) { // read all the files except which start with . or ..
			ofile << dirstr->d_name << endl;
			//			cout<< dirstr->d_name<<endl;
		}

	}
}
void ReadDir(const string &dname, vector<string> &lname) {

	DIR *fs = opendir(dname.c_str()); // Give Path directory to read
	if (fs == NULL) {
		cout << "Could not open dir \n" << endl;
		return;
	}

	struct dirent *dirstr;

	while ((dirstr = readdir(fs)) != NULL) {

		if (strncmp(dirstr->d_name, ".", 1) != 0 && strncmp(dirstr->d_name,
				"..", 2) != 0) { // read all the files except which start with . or ..
			lname.push_back(dirstr->d_name);

		}

	}
}
void ReadImage(UCHAR *mat, const Image & image, UINT offset) {
	int m = image.rows(), n = image.columns();
	ColorGray mycolor;
	unsigned int count = 0;
	for (UINT i = offset; i < m - offset; i++) {

		for (UINT j = offset; j < n - offset; j++) {
			mycolor = image.pixelColor(j, i);
			mat[count++] = (unsigned char) (mycolor.shade() * 255);
		}
	}
} //`Magick++-config --cppflags --cxxflags --ldflags --libs`
void ReadImage(const Image &image, vector<REAL> &pixval, REAL gnorm) {
	UINT m = image.rows(), n = image.columns();
	ColorGray mycolor;
	UINT count = 0;
	for (UINT i = 0; i < m; ++i) {
		for (UINT j = 0; j < n; ++j) {
			mycolor = image.pixelColor(j, i);
			pixval[count++] = gnorm == 1 ? mycolor.shade() : pow(
					(double) mycolor.shade(), (double) gnorm);
		}
	}
}
void Scale(Image &simage, UINT rows, UINT cols, bool aspect) {
	int w = (UINT) (2 * round(1 / simage.columns() * cols) + 1), // width of gaussian mask
			h = (UINT) (2 * round(1 / simage.rows() * rows) + 1); // height of gaussian mask
	int size = w < h ? w : h;
	float sigma = (size / 2 - 1) * 0.3 + 0.8;
	Geometry gobj = Geometry(cols, rows, 0, 0);
	gobj.aspect(aspect);
	simage.gaussianBlur(size, sigma);
	simage.scale(gobj);
}
void Scale(Image &simage, REAL scale) {
	Scale(simage, (UINT) (simage.rows() / scale),
			(UINT) (simage.columns() / scale), true);
	//	Scale(simage,simage.rows()/scale,simage.columns()/scale,true);
	//	int w = (UINT)(2*round(1/simage.columns()*cols) + 1),     // width of gaussian mask
	//	h = (UINT)(2*round(1/simage.rows()*rows ) + 1);       // height of gaussian mask
	//	int size=w < h ? w:h;
	//	float sigma = (size/2 - 1)*0.3 + 0.8 ;
	//	Geometry gobj=Geometry(cols ,rows, 0 , 0);
	//	gobj.aspect(aspect);
	//	simage.gaussianBlur(size,sigma);
	//	simage.scale(gobj);
}
string int2str(int input) {

	ostringstream oss;
	oss << input;
	return oss.str();
}
void CopyConstBoundaries(Image& input, UINT size) { // Copy the Constant Boundaries as Image border
	UINT icols = input.columns(), irows = input.rows(), ocols = icols + 2
			* size, orows = irows + 2 * size;
	Image output(Geometry(ocols, orows, 0, 0), "Black");
	output.modifyImage();
	output.type(input.type());
	PixelPacket *ocache = output.getPixels(0, 0, ocols, orows), *tcache1,
			*tcache2;
	PixelPacket *icache = input.getPixels(0, 0, icols, irows);

	// copy the center image
	for (UINT i = 0; i < irows; ++i)
		for (UINT j = 0; j < icols; ++j)
			*(ocache + (i + size) * ocols + (j + size)) = *(icache + i * icols
					+ j);

	// Copy to Right & Left
	for (UINT i = 0; i < irows; ++i) {
		tcache1 = ocache + (i + size) * ocols + size + icols;
		tcache2 = ocache + (i + size) * ocols;
		const PixelPacket & ref1 = *(icache + (i + 1) * icols - 1), &ref2 =
				*(icache + i * icols);
		for (int j = 0; j < size; ++j) {
			*(tcache1 + j) = ref1;
			*(tcache2 + j) = ref2;
		}
	}

	PixelPacket *icache1 = (icache + (irows - 1) * icols), *icache2 = icache;
	// Copy to above & Below
	for (UINT i = 0; i < size; ++i) {
		tcache1 = ocache + (i + irows + size) * ocols + size; // below
		tcache2 = ocache + i * ocols + size;

		for (int j = 0; j < icols; ++j) {
			*(tcache1 + j) = *(icache1 + j);
			*(tcache2 + j) = *(icache2 + j);
		}
	}
	// Copy at the Four corners the four pixels....
	// four pixel values
	const PixelPacket & pix1 = *icache, &pix2 = *(icache + icols), &pix3 =
			*(icache + (irows - 1) * icols), &pix4 = *(icache + irows * icols
			- 1);

	for (UINT i = 0; i < size; ++i) {
		icache1 = (ocache + i * ocols);
		icache2 = (ocache + i * ocols + icols + size);
		tcache1 = (ocache + (irows + size + i) * ocols);
		tcache2 = (ocache + (irows + size + i) * ocols + icols + size);
		for (int j = 0; j < size; ++j) {
			*(icache1 + j) = pix1;
			*(icache2 + j) = pix2;
			*(tcache1 + j) = pix3;
			*(tcache2 + j) = pix4;
		}
	}
	output.syncPixels();
	input = output;
}

void CopyBoundaries(Image& input, Image& output, char dir) {
	output.modifyImage();
	output.type(TrueColorType);
	PixelPacket *ocache = output.getPixels(0, 0, output.columns(),
			output.rows());
	const PixelPacket *icache = input.getConstPixels(0, 0, input.columns(),
			input.rows());
	int dx = output.columns() - input.columns();
	int dy = output.rows() - input.rows();
	switch (dir) {
	case 'r':
		for (UINT i = 0; i < input.rows(); ++i)
			for (UINT j = 0; j < input.columns(); ++j)
				*(ocache + i * output.columns() + j) = *(icache + i
						* input.columns() + j);
		for (UINT i = 0; i < input.rows(); ++i)
			for (int j = 0; j < dx; ++j)
				*(ocache + input.columns() + i * output.columns() + j)
						= *(icache + input.columns() + i * input.columns() - 1);
		break;
	case 'l':
		for (UINT i = 0; i < input.rows(); ++i)
			for (UINT j = 0; j < input.columns(); ++j)
				*(ocache + i * output.columns() + j + dx) = *(icache + i
						* input.columns() + j);
		for (UINT i = 0; i < input.rows(); ++i)
			for (int j = 0; j < dx; ++j)
				*(ocache + i * output.columns() + dx - j - 1) = *(icache + i
						* input.columns());
		break;
	case 'd':
		for (UINT i = 0; i < input.rows(); ++i)
			for (UINT j = 0; j < input.columns(); ++j)
				*(ocache + i * output.columns() + j) = *(icache + i
						* input.columns() + j);
		for (int i = 0; i < dy; ++i)
			for (UINT j = 0; j < input.columns(); ++j)
				*(ocache + (input.rows() + i) * output.columns() + j)
						= *(icache + (input.rows() - 1) * input.columns() + j);
		break;
	case 'u':
		for (UINT i = 0; i < input.rows(); ++i)
			for (UINT j = 0; j < input.columns(); ++j)
				*(ocache + (i + dy) * output.columns() + j) = *(icache + i
						* input.columns() + j);
		for (int i = 0; i < dy; ++i)
			for (UINT j = 0; j < input.columns(); ++j)
				*(ocache + (dy - i - 1) * output.columns() + j) = *(icache + j);
		break;
	case 't': // copy top in the top left corner
		for (UINT i = 0; i < input.rows(); ++i) // copy the image in the lower half
			for (UINT j = 0; j < input.columns(); ++j)
				*(ocache + (i + dy) * output.columns() + j + dx) = *(icache + i
						* input.columns() + j);
		for (int i = 0; i < dy; ++i)
			for (UINT j = 0; j < input.columns(); ++j)
				*(ocache + (dy - i - 1) * output.columns() + j + dx) = *(icache
						+ j);
		for (UINT i = 0; i < output.rows(); ++i)
			for (int j = 0; j < dx; ++j)
				*(ocache + i * output.columns() + dx - j - 1) = *(ocache + i
						* output.columns() + dx);
		break;
	case 'b': // copy bottom in the bottom right corner
		for (UINT i = 0; i < input.rows(); ++i)
			for (UINT j = 0; j < input.columns(); ++j)
				*(ocache + i * output.columns() + j) = *(icache + i
						* input.columns() + j);
		for (int i = 0; i < dy; ++i)
			for (UINT j = 0; j < input.columns(); ++j)
				*(ocache + (input.rows() + i) * output.columns() + j)
						= *(icache + (input.rows() - 1) * input.columns() + j);
		for (UINT i = 0; i < output.rows(); ++i)
			for (int j = 0; j < dx; ++j)
				*(ocache + input.columns() + i * output.columns() + j)
						= *(ocache + input.columns() + i * output.columns() - 1);
		break;
	case 'a': // copy above in the top right corner
		for (UINT i = 0; i < input.rows(); ++i)
			for (UINT j = 0; j < input.columns(); ++j)
				*(ocache + (i + dy) * output.columns() + j) = *(icache + i
						* input.columns() + j);
		for (int i = 0; i < dy; ++i)
			for (UINT j = 0; j < input.columns(); ++j)
				*(ocache + (dy - i - 1) * output.columns() + j) = *(icache + j);
		for (UINT i = 0; i < output.rows(); ++i)
			for (int j = 0; j < dx; ++j)
				*(ocache + input.columns() + i * output.columns() + j)
						= *(ocache + input.columns() + i * output.columns() - 1);
		break;
	case 's': // copy in the left bottom corner
		for (UINT i = 0; i < input.rows(); ++i) // copy the image in the upper right half
			for (UINT j = 0; j < input.columns(); ++j)
				*(ocache + i * output.columns() + j + dx) = *(icache + i
						* input.columns() + j);

		for (int i = 0; i < dy; ++i)
			for (UINT j = 0; j < input.columns(); ++j)
				*(ocache + (input.rows() + i) * output.columns() + j)
						= *(icache + (input.rows() - 1) * input.columns() + j);

		for (UINT i = 0; i < output.rows(); ++i)
			for (int j = 0; j < dx; ++j)
				*(ocache + i * output.columns() + dx - j - 1) = *(ocache + i
						* output.columns() + dx);
		break;
	default:
		cout << " Invalid Argument : Returning .... " << endl;
		return;
	}
	output.syncPixels();
	//	output.write("temp.png");
}
void Padding(Image& input, Image& output, UINT size, char dir) {
	output.modifyImage();
	output.type(TrueColorType);
	PixelPacket *ocache = output.getPixels(0, 0, output.columns(),
			output.rows());
	const PixelPacket *icache = input.getConstPixels(0, 0, input.columns(),
			input.rows());
	switch (dir) {
	case 'r':
		for (UINT i = 0; i < input.rows(); ++i)
			for (UINT j = 0; j < input.columns(); ++j)
				*(ocache + i * output.columns() + j) = *(icache + i
						* input.columns() + j);
		break;
	case 'l':
		for (UINT i = 0; i < input.rows(); ++i)
			for (UINT j = 0; j < input.columns(); ++j)
				*(ocache + i * output.columns() + j + size) = *(icache + i
						* input.columns() + j);
		break;
	case 'd':
		for (UINT i = 0; i < input.rows(); ++i)
			for (UINT j = 0; j < input.columns(); ++j)
				*(ocache + i * output.columns() + j) = *(icache + i
						* input.columns() + j);
		break;
	case 'u':
		for (UINT i = 0; i < input.rows(); ++i)
			for (UINT j = 0; j < input.columns(); ++j)
				*(ocache + (i + size) * output.columns() + j) = *(icache + i
						* input.columns() + j);
		break;
	default:
		cout << " Invalid Argument : Returning .... " << endl;
		return;
	}
	output.syncPixels();
}

void Reflect(Image& input, Image& output, UINT size, char dir) {
	output.modifyImage();
	output.type(TrueColorType);
	PixelPacket *ocache = output.getPixels(0, 0, output.columns(),
			output.rows());
	const PixelPacket *icache = input.getConstPixels(0, 0, input.columns(),
			input.rows());
	switch (dir) {
	case 'r':
		for (UINT i = 0; i < input.rows(); ++i)
			for (UINT j = 0; j < input.columns(); ++j)
				*(ocache + i * output.columns() + j) = *(icache + i
						* input.columns() + j);
		for (UINT i = 0; i < input.rows(); ++i)
			for (UINT j = 0; j < size; ++j)
				*(ocache + input.columns() + i * output.columns() + j)
						= *(icache + input.columns() + i * input.columns() - j
								- 1);
		break;
	case 'l':
		for (UINT i = 0; i < input.rows(); ++i)
			for (UINT j = 0; j < input.columns(); ++j)
				*(ocache + i * output.columns() + j + size) = *(icache + i
						* input.columns() + j);
		for (UINT i = 0; i < input.rows(); ++i)
			for (UINT j = 0; j < size; ++j)
				*(ocache + i * output.columns() + size - j - 1) = *(icache + i
						* input.columns() + j);
		break;
	case 'd':
		for (UINT i = 0; i < input.rows(); ++i)
			for (UINT j = 0; j < input.columns(); ++j)
				*(ocache + i * output.columns() + j) = *(icache + i
						* input.columns() + j);
		for (UINT i = 0; i < size; ++i)
			for (UINT j = 0; j < input.columns(); ++j)
				*(ocache + (input.rows() + i) * output.columns() + j)
						= *(icache + (input.rows() - i - 1) * input.columns()
								+ j);
		break;
	case 'u':
		for (UINT i = 0; i < input.rows(); ++i)
			for (UINT j = 0; j < input.columns(); ++j)
				*(ocache + (i + size) * output.columns() + j) = *(icache + i
						* input.columns() + j);
		for (UINT i = 0; i < size; ++i)
			for (UINT j = 0; j < input.columns(); ++j)
				*(ocache + (size - i - 1) * output.columns() + j) = *(icache
						+ i * input.columns() + j);
		break;
	default:
		cout << " Invalid Argument : Returning .... " << endl;
		return;
	}
	output.syncPixels();
}
void ReflectImage(Image& output, UINT size, char dir) {
	Image input = output;
	output.modifyImage();
	//	output.type(TrueColorType);
	PixelPacket *ocache;
	const PixelPacket *icache = input.getConstPixels(0, 0, input.columns(),
			input.rows());
	switch (dir) {
	case 'r':
		output.size(Geometry(input.columns() + size, input.rows()));
		ocache = output.getPixels(0, 0, output.columns(), output.rows());
		for (UINT i = 0; i < input.rows(); ++i)
			for (UINT j = 0; j < input.columns(); ++j)
				*(ocache + i * output.columns() + j) = *(icache + i
						* input.columns() + j);
		for (UINT i = 0; i < input.rows(); ++i)
			for (UINT j = 0; j < size; ++j)
				*(ocache + input.columns() + i * output.columns() + j)
						= *(icache + input.columns() + i * input.columns() - 1);
		break;
	case 'l':
		output.size(Geometry(input.columns() + size, input.rows()));
		ocache = output.getPixels(0, 0, output.columns(), output.rows());
		for (UINT i = 0; i < input.rows(); ++i)
			for (UINT j = 0; j < input.columns(); ++j)
				*(ocache + i * output.columns() + j + size) = *(icache + i
						* input.columns() + j);
		for (UINT i = 0; i < input.rows(); ++i)
			for (UINT j = 0; j < size; ++j)
				*(ocache + i * output.columns() + size - j - 1) = *(icache + i
						* input.columns());
		break;
	case 'd':
		output.size(Geometry(input.columns(), input.rows() + size));
		ocache = output.getPixels(0, 0, output.columns(), output.rows());
		for (UINT i = 0; i < input.rows(); ++i)
			for (UINT j = 0; j < input.columns(); ++j)
				*(ocache + i * output.columns() + j) = *(icache + i
						* input.columns() + j);
		for (UINT i = 0; i < size; ++i)
			for (UINT j = 0; j < input.columns(); ++j)
				*(ocache + (input.rows() + i) * output.columns() + j)
						= *(icache + (input.rows() - 1) * input.columns() + j);
		break;
	case 'u':
		output.size(Geometry(input.columns(), input.rows() + size));
		ocache = output.getPixels(0, 0, output.columns(), output.rows());
		for (UINT i = 0; i < input.rows(); ++i)
			for (UINT j = 0; j < input.columns(); ++j)
				*(ocache + (i + size) * output.columns() + j) = *(icache + i
						* input.columns() + j);
		for (UINT i = 0; i < size; ++i)
			for (UINT j = 0; j < input.columns(); ++j)
				*(ocache + (size - i - 1) * output.columns() + j) = *(icache
						+ j);
		break;

	default:
		cout << " Invalid Argument : Returning .... " << endl;
		return;
	}
	output.syncPixels();
}
void BilinearInterpolation(vector<REAL>& input, UINT cols, UINT rows,
		REAL scale, vector<float>& output)
// Performs the binlinear interpolation on the given input 2D image...
{
	UINT nrows = (UINT) (rows * scale), ncols = (UINT) (cols * scale), x0, y0,
			x1, y1;
	output.resize(nrows * ncols, 0); //
	REAL ii, ij, dx, dy;

	for (UINT i = 0; i < nrows; ++i) {
		ii = i / scale;
		dy = ii - floor(ii);
		y0 = (UINT) floor(ii);
		y1 = (UINT) ceil(ii); // y0 = ii, y1 = y0+1;
		// so copy the original pixels at boundary...
		y1 = y1 > rows - 1 ? rows - 1 : y1;

		for (UINT j = 0; j < ncols; ++j) {
			ij = j / scale;
			dx = ij - floor(ij);
			x0 = (UINT) floor(ij);
			x1 = (UINT) ceil(ij); // x0 = ij, x1 = x0+1;
			x1 = x1 > cols - 1 ? cols - 1 : x1;
			output[j + i * ncols] = (1 - dx) * (1 - dy) * input[x0 + y0 * cols]
					+ (1 - dx) * dy * input[x0 + y1 * cols] + (1 - dy) * dx
					* input[x1 + y0 * cols] + dx * dy * input[x1 + y1 * cols];
		}
	}
}
void BilinearInterpolation(vector<REAL>& input, UINT cols, UINT rows, UINT dim,
		REAL scale, vector<REAL>& output)
// Performs the binlinear interpolation on the given input 2D image...
{
	UINT nrows = (UINT) (rows * scale), ncols = (UINT) (cols * scale), x0, y0,
			x1, y1;
	output.resize(nrows * ncols * dim, 0); //
	REAL ii, ij, dx, dy;

	for (UINT k = 0; k < dim; ++k) {
		for (UINT i = 0; i < nrows; ++i) {
			ii = i / scale;
			dy = ii - floor(ii);
			y0 = (UINT) floor(ii);
			y1 = (UINT) ceil(ii); // y0 = ii, y1 = y0+1;
			// so copy the original pixels at boundary...
			y1 = y1 > rows - 1 ? rows - 1 : y1;

			for (UINT j = 0; j < ncols; ++j) {
				ij = j / scale;
				dx = ij - floor(ij);
				x0 = (UINT) floor(ij);
				x1 = (UINT) ceil(ij); // x0 = ij, x1 = x0+1;
				x1 = x1 > cols - 1 ? cols - 1 : x1;
				output[dim * j + k + i * ncols * dim] = (1 - dx) * (1 - dy)
						* input[dim * x0 + k + y0 * cols * dim] + (1 - dx) * dy
						* input[dim * x0 + k + y1 * cols * dim] + (1 - dy) * dx
						* input[dim * x1 + k + y0 * cols * dim] + dx * dy
						* input[dim * x1 + k + y1 * cols];
			}
		}
	}
}

void ReadFile(const string &nfile, vector<string>&vstr) {
	ifstream ifile(nfile.c_str());
	if (!ifile) {
		cout << " Couldn't Open the Input File" << nfile;
		exit(-1);
	}
	string str;
	while (ifile >> str && ifile)
		vstr.push_back(str);
	ifile.close();

}
void Compare(const string& nfile1, const string& nfile2, vector<string>&diff) {
	vector<string> vstr1, vstr2;
	// Read File
	ReadFile(nfile1, vstr1);
	ReadFile(nfile2, vstr2);
	ofstream ofile("diff.txt");
	if (vstr2.size() > vstr1.size()) {
		for (UINT i = 0; i < vstr2.size(); ++i) {
			if (find(vstr1.begin(), vstr1.end(), vstr2[i]) == vstr1.end()) {
				ofile << vstr2[i] << endl;
				diff.push_back(vstr2[i]);
			}
		}
	} else {
		for (UINT i = 0; i < vstr1.size(); ++i) {
			if (find(vstr2.begin(), vstr2.end(), vstr1[i]) == vstr2.end()) {
				ofile << vstr1[i] << endl;
				diff.push_back(vstr1[i]);
			}
		}

	}

}
//template class<T>
//T FindMaximum(vector<T>& ivec,UINT &index)
//{
//	index=0;
//	T maxval=ivec[0];
//	UINT k=0;
//	for(vector<T>::iterator iter=ivec.begin(); iter!=ivec.end();++iter,++k)
//			if( *iter > maxval)
//			{
//				index = k;
//				maxval = *iter;
//			}
//	return maxval;
//}
//void Convolve(Image &image,vector<float> &filter,char dir)
//{
//	image.modifyImage();
//	PixelPacket *cache=image.getPixels(0,0,input.columns(),input.rows());
//	if(dir == 'h')
//	{
//		Image oimage(Geometry(input.columns()+filter.size()/2,input.rows()),"black");
//		for(int i=filter.size()/2; i < image.rows();++i)
//			for(int j=0; j < image.columns();++j)
//			{
//				if(i < filter.size()/2 || i > image)
//			}
//	}
//	else
//	{
//
//	}
//	image.syncPixels();
//}

void DrawRectangle(Image &image, int startx, int starty, int width, int height,
		const string &col) {

	PixelPacket *pc = image.getPixels(0, 0, image.columns(), image.rows());
	int columns = image.columns();
	int rows = image.rows();
	Color cobj(col.c_str());
	int endx = MIN(width+startx,image.columns()), endy =
			MIN(height+starty,image.rows());
	//	startx = (width + startx) > columns-1 ? (startx - (width+startx - columns-2)) : startx;
	//	starty = (height + starty) > rows-1 ? (starty - (height+starty - rows-2)) : starty;
	for (int i = startx; i < endx; i++) {
		*(pc + starty * columns + i) = cobj;
		*(pc + (starty + 1) * columns + i) = cobj;
		*(pc + (starty + height - 1) * columns + i) = cobj;
		*(pc + (starty + height) * columns + i + 1) = cobj;
	}
	for (int i = starty; i < endy; i++) {
		*(pc + i * columns + startx) = cobj;
		*(pc + i * columns + startx + 1) = cobj;
		*(pc + i * columns + (startx + width - 1)) = cobj;
		*(pc + i * columns + (startx + width)) = cobj;
	}
	image.syncPixels();
}

void CropImage(Image &image, int xmin, int ymin, int width, int height) {
	if (xmin < 0) {
		ReflectImage(image, ABS(xmin), 'l');
		xmin = 0;
	}
	if (ymin < 0) {
		ReflectImage(image, ABS(ymin), 'u');
		ymin = 0;
	}

	if (xmin + width >= image.columns())
		ReflectImage(image, ABS( xmin+width - image.columns()+1), 'r');

	if (ymin + height >= image.rows())
		ReflectImage(image, ABS(ymin+height - image.rows()+1), 'd');

	image.crop(Geometry(width, height, xmin, ymin));
}
void CropImage2(Image &image, int xmin, int ymin, int xmax, int ymax) {
	if (xmin < 0) {
		ReflectImage(image, ABS(xmin), 'l');
		xmin = 0;
	}
	if (ymin < 0) {
		ReflectImage(image, ABS(ymin), 'u');
		ymin = 0;
	}

	if (xmax >= image.columns())
		ReflectImage(image, ABS( xmax - image.columns()+1), 'r');

	if (ymax >= image.rows())
		ReflectImage(image, ABS(ymax - image.rows()+1), 'd');

	image.crop(Geometry(xmax - xmin, ymax - ymin, xmin, ymin));
}
unsigned long getActiveMemory() {
	std::string token;
	std::ifstream file("/proc/meminfo");
	while (file >> token) {
		if (token == "Active:") {
			unsigned long mem;
			if (file >> mem) {
				return mem;
			} else {
				return 0;
			}
		}
		// ignore rest of the line
		file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}
	return 0; // nothing found
}

UINT GetNumberLines(const string & fname) {
	ifstream ifile(fname.c_str());

	if (!ifile) {
		cout << "Error Couldn't Open the file for Reading ... " << fname
				<< endl;
		exit(EXIT_FAILURE);
	}
	UINT lcount = 0;
	string str;
	while (ifile >> str && ifile)
		++lcount;

	return lcount;
}
void ResizeImage(Image & image, UINT width, UINT height, FilterType ftype,
		bool exact) {
	if (image.rows() <= width && image.columns() <= height)
		return;
	Geometry gobj = Geometry(width, height, 0, 0);
	gobj.aspect(exact);
	if (ftype == Bilinear)
		image.filterType(TriangleFilter);
	else if (ftype == Gaussian)
		image.filterType(GaussianFilter);
	image.resize(gobj);
}
void fwrite_u8_pnm(const char* file, int nx, int ny, int nc,
		const u8 im[/*0:nx-1@dx,0:ny-1@dy,0:nc-1@dc*/], int dx, int dy, int dc)

{
	FILE* out = (file) ? fopen(file, "w") : stdout;
	if (!out || (nc != 1 && nc != 3) || fprintf(out, "P%d\n%d %d\n%d\n",
			(nc == 3) ? 6 : 5, nx, ny, 255) <= 0) {
		fail: cout << "Can't write PNM image file :" << file << endl;
		return;
	}
	if (dc == 1 && dx == nc && dy == nc * nx) {
		if (fwrite(im, sizeof(u8), nc * nx * ny, out) != nc * nx * ny)
			goto fail;
	} else {
		int x, y, c;
		for (y = 0; y < ny; y++) {
			for (x = 0; x < nx; x++) {
				for (c = 0; c < nc; c++) {
					if (putc(im[c * dc + x * dx + y * dy], out) == EOF)
						goto fail;
				}
			}
		}
	}
	if (out && file)
		fclose(out);
}

void ParseListFile(const string &nfile, vector<string>&ntimages) {
	string str;
	ifstream ifile(nfile.c_str(), ios::in);
	if (!ifile) {
		cout << " Error: Couldn't Parse the List file " << endl;
		exit(EXIT_FAILURE);
	}
	while (ifile >> str && ifile)
		ntimages.push_back(str);
	ifile.close();
}
void Dec2Bin(UINT val, string& str) { // Considering only 8 bit integers

	//	str.resize(8);
	//	fill(str.begin(),str.end(),'0');
	//	UINT k=7;
	//	while(val)
	//	{
	//		str[k--] = val%2+'0';
	//		val/=2;
	//	}
	int i = val;
	str.resize(sizeof(UINT), '0');
	char *p = &str[0];
	while (i > 0) {
		(i & 0x1) ? (*p++ = '1') : (*p++ = '0');
		i >>= 1;
	}
}
UINT Bin2Dec(const string&str) {
	UINT val = 0, count = 0;

	for (UINT k = str.length(); k >= 0; --k, ++count)
		val += (UINT) pow(2.0, (double) count) * (str[k] - '0');
	return val;
}

string SplitFilename(const string& str) {
	int found;
	//  cout << "Splitting: " << str << endl;
	found = str.find_last_of(".");
	return str.substr(0, found);
}
string SplitFilenameFromDir(const string& str) {
	size_t found;
	//	cout << str << flush;
	found = str.find_last_of("/\\");
	//  cout << " folder: " << str.substr(0,found) << endl;
	return str.substr(found + 1);
}
void FindCrossValidationIndeces(UINT kfold, UINT tnexamples,
		vector<vector<UINT> >& trainidx, vector<vector<UINT> >& testidx) {
	//OLD
	//	UINT trnsize = (tnexamples / kfold) * (kfold - 1), tstsize = tnexamples	- trnsize;
	UINT tstsize = tnexamples / kfold, trnsize = tnexamples - tstsize;
	trainidx.resize(kfold, vector<UINT> (trnsize, 0));
	testidx.resize(kfold, vector<UINT> (tstsize, 0));
	vector<UINT> tmpvec(tnexamples, 0);
	for (UINT i = 0; i < tnexamples; ++i)
		tmpvec[i] = i;

	RandomPermute(tmpvec); // Randomly Permute the indeces...
	UINT tstcount = 0;
	for (UINT i = 0; i < kfold; ++i) {
		// indeces used for testing
		for (UINT j = 0; j < tstsize; ++j) {
			testidx[i][j] = tmpvec[(j + tstcount) % tnexamples];
		}
		tstcount += tstsize;
		for (UINT j = 0; j < trnsize; ++j) {
			trainidx[i][j] = tmpvec[(j + tstcount) % tnexamples];
		}

	}

	for (UINT i = 0; i < kfold; ++i) {
		sort(testidx[i].begin(), testidx[i].end());
		sort(trainidx[i].begin(), trainidx[i].end());
	}

}

void ProcessMemoryUsage(double& vm_usage, double& resident_set) {
	using std::ios_base;
	using std::ifstream;
	using std::string;

	vm_usage = 0.0;
	resident_set = 0.0;

	// 'file' stat seems to give the most reliable results
	//
	ifstream stat_stream("/proc/self/stat", ios_base::in);

	// dummy vars for leading entries in stat that we don't care about
	//
	string pid, comm, state, ppid, pgrp, session, tty_nr;
	string tpgid, flags, minflt, cminflt, majflt, cmajflt;
	string utime, stime, cutime, cstime, priority, nice;
	string O, itrealvalue, starttime;

	// the two fields we want
	//
	unsigned long vsize;
	long rss;

	stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
			>> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
			>> utime >> stime >> cutime >> cstime >> priority >> nice >> O
			>> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

	long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
	vm_usage = vsize / 1024.0;
	resident_set = rss * page_size_kb;
}
void DisplayMemoryUsage() {

	double vm_usage, resident_set;
	ProcessMemoryUsage(vm_usage, resident_set);
	cout << "\n Resident =" << resident_set << " Virtual =" << vm_usage
			<< flush;
}

bool CanReadFile(const string &fname) {
	ifstream ifile(fname.c_str(), ios::in);
	if (!ifile || ifile.eof())
		return false;
	ifile.close();
	return true;
}
void Split(int *start, int *end, int *&ls, int *&le, int *&rs, int *&re,
		int threshold)
/*Will split the vector in to Two portions*/
{

	int *pf = start;/*forward*/
	int *pr = end - 1;/*reverse*/
	int iter = 1;
	while (pf < pr) {
		while (*pf < threshold && pf < pr)
			++pf;
		while (*pr >= threshold && pf < pr)
			--pr;
		// swap the data
		if (pf < pr) {
			int data = *pf;
			*pf = *pr;
			*pr = data;
		}

	}
	ls = start;
	le = pf;
	rs = pr;
	re = end;
}
void MyGaussianFilter(REAL sigma, vector<REAL> &filter) {
	REAL xmin = ceil(-3.0 * sigma - 0.5), xmax = floor(3 * sigma + 0.5);
	UINT len = xmax - xmin + 1;
	filter.resize(len, 0);

	REAL sum = 0, x = xmin, sigmas = sigma * sigma;
	for (vector<REAL>::iterator iter = filter.begin(); iter != filter.end(); ++iter, ++x) {
		*iter = exp(-(x * x) / (2 * sigmas));
		sum += *iter;
	}

	for (vector<REAL>::iterator iter = filter.begin(); iter != filter.end(); ++iter) {
		*iter /= sum;
	}
}
void PadImage(Image& out, UINT dim[]) {
	// dim specifies the dimension in the top,right,bottom and left dimensions to be padded.
	Image input = out;

	//	output.type(TrueColorType);
	//	PixelPacket *ocache;
	UINT icols = input.columns(), irows = input.rows();
	const PixelPacket *icache = input.getConstPixels(0, 0, icols, irows);

	UINT ocols = icols + dim[1] + dim[3], orows = irows + dim[0] + dim[2];
	Image output(Geometry(ocols, orows), "White");
	output.type(TrueColorType);
	output.modifyImage();
	cout << output.rows() << " " << output.columns() << endl << flush;
	//	output.display();
	//	input.display();
	PixelPacket *ocache = output.getPixels(0, 0, ocols, orows);

	//  First Copy the Centrel Pixel Values...
	for (UINT i = 0; i < irows; ++i)
		for (UINT j = 0; j < icols; ++j)
			*(ocache + (i + dim[0]) * ocols + (j + dim[3])) = *(icache + i
					* icols + j);

	// top
	if (dim[0] > 0)
		for (UINT i = 0; i < dim[0]; ++i)
			for (UINT j = 0; j < icols; ++j)
				*(ocache + i * ocols + (j + dim[3])) = *(icache + j);
	//bottom
	if (dim[2] > 0)
		for (UINT i = 0; i < dim[2]; ++i)
			for (UINT j = 0; j < icols; ++j)
				*(ocache + (irows + i + dim[0]) * ocols + (j + dim[3]))
						= *(icache + (irows - 1) * icols + j);
	// right
	if (dim[1] > 0)
		for (UINT i = 0; i < orows; ++i)
			for (UINT j = 0; j < dim[1]; ++j)
				*(ocache + (icols + dim[3]) + i * ocols + j) = *(ocache + icols
						+ dim[3] + i * ocols - 1);
	// left
	if (dim[3] > 0)
		for (UINT i = 0; i < orows; ++i)
			for (UINT j = 0; j < dim[3]; ++j)
				*(ocache + i * ocols + j) = *(ocache + i * ocols + dim[3]);

	output.syncPixels();
	out = output;
}

void PrintIdx(vector<UINT>&trainidx, vector<UINT>&testidx, UINT ptrainoffset,
		UINT ptestoffset) {
	cout << "\n************* Train Idx = ************* " << endl;
	UINT count = 0;
	for (vector<UINT>::iterator iter = trainidx.begin(); iter != trainidx.end(); ++iter, ++count) {
		cout << *iter << " ";
		if (ptrainoffset != 0 && ptrainoffset == count + 1)
			cout << ",";
	}
	cout << endl;
	cout << "§§§§§§§§§§§§§ Test Idx = §§§§§§§§§§§§§ " << endl;
	count = 0;
	for (vector<UINT>::iterator iter = testidx.begin(); iter != testidx.end(); ++iter, ++count) {
		cout << *iter << " ";
		if (ptestoffset != 0 && ptestoffset == count + 1)
			cout << ",";
	}
}
void PrintErrorExit(string msg) {
	cout << endl << msg << " So Quittttttttting " << endl;
	exit(EXIT_SUCCESS);
}

void GetCrossValidationIndeces(UINT kfold, UINT npos, UINT nneg,
		vector<vector<UINT> >&trainidx, vector<vector<UINT> >&testidx) {

	vector<vector<UINT> > ptrainidx, ntrainidx, ptestidx, ntestidx;
	FindCrossValidationIndeces(kfold, npos, ptrainidx, ptestidx);
	FindCrossValidationIndeces(kfold, nneg, ntrainidx, ntestidx);
	WriteMatrix2File("postrainIdx.dat", ptrainidx);
	WriteMatrix2File("negtrainIdx.dat", ntrainidx);
	WriteMatrix2File("postestIdx.dat", ptestidx);
	WriteMatrix2File("negtestIdx.dat", ntestidx);

	UINT ttrain = ptrainidx[0].size() + ntrainidx[0].size(), ttest =
			ptestidx[0].size() + ntestidx[0].size();
	trainidx.resize(kfold, vector<UINT> (ttrain, 0));
	testidx.resize(kfold, vector<UINT> (ttest, 0));
#ifdef PRINT_INDECES
	cout << " Npositves = " << npos << endl << " NNegatives = " << nneg << endl
	<< " KFole = " << kfold << endl;
	cout << " Total Training Size=" << ttrain << endl << " Total Test Size="
	<< ttest << endl;
	// Merge the two training & Validation indeces....
	cout << " \n #################### -Ve Set ##################" << endl;
#endif
	for (UINT i = 0; i < kfold; ++i) {
		copy(ptrainidx[i].begin(), ptrainidx[i].end(), trainidx[i].begin());
		copy(ptestidx[i].begin(), ptestidx[i].end(), testidx[i].begin());
		UINT trnoffset = ptrainidx[i].size(), tstoffset = ptestidx[i].size();
		UINT count = trnoffset;
		for (vector<UINT>::iterator iter = ntrainidx[i].begin(); iter
				!= ntrainidx[i].end(); ++iter, ++count)
			trainidx[i][count] = *iter + npos;

		count = tstoffset;
		for (vector<UINT>::iterator iter = ntestidx[i].begin(); iter
				!= ntestidx[i].end(); ++iter, ++count)
			testidx[i][count] = *iter + npos;
#ifdef PRINT_INDECES
		PrintIdx(ntrainidx[i], ntestidx[i]);
	}
	cout << " \n #################### +Ve Set ##################" << endl;
	for (UINT i = 0; i < kfold; ++i) {
		PrintIdx(ptrainidx[i], ptestidx[i]);
	}
	cout << " \n #################### Complete Set ##################" << endl;
	for (UINT i = 0; i < kfold; ++i) {
		PrintIdx(trainidx[i], testidx[i], ptrainidx[i].size(),
				ptestidx[i].size());
	}
#else
	}
#endif

}
string strtokenize(const string& str, char ch, bool linstance) {

	int pos;
	if (!linstance)
		pos = str.find(ch); // position of "live" in str
	else
		pos = str.find_last_of(ch);
	return str.substr(0, pos); // get from "live" to the end
}
void AppendPath(const string &pathdir, vector<string>&fnames) {
	for (vector<string>::iterator iter = fnames.begin(); iter != fnames.end(); ++iter) {
		string t = pathdir;
		*iter = t + *iter;
	}
}
void WriteStrings2File(const string &fname, vector<string> &feat) {
	ofstream ofile(fname.c_str(), ios::out);
	if (!ofile) {
		cout << " Couldn't Create the File " << fname << endl;
		return;
	}
	for (vector<string>::iterator iter = feat.begin(); iter != feat.end(); ++iter)
		ofile << *iter << endl;
	ofile.close();
}
void validate(boost::any& v, const std::vector<std::string>& values,
		Channel* target_type, int) {
	const std::string& s =
			boost::program_options::validators::get_single_string(values);
	v = boost::any(Channel(boost::lexical_cast<int>(s)));
}

#ifdef WITH_BLITZ
void Write2File(const string & ofname, const Array2DReal & flowu) {
	//	string ofname = "flow.txt";
	std::ofstream ofile(ofname.c_str());
	blitz::TinyVector<int, 2> extent = flowu.extent();
	for (int i = 0; i < extent[0]; ++i) {
		for (int j = 0; j < extent[1]; ++j) {
			ofile << flowu(i, j) << " ";
		}
		ofile << endl;
	}
	ofile.close();
}
Array2DReal ReadData(const string &fname) {
	ifstream ifile(fname.c_str(), ios::in | ios::binary);
	int m, n;
	//	[m,n]=size(data);
	//	fwrite(fid,m,'integer*4'); % Number of Examples
	//	fwrite(fid,n,'integer*4');% Dim of features
	//	fwrite(fid,data','real*4'); % write the data

	ifile.read((char*) &m, sizeof(int));
	ifile.read((char*) &n, sizeof(int));
	REAL *data = new REAL[m * n];
	ifile.read((char*) data, sizeof(REAL) * (m * n));

	blitz::Array<REAL, 2> matrix(data, blitz::shape(m, n),blitz::deleteDataWhenDone);

	return matrix;
}
#endif
#if __defMAIN
int main()
{
	Image image("./tmp/frame0.png");
	//	image.type(Graycolor);
}
#endif
