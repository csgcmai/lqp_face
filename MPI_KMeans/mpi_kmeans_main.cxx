#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <ext/numeric>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "mpi_kmeans.h"

namespace po = boost::program_options;

static uInt count_lines(const std::string& filename) {
	std::ifstream in(filename.c_str());
	if (in.fail()) {
		std::cerr << "count_lines, failed to open \"" << filename << "\"."
				<< std::endl;
		exit(EXIT_FAILURE);
	}
	uInt lines = 0;
	std::string line;
	while (in.eof() == false) {
		std::getline(in, line);
		if (line.size() == 0)
			continue;
		lines += 1;
	}
	in.close();
	return (lines);
}

static void write_cluster_centers(const std::string& output_filename,
		const std::vector<std::vector<PREC> >& data_CX) {

	std::cout << "Writing cluster centers to \"" << output_filename << "\""
			<< std::endl;

	std::ofstream wout(output_filename.c_str());
	if (wout.fail()) {
		std::cerr << "Failed to open \"" << output_filename
				<< "\" for writing." << std::endl;
		exit(EXIT_FAILURE);
	}

	wout << std::setprecision(12);
	for (uInt m = 0; m < data_CX.size(); ++m) {
		for (uInt n = 0; n < data_CX[m].size(); ++n) {
			wout << (n == 0 ? "" : " ") << data_CX[m][n];
		}
		wout << std::endl;
	}
	wout.close();

	return;
}

static void write_cluster_centers(const std::string& output_filename,
		PREC *data_CX, uInt nof_clusters, uInt dims) {

	std::cout << "Writing cluster centers to \"" << output_filename << "\""
			<< std::endl;

	std::ofstream wout(output_filename.c_str());
	if (wout.fail()) {
		std::cerr << "Failed to open \"" << output_filename
				<< "\" for writing." << std::endl;
		exit(EXIT_FAILURE);
	}

	wout << std::setprecision(12);
	uInt cntr = 0;
	for (uInt m = 0; m < nof_clusters; ++m) {
		for (uInt n = 0; n < dims; ++n) {
			wout << (n == 0 ? "" : " ") << data_CX[cntr];
			cntr += 1;
		}
		wout << std::endl;
	}
	wout.close();

	return;
}

static int read_problem_data(const std::string& train_filename, std::vector<
		std::vector<PREC> >& data_X) {
	data_X.clear();

	std::ifstream in(train_filename.c_str());
	if (in.fail()) {
		std::cerr << "Failed to open file \"" << train_filename
				<< "\" for reading." << std::endl;
		std::cerr << "Try mpi_assign --help" << std::endl;
		exit(EXIT_FAILURE);
	}

	std::string line;
	uInt ndims = 0;
	while (in.eof() == false) {
		std::getline(in, line);
		if (line.size() == 0)
			continue; // skip over empty lines

		// remove trailing whitespaces 
		line.erase(line.find_last_not_of(" ") + 1);

		std::vector<PREC> current_data;
		std::istringstream is(line);
		while (is.eof() == false) {
			PREC value;
			is >> value;
			current_data.push_back(value);
		}

		// Ensure the same number of dimensions for each point
		if (ndims == 0)
			ndims = current_data.size();
		assert(ndims == current_data.size());
		data_X.push_back(current_data);
	}
	in.close();

	if (data_X.size() == 0) {
		std::cerr << "No points read from file \"" << train_filename << "\""
				<< std::endl;
		std::cerr << "Try mpi_assign --help" << std::endl;
		exit(EXIT_FAILURE);
	}

	return (data_X.size());

}
static void Write2File(const std::string &fname, uInt *pcount, uInt size) {
	std::ofstream ofile(fname.c_str(), std::ios::out);

	if (!ofile) {
		std::cout << "\n Writing Assignments to File Failed " << fname << "\n";
	} else {
		uInt i = 0;
		for (; i < size && ofile; ++i)
			ofile << pcount[i] << std::endl;
		assert(i == size);
		ofile.close();
	}

}
static void ReadPointCount(const std::string &fname, uInt *pcount, uInt size) {
	std::ifstream ifile(fname.c_str(), std::ios::in);

	if (!ifile) {
		std::cout
				<< "\n No Count File Found So BY default 1 count for each data point \n";
	} else {
		uInt i = 0;
		for (; i < size && ifile; ++i)
			ifile >> pcount[i];
		assert(i == size);
		ifile.close();
	}
}
static int read_problem_data(const std::string& train_filename, PREC *data_X) {

	std::ifstream in(train_filename.c_str());
	if (in.fail()) {
		std::cerr << "Failed to open file \"" << train_filename
				<< "\" for reading." << std::endl;
		std::cerr << "Try mpi_assign --help" << std::endl;
		exit(EXIT_FAILURE);
	}

	uInt nof_points = count_lines(train_filename);
	if (nof_points == 0) {
		std::cerr << "No points read from file \"" << train_filename << "\""
				<< std::endl;
		std::cerr << "Try mpi_assign --help" << std::endl;
		exit(EXIT_FAILURE);
	}

	std::string line;
	uInt ndims = 0;
	uInt cntr = 0;
	while (in.eof() == false) {
		std::getline(in, line);
		if (line.size() == 0)
			continue; // skip over empty lines

		// remove trailing whitespaces 
		line.erase(line.find_last_not_of(" ") + 1);

		std::vector<PREC> current_data;
		std::istringstream is(line);
		while (is.eof() == false) {
			PREC value;
			is >> value;
			current_data.push_back(value);
		}

		// Ensure the same number of dimensions for each point
		if (ndims == 0)
			ndims = current_data.size();
		assert(ndims == current_data.size());
		if (data_X == NULL)
			data_X = (PREC *) malloc(nof_points * ndims * sizeof(PREC));

		for (uInt i = 0; i < ndims; ++i) {
			data_X[cntr] = current_data[i];
			cntr += 1;
		}
	}
	in.close();

	return (nof_points);

}

int main(int argc, char* argv[]) {

	std::string train_filename;
	std::string output_filename;
	int nof_clusters;
	int nof_restarts;
	int maxiter;

	// Set Program options
	po::options_description generic("Generic Options");
	generic.add_options()("help", "Produce help message")("verbose",
			"Verbose output");

	po::options_description input_options("Input/Output Options");
	input_options.add_options()("data",
			po::value<std::string>(&train_filename)->default_value("data.txt"),
			"Training file, one datum per line")("output", po::value<
			std::string>(&output_filename)->default_value("output.txt"),
			"Output file, one cluster center per line");

	po::options_description kmeans_options("K-Means Options");
	kmeans_options.add_options()("k",
			po::value<int>(&nof_clusters)->default_value(100),
			"Number of clusters to generate")("restarts", po::value<int>(
			&nof_restarts)->default_value(0),
			"Number of K-Means restarts. (0: single run)")("maxiter",
			po::value<int>(&maxiter)->default_value(0),
			"Maximum number of K-Means iterations. (0: infinity)");

	po::options_description all_options;
	all_options.add(generic).add(input_options).add(kmeans_options);
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(all_options).run(),
			vm);
	po::notify(vm);

	bool verbose = vm.count("verbose");

	if (vm.count("help")) {
		std::cerr << "K-Means clustering" << std::endl;
		std::cerr << all_options << std::endl;
		std::cerr << std::endl;
		std::cerr << "Example:" << std::endl;
		std::cerr
				<< "  mpi_kmeans --k 2 --data example.txt --output clusters.txt"
				<< std::endl;
		exit(EXIT_SUCCESS);
	}

	// read in the problem
	std::cout << "Training file: " << train_filename << std::endl;
	std::vector<std::vector<PREC> > data_X; // so far kmeans does not support std::<vector>
	int nof_points = read_problem_data(train_filename, data_X);
	assert(nof_points>0);

	uInt dims = data_X[0].size(), *pcount = new uInt[nof_points];
	std::fill(pcount, pcount + nof_points, 1); // by default all points contain a 1 count...
	assert(dims>0);
	ReadPointCount("example_count.txt", pcount, nof_points);
	// convert points to PREC*
	//PREC *X = (PREC *) malloc(nof_points * dims * sizeof(PREC));
	PREC *X = new PREC[nof_points * dims];
	uInt cntr = 0;
	for (uInt m = 0; m < data_X.size(); ++m) {
		for (uInt n = 0; n < data_X[m].size(); ++n) {
			X[cntr] = data_X[m][n];
			cntr += 1;
		}
		data_X[m].clear();
	}
	data_X.clear();

	// start K-Means
	std::cout << "Starting Kmeans ..." << std::endl;
	std::cout << " ... with " << nof_points << " training points " << std::endl;
	std::cout << " ... for " << nof_clusters << " clusters " << std::endl;

	uInt *assignment = new uInt[nof_points];//(uInt *) malloc(nof_points* sizeof(uInt));
	std::fill(assignment, assignment + nof_points, 0);
	//	PREC *CX = (PREC *) calloc(nof_clusters * dims, sizeof(PREC));
	PREC *CX = new PREC[nof_clusters * dims];
	std::fill(CX, CX + nof_clusters * dims, 0);
	PREC sse = kmeans(CX, X, assignment, dims, nof_points, nof_clusters,
			maxiter, nof_restarts, pcount);
	Write2File("bestAssignement.txt", assignment, nof_points);
	//	free(X);
	delete[] X;
	assert(CX);

	std::cout << "Done!" << std::endl;
	std::cout << "Sum of Squared Error : " << sse << std::endl;

	// write the clusters
	// write_cluster_centers(output_filename,data_X);
	write_cluster_centers(output_filename, CX, nof_clusters, dims);
	delete[] CX;
	// done
	exit(EXIT_SUCCESS);
}
