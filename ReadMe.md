Face Recognition Using Local Quantized Patterns
=================================================

Background
==========

Project webpage: http://sites.google.com/site/sibtulhussain/research

This is an implementation of our face verification system [2] based on
Local Quantized Pattern features [1, 3]. The implementation is a
replica of the MATLAB learning code used in [2], with feature
computation code from [1].

The distribution contains routines for computing local pattern
features (LBP, LTP, LQP, etc.) for unsupervised learning, and also
code for face verification model training on the "Labelled Faces in
the Wild" aligned datasets [4, 5]. The feature computation is written in C++
while the learning code is written in python. We also provide examples
of high-level bash and python scripts for running the feature
computation.

The research was supported by the Higher Education Commission (HEC) of
Pakistan, the European Commission research project CLASS, and the
French ANR project ANR-08-SECU-008-01/SCARFACE.

For questions concerning the distribution, please contact Sibt ul
Hussain at sibt.ul.hussain@gmail.com.

References
==========

[1]. S. Hussain and B. Triggs, "Visual Recognition using Local Quantized
     Patterns", In 12th European Conference on Computer Vision (ECCV), 2012.  
[2]. S. Hussain, T. Napoleon and F. Jurie, "Face Recognition using Local
     Quantized Patterns", In 23rd British Machine Vision Conference (BMVC), 2012.  
[3]. S. Hussain, "Machine Learning Methods for Visual Object Detection", PhD Thesis,
     University of Grenoble, 2011. Available as ISBN 978-3841790682 from Editions 
     Universitaires Europeennes, 2012.
[4]. Gary B. Huang and Manu Ramesh and Tamara Berg and Erik Learned-Miller,
     "Labeled Faces in the Wild: A Database for Studying Face Recognition
     in Unconstrained Environments", Technical report, University of
     Massachusetts, 2007. [Available here](http://vis-www.cs.umass.edu/lfw/)
[5]. Lior Wolf, Tal Hassner, and Yaniv Taigman, "Effective Face Recognition by
     Combining Multiple Descriptors and Learned Background Statistics", IEEE
     Trans. on Pattern Analysis and Machine Intelligence (TPAMI), 33(10), Oct. 2011    

License
=======

This is unsupported software released "as is" under the BSD
license.  For details see license.lic.  If you use this software for a
publication, please cite references [1] and [2] above.
    
Environment
===========

The software was tested mainly on 64 bit Ubuntu 12.04.03 but it should
compile and run on other recent Linux workstation distributions
provided that the listed dependencies are installed. For example
Fedora 19 x86_64 also works. MS Windows and Macintosh are not
currently supported.

Running the complete LFW example requires a workstation with 8Gb of
memory. 4Gb suffices if only the features are computed.

The code was tested with the following libraries on Ubuntu 12.04.03:

g++ version 4.6.3,
[ImageMagick](http://www.imagemagick.org/script/install-source.php) version 6.8.1,
[Boost](http://www.boost.org/) version 1.46,
[Eigen](http://eigen.tuxfamily.org) version 2.0 (included in this distribution),
[MPI_KMEANS](http://mloss.org/software/view/48/) version 1.5 (included in this distribution),
[IPython](http://ipython.org) version 0.12.1 (optional, for running ipynb notebooks)

There may be compatibility issues for other versions of these libraries.

Compiling the Basic Feature Code
================================

1. To compile the basic feature computation binary you will need to
install ImageMagick and C++ boost:

- Install the ImageMagick development libraries. The easiest approach
is to install the standard binary packages for your distribution (e.g.
on ubuntu : sudo apt-get install libmagick++-dev). However note that
by default these libraries are compiled with OpenMP enabled which can
cause multi-processor machines to thrash making the overall feature
computation slow. If this causes problems, the recommended fix is to
download the source code
(http://www.imagemagick.org/script/download.php) and compile it
locally without openmp support i.e. ./configure --disable-openmp

- Install the C++ boost libraries, e.g. on ubuntu: sudo apt-get
install libboost-all-dev

2. Run make in the features directory:
 cd features
 make all 
 (or "make -j all" to compile the files in parallel in a
 multi-processor setup). See the Makefile for other compilation
 options.

3. A successful compilation generates a single executable mainFeatures
in the build directory. "mainFeatures --help" lists all of its
available options.

Running the Basic Feature Code
==============================

To compute features with the default parameter settings, run
computeFeatures.sh with the desired feature type (LQP, LBP or LTP) and
image lists:

bash computeFeatures.sh lqp target-images.txt output-dir codebook-images.txt
	
Here, LQP features will be computed for all of the images listed in
'target-images.txt' and stored under directory 'output-dir', using a
codebook learned from the images listed in 'codebook-images.txt'
(optional - needed only for LQP). See computeFeatures.sh for further
help.

Note that the features of each input image are stored in a separate
binary file as a linear array. You can read these files from Matlab
(readFeatureFile in the matlab directory) or python code (function
read_feature_file in face-rec/code/trainAndEvaluateLFW.py).


Installing the Python Code
==========================

This is needed for the demo and face verification experiments and also
for running feature computation under Python.

Besides Python itself, you will need to install the following
scientific packages: numpy, scipy and matplotlib. (On Ubuntu: sudo
apt-get install python-scipy python-numpy python-matplotlib; note that
Fedora needs both python-matplotlib and python-matplotlib-tk
RPM's). To use the optional ipython (*.ipynb) notebooks you will also
need to install ipython.

Using the Python LFW Code
=========================

The complete python code for our face-verification algorithm [2] can
be found in the face-rec directory. For further information and
parameter settings, see 'face-rec/config.py'. 

To run this you will need to download an aligned version of the
"Labeled Faces in the Wild" face dataset. The current parameter
settings were optimized for the old LFWA alignment
(http://www.openu.ac.il/home/hassner/data/lfwa/). The code also runs
well with more recent alignments such as the deep-funnelled one,
however the parameter settings would need to be tweaked to get the
optimal results on these.  

Extract the aligned LFW dataset somewhere on your hard disk, and also
make an empty directory for feature and model output.

1. Demo Run
-----------

A complete run on LFW takes several hours. First we suggest that you
test your installation with a short demo.  'cd' to face-rec/code and
call:

python demo.py output-dir lfwa-path 2

where 'output-dir' is the desired output directory, 'lfwa-path' is the
path to the LFWA input and N is the number of tests to run (default
3, using LTP features).

Running the complete demo program requires 3Gb of memory and takes
around 1.5 hours on an i7 machine. Running only two experiments takes
around 8 minutes. See 'demo.py' for further details, options and
possible configurations for other feature types.

2. Full LFW Run
----------------

Update the variables 'idir' and 'odir' in 'face-rec/config.py' to
point to the desired LFWA input and output directories. Then 'cd' to
the face-rec/code directory and call:

python run.py --configfile=../config.py

On successful completion the computed features for each specified
feature-type will be in a directory 'odir/features/feature-type', with
the learned models and parameter settings in
'odir/features/feature-type/data'.  Learned models are numpy binary
files that can be loaded into python using the numpy load function.

To compute only the features, either update config.py or call run.py
with the desired command line arguments -- call 'run.py --help' to
list the available options.

3. Visualizing the Results
--------------------------

When the experiments finish you can visualize some of the results by
running the demo notebook:

ipython notebook demo.ipynb

You can examine some cached results without running the experiments by
either opening 'demo-nb.pdf' or running the corresponding demo
notebook. We also provide excerpts from our original experimental
notebook in 'lqp-face-recog-release-nb.pdf' and
'lqp-face-recog-release.ipynb'.

Acknowledgements
================

The feature computation code was heavily influenced (among others) by
the public releases of the Felzenszwalb et.al. "Discriminatively
Trained Deformable Part Models" code and the MVG Osolo "LBP" code. We
also acknowledge the public release of Eigen and MPI_KMEANS packages.

We would also like to thank to Alexis Mignon for releasing his python
code for PCA computation and Thibault Napoleon for his valuable
feedback.
