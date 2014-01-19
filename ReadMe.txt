Local Quantized Patterns Based Face Recognition
=================================================

Information
===========

Project webpage: http://sites.google.com/site/sibtulhussain/research/.

This is an implementation of our face verification system [2] based on Local
Quantized Pattern features [1, 3]. The current implementation is the replica of
MATLAB learning code used in [2], however the feature computation code is the original
one [1].

The current distribution contains code for computing local pattern features
(LBP, LTP, LQP, etc.) and unsupervised learning, as well as models trained on
Labelled Faces in Wild datasets. Feature computation code is implemented in C++,
while the learning code is written in python. We also provide scripts (both in
python and bash) for computing features only. The software was tested on
Ubuntu 12.04.03 with the libraries mentioned as follows. 

[ImageMagick](http://www.imagemagick.org/script/install-source.php) version 6.8.1,
g++ version 4.6.3,
[Boost](http://www.boost.org/) version 1.46,
[Eigen](http://eigen.tuxfamily.org) version 2.0 (included in this distribution),
[MPI_KMEANS](http://mloss.org/software/view/48/) version 1.5 (included in this distribution),
[IPython](http://ipython.org) version 0.12.1 (optional for notebooks)

There may be compatibility issues with other versions of libraries.

For questions concerning the code please contact Sibt ul Hussain at
*sibt dot ul dot Hussain at gmail dot com*.

This project has been supported by grants from Higher Education Commission (HEC)
of Pakistan, European Union Project CLASS, and French ANR i.e.
ANR-08-SECU-008-01/SCARFACE.

License
=========
For license terms please see license.lic. 

How to Cite
===========
When citing our system, please cite references [1] and [2].

References:
===========
[1]. S. Hussain and B. Triggs, "Visual Recognition using Local Quantized
     Patterns", In 12th European Conference on Computer Vision (ECCV), 2012.  
[2]. S. Hussain, T. Napoleon and F. Jurie, "Face Recognition using Local
     Quantized Patterns", In 23rd British Machine Vision Conference (BMVC), 2012.  
[3]. S. Hussain, "Machine Learning Methods for Visual Object Detection", PhD Thesis
     University of Grenoble, 2011.
    
Software Requirements:
======================
Although we have successfully tested this system on Ubuntu 12.04 however it
should compile and run on most of the recent Linux distributions. To compile the
code successfully you will need to install ImageMagick development files and C++
boost libraries. Details of the libraries installation process are given below.

1. Install the ImageMagick development libraries. Easiest but not recommended
method is to download the already available compiled sources for your os (e.g.
on ubuntu : sudo apt-get install libmagick++-dev). However, it is not
recommended because by default, these libraries are compiled with openmp flag
and lead to too much threshing in multi-processor setup and make the overall
computation slow. The recommended method is to download the source code
(http://www.imagemagick.org/script/download.php) and compile it locally without
openmp support i.e. ./configure --disable-openmp


2. Install the C++ boost libraries for your OS. e.g. on Ubuntu boost libraries
   can be installed as follows: sudo apt-get install libboost-all-dev


Hardware Requirements:
======================
Running the complete system requires a machine with 8Gb memory (4Gb, if only
features are computed). 


Compiling
============
1. To compile, call the make in features directory, i.e.
 cd features
 make all 
  or
 make -j all (to compile the files in parallel in a multi-processor setup)
 Please see Makefile for other compilation options. 

2. On successful completion executable file named mainFeatures will be generated
and can be located in the build directory. Help on all the available options can
be obtained by calling mainFeatures with help flag i.e. mainFeatures --help

Basic Usage
=============
1. Bash (For computing features only):
--------------------------------------

To compute the features (with default parameters) simply call computeFeatures.sh
providing features type (LQP, LBP or LTP) and files containing list of images
for training and codebook learning (for LQP only) as arguments e.g.

bash computeFeatures.sh lqp ~/dataset/lfw-list.txt ~/experiments/data/
~/dataset/lfw-list-view1.txt
	
Here LQP features will be computed from all the images present in the list file
'lfw-list.txt' using the codebook learned using images present in the list file
'lfw-list-view1.txt'; the computed features will be stored in output directory
'~/experiments/data/'. See the file computeFeatures.sh for further help.

Please note that computed features for each input image are stored linearly in a
separate binary file. You can read this binary file using the provided Matlab
(readFeatureFile in matlab directory) or python code (function read_feature_file
in face-rec/code/trainAndEvaluateLFW.py).

2. Python (Both for features computation and face-verification)
---------------------------------------------------------------
**2.1. Requirements:**
For Python code to work, you will need to install following scientific packages:
numpy, scipy and matplotlib (On Ubuntu: sudo apt-get install python-scipy
python-numpy python-matplotlib; note that Fedora needs both python-matplotlib
and python-matplotlib-tk RPM's).

**2.2. Usage:**
We provide the complete python code for our face-verification algorithm [2], the
code can be found in face-rec directory. Further information can be obtained by
reading the configuration file ('config.py') in the face-rec directory. To use
the code you will need to download any of the LFW aligned version --- although
our code have been extensively tested with following version
(http://www.openu.ac.il/home/hassner/data/lfwa/), however it also works out of
box for the recently released deeply-funnelled version --- and extract it
somewhere on your hard disk. Next, you will need to update the first two
variables (idir and odir) in 'config.py' by pointing them to the respective
path. Once done simply 'cd' to code directory and execute the 'run.py' script by
giving it path of configuration file. i.e.

python run.py --configfile=../config.py

On successful completion computed features can be located in directory
'odir/features/feature-type', whereas learned models on different feature types
with different parameters can be located in 'odir/features/feature-type/data'.
Learned models are numpy binary files and can be loaded into python by calling
load function of numpy. 

For computing features only, you can either configure config.py or call run.py
with command line arguments. Call run.py with '--help' flag to see all the
available command line options.

**2.3. Demo:** 
However before running a thorough set of experiments we recommend you to run the
demo script (face-rec/code/demo.py) providing paths to output directory and
input directory (LFWa folder) and number of experiments --- this argument for
number of experiments is optional, by default demo script runs three experiments
using LTP features. E.g. to run demo.py with a set of 2 experiments we call it
as follows: 

python demo.py /tmp/experiments /data/lfwa 2

This script will hopefully enable you to find out potential problems (if any)
with your installation. Note that running the complete demo program requires a
system with 3Gb memory and takes around 1.5 hours on an i7 machine. Running only
two experiments takes around 8 minutes. See 'demo.py' for further details,
options and possible configurations for other feature types. 

Once done with the experiments, you can run the demo notebook (demo.ipynb) to
visualize the results.  

You can visualize the demo results directly without running the experiments by
either opening the demo-nb.pdf or running the provided demo notebook. We also
provide some excerpts from our experiments notebook (see lqp-face-recog-release-nb.pdf
or lqp-face-recog-release.ipynb). 

Acknowledgements:
=================
Our features computation code draws heavy influences (among others) from the
public release of Felzenszwalb et.al "Discriminatively Trained Deformable Part
Models" and MVG Osolo "LBP" code. We also acknowledge the public release of
Eigen and MPI_KMEANS packages.

We are athankful to Alexis Mignon for releasing his python code for computing PCA and
Thibault Napoleon for his valuable feedback.
