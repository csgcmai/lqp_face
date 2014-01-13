Information
===========

Project webpage: .

This is an implementation of our face recognition system [2]  based on
Local Quantized Pattern features [1]. The current implementation is the 
replica of learning code used in [2], however the feature computation code is original one [1, 2].

The current distribution contains code for computing local pattern features (LBP, LTP, LQP, etc.) 
and unsupervised learning, as well as models trained on Labelled Faces in Wild datasets.

Feature computation code is implemented in C++, while the learning code is written in python.
We also provide scripts (both in both python and bash) for computing features only.  The software was tested on Ubuntu 12.04.03 with the libraries mentioned as follows. 

ImageMagick: Version 6.8.1 Q16
g++: Version 4.6.3
Boost: Version 1.46

There may be compatibility issues with other versions of libraries.

For questions concerning the code please contact Sibt ul Hussain at
<sibt Dot ul Dot Hussain AT gmail DOT com>.

This project has been supported by the Grants from Higher Education Commission (HEC) of Pakistan, European Union Project CLASS, and French ANR i.e. ANR-08-SECU-008-01/SCARFACE.

How to Cite
===========
When citing our system, please cite reference [1, 2].



References:
===========
[1]. S. Hussain and B. Triggs, "Visual Recognition using Local Quantized Patterns",  In 12th European Conference on Computer Vision (ECCV), 2012.  
[2]. S. Hussain, T. Napoleon and F. Jurie, "Face Recognition using Local Quantized Patterns", In 23rd British Machine Vision Conference (BMVC), 2012.


Requirements:
==============
To compile the code successfully you will need to install ImageMagick development files and C++ boost libraries. Details of the complete process are given below.

1. Install the ImageMagick development libraries:

Easiest but not recommended method is to download the already available compiled sources for your os (e.g. on ubuntu : sudo apt-get install libmagick++-dev). However, it is not recommended because by default, these libraries are compiled with openmp flag and lead to too much thrashing in multi-processor setup and make the overall computation slow. The recommended method is to download the source code (http://www.imagemagick.org/script/download.php) and compile it locally without openmp support i.e. ./configure --disable-openmp


2. Install the C++ boost libraries for your OS. e.g. on Ubuntu boost libraries can be installed as follows: sudo apt-get install libboost-all-dev




Compiling
============
1. To compile, call the make in features directory, i.e.
 cd features
 make all 
  or
 make -j all (to compile the files in parallel in multi-processor setup)

2. On successful completion executable file named mainFeatures will be generated and can be located in the build directory. Help on all the  available options can be obtained by calling mainFeatures with help flag i.e. mainFeatures --help

Basic Usage
=============
1. Bash Solution (For Computing Features only):
--------------------------------------------
To compute the features (with default parameters) simply call computeFeatures.sh with features type and image files containing list of images for training and codebook construction  as arguments e.g.

To compute LQP/LBP/LTP features only.

	bash computeFeatures.sh lqp ~/dataset/lfw-list.txt ~/experiments/data/ ~/dataset/lfw-list-view1.txt
	
Where file.txt contain list of images. See the file computeFeatures.sh for further help.

1. Python Solution 
---------------
We provide the complete code for our face-verification algorithm [2], the code can be found in face-rec directory. Further information can be obtained by reading the configuration  file ('config.py') in the face-rec. To use the code you will need to download any of the LFW aligned version --- although our code have been extensively tested with following version (http://www.openu.ac.il/home/hassner/data/lfwa/) it also works out of box for the recently released deeply-funnelled version --- and extract it somewhere on your hard disk. Next, you will need to update the first two variables (idir and odir) in 'config.py' by pointing them to the respective path. Once done simply 'cd' to code directory and execute the 'run.py' script by giving it path of configuration file. i.e.

python run.py --configfile=../config.py

However note that for Python code to work, you will need to install following scientific packages numpy, scipy and matplotlib (On Ubuntu: sudo apt-get install python-scipy python-numpy python-matplotlib; note that Fedora needs both python-matplotlib and python-matplotlib-tk RPM's).





Acknowledgements:
--------------------





