Alignment effects on Performance:

Requirements:
---------------
To compile the code successfully you will need to install ImageMagick development files and C++ boost libraries. The details of the complete process are given below.

1. Install the ImageMagick development libraries:

Easiest but not recommended method is to download the already available compiled sources for your os (e.g. on ubuntu : sudo apt-get install libmagick++-dev). However, it is not recommended because by default, these libraries are compiled with openmp flag and lead to too much thrashing in multi-processor setup and make the overall computation slow. The recommended method is to download the source code (http://www.imagemagick.org/script/download.php) and compile it locally without openmp support i.e. ./configure --disable-openmp


2. Install the C++ boost libraries for your os. e.g. on Ubuntu boost libraries can be installed as follows: sudo apt-get install libboost-all-dev

Tested on Ubuntu 12.04
------------------------
ImageMagick: Version 6.8.1 Q16
g++: Version 4.6.3
Boost: Version 1.46


Compiling
-----------------
1. To compile, call the make in features directory, i.e.
 cd features
 make all 
  or
 make -j all (to compile the files in parallel in multi-processor setup)

2. On successful completion exe file named mainFeatures will be generated and can be located in the build directory. All the different available options documentation can be read by calling mainFeatures with help flag i.e. mainFeatures --help

Usage
-----------------
-----------------

1. Python Solution 
---------------
We provide the complete code (replica) for our face-verification algorithm [2], the code can be found in lfw directory. Further information can be obtained by reading the configuration  file ('config.py') in the lfw. To use the code you will need to download any of the LFW aligned version --- although our code have been extensively tested with following version (http://www.openu.ac.il/home/hassner/data/lfwa/) it also works out of box for the recently released deeply-funnelled version --- and extract it somewhere on your hard disk. Next, you will need to update the first two variables (idir and odir) in 'config.py' by pointing them to the respective path. Once done simply 'cd' to code directory and execute the 'run.py' script by giving it path of configuration file. i.e.

python run.py --configfile=../config.py

However note that for Python code to work, you will need to install following scientific packages numpy, scipy and matplotlib (On Ubuntu: sudo apt-get install python-scipy python-numpy python-matplotlib; note that Fedora needs both python-matplotlib and python-matplotlib-tk RPM's).


2. Bash Solution (For Computing Features only):
--------------------------------------------
To compute the features (with default parameters) simply call computeFeatures.sh with the type of features and files containing list of images for training and codebook construction e.g.

To compute LQP/LBP/LTP features only.

	bash computeFeatures.sh lqp ~/dataset/lfw-list.txt ~/experiments/data/ ~/dataset/lfw-list-view1.txt
	
Where file.txt contain list of images. See the file computeFeatures.sh for further help.


Acknowledgements:
--------------------


References:
--------------
If you use our code, please cite the following papers.

1. S. Hussain and B. Triggs, "Visual Recognition using Local Quantized Patterns",  In 12th European Conference on Computer Vision (ECCV), 2012.  
2. S. Hussain, T. Napoleon and F. Jurie, "Face Recognition using Local Quantized Patterns", In 23rd British Machine Vision Conference (BMVC), 2012.



