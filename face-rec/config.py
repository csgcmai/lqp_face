#configuration file with default options,
# There are four main sections, General, Features, LQP and Learning
# you can comment any of the [Features,LQP] and Learning section
# according to your requirement.
[General] 
#general options
idir=/home/hussain/datasets/LFW/lfwa #path of lfwa folder
odir=/scratch/testing/experiments/ #path where cropped_images, learned model and computed features will be stored
dataset=LFW # name of dataset to use; it can be either LFW or FERET[currently not supported]
width=80 # width of cropped images
height=150 #height of cropped images
padding=10 # (same as cellsize) use a padding of one cell on each side this value must be same as the option cell-size has in the features section
xoffset=1 # offsets to be added to the crop window placed over the original aligned images
yoffset=-4
cbdataset=train-val#complete # use only with LQP Features. This is used to choose subset of dataset for codebook learning e.g. in case of LFW it can be either view1 training validation ('train-val') subset or complete view1 set('complete') 
ftype=LQP # Feature types. Choice can be LBP, LTP, LBP+LTP or LQP
usergb=True # if color images, use  color information during feature computations.

[Features]
#options for feature computation
listfile=""  # a list file containing list of cropped images to  compute features
cellsize=10 # cellsize for the histogram grid
tol=5 #[5,7] # tolerance values used for LTP or LQP features (can be a list) e.g. -t='[5,7]'
  
[LQP] #LQP Options
lqptype=2 #  LQP Type can be  Circular(2.) or Hor+Ver+Diag+ADiag (9.)
lqpsize=7 #radius of LQP Disk or other types (can be a list), e.g. -s='[5,7]'
coding=4 #coding (encoding) type can be :0. Binary, 1. Ternary, 4. Split-Ternary
cbsize=150 #  Codebook size (number of visual words) used for LQP   computation (can be a list), e.g. -w='[100,150]'
cbfile="" #  a list file containing list of images for learning the codebook
                   
[Learning] 
#options for model learning
view=complete #view2 #complete#  Choice of the dataset, options cans be view1: used for
#                        parameter tuning purposes view2: used only for model
#                        evaluation complete: a model parameters will be first
#                        tuned on view1 and results will be reported on view2
ttype=with-pca #   Choice of Training with or without PCA (for feature
#                        evaluation) Available options are with-pca or without-
#                        pca
featdir=/tmp_data/lfw/features/LTP-tolerance=5 #    directory path where computed features have been   stored
dist=cosine #chi-square #             Distance metric for comparing features. Available options are cosine, chi-square and L2
pcadim=[100,  200,  300,  400,  500,  600,  700,  800,  900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000] #          Number of PCA components [you can pass a list i.e.
#                        -p=[100, 200, 300] in which case list will beused
#                       for picking best performing number on view1 ]