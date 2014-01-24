"""
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
"""
""" This file contain utility functions for image cropping, finding mean images, etc."""
from struct import *
import argparse
import numpy as np
import cPickle as cp
import time
from scipy.io import loadmat
from scipy.io import savemat
import os
import sys
import PCA as PCAAlgo
from PIL import Image 
from shutil import *
# from ImageOps import crop
# from numpy.f2py.auxfuncs import isarray
  


def str2bool(v):
  return v.lower() in ("yes", "true", "t", "1")

def rgb2gray(rgb):
    '''converts an rgb image to a grayscale one'''
    return np.dot(rgb[..., :3], [0.299, 0.587, 0.144])

def flattened_images_subdir(idir, odir, usergb):
    ''' copy the images from all the subdirectories into a single flattened directory in odir '''
    sdirs = get_dirs_withpath(idir)
    limages = get_imlist(idir)
    if usergb:
        rgbsuffix = '-rgb'
    else:
        rgbsuffix = ''
    if len(sdirs) != 0 and len(limages) == 0:
        odir = os.path.join(odir, 'lfwa-flattened' + rgbsuffix) ;
        if os.path.exists(odir):
            return odir, rgbsuffix
        print "Writing Flattened images from %s to %s :" % (idir, odir)
        create_directory(odir.strip())
        for d in sdirs:
            limages = os.listdir(d)
            for i in limages:
                if os.path.splitext(i)[1] == '.jpg' or os.path.splitext(i)[1] == '.png':
                    simg = os.path.join(d, i)
                    dimg = os.path.join(odir, i)
                    im = Image.open(simg)
                    #TODO: Replace this portion with png.py module...
                    if im.mode == "RGB" and not usergb:
                        im.convert('L').save(dimg);
                    else:
                        im.save(dimg);
        return odir, rgbsuffix  # return the newly created flattened directory as input directory
    return idir, rgbsuffix
def concatenate_data(data):
    ''' concatenate data sitting at individual indeces'''
    x = data[0]
    for a in data[1:]:
        x = np.concatenate((x, a))
    return x
def generate_data(view):
    ''' generate data for training and testing
        Note that view2, 9 splits are randomly permuted.
    '''
    view = view.lower()
    check_values(view, ["view1", "view2"])
#     datadir = os.path.join, "data")
    datadir = os.path.join(os.path.pardir, "data")
    ifname = os.path.join(datadir , view + '-raw.npy')
    ofname = os.path.join(datadir , view + '-data.npy')
    if view == "view1":
        [inames, ilabels, test, testlabels] = np.load(ifname)
        npos = inames.shape[0] / 2;
        nexamp = npos / 2;
      
        
        pcaTrain = np.concatenate((inames[0:nexamp, :], inames[npos:npos + nexamp, :]));
        pcalabel = np.concatenate((ilabels[0:nexamp], ilabels[npos:npos + nexamp]));
        
        thTrain = np.concatenate((inames[nexamp:npos, :], inames[npos + nexamp:, :]));
        thlabel = np.concatenate((ilabels[nexamp:npos], ilabels[npos + nexamp:]));
        
        
    else:
        [inames, ilabels] = np.load(ifname)
        nfolds = 10
        test = []
        testlabels = []
        pcaTrain = []
        pcalabel = []
        thTrain = []
        thlabel = []
        for k in range(nfolds):
            test.append(inames[0, k])
            testlabels.append(ilabels[0, k])    
            trnidx = np.setdiff1d(range(10), [k])
            np.random.shuffle(trnidx)
            pcaTrain.append(concatenate_data((inames[0, trnidx[0:4]])))
            pcalabel.append(concatenate_data((ilabels[0, trnidx[0:4]])))
            
            thTrain.append(concatenate_data((inames[0, trnidx[4:]])))
            thlabel.append(concatenate_data((ilabels[0, trnidx[4:]])))
    np.save(ofname, [pcaTrain, pcalabel, thTrain, thlabel, test, testlabels])
    return pcaTrain, pcalabel, thTrain, thlabel, test, testlabels
def check_data():
    ''' function checks for the files that store information for view1 and  view2
        datasets
        if these file do not exist it will generate the required ones  and write them files into 
        ../data/ directory '''
    datadir = os.path.join(os.path.pardir, "data")
    v1 = datadir + os.path.sep + 'view1.npy'
    v2 = datadir + os.path.sep + 'view2.npy'
    if not os.path.exists(v1):
        pcaTrain, pcalabel, thTrain, thlabel, test, testlabels = generate_data('view1')
        np.save(v1, [pcaTrain, pcalabel, thTrain, thlabel, test, testlabels])
    if not os.path.exists(v2):
        pcaTrain, pcalabel, thTrain, thlabel, test, testlabels = generate_data('view2')
        np.save(v2, [pcaTrain, pcalabel, thTrain, thlabel, test, testlabels])
        
def check_values(r, coptions):
        ''' Validate given input values for the valid values '''
        if r not in coptions:
            raise ValueError("Wrong option" + r + ", choose among the following options: " + ' '.join(coptions))
def check_islist(val):
    ''' convert a scalar to a list, if the input is already a list returns it as it is '''  
    if type(val) == list:
            return val
    return [val]
def get_list_values(r):
    ''' Return list of values by tokenizing input string list'''
    if isinstance(r, list) or type(r) == int or type(r) == float:
        return r    
    else:
        if r.find('[') != -1:
            if r.find(',') != -1:
                values = map(int, r.strip('[').strip(']').split(','))
            else:
                values = map(int, r.strip('[').strip(']').split())
        else:
            values = int(r)
        return values
def write_list(fname, imlist, odir=''):
    ''' Write list of images into a list file ''' 
    ofile = open(fname, 'w');                
    
    for im in imlist:
#         print im
        if len(odir) != 0:
            im = os.path.join(odir, os.path.basename(im))
        ofile.write(im + '\n');
    ofile.close()
def read_txt_file(fname):
    """ read a txt file, each line is read as a separate element in the returned list"""
    return open(fname).read().splitlines()
def get_number_lines(file):
    """
        return numbers of lines in a file...
       
    """
    return len(read_txt_file(file))
    
def get_dirs(path):   
    ''' return path of (sub-)directories present at the given path '''
    check_path(path)
    dirs = [] 
    for n in os.listdir(path):
        dpath = os.path.join(path, n); 
        if os.path.isdir(dpath):
            dirs.append(dpath)    
    return (dirs if len(dirs) > 0 else None)
def get_imlist(path):
    """Get list of images (supports either .jpg or .png images) available at a given path """
    check_path(path)
    imlist = []
    for filename in os.listdir(path):
        if os.path.splitext(filename)[1] == '.jpg' or os.path.splitext(filename)[1] == '.png':
            imlist.append(os.path.join(path, filename))
    return imlist
# 
# def get_dirs(path='.'):
#     list = []
#     if path:
#         if os.path.isdir(path):
#             try:
#                 names = os.listdir(path)
#             except os.error:
#                 return []
#         names.sort()
#         for name in names:
#             if os.path.isdir(os.path.join(path, name)):
#                 list.append(name)
#         return list
#     
def get_dirs_withpath(path='.'):
    list = []
    names = []
    if os.path.isdir(path):
        try:
            names = os.listdir(path)
        except os.error:
            return names
    names.sort()
    for name in names:
        if os.path.isdir(os.path.join(path, name)) and name != '.svn' and name != '.git':
            list.append(os.path.join(path, name))
    return list

def create_directory(dir):
    ''' Create a new directory recursively
        Note that directory name must note contain any spaces, otherwise it will generate errors'''
    if os.path.isdir(dir):
        return

    parent, base = os.path.split(dir)
    create_directory(parent)
    os.mkdir(dir, 0777)

def check_path(fname, message=""):
    ''' Function check for the validity of path (existence of file or directory path),
     if not found raise exception'''
    if len(message) == 0:
        message = "path " + fname + " Not found"
    if not os.path.exists(fname):
            print message
            raise ValueError(message)
def get_cropped_images_dir_name(options):
    ''' Return the directory name for storing cropped images '''
    dirname = "cropped_%s_images_%dX%d_xoffset=%d_yoffset=%d%s" % (options.dataset, options.width,
                                                                   options.height, options.xoffset,
                                                                   options.yoffset, options.rgbsuffix)
    return os.path.join(options.odir, dirname)
def crop_images_with_offsets(options):
    """ crop images available at the given path 
        option namespace has following members 
        ipath = input which contains the set of images to be cropped
        width = width of cropped image.    
        height = height of cropped image.
        opath = path where cropped images will be saved...
        padding= number of pixels to add on each side of cropped image
        dataset= if provided perform specific preprocessing for a given dataset 
    """     
    assert(options.width > 0 and options.height > 0 and  options.padding >= 0)

    cidir = get_cropped_images_dir_name(options)    
    if not os.path.exists(cidir):
        print "Making directory ", cidir
        create_directory(cidir)
        
    if options.idir == cidir:
        raise ('Error : Input Path: (' + options.idir + ') == Output Path (' + cidir + ')')
            
    imlist = get_imlist(options.idir)
    lfname = options.odir + "%s-list-%dX%d_xoffset=%d_yoffset=%d.txt" % (options.dataset, options.width, options.height, options.xoffset, options.yoffset) 
    print "Writing the list of images in following file ", lfname 
    
    ofile = open(lfname, 'w');
#    return 
    for count, image in enumerate(imlist):
        im = Image.open(image)        
        if im.mode == 'RGB' and not options.usergb:
            im = im.convert('L');  # convert to grayscale
        if options.dataset == 'FERET':
            # Rescale the image... 
            print ' Rescaling the images for FERET Dataset'
            im.resize((im[0] * 0.7, im[1] * 0.7));
        im = np.array(im);
        x = np.floor((im.shape[1] - (options.width + 2 * options.padding)) / 2) - 1;
        y = np.floor ((im.shape[0] - (options.height + 2 * options.padding)) / 2) - 1;
        
        x = x + options.xoffset
        y = y + options.yoffset;
        print 'Processing Image Number %d/%d, Cropping rectangle location x = %d, y=%d' % (count + 1, len(imlist), x, y)
        
        
        if not options.usergb:
            cim = (im [y:y + options.height + 2 * options.padding, x:x + options.width + 2 * options.padding ]);
        else:
            cim = (im [y:y + options.height + 2 * options.padding, x:x + options.width + 2 * options.padding, :]);
        
        
        cr_image = Image.fromarray(cim)
        
        # In case you encounter errors fromarray function please uncomment
        # following two lines of code for saving 
        # np array in .png format using png package (included)
        # and comment out the above line "cr_image=Image.fromarray(cim)"
        
        # Note: This only works for gray-scale images... 
        
#         import png
#         cr_image = png.fromarray(cim, 'L') 
        
        imname = os.path.join(cidir, os.path.basename(image).rsplit('.')[0] + '.png')        
        cr_image.save(imname)
        
        ofile.write(imname + "\n")
    ofile.close()

      
def mean_image(dir):
    """ Computes Mean of the images present at the given path """ 
    images = get_imlist(dir);
    im = np.zeros(np.array(Image.open(images[0])).shape);
    for image in images:
        im = im + np.array(Image.open(image));  
    im = im / len(images);
    return im;

def train_pca(data, ncomp):
    """ 
        Train Whitened PCA on the given data
        data = m x n numpy array 
            where m is number of examples, n number of dimensions of each example
        ncomp = number of components to retain.  
    """
    print " Learning PCA ...."
    st = time.time()
    ncomp = np.min((np.min(data.shape), ncomp))
    use_exterior = False
    if data.shape[0] < data.shape[1]:
        use_exterior = True  # use exterior product for fast computation of PCA, see PCA.py for details
    
    PCA = PCAAlgo.PCA(ncomp, extern=use_exterior);
    PCA.fit(data);
    pdata = PCA.transform(data, whiten=True);
    et = time.time() - st;
    print " Time Taken by Learning PCA = %0.1f sec" % et
    return PCA, pdata
def crop_images(options):
    ''' crop images using given arguments and options'''
    print '''#############################################################'''
    print '''#         Summary of  Options For Image Cropping            #'''
    print '''#############################################################'''
    print 'Dataset = %s' % (options.dataset)
    print 'Input directory =  %r' % options.idir
    print 'Output directory for cropped images  =  %r' % options.odir
    print 'Paddding Size =  %r' % options.padding
    print 'Cropped Image Width(+2*padding) =  %d' % (options.width + 2 * options.padding)
    print 'Cropped Image Height(+2*padding) =  %d' % (options.height + 2 * options.padding)
    print 'Xoffset = %r, Yoffset = %r' % (options.xoffset, options.yoffset)
    print 'Working on ' + 'Gray-Scale Images' if not options.usergb else "RGB Images" 
    print '''#############################################################'''

    check_values(options.dataset, ["FERET", "LFW"])
    if not os.path.exists(options.idir):
        raise ValueError("Image directory %s does not exist" % options.idir)
    
    print options
    print "Cropping Images found at following directory ", options.idir
    crop_images_with_offsets(options)
def crop_images_args(argv=None):
    ''' Parse command line arguments and do image cropping''' 
    if argv == None:
        argv = sys.argv
    print 'ARGV:', argv[1:]
    parser = argparse.ArgumentParser(description=" Image Cropping Parameters...", version=1.0)
    ggroup = parser.add_argument_group('General Options')  # general Group
    ggroup.add_argument('-i', '--idir',
                        dest='idir',
                        help="directory path where images are stored")
    ggroup.add_argument('-o', '--odir',
                        dest='odir', default="/tmp/",
                        help="directory path where cropped images will be stored")                   
    ggroup.add_argument('-w', '--width',
                        dest='width',
                        help="Width of cropped image", default=80)
    ggroup.add_argument('-H', '--height',
                        dest='height',
                        help="height of cropped image", default=150)
    ggroup.add_argument('-p', '--padding',
                        dest='padding',
                        help="Size of Padding to be added around the cropped image (Must be of the cell size)", default=10)
    ggroup.add_argument('-x', '--x-offset',
                        dest='xoffset',
                        help=('''(+,-) offset position of cropping rectangle from the center image '''),
                        default=1, type=int)
    ggroup.add_argument('-y', '--y-offset',
                        dest='yoffset',
                        help=('''(+,-) offset position of cropping rectangle from the center image '''),
                        default=-4, type=int)
    ggroup.add_argument('-D', '--dataset',
                        dest='dataset',
                        help=('''Images belong to which dataset, LFW or FERET'''),
                        default="LFW")
    ggroup.add_argument('--usergb', dest="usergb",
                            help="Use color information during feature computations ",default="False")
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)         
    options = parser.parse_args()
    # ftype, featdir, listfile
    options.usergb = str2bool(options.usergb)
    crop_images(options)
if __name__ == "__main__":
    sys.exit(crop_images_args())
    
#   crop_images('/home/hussain/dataset/', 80 , 150, '/tmp/', 10, 'LFW');
# usage='python util.py -i /home/hussain/datasets/LFW/lfwa/ -o /tmp/'    
    
#    if options.ttype == "with-pca":
#        pcadim = options.pcadim
#        pcadim = map(int, pcadim.strip('[').strip(']').split(','))
#        if not isinstance(pcadim, list):
#            pcadim = [pcadim]
#        print "Training and Evaluating a PCA Model with following given dimensions  = ", pcadim
#        LFW_PCA(options.view, options.featdir, options.ftype, pcadim)
#    else:
#        print "Training and Evaluating a Model without PCA  "
#        LFW_WPCA(options.view, options.featdir, options.ftype)    
