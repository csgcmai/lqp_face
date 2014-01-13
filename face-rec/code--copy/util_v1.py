"""
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
"""
""" This file contain functions to read LFW dataset format and generates training and test data for folds"""
from struct import *
import numpy as np
import cPickle as cp
import time
from scipy.io import loadmat
from scipy.io import savemat
import os
import sys
# import PCA as pca
from PIL import Image 
# from twisted.words.test.test_msn import printError

def get_imlist(path):
    """Get list of images (supports either .jpg or .png images) available at a given path """
    imlist = []
    for filename in os.listdir(path):
        if os.path.splitext(filename)[1] == '.jpg' or os.path.splitext(filename)[1] == '.png':
            imlist.append(os.path.join(path, filename))
    return imlist

def crop_images(ipath, width, height, opath, padding=0 , dataset=''):
    """ crop images available at the given path 
        ipath = input which contains the set of images to be cropped
        width = width of cropped image.
        height = height of cropped image.
        opath = path where cropped images will be saved...
        padding= number of pixels to add on each side of cropped image
        dataset= if provided perform specific preprocessing for a given dataset 
    """     
    assert(width > 0 and height > 0 and  padding >= 0)
    
    if not os.path.exists(opath):
        os.mkdir(opath);
        
    if ipath == opath:
        print 'Error : Input Path: (', ipath, ') == Output Path (', opath, ')'
        sys.exit(1)
            
    imlist = get_imlist(ipath)
    
    for image in imlist:
        im = Image.open(image)        
        if im.mode == 'RGB':
            im.convert('L');  # convert to grayscale
        if dataset == 'FERET':
            # Rescale the image... 
            print ' Rescaling the images for FERET Dataset'
            im.resize((im[0] * 0.7, im[1] * 0.7));
        im = np.array(im);
        x = np.floor((im.shape[1] - (width + 2 * padding)) / 2) - 1;
        y = np.floor ((im.shape[0] - (height + 2 * padding)) / 2) - 1;
        cr_image = Image.fromarray(im [y:y + height + 2 * padding, x:x + width + 2 * padding ]);        
        cr_image.save(os.path.join(opath, os.path.basename(image)))
        
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
    
    pca = pca.PCA(ncomp, extern=use_exterior);
    pca.Fit(data);
    pdata = pca.transform(X, whiten=True);
    et = time.time() - st;
    print " Time Taken by Learning PCA = %0.1f sec" % et
    return pca, pdata

def compute_distance(lfeat, rfeat, dtype, flfeat, frfeat):
    """ 
        compute distance between given matrices for pair of images i.e. lfeat and rfeat 
        lfeat = m x n matrix for left item of the pair, where m is number of examples and n is number of features...
        if flipped features are given then minimum between among the four possible options
        i.e min( min(lfeat,rfeat) , min(lfeat,frfeat), min(flfeat,rfeat), min(flfeat,frfeat) )
        is returned
        dtype = type of distance metric 
        """
    if dtype.lower() == "cosine":
        norm1 = np.sqrt (np.sum (np.power(lfeat, 2), 1));
        norm2 = np.sqrt (np.sum (np.power(rfeat, 2), 2));
        dis1 = (np.sum (lfeat * rfeat, 1) / (norm1 * norm2));
    
        norm3 = np.sqrt (np.sum (np.power(frfeat, 2), 2));            
        dis2 = np.sum (lfeat * frfeat, 1) / (norm1 * norm3);
        
        norm4 = np.sqrt (np.sum (np.power(flfeat, 2), 2));
        dis3 = np.sum (rfeat * flfeat, 1) / (norm2 * norm4);
        
        dis4 = np.sum (flfeat * frfeat, 1) / (norm3 * norm4);
        
        dist = -np.max(np.max(np.max(dis1, dis2), dis3), dis4);
        
    elif dtype.lower() == 'l2':
        dist = np.sum((np.power(lfeat - rfeat), 2), 2);
    else:
        print " Error: Invalid option for computing distance"
        sys.exit()
        
def find_optimal_threshold(distance, labels, method):
    """
        Find the optimal threshold for the decision process
        distance = computed distances betweeen image pairs
        labels = class labels
        method = method used for finding the threshold 
    """
    npos = np.sum(labels == 1);  # number of positive examples;
    sdist = np.sort(dist);
    if method.upper() == "EER":
         threshold = sdist(npos);  # number of -ve's =  number of +ve's
    else:
        print " Error: Invalid method type for the configuration of Threhsold"
        sys.exit()
    return threshold

def read_LFW(fname): 
    """ read LFW dataset annotation information with names and labels... """
#    dir = os.getcwd()
# os.chdir(dirname)
    if not os.path.exists(fname):
        print 'LFW File :', fname, 'does not exist'
        exit()
    lines = open(fname).readlines();
    print lines
    for l in lines:
        print l

def read_feature_file(fname):
    """
        Reads the feature from the given fname:
    """
    fin = open(fname, 'rb');
    buffer = fin.read();
    tvar = list(unpack('3i', buffer[0:12]));
    nrows = tvar[1];   
    colsize = tvar[2];
    features = np.array(unpack('*f', buffer[12:]));
    assert(nrows * colsize == features.len)
    fin.close()
    return features, nrows, colsize

def read_descriptors(inames, featdir, ftype, dstr=" ------ Descriptors ---------"):
    """
     Read Descriptors for the given list of images and return 
    function [litem,ritem,flitem,fritem,nd]=readDescriptors(inames,featdir,ftype,dstr,pprocess,nd)
    read the descriptors written by mfiles
    iname: name of images
    featdir:feature directory containing computed features.
    ftype: [LBP, LTP, LQP]
    dstr: display string 
    Output=
    litem: features for image pair left item;
    ritem: features for image pair right item;
    litem: features for flipped-image pair left item;
    ritem: features for flipped-image pair right item;
    """
#    cellsize = 10;
#    cbsize = 300;  # 2 x 150    
    ftype = ftype.upper();
    nimages = len(inames);
    feat, m, n = read_feature_file(os.path.join(featdir, inames[0][0], '.jpg-', ftype));
#    n = n / cbsize; ncells = m * n;
        
    litem = zeros(nimages, feat.size, float32);
    flitem = zeros(nimages, feat.size, float32);
    ritem = zeros(nimages, feat.size, float32);
    fritem = zeros(nimages, feat.size, float32);
    st = time.time();
    
    print 'Starting Reading Descriptors for %s' % dstr;

    for k, images in enumerate(imnames):
        litem[k, :] = read_feature_file(featdir, os.path.join(imnames[k][0], '.jpg-', ftype))
        ritem[k, :] = read_feature_file(featdir, os.path.join(imnames[k][1], '.jpg-', ftype));
        
        flitem[k, :] = read_feature_file(featdir, os.path.join(imnames[k][0], '.jpg-', ftype, '-f'));
        fritem[k, :] = read_feature_file(featdir, os.path.join(imnames[k][1], '.jpg-', ftype, '-f'));    
    
    if k % 200 == 0 :
        print ' Reading %d/%d' % (k, len(imnames)) 
    et = time.time() - st;
    
    print 'Finished Reading Descriptors for %s \n Total Time Taken by Reading Operation =%f\n' % (dstr, et);
    return litem, ritem, flitem, fritem
    
def read_compute_distance(pca, images, featdir, ftype, dtype):
    """
        Read features and then Project them to learn PCA components
        and finally compute distance between the projected descriptors...        
        pca = pca class used to project features...
        images = name of images to read
        featdir= feature directory
        ftype = [LBP, LTP, LQP]
        dtype = distance type [cosine, L2]
        
        output:
        dist= computed distance between pair of projected descriptors...
    """
    litem, ritem, flitem, fritem = read_descriptors(images, featdir, ftype)
        # project the features
    litem = pca.transform(litem, whiten=True);
    ritem = pca.transform(ritem, whiten=True);
    flitem = pca.transform(flitem, whiten=True);
    fritem = pca.transform(fritem, whiten=True);
    # find distance between pair of iamges...
    dist = compute_distance(lfeat, rfeat, dtype, flfeat, frfeat)
    
    return dist
def compute_accuracy(dist, labels, threshold):
    """
         Compute the average accuracy over the given set of images.
         dist: computed distance between pair of images.
         labels: true class labels
         threshold: decision threshold.
    """
    trueclass = np.sum(dist <= threshold and labels == 1) + np.sum(dist > threshold and laels != 1)
    accuracy = trueclass / len(labels)
    return accuracy
def train_test(data, featdir, ftype, pcadim, filesuffix):
    """
    Train the PCA and Evaluate on the test set..
    data = list of images name for training and testing.
    ftype = feature type
    pcadim= number of pca-dimensions...
    filesuffix= used for storing temporary data...
    """
    pcaTrain, pcalabel, thTrain, thlabel, test, testlabels = data
     
    if isinstance(pcadim, list):
        ncomp = max(pcadim)
    else:
        ncomp = pcadim
     # read features for the pca learning...
     
    pcafile = "pca-" + str(ncomp) + "-" + filesuffix + ".npy"
    if not os.path.exists(path): 
        litem, ritem, flitem, fritem = read_descriptors(pcaTrain, featdir, ftype)
        pcafeat = np.concatenate((litem, flitem))
        pca, pfeat = train_pca(pcafeat, ncomp)
        np.save(pcafile, [pca, pfeat]);
    else:
        pca, pfeat = np.load(pcafile)
    
    count = 1;
    accuracy = []
    for k in pcadim:
        suffix = filesuffix + "-ncomp=" + str(k)
        # read and project features and compute the distance between pair of images...
        dist = read_compute_distance(pca, thTrain, featdir, ftype, "cosine")
        thershold = find_optimal_threshold(dist, thlabel, "EER")
        dist = read_compute_distance(pca, test, featdir, ftype, "cosine")
        accuracy.append(compute_accuracy(dist, testlabels, threshold))
        print "For PCA Dim = %d, Classification accuracy(on view = %s)= %f" % (k, view, accuracy[-1])
    return accuracy         
     

def LFW(view, featdir, ftype='LQP', chunk=1800):
    """
        Read the features for the given view and train & evaluate the PCA model for the given set
        of examples:
        view= data from view1 or view2 to use (see LFW documentation)
        featdir= directory path of features.  
        ftype =  any of [LBP, LTP, LQP] features to use
        chunk =  used only during the cross-validation to find the best number of PCA components. 
     """
    view = view.lower()
    ftype = ftype.upper()
    ncomp = 2000;  # number of pca components
    
    if view not in ["view1", "view2"]:
        print "Error: Provide the Dataset view (either view1 or view2)to use for the training "
        sys.exit(1)
    if ftype not in ["LBP", "LTP", "LQP"]:
        print "Error: Wrong Feature Type, choose any one of the following features [LBP, LTP, LQP] "
        sys.exit(1)
    
    
    pcadim = range(chunk, ncomp, chunk)
    pcafile = "./data/" + view + "-pca-" + ftype + ".npy"    
   
    if os.path.exists(pcafile):
        pca, trnfeat, trnlabels, testfeat, testlabels, accuracy = load(pcafile)
    else:
        data = np.load(view + '.npy')
        if view == "view1":
            pca, trnfeat, trnlabels, testfeat, testlabels, accuracy = train_test(data, featdir, ftype, pcadim)
            np.save(pcafile, [pca, trnfeat, trnlabels, testfeat, testlabels, accuracy])
        else:
            accuracy = []
            pca = []
            trnfeat = []
            trnlabels = []
            testfeat = []
            testlabels = []
            for k in range(10):  # 10 fold
                # extract training data
                tdata = []    
                for t in range(6):
                    tdata.append(data[t][0][k])
 
                #                  
                ipca, itrnfeat, itrnlabels, itestfeat, itestlabels, iaccuracy = train_test(data, featdir, ftype, pcadim)
                accuracy.append(iaccuracy)
                pca.append(ipca)
                
                trnfeat.append(itrnfeat)
                trnlabels.append(itrnlabels)
                
                testfeat.append(itestfeat)
                testlabels.append(itestlabels)
            np.save(pcafile, [pca, trnfeat, trnlabels, testfeat, testlabels, accuracy])

if __name__ == "__main__":
    crop_images('/home/hussain/dataset/', 80 , 150, '/tmp/', 10, 'LFW');
