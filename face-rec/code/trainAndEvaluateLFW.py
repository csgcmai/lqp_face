"""
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
"""
''' This file contains functions for the training and evaluation of models with different options on LFW dataset'''
from util import *
import math as math
import os
from matplotlib.backends.backend_tkagg import round
def compute_distance(lfeat, rfeat, dtype, flfeat, frfeat):
    """ 
        compute distance between given matrices for pair of images i.e. lfeat and rfeat 
        lfeat = m x n matrix for left item of the pair, where m is number of examples and n is number of features...
        if flipped features are given then minimum between among the four possible options
        i.e min( min(lfeat,rfeat) , min(lfeat,frfeat), min(flfeat,rfeat), min(flfeat,frfeat) )
        is returned
        dtype = type of distance metric 
        
        """
        
    print 'Computing pair-wise similarity using %s Similarity Metric' % dtype    
    if dtype.lower() == "cosine":
        norm1 = np.sqrt (np.sum (np.power(lfeat, 2), 1));
        norm2 = np.sqrt (np.sum (np.power(rfeat, 2), 1));
        dist1 = (np.sum (lfeat * rfeat, 1) / (norm1 * norm2));
    
        norm3 = np.sqrt (np.sum (np.power(frfeat, 2), 1));            
        dist2 = np.sum (lfeat * frfeat, 1) / (norm1 * norm3);
        
        norm4 = np.sqrt (np.sum (np.power(flfeat, 2), 1));
        dist3 = np.sum (rfeat * flfeat, 1) / (norm2 * norm4);
        
        dist4 = np.sum (flfeat * frfeat, 1) / (norm3 * norm4);
        
        dist = -np.maximum(np.maximum(np.maximum(dist1, dist2), dist3), dist4);
        
    elif dtype.lower() == 'l2':
        dist1 = np.sum(np.power(lfeat - rfeat, 2), 1);
        dist2 = np.sum(np.power(lfeat - frfeat, 2), 1);
        dist3 = np.sum(np.power(flfeat - rfeat, 2), 1);
        dist4 = np.sum(np.power(flfeat - frfeat, 2), 1);
        dist = np.minimum(np.minimum(np.minimum(dist1, dist2), dist3), dist4);
    elif dtype.lower() == 'chi-square':
        lfeat = lfeat * lfeat;  # because we are using feature normalized via L1Sqrt method
        flfeat = flfeat * flfeat;  # because we are using feature normalized via L1Sqrt method
        rfeat = rfeat * rfeat;  # because we are using feature normalized via L1Sqrt method
        frfeat = frfeat * frfeat;  # because we are using feature normalized via L1Sqrt method
        
        dist1 = 0.5 * np.sum(np.power(lfeat - rfeat, 2) / (lfeat + rfeat + np.spacing(1)), 1)
        dist2 = 0.5 * np.sum(np.power(lfeat - frfeat, 2) / (lfeat + frfeat + np.spacing(1)), 1)
        dist3 = 0.5 * np.sum(np.power(flfeat - rfeat, 2) / (flfeat + rfeat + np.spacing(1)), 1)
        dist4 = 0.5 * np.sum(np.power(flfeat - frfeat, 2) / (flfeat + frfeat + np.spacing(1)), 1)
        dist = np.minimum(np.minimum(np.minimum(dist1, dist2), dist3), dist4);
    else:
        raise ValueError("Invalid option for computing distance" + dtype)
    return dist
        
def find_optimal_threshold(distance, labels, method):
    """
        Find the optimal threshold for the decision process
        distance = computed distances betweeen image pairs
        labels = class labels
        method = method used for finding the threshold 
    """
    # redundant because number of +ve's == -ve's
    npos = np.sum(labels == 1);  # number of positive examples;
    sdist = np.sort(distance);
    if method.upper() == "EER":
         threshold = sdist[npos];  # number of -ve's =  number of +ve's
    else:
        raise ValueError(" Invalid method %s type for the configuration of Threshold" % method)
    return threshold

def read_LFW(fname): 
    """ read LFW dataset annotation information with names and labels... """
#    dir = os.getcwd()
# os.chdir(dirname)
    if not os.path.exists(fname):
        raise ValueError('LFW File :' + fname + 'does not exist')
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
    features = np.array(unpack(str(nrows * colsize) + 'f', buffer[12:]));
    assert(nrows * colsize == features.size)
    fin.close()
    return features, nrows, colsize

def read_descriptors(inames, featdir, ftype, dstr=""):
    """
     Read Descriptors for the given list of images and return 
    function [litem,ritem,flitem,fritem,nd]=readDescriptors(inames,featdir,ftype,dstr,pprocess,nd)
    read the descriptors written by cpp code
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
    feat, m, n = read_feature_file(os.path.join(featdir, inames[0][0][0] + '.png-' + ftype));
#    n = n / cbsize; ncells = m * n;
        
    litem = np.zeros((nimages, feat.size), np.float32);
    flitem = np.zeros((nimages, feat.size), np.float32);
    ritem = np.zeros((nimages, feat.size), np.float32);
    fritem = np.zeros((nimages, feat.size), np.float32);
    st = time.time();
    
    print 'Starting Reading Descriptors %s' % dstr;

    for k, images in enumerate(inames):
        lename = os.path.join(featdir, inames[k][0][0] + '.png-' + ftype);
        rename = os.path.join(featdir, inames[k][1][0] + '.png-' + ftype);
        
        litem[k, :] = read_feature_file(lename)[0];
        ritem[k, :] = read_feature_file(rename)[0];
        
        flitem[k, :] = read_feature_file(lename + '-f')[0];
        fritem[k, :] = read_feature_file(rename + '-f')[0];    
    
        if k % 200 == 0 :
            print ' Reading %d/%d' % (k, len(inames)) 
    et = time.time() - st;
    
    print 'Finished Reading Descriptors %s \n Total Time Taken by Reading Operation =%f\n' % (dstr, et);
    return litem, ritem, flitem, fritem
    
def read_compute_distance_pca(PCA, images, featdir, ftype, dtype, pcadim):
    """
        Read features and then Project them to learn PCA components
        and finally compute distance between the projected descriptors...        
        PCA = PCA class used to project features...
        images = name of images to read
        featdir= feature directory
        ftype = [LBP, LTP, LQP]
        dtype = distance type [cosine, L2]
        pcadim= list of PCA dimensions used during validation purposes....
        output:
        dist= computed distance between pair of projected descriptors...
    """
    litem, ritem, flitem, fritem = read_descriptors(images, featdir, ftype)
#    if not isinstance(pcadim, list):
#            pcadim = [pcadim]
    
    # project the features
    dist = np.zeros((litem.shape[0], len(pcadim)))
    for count, k in enumerate(pcadim):
        plitem = PCA.transform(litem, whiten=True, ncomp=k);
        pritem = PCA.transform(ritem, whiten=True, ncomp=k);
        pflitem = PCA.transform(flitem, whiten=True, ncomp=k);
        pfritem = PCA.transform(fritem, whiten=True, ncomp=k);
        dist[:, count] = compute_distance(plitem, pritem, dtype, pflitem, pfritem)
    return dist
def read_compute_distance(images, featdir, ftype, dtype):
    """
        Read features and then Project them to learn PCA components
        and finally compute distance between the projected descriptors...        
        PCA = PCA class used to project features...
        images = name of images to read
        featdir= feature directory
        ftype = [LBP, LTP, LQP]
        dtype = distance type [cosine, L2]
        
        output:
        dist= computed distance between pair of projected descriptors...
    """
    chunk = min(1000, len(images));
    nchunks = math.ceil(len(images) / float(chunk))
    print " Reading features for %d images in %d chunks of  %d images" % (len(images), nchunks, chunk)
    dist = np.zeros(len(images)) 
    for i, s in enumerate(range(0, len(images), chunk)):
        print '------------->  Reading Chunk # %d/%d' % (i + 1, nchunks)
        limit = min(len(images), s + chunk)
        litem, ritem, flitem, fritem = read_descriptors(images[s:limit], featdir, ftype)
        # find distance between pair of images...
        s1 = compute_distance(litem, ritem, dtype, flitem, fritem)
        dist[s:limit] = s1     
    return dist
def compute_accuracy(dist, labels, threshold):
    """
         Compute the average accuracy over the given set of images.
         dist: computed distance between pair of images.
         labels: true class labels
         threshold: decision threshold.
    """
    trueclass = np.sum(np.logical_and((dist <= threshold).reshape(dist.shape[0], 1), labels == 1)
           ) + np.sum(np.logical_and((dist > threshold).reshape(dist.shape[0], 1), labels != 1)) 
    accuracy = float(trueclass) / len(labels)
    return accuracy
def train_test_pca(data, featdir, ftype, pcadim, filesuffix, distance):
    """
    Train the PCA and Evaluate on the test set..
    data = list of images name for training and testing.
    ftype = feature type
    pcadim= number of PCA-dimensions...
    filesuffix= used for storing temporary data...
    """
    pcaTrain, pcalabel, thTrain, thlabel, test, testlabels = data
     
    if isinstance(pcadim, list):  # and len(pcadim) > 1:
        ncomp = max(pcadim)
    else:
        ncomp = pcadim
        pcadim = [pcadim]
     # read features for the PCA learning...
    ddir = os.path.join(featdir, "data");
    if not os.path.exists(ddir):
        print " Creating following Data directory : ", ddir
        os.mkdir(ddir)
    pcafile = os.path.join(ddir, "PCA-") + str(ncomp) + "-" + distance + "-" + filesuffix + ".npy"
    if not os.path.exists(pcafile): 
        litem, ritem, flitem, fritem = read_descriptors(pcaTrain, featdir, ftype)
        pcafeat = np.concatenate((litem, ritem))
        PCA, pfeat = train_pca(pcafeat, ncomp)
        np.save(pcafile, [PCA, pfeat]);
    else:
        PCA, pfeat = np.load(pcafile)
    
    # read and project features and compute the distance between pair of images for both training and test set...
    trndist = read_compute_distance_pca(PCA, thTrain, featdir, ftype, distance, pcadim)
    testdist = read_compute_distance_pca(PCA, test, featdir, ftype, distance, pcadim)
    
    accuracy = []
    for count, k in enumerate(pcadim):
#        suffix = filesuffix + "-ncomp=" + str(k)   
        threshold = find_optimal_threshold(trndist[:, count], thlabel, "EER")
        accuracy.append([k, compute_accuracy(testdist[:, count], testlabels, threshold)])
        print "For PCA Dim = %d, Classification accuracy(on view = %s)= %f" % (k, filesuffix, accuracy[-1][1])
    return  accuracy, PCA         
     
def train_test(data, featdir, ftype, filesuffix, setsuffix="", distance='cosine'):
    """
    Simply set the threshold and Evaluate it on the test set..
    data = list of images name for training and testing.
    ftype = feature type
    filesuffix= used for storing temporary data...
    """
    pcaTrain, pcalabel, thTrain, thlabel, test, testlabels = data        
   
    # concatenate the features...
    thTrain = np.concatenate((pcaTrain, thTrain))
    dist = read_compute_distance(thTrain, featdir, ftype, distance)
    thlabel = np.concatenate((pcalabel, thlabel))
    threshold = find_optimal_threshold(dist, thlabel, "EER")
    dist = read_compute_distance(test, featdir, ftype, distance)
    accuracy = (compute_accuracy(dist, testlabels, threshold))
    print "For Feature Type = %s, Classification accuracy(on view = %s%s)= %f" % (ftype, filesuffix, setsuffix, accuracy)
    return accuracy         
     
def LFW_PCA(view, featdir, ftype, pcadim, distance):
    """
        Read the features for the given view and train & evaluate the PCA model for the given set
        of examples:
        view= data from view1 or view2 to use (see LFW documentation)
        featdir= directory path of features.  
        ftype =  any of [LBP, LTP, LQP] features to use
        pcadim  = Number of Principal components to learn during training, can be a list which is
                 used only during the cross-validation to find the best number of PCA components. 
     """
    view = view.lower()
    ftype = ftype.upper()    
    
    if view not in ["view1", "view2"]:
        raise ValueError(' Provide the Dataset view (either view1 or view2)to use for the training ')
    if ftype not in ["LBP", "LTP", "LQP"]:
        raise ValueError('Error: Wrong Feature Type, choose any one of the following features [LBP, LTP, LQP] ')
    if isinstance(pcadim, list) :
        mpcadim = max(pcadim)
    else:
        mpcadim = pcadim;
    
   
    ddir = os.path.join(featdir, "data");
    if not os.path.exists(ddir):
        print " Creating following Data directory : ", ddir
        os.mkdir(ddir)
    
    pcafile = os.path.join(ddir, view) + "-" + ftype + "-PCA-" + str(mpcadim) + "-" + distance + "-results.npy"    
   
    
    if os.path.exists(pcafile):
        accuracy = np.load(pcafile)
        if view == "view1":
            print (" \n ----------------------- \n"
            "Classification Accuracy of %s's on %s \n"
            "-----------------------\n" % (ftype, view))
            print " PCA Dimension \t Accuracy "
            for val in accuracy:
                print "%f\t%f" % tuple(val)
            idx = np.argmax(accuracy[:, 1])
            return accuracy[idx, 0], accuracy[idx, 1]
        else:
            print (" \n ----------------------- \n" 
            "Classification Accuracy of %s's on %s \n-----------------------\n" % (ftype, view))
            print "#        PCA-Dim        Accuracy"
            mean = 0;
            for count, val in enumerate(accuracy):
                mean = mean + val[1]
                print str(count + 1) + ".\t %f\t%f" % tuple(val)
#             mean = mean / (count + 1)
            print " Mean =", (mean / (count + 1) * 100)
            return val[0], mean
#           print "For Feature Type = %s, Classification accuracy(on view = %s)= %f" % (ftype, view, sum(accuracy) / len(accuracy))
    else:
        accuracy = []
        data = np.load(os.path.join(os.pardir, os.path.join("data" , view + '.npy')))
        if view == "view1":
            fsuffix = 'view1-' + ftype 
            accuracy, ignore = train_test_pca(data, featdir, ftype, pcadim, fsuffix, distance)
            np.save(pcafile, accuracy)
            maxa = -np.inf
            bdim = 0
            for a in accuracy:
                if a[1] > maxa:
                    maxa = a[1]
                    bdim = a[0]                  
            return  bdim, maxa
        else:  # check that PCA has been already learned on the view1 data...
           

            for k in range(10):  # 10 fold
                # extract training data
                tdata = []    
                for t in range(6):
                    tdata.append(data[t][k])
 
                #
                fsuffix = 'view2-Set' + str(k + 1) + ftype;                  
                iaccuracy, ignore = train_test_pca(tdata, featdir, ftype, mpcadim, fsuffix, distance)
                accuracy.append(iaccuracy[0])
#                PCA.append(ipca)
#           print " Mean Accuracy on View2 = ", sum(accuracy) / len(accuracy)
            print (''' \n ----------------------- \n''' 
            '''Classification Accuracy of %s's on %s \n
             -----------------------\n''' % (ftype, view))
            print " Subset No \t PCA Dim \t Accuracy "
            mean = 0;
            for count, val in enumerate(accuracy):
                mean = mean + val[1]
                print str(count + 1) + ".\t %f\t%f" % tuple(val)
#             mean = mean / (count + 1)
            print " Mean =", (mean / (count + 1) * 100) 
            
        np.save(pcafile, accuracy)
        return mpcadim, mean

def LFW_WPCA(view, featdir, ftype='LQP', distance='cosine'):
    """
        Read the features for the given view and train & evaluate the simple (without PCA) model for the given set
        of examples:
        view= data from view1 or view2 to use (see LFW documentation)
        featdir= directory path of features.  
        ftype =  any of [LBP, LTP, LQP] features to use         
     """
    view = view.lower()
    ftype = ftype.upper()
    
    if view not in ["view1", "view2"]:
        raise ValueError(" Provide the Dataset view (either view1 or view2) to use for the training ")
    if ftype not in ["LBP", "LTP", "LQP"]:
        raise ValueError(" Wrong Feature Type, choose any one of the following features [LBP, LTP, LQP] ")
    ddir = os.path.join(featdir, "data");
    if not os.path.exists(ddir):
        print " Creating following Data directory : ", ddir
        os.mkdir(ddir)
        
    rfile = os.path.join(ddir, view) + "-" + ftype + "-" + distance + "-simple-threshold-model-results.npy"
    
    if os.path.exists(rfile):
        accuracy = np.load(rfile)
        if view == "view1":
            print "For Feature Type = %s, Classification accuracy(on view = %s)= %f" % (ftype, view, accuracy)
        else:
            print (" \n ----------------------- \n" 
            "For Simple Model(Without PCA) = Classification Accuracy of %s's on %s \n-----------------------\n" % (ftype, view))
            print "Set#    Accuracy"
            mean = 0;
            for count, val in enumerate(accuracy):
                mean = mean + val
                print str(count + 1) + ".\t %f" % val
            print " Mean =", (mean / (count + 1) * 100)
        
    else:
        data = np.load(os.path.join(os.pardir, os.path.join("data" , view + '.npy')))
        if view == "view1":
            accuracy = train_test(data, featdir, ftype, 'view1', distance=distance)
            np.save(rfile, accuracy)
        else:  # check that
            accuracy = []
            
            for k in range(10):  # 10 fold
                # extract  data
                # for each iteration data is extracted as [data[0][0][0],data[1][0][0],data[2][0][0],data[3][0][0],data[4][0][0],data[5][0][0]]           
                tdata = []    
                for t in range(6):
                    #  tdata.append(data[t][0][k])
                    tdata.append(data[t][k])
 
                #                  
                iaccuracy = train_test(tdata, featdir, ftype, 'view2', '-Set ' + str(k), distance)
                accuracy.append(iaccuracy)        
            
            print (" \n ----------------------- \n" 
            "For Simple Model(Without PCA) = Classification Accuracy of %s's on %s \n-----------------------\n" % (ftype, view))
            print "Set#    Accuracy"
            mean = 0;
            for count, val in enumerate(accuracy):
                mean = mean + val
                print str(count + 1) + ".\t %f" % val
            print " Mean =", (mean / (count + 1) * 100)
            
            np.save(rfile, accuracy)
def learn_model(options):
    ''' Learn PCA or Without PCA model on LFW dataset '''
    
    check_values(options.view, ["view1", "view2", "complete"])
    check_values(options.ftype, ["LBP", "LTP", "LQP"])
    check_values(options.ttype, ["with-pca", "without-pca"])
    check_values(options.dist.lower(), ["cosine", "l2", "chi-square"])    
    if not os.path.exists(options.featdir):
        raise ValueError("Feature directory (%s) does not exist." + 
         "Please run the feature computation before running the training/evaluation code" % options.featdir)
    

    print '''#############################################################'''
    print '''#                Summary of Given Options                   #'''
    print '''#############################################################'''
    print 'Running on   %s dataset %s ' % (options.view,
                                           (" (View 1 is used for parameter tuning and view 2 for evaluation" if options.view == "complete" else ""))
    print 'Trianing a model  %s on %s features found in the directory %s ' % (options.ttype, options.ftype, options.featdir)
    if options.ttype == 'with-pca':
        print  'PCA model is trained using following components ', options.pcadim
    print 'Features are matched using =  %s' % ("cosine similarity" if options.dist == "cosine" else options.dist + " distance")
          
    print 'Learned Model will be stored in directory: %s ' % os.path.join(options.featdir, 'data')
    print '''#############################################################'''
    
    
    if options.ttype == "with-pca":
        pcadim = get_list_values(options.pcadim)
        pcadim = check_islist(pcadim)
        
        if options.view in ["view1", "view2"]:
            print "Training and Evaluating a PCA Model with following given dimensions  = ", pcadim
            fpcadim, accuracy = LFW_PCA(options.view, options.featdir, options.ftype, pcadim, options.dist)
        else:
            print " Doing a complete trainig cycle \n ----------> First: Doing Parameter tuning on view1 <---------- "
            
            fpcadim, accuracy = LFW_PCA("view1", options.featdir, options.ftype, pcadim, options.dist)
#            if not isinstance(fpcadim, list):
#                fpcadim = [fpcadim]
            print " ----------> Second : Evaluating the model with PCA dimensions =%d  on view2 <---------- " % fpcadim
            fpcadim, accuracy = LFW_PCA("view2", options.featdir, options.ftype, int(fpcadim), options.dist)      
    else:
        if options.view == "complete":
            raise ValueError(' Without-PCA: Model can be only used with view1 or view2 option')
        print "Training and Evaluating a Model without PCA  "
        LFW_WPCA(options.view, options.featdir, options.ftype, options.dist)    
    
    
def main(argv=None):
    ''' Parse command line/Configuration file  arguments '''
    if argv is None:
        argv = sys.argv
    parser = argparse.ArgumentParser(description='''Face Verification on LFW''', version=1.0)
    ggroup = parser.add_argument_group('General Options')  # general Group
    ggroup.add_argument('-V', '--view',
                        dest='view',
                        help="Choice of the dataset, options cans be \n"
                        + " view1: used for parameter tuning  purposes\n"
                        + " view2: used only for model evaluation \n"
                        + " complete: a model parameters will be first tuned on view1  "
                        + "  and results will be reported on view2  \n",
                        default='view1')
    ggroup.add_argument('-f', '--features', dest="ftype",
                    default="LQP", help="Choose among LBP, LTP, LBP+LTP or LQP feature type")
    ggroup.add_argument('-t', '--train-type', dest="ttype",
                        default="with-pca", help="Choice of Training with or without PCA (for feature evaluation)\n"
                         + "Available options are with-pca or without-pca")
    ggroup.add_argument('-d', '--dir',
                        dest='featdir', default="/tmp/",
                        help="directory path where computed features have been stored")                   
    ggroup.add_argument('-D', '--distance',
                        dest='dist', default="cosine",
                        help="Distance metric for comparing features. Available ptions are cosine, chi-square and L2")
    ggroup.add_argument('-p', '--pca-dim',
                        dest='pcadim',
                        help=('''Number of PCA components\n'''
                        '''[you can pass a list i.e. -p='[100, 200, 300]' in which case list will be'''
                         '''used for picking best performing number on view1 ]'''), default='[1600]')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)          

    try:       
        options = parser.parse_args()
       
        options.dist = options.dist.lower()
        
#        print ''' Now Building Model........... '''
    
    except ValueError:
                print 'Please provide valid options, for help please use following command: python trainAndEvaluateLFW.py -h'
    learn_model(options)


if __name__ == "__main__":
    sys.exit(main())
#   crop_images('/home/hussain/dataset/', 80 , 150, '/tmp/', 10, 'LFW');    
#   LFW_WPCA('view2', '/home/hussain/research/code-release/new-data-xoffset=1-yoffset=-4/lbp-norm-1/', 'LBP')
#    try:
#    LFW_PCA('view2', '/home/hussain/research/code-release/new-data-xoffset=1-yoffset=-4/LTP-tolerance=5', 'LTP', 1000, 1000)
    
    
#    except Exception as ex:
#        print ex
#        raw_input()

#    ''' e.g. python trainAndEvaluateLFW.py -t with-pca -f LTP --view view1'''
#                                     ''' -d /tmp/LTP-tolerance=5 -p="[200,500]"\n
#                                     will train and evaluate LTP features based model on LFW dataset'''
#    print 'ARGV:', sys.argv[1:]
    
