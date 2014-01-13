"""
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
"""
"""This file contain wrapper class for CPP code for feature computation with different arguments
   Can be used as a module, or with command line options... see main 
"""
from struct import *
import numpy as np
import cPickle as cp
from scipy.io import loadmat
from scipy.io import savemat
import datetime, time
import os, sys, socket, string
from shutil import *
import multiprocessing, argparse
import util as u
# from sepolgen.yacc import compute_lookback_includes
class Features:
    """
        Main Wrapper Class used to compute different types of features by calling CPP Code:
        Currently Supported features are LQP, LBP, LTP and LBP+LTP
    """
    def __init__(self, options, exefile='mainFeatures'):    
        """
            Initialize the different parameters for calling CPP program and provides
            main interface function to compute features...
            ftype= features type among [LBP, LTP, LBP+LTP]
            featdir= directory path for storing feature files.
            listfile= list file containing absolute path of images for computing features... 
        """
        self.display(options)
        ftype = options.ftype.upper()
        
        if not ftype in ["LQP", "LBP", "LTP", "LBP+LTP"]:
            print 'Error: Any one of the following Feature Types should be chosen  "LQP","LBP","LTP","LBP+LTP", not %s' % ftype
            raise ValueError
        
        self.ftype = options.ftype;
        if not os.path.exists(options.featdir):
            print 'Making Directory =', options.featdir
#             os.mkdir(options.featdir)
            u.create_directory(options.featdir)
        self.featdir = options.featdir;
        
        if not os.path.exists(options.listfile):
            print "Error: Given file(%s) for list of images does not exist " % options.listfile
            raise ValueError
        self.listfile = options.listfile
        self.rdir = os.getcwd()  # root directory
        self.exe = os.path.normpath(os.path.join(self.rdir, os.path.join(os.pardir, os.pardir), 'build', exefile))
        
        values = {'exefile':self.exe, 'color_channel':4,
                   'add_rgb':1, 'norm_sep':0, 'cell_size':options.cellsize,
                   'lbp_stride':options.cellsize, 'width':options.width, 'height':options.height,
                   'listfile':options.listfile, 'cbfile':options.cbfile }
        cmd = string.Template("${exefile} --Color-Channel=${color_channel} "
                              " --Add-RGBLBP=${add_rgb}  --Norm-Sep=${norm_sep} "
                              " --Cell-Size=${cell_size} --LBP-Stride=${lbp_stride} "
                              " --Win-Width=${width} --Win-Height=${height} --Validation-File=${cbfile}"
                              " --TrainingFile=${listfile}")
        cmd = cmd.safe_substitute(values)
       
        self.dirs = []  # contain the name of feature directories...
        cwd = os.getcwd()
        if self.ftype == "LQP":
            self.compute_lqp(cmd, options)
        else:
            self.compute_localpattern(cmd, options)
        os.chdir(cwd)
    def get_dir_names(self):
        return self.dirs
    def display(self, options):
        ''' display the options used for feature computations ''' 
        print '''#############################################################'''
        print '''#                Summary of Options for Computing Features  #'''
        print '''#############################################################'''
        print 'Features =  %r' % options.ftype
        print 'List File =  %r' % options.listfile
        print 'Codebook File (list of images used for LQP codebook learning) =  %r' % options.cbfile
        print 'Directory for Feature Storage =  %r' % options.featdir
        
        print 'Image Width =  %r' % options.width
        print 'Image Height =  %r' % options.height
        print 'Cell Size =  %r' % options.cellsize
        if options.ftype.upper() in ['LQP', 'LTP', 'LBP+LTP']:
            print 'Tolerance =  %r' % options.tol
        
        if options.ftype.upper() == 'LQP':
            print 'LQP Type =  %r' % ("Circular" if options.lqptype == 2 else "Hor+Ver+Diag+ADiag") 
            print 'LQP Size =  %r' % options.lqpsize
            print 'Encoding=  %r' % ['Binary', 'Ternary', '', '', 'Split-Ternary'][options.coding]
            print 'Code Book Size (Number of Visual Words) =  %r' % options.cbsize
        self.ftype = options.ftype.upper()
        print '''#############################################################'''
        
        print ''' Now Extracting Features........... '''
    def check_features_computed(self):
        """
            Functions check whether features have already been computed or not...
            
        """
        nimages = u.get_number_lines(self.listfile)
        count = 0;
        for filename in os.listdir(os.getcwd()):
            if os.path.splitext(filename)[1] == '.jpg-' + self.ftype or os.path.splitext(filename)[1] == '.png-' + self.ftype:
                count = count + 1;
        return nimages == count
    def compute_features(self, cmd):
        """
            Function calls the CPP exe to compute features and logs the output information into a log file
            cmd= calling command with set of arguments.
            On successful completion produces descriptors for given list of images... 
        """
        if self.check_features_computed():
            print "Features have already been computed.... "
            return 
        logfile = os.getcwd() + os.path.sep + "record.log"
        fid = open(logfile, 'wt');
        fid.write("Running on " + socket.gethostname())
        fid.write("\nIn the Dir = " + os.getcwd())
        today = datetime.date.today()
        date = today.strftime('%a %b %d %Y')
        fid.write("\nOn Following Date=" + date)
        print "Appended Command = ", cmd
        fid.write("\nAppended Command = " + cmd)
        fid.close();
    
        # Different Program option values, as dictionary items used during 
        
        cmd = cmd + " >> " + logfile
        st = time.time()
        print " Computing Features ... \n", cmd
        os.system(cmd)
        tt = time.time() - st
        print " Finished Computing Features, Total Time Taken = %f sec (%f min)" % (tt, tt / 60)

    def compute_lqp(self, addcmd, options):
        """
            Compute LQP Features, 
            can be configured to run with different set of LQP types, coding types, tolerance, codebook size etc.
        """  
        values = { 'feature_type':14, 'cb_distance':0,
                 'ncluster_rounds':10,
                 'code_prune_count':10,
                 'lqp_stride':1, 'norm':1, 'patch_type':options.lqptype, 'code_type':options.coding}
        cmd = string.Template(" --FeatureType=${feature_type}  --PatchType=${patch_type} "\
                            " --CodeBook-DMetric=${cb_distance} --ClusteringRounds=${ncluster_rounds} "\
                            "--Patch-PruneCount=${code_prune_count} --Patch-PatchStride=${lqp_stride} "\
                            " --Patch-CodingType=${code_type} --Normalization=${norm} ")
        cmd = addcmd + cmd.safe_substitute(values)
         
        lqpsize = u.check_islist(options.lqpsize)
        cbsize = u.check_islist(options.cbsize)
        tolerance = u.check_islist(options.tol)
        featdir = options.featdir
                
        count = 1
        for tol in tolerance:
            for lqps in lqpsize:
                for cbs in cbsize:
                    os.chdir(featdir)
                    dirname = 'lqp-size=' + str(lqps) + '-codebooksize=' + str(cbs) + '-tolerance=' + str(tol)
                    dirname = os.path.join(featdir, dirname)
                    if not os.path.exists(dirname):
                        print 'Making Directory =', dirname
                        u.create_directory(dirname)
                    self.dirs.append(dirname)
                    os.chdir(dirname)              
                    values = {'tol':tol, 'cbsize':cbs, 'patch_size':lqps}
                    addcmd = string.Template(" --LTP-Tolerance=${tol} --CodeBookSize=${cbsize}  --Patch-PatchSize=${patch_size} ")
                    addcmd = cmd + addcmd.safe_substitute(values)
                    print  "Running Process Number = ", count
                    count = count + 1
                    self.compute_features(addcmd)
#                    p = multiprocessing.Process(target=self.compute_features, args=(addcmd,))
#                    p.daemon = True  # Make it a daemon process
#                    p.start()
    
    def compute_localpattern(self, addcmd, options):
        """
            Function is used to compute Local Pattern (LBP, LTP, Features:
            
        """
        features = {"LBP":1, "LTP":2, "LBP+LTP":3}
        values = {'feature_type': features[self.ftype], 'norm':1}
        cmd = string.Template(" --FeatureType=${feature_type}  --Normalization=${norm} ")
        cmd = addcmd + cmd.safe_substitute(values)
        tolerance = options.tol  # used only for ltp and lbp+ltp features, lbp features ignores this options
        
        count = 1
        if self.ftype.upper() != "LBP":
            for tol in tolerance:
                os.chdir(self.featdir)
                dirname = self.ftype + '-tolerance=' + str(tol)
                dirname = os.path.join(self.featdir, dirname)
                if not os.path.exists(dirname):
                    print 'Making Directory =', dirname
                    u.create_directory(dirname)
                self.dirs.append(dirname)
                os.chdir(dirname)              
                values = {'tol':tol}
                addcmd = string.Template(" --LTP-Tolerance=${tol} ")
                addcmd = cmd + addcmd.safe_substitute(values)
                print  "Running Process Number = ", count
                count = count + 1
                self.compute_features(addcmd)
        else:
                os.chdir(self.featdir)
                dirname = self.ftype
                dirname = os.path.join(self.featdir, dirname)
                if not os.path.exists(dirname):
                    print 'Making Directory =', dirname
                    u.create_directory(dirname)
                self.dirs.append(dirname)
                os.chdir(dirname)              
                self.compute_features(cmd)
#            p = multiprocessing.Process(target=self.compute_features, args=(addcmd,))
#            p.daemon = True  # Make it a daemon process
#            p.start()
def main(argv=None):
    """ Parse Command line arguments and call the computation function... """
    # Do argv default this way, as doing it in the functional
    # declaration sets it at compile time.
    if argv is None:
        argv = sys.argv
    parser = argparse.ArgumentParser(description=" Feature Computation ...", version=1.0)
    ggroup = parser.add_argument_group('General Options')  # general Group
    ggroup.add_argument('-f', '--features', dest="ftype",
                        default="LQP", help="Choose among LBP, LTP, LBP+LTP or LQP feature type")
    ggroup.add_argument('-l', '--listfile',
                        dest='listfile',
                        help="a list file containing list of cropped images for which  features will be computed ")
    ggroup.add_argument('-d', '--dir',
                        dest='featdir',
                        help="directory path for storing feature files")                   
    ggroup.add_argument('--width',
                        dest='width',
                        help="Width of Cropped Images ", default=80, type=int)
    ggroup.add_argument('--height',
                        dest='height',
                        help="Height of Cropped Images ", default=150, type=int)
    ggroup.add_argument('-c', '--cellsize',
                        dest='cellsize',
                        help="cellsize of histogram grid ", default=10, type=int)
    ggroup.add_argument('-t', '--tolerance',
                        dest='tol',
                        help="tolerance values used for LTP or LQP features (can be a list) e.g. -t='[5,7]'",
                        default=[5])
    
    lqpgroup = parser.add_argument_group('LQP Options')  # LQP options
    lqpgroup.add_argument('--lqptype',
                        dest='lqptype',
                        help="LQP Type 2. Circular, 9. Hor+Ver+Diag+ADiag",
                        default=2, type=int)
    lqpgroup.add_argument('-s', '--lqpsize',
                        dest='lqpsize',
                        help="LQP size radius of LQP Disk or other types (can be a list), e.g. -s='[5,7]'",
                        default=[5])
    lqpgroup.add_argument('-e', '--encoding',
                        dest='coding',
                        help="0. Binary, 1. Ternary, 4. Split-Ternary",
                        default=4, type=int) 
    lqpgroup.add_argument('-w', '--words',
                        dest='cbsize',
                        help="Codebook size (number of visual words) used for LQP computation (can be a list), e.g. -w='[100,150]'",
                        default=[150])
    lqpgroup.add_argument('--cbfile',
                        dest='cbfile',
                        help="a list file containing list of images for learning the codebook", default="")
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)             
    try:
        options = parser.parse_args()
#        print options
        u.check_values(options.ftype, ["LBP", "LTP", "LQP", "LBP+LTP"])
        options.cbsize = u.get_list_values(options.cbsize)
        options.lqpsize = u.get_list_values(options.lqpsize)
        options.tol = u.get_list_values(options.tol)       
        feat = Features(options)    
    except ValueError:
            print 'Please provide valid options, for help please use following command: python Features.py -h'
    
if __name__ == "__main__":
    import run as r
    parser = r.ArgumentParser("FeaturesComputation", 1.0)
    sys.exit(main())
#    print 'ARGV:', sys.argv[1:]
   

    # ftype, featdir, listfile
    
# python Features.py --help
# ARGV: ['--help']
# usage: Features.py [-h] [-v] [-f FTYPE] [-l LISTFILE] [-d FEATDIR]
#                          [--width WIDTH] [--height HEIGHT] [-c CELLSIZE]
#                          [-t TOL] [--type LQPTYPE] [-s LQPSIZE] [-w CBSIZE]
#                          [-e CODING]
#
# Feature Computation ...
#
# optional arguments:
#  -h, --help            show this help message and exit
#  -v, --version         show program's version number and exit
#
# General Options:
#  -f FTYPE, --features FTYPE
#                        Choose among LBP, LTP, LBP+LTP or LQP feature type
#  -l LISTFILE, --listfile LISTFILE
#                        list file containing list of cropped images to compute
#                        features
#  -d FEATDIR, --dir FEATDIR
#                        directory path for storing feature files
#  --width WIDTH         Width of Cropped Images
#  --height HEIGHT       Height of Cropped Images
#  -c CELLSIZE, --cellsize CELLSIZE
#                        cellsize of histogram grid
#  -t TOL, --tolerance TOL
#                        tolerance values used for LTP or LQP features
#
# LQP Options:
#  --type LQPTYPE        LQP Type 2. Circular, 9. Hor+Ver+Diag+ADiag
#  -s LQPSIZE, --size LQPSIZE
#                        LQP size radius of LQP Disk or other types
#  -w CBSIZE, --words CBSIZE
#                        Codebook size (number of visual words) used for LQP
#                        computation
#  -e CODING, --encoding CODING
#                        0. Binary, 1. Ternary, 4. Split-Ternary
# python Features.py -f LQP -l /home/hussain/datasets/LFW/lfw-list-xoffset=1-yoffset=-4.txt  -d /home/hussain/research/code-release/new-data-xoffset=1-yoffset=-4/ --width=80 --height=150 -c=10 -t=5 --type= 2 -s=5 -w=150 -e=4
