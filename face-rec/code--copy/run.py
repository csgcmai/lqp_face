"""
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
"""
''' A simple class to process commandline and configuration file options to run code...
 ''' 
from struct import *
import numpy as np
import time
import sys
import os
import util as u
import ConfigParser
import Features as feat
import trainAndEvaluateLFW as tlfw
import argparse
from pylab import *
class CommentlessFile(file):
    ''' A class for removing comments from the config file '''
    def readline(self):
        line = super(CommentlessFile, self).readline()
        if line:
            line = line.split(';', 1)[0].strip()
            line = line.split('#', 1)[0].strip()
            return line + '\n'
        else:
            return ''
def prepare_list(idir, odir, dataset='train-val'):
    ''' Function generates list of images and write them to two different files
        
        imglistfile = contain list of all the images present in the dataset
        cblistlist = contain list of images present in view1 dataset...
        dataset = Which data to use from view1 for learning codebook, complete or only train-val
     '''
    u.check_values(dataset, ['train-val', 'complete'])
    imglist = u.get_imlist(idir)
    imglistfile = os.path.join(odir, "lfw-list-images.txt");
    cblistfile = os.path.join(odir, "lfw-view1.txt");
    u.write_list(imglistfile, imglist, idir)
    # load the view1 data
    pcaTrain, pcalabel, thTrain, thlabel, test, testlabels = np.load(os.path.join(os.pardir, os.path.join("data" , 'view1.npy')))
    if dataset == 'complete':
        imglist = np.concatenate((np.concatenate((pcaTrain, thTrain)), test));
    else:
        imglist = np.concatenate((pcaTrain, thTrain));
#     print len(imglist)
#    imglist=[i[0] for i in imglist.flatten()]
    imglist = list(set([i[0] for i in imglist.flatten()]))
    imglist = [i + ".jpg" for i in imglist]
    u.write_list(cblistfile, imglist, idir)
    
    return imglistfile, cblistfile
class ArgumentParser:
    ''' A generic class to parse command-line arguments for different sub-modules  ''' 
    def __init__(self, v=1.0):
        ''' Initializes Argument Parser class to 
        parse a particular class type arguments, options can be: [features,learning,complete] '''
        
        self.args = sys.argv[1:]
        u.check_data()  # first check that data with split pairs information exist
        config_parser = argparse.ArgumentParser(
                                              description="Local Pattern Based Face Verification", add_help=False, version=v)
        pgroup = config_parser.add_argument_group("Required Arguments")
        pgroup.add_argument("-c", "--configfile", dest="configfile",
                                 help="Please specify a config file,"
                                 + " if absent specify the module type (see below) and command line arguments with provided (or default) values will be used ")
        pgroup.add_argument("-m", "--module", dest="module",
                                 help="Please specify your choice for the system stage, options are:"
                                 + " features --> only features are computed \n"
                                 + " learning --> model is learned over computed features "
                                 + " complete --> a complete system i.e. features+learning ")
                   
        
        args, self.args = config_parser.parse_known_args(self.args)
        self.config = None


        self.type = "complete"  
        if args.configfile:
            self.config = ConfigParser.SafeConfigParser()
            self.config.readfp(CommentlessFile(args.configfile))
#            if not "Learning" in  self.config.sections() and "Features" in self.config.sections():
#                show_learning_help=False;
#        elif:
            if "Features" in self.config.sections() and "Learning" in self.config.sections():
                self.type = "complete"
            else : 
                self.type = "Features" if "Features" in self.config.sections() else "Learning"
            self.type = self.type.lower()
        elif not args.module is  None:   
            self.type = args.module.lower()
            u.check_values(self.type, ["complete", "general", "features", "learning"])
#        else:
#            print 'Please provide a module option OR a config file'
#            print config_parser.print_help()
#            sys.exit(1)
                
        self.gparser = argparse.ArgumentParser(
                                               # Inherit options from config_parser
                                                parents=[config_parser],
                                                # print script description with -h/--help
                                                description=__doc__,
                                                # Don't mess with format of description
                                                # formatter_class=argparse.ArgumentDefaultsHelpFormatter
                                                formatter_class=argparse.RawDescriptionHelpFormatter, add_help=False,
                                              )

        # parse arguments...
        self.options = self.parse_general()
        if self.type == "features":
            self.foptions = self.parse_features()
        elif self.type == "learning":
            self.loptions = self.parse_learner()
        else:
            self.foptions = self.parse_features()
            self.loptions = self.parse_learner()
        self.execute()
        
    def execute(self):
        ''' Function calls different modules to perform different tasks '''
        # generic options compute cropped images, generate a list of cropped images and list of codebook
        options = self.options
        cidir = u.get_cropped_images_dir_name(options)
        if not os.path.exists(cidir) or len(os.listdir(cidir)) < len(os.listdir(options.idir)):
            print "--------- Cropping Images ------------"
            u.crop_images(options)
            mim = u.mean_image(cidir)
            mimname = os.path.join(os.path.normpath(os.path.join(cidir, os.path.pardir)), 'meanim.png')
            imsave(mimname, mim, cmap=cm.gray);
        imglistfile, cblistfile = prepare_list(cidir, options.odir, options.cbdataset)
        if self.type in [ "features", "complete"]:
            print "--------- Computing  Features ------------"
            if len(self.foptions.cbfile) == 0:
                self.foptions.cbfile = cblistfile
            if len(self.foptions.listfile) == 0:
                self.foptions.listfile = imglistfile
            features = feat.Features(self.foptions)
            featdir = features.get_dir_names()
        if self.type in [ "learning", "complete"]:
            print "--------- Learning Model ------------"
            if self.type == "learning":
                featdir = self.loptions.featdir
#             else:
#                 featdir = u.get_dirs(self.foptions.featdir)
            if  type(featdir) != list:
                featdir = [featdir]
            if featdir == None or len(featdir) == 0:
                raise ValueError('Feature directory is empty Please provide valid  path for features directory ')
            
            for fd in featdir:
                self.loptions.featdir = fd
                tlfw.learn_model(self.loptions)
    def parse_general(self):
        ''' Parse general options '''
        defaults = {
                   "idir" : "",
                   "odir" : "",
                   "dataset":"LFW",
                   "width":80,
                   "height":150,
                   "padding":10,
                   "xoffset":1,
                   "yoffset":-4,
                   "cbdataset":"train-val",
                   "ftype":"LQP",
                   "usergb":True
                  }
        if not self.config is None:
            defaults = dict(self.config.items("General"))
        self.gparser.set_defaults(**defaults)

        ggroup = self.gparser.add_argument_group('General Options')  # general Group
        ggroup.add_argument('-i', '--idir',
                        dest='idir',
                        help="directory path where dataset (e.g. LFW-a) images are stored")
        ggroup.add_argument('-o', '--odir',
                        dest='odir',
                        help="directory path where learned models and cropped images will be stored")
        ggroup.add_argument('--dataset',
                            dest='dataset',
                            help=('''Images belong to which dataset, LFW or FERET'''))
        ggroup.add_argument('--width',
                            dest='width',
                            help="Width of cropped image", type=int)
        ggroup.add_argument('--height',
                            dest='height',
                            help="height of cropped image", type=int)
        ggroup.add_argument('--padding',
                            dest='padding',
                            help="Size of Padding to be added around the cropped image (Must be of the cell size)", type=int)
        ggroup.add_argument('-x', '--x-offset',
                            dest='xoffset',
                            help=('''(+,-) offset position of cropping rectangle from the center image '''),
                            type=int)
        ggroup.add_argument('-y', '--y-offset',
                            dest='yoffset',
                            help=('''(+,-) offset position of cropping rectangle from the center image '''),
                            default=-4, type=int)
        ggroup.add_argument('--features', dest="ftype",
                            help="Choose among LBP, LTP, LBP+LTP or LQP feature type")
        ggroup.add_argument('--codebook-dataset',
                        dest='cbdataset',
                        help="(For LQP) Which portion of view1 to use for codebook learning, option can be either train-val or compelete")
        ggroup.add_argument('--usergb', dest="usergb",
                            help="Use color information during feature computations (by default = True) ")
       
        try:
            options, self.args = self.gparser.parse_known_args(self.args)
            u.check_values(options.dataset, ["FERET", "LFW"])
            options.idir = u.flattened_images_subdir(options.idir, options.odir, options.usergb)
            print options
            return options
        except ValueError:
                print 'Please provide valid options, for help please use following command: python run.py -h'

    def parse_features(self):
        ''' Function parses the feature subsection (commandline) arguments '''
#        if argv is None:
#            argv = sys.argv
#    #these are default values used if no config file or command line arguments are given...
        defaults = {
#                    "ftype" : "LQP",
                   "listfile" : "",
#                   "featdir":"",
#                   "width":80,
#                   "height":150,
                   "cellsize":10,
                   "tol":5,
                   "lqptype":2,
                   "lqpsize":5,
                   "coding":4,
                   "cbsize":150,
                   "cbfile":""
                  }
        
        if not self.config is None:
            defaults = dict(self.config.items("Features") + self.config.items("LQP"))
        
        self.fparser = argparse.ArgumentParser(
                                                parents=[self.gparser],
                                                description=__doc__,
                                                formatter_class=argparse.RawDescriptionHelpFormatter,
                                                add_help=False if self.type == 'complete' else True 
                                              )            
 
        self.fparser.set_defaults(**defaults)
         
              
        ggroup = self.fparser.add_argument_group('Features Configuration Options')  # general Group
        
        ggroup.add_argument('--listfile',
                            dest='listfile',
                            help="a list file containing list of cropped images for which  features will be computed ")
        ggroup.add_argument('--cellsize',
                            dest='cellsize',
                            help="cellsize of histogram grid ", type=int)
        ggroup.add_argument('--tolerance',
                            dest='tol',
                            help="tolerance values used for LTP or LQP features (can be a list) e.g. --tolerance='[5,7]'")
        
        
        lqpgroup = self.fparser.add_argument_group('LQP Options')  # LQP options
        lqpgroup.add_argument('--lqptype',
                            dest='lqptype',
                            help="LQP Type 2. Circular, 9. Hor+Ver+Diag+ADiag",
                            type=int)
        lqpgroup.add_argument('--lqpsize',
                            dest='lqpsize',
                            help="LQP size radius of LQP Disk or other types (can be a list), e.g. -lqpsize='[5,7]'")
        lqpgroup.add_argument('--coding',
                            dest='coding',
                            help="0. Binary, 1. Ternary, 4. Split-Ternary",
                            type=int) 
        lqpgroup.add_argument('--words',
                            dest='cbsize',
                            help="Codebook size (number of visual words) used for LQP computation (can be a list), e.g. -w='[100,150]'")
        lqpgroup.add_argument('--cbfile',
                            dest='cbfile',
                            help="a list file containing list of images for learning the codebook")      
                   
        try:
            options, self.args = self.fparser.parse_known_args(self.args)
            print options
            u.check_values(options.ftype, ["LBP", "LTP", "LQP", "LBP+LTP"]) 
            options.cbsize = u.get_list_values(options.cbsize)
            options.lqpsize = u.get_list_values(options.lqpsize)
            options.tol = u.get_list_values(options.tol)   
            options.featdir = os.path.join(self.options.odir, "features");
            return options
        except ValueError:
                print 'Please provide valid options, for help please use following command: python Features.py -h'
        
    def parse_learner(self):
        ''' Function parses the learning subsection (commandline) arguments '''
        defaults = {
                   "view" : "complete",
                   "ttype" : "with-pca",
                   "dist":"cosine",
                   "pcadim":'[100,200,500,750,1000,1200,1400,1500,1600,1700,1800,1900,2000]'
                  }
        
        self.lparser = argparse.ArgumentParser(# Inherit options from config_parser
                                                    parents=[self.fparser if self.type == "complete" else self.gparser],
                                                    # print script description with -h/--help
                                                    description=__doc__,
                                                    # Don't mess with format of description
                                                    formatter_class=argparse.RawDescriptionHelpFormatter)
        if not self.config is None:
            defaults = dict(self.config.items("Learning"))
 
        self.lparser.set_defaults(**defaults)
         
              
        
        ggroup = self.lparser.add_argument_group('Learning Options')  # general Group
        ggroup.add_argument('-V', '--view',
                            dest='view',
                            help="Choice of the dataset, options cans be \n"
                            + " view1: used for parameter tuning  purposes\n"
                            + " view2: used only for model evaluation \n"
                            + " complete: a model parameters will be first tuned on view1  "
                            + "  and results will be reported on view2  \n")
        if self.type == "learning":
            ggroup.add_argument('--featdir',
                        dest='featdir',
                        help="directory path where computed features have been stored")
        
        ggroup.add_argument('--train-type', dest="ttype",
                            help="Choice of Training with or without PCA (for feature evaluation)\n"
                             + "Available options are with-pca or without-pca")                      
        ggroup.add_argument('--distance',
                            dest='dist',
                            help="Distance metric for comparing features. Available options are cosine, chi-square and L2")
        ggroup.add_argument('--pca-dim',
                            dest='pcadim',
                            help=('''Number of PCA components\n'''
                            '''[you can pass a list i.e. -p='[100, 200, 300]' in which case list will be '''
                             '''used for picking best performing number on view1 ]'''))        
        try:
            options, self.args = self.lparser.parse_known_args(self.args)
#            options = self.parser.parse_known_args()
            print options 
            return options
        except ValueError:
                print 'Please provide valid options, for help please use following command: python run.py -h'
    
#    u.write_list(cblistfile, imglist)
#==============================================================================
def main(argv=None):
    ''' Function calls different scripts to learn the complete model on LQP features...'''
    ArgumentParser(1.0)
    
if __name__ == "__main__":
#     prepare_list('/home/hussain/datasets/LFW/lfwa', '/tmp')
     print 'ARGV:', sys.argv[1:]
     sys.exit(main())
##==============================================================================
