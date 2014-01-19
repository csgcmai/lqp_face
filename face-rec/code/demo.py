#!/bin/python
"""
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
"""


# Test script for testing the features...
# This script first generates the configuration files and then runs experiments
import run as r
import os, sys
from ConfigParser import SafeConfigParser
def getConfigParameters():
    # print os.getcwd()
    parser = SafeConfigParser()
    parser.read('face-rec/config_LTP_view1_wpca.py')
    odir = parser.get('General', 'odir').split('=', 1)[0]
    ftype = parser.get('General', 'ftype').split('=', 1)[0]
    return odir, ftype
    
def parseConfigFile(idir='', odir='./'):
    '''
        Generate configuration files for the test set.
        These configuration files are for LTP features.
    '''
    pwd = os.getcwd();
    os.chdir('../');
    parser = SafeConfigParser()
    parser.read('config.py')
    # Configuration files for three different expeirments
    # 1. LTP Features with Chi-square similarity (without PCA learning) on view1
    # 2. LTP Features with Chi-square similarity (without PCA learning) on view2.
    # 3. LTP Features with cosine similarity on view1 (for cross-validation) and view2 (testing).
    cfnames = ['config_LTP_view1_wpca.py', 'config_LTP_view2_wpca.py', 'config_LTP_pca.py']
    if len(idir) > 0:
        parser.set('General', 'idir', idir)
    
    parser.set('General', 'odir', odir)
    parser.set('General', 'ftype', 'LTP')
    
    for fname in cfnames:
        wf = open(fname, 'w');
        if fname == 'config_LTP_view1_wpca.py':
            parser.set('Learning', 'view', 'view1');
            parser.set('Learning', 'ttype', 'without-pca')
            parser.set('Learning', 'dist', 'chi-square')
        elif fname == 'config_LTP_view2_wpca.py':
            parser.set('Learning', 'view', 'view2');
            parser.set('Learning', 'ttype', 'without-pca')
            parser.set('Learning', 'dist', 'chi-square')          
        else:
            parser.set('Learning', 'view', 'complete')
            parser.set('Learning', 'ttype', 'with-pca')
            parser.set('Learning', 'dist', 'cosine')
        parser.write(wf)
        wf.close()
    os.chdir(pwd)
    return cfnames
    
    for section_name in parser.sections():
        print 'Section:', section_name
        print ' Options:', parser.options(section_name)
        for name, value in parser.items(section_name):
            print ' % s = % s' % (name, value)
        print
def main():
    ''' main testing function... '''
#    print len(sys.argv)
    if not (len(sys.argv) >= 2 and len(sys.argv) <= 4):
        print "Invalid Number of Arguments: Please provide path to output directories [input directory path optional], i.e."
        print "python demo.py path_to_odir path_to_idir[optional] number_of_experiments[optional]"
        return
    else:
        odir = sys.argv[1]
        print 'Output Directory = ', odir
        idir = ''
        nexp = 3
        if len(sys.argv) == 3:
            idir = sys.argv[2]
            print 'Input Directory = ', idir
        if len(sys.argv) == 4:
            nexp = int(sys.argv[3])
        
        print 'Number of Experiments = ', nexp
            
        if idir != '' and not os.path.exists(idir):
            raise ValueError("Input Image directory %s does not exist" % idir)
#         print "Idir= %s, Odir=%s" % (idir, odir)
        # cfnames = genConfigFiles(idir, odir)
        cfnames = parseConfigFile(idir, odir)
        for i  in range(0, min(nexp, len(cfnames))):
            print "Running Test case Number %d" % (i + 1)
            sys.argv = ['', '--configfile=' + os.path.join(os.path.pardir, cfnames[i])];
            print sys.argv
            # now call the running routine...
            r.main()   
          
     
if __name__ == "__main__":
    #	print 'ARGV:', sys.argv[1:]
    sys.exit(main())       
