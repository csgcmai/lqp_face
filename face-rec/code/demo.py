#!/bin/python
# Test script for testing the features...
# This script first generates the configuration files and then perform the learning of data
import run as r
import os, sys
from ConfigParser import SafeConfigParser
def getOutputDir():
    parser = SafeConfigParser()
    parser.read('../config_LTP_view1_wpca.py')
    return parser.get('General', 'odir').split('=', 1)[0]
    
def parseConfigFile(idir='', odir='./'):
    '''
        Generate configuration files for the test set.
        These configuration files are for LTP features.
    '''
    pwd = os.getcwd();
    os.chdir('../');
    parser = SafeConfigParser()
    parser.read('config.py')
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
# 
# def genConfigFiles(idir='./', odir='./'):
#     '''
#         Generate configuration files for the test set.
#         These configuration files are for LTP features.
#     '''
#     pwd = os.getcwd();
#     os.chdir('../');
#     parseFile()
#     f = open('config.py', 'r');
#     
#     cfnames = ['config_LTP_view1_wpca.py', 'config_LTP_view2_wpca.py', 'config_LTP_pca.py']
#     for fname in cfnames:
#         wf = open(fname, 'w');
#         for l in f:
#             # print l, len(l)
#             l = l.lstrip().rstrip()
#             if len(l) > 0 and  l.lstrip()[0] != '#':
#                 sstr = l.split('=')
#                 if sstr[0] == 'idir' and len(idir) > 0:
#                     print 'Input Directory = ', sstr[1]
#                     l = 'idir=' + idir;
#                 elif sstr[0] == 'odir':
#                     print 'Output Directory = ', sstr[1]
#                     l = 'odir=' + odir;                    
#                 elif sstr[0] == 'ftype':
#                     print 'Feature Type = ', sstr[1]
#                     l = 'ftype=LTP';
#                 elif sstr[0] == 'view':
#                     
#                     print 'View = ', sstr[1]
#                     if fname == 'config_LTP_view1_wpca.py':
#                         l = 'view=view1';
#                     elif fname == 'config_LTP_view2_wpca.py':
#                         l = 'view=view2';
#                     else:
#                         l = 'view=complete';
#                     
#                 elif sstr[0] == 'dist':
#                     print 'Distance = ', sstr[1]            
#                     if fname == 'config_LTP_pca.py':
#                         l = 'dist=cosine'
#                     else:
#                         l = 'dist=chi-square'
#                 elif sstr[0] == 'ttype':
#                     print 'Distance = ', sstr[1]            
#                     if fname == 'config_LTP_pca.py':
#                         l = 'ttype=with-pca'
#                     else:
#                         l = 'ttype=without-pca'                    
#             wf.write(l + '\n')
#         wf.close()
#         f.seek(0)
#         os.chdir(pwd)
#     return cfnames
def main():
    ''' main testing function... '''
    print len(sys.argv)
    if len(sys.argv) != 2 and len(sys.argv) != 3:
        print "Invalid Number of Arguments: Please provide path to output directories [input directory path optional], i.e."
        print "python demo.py path_to_odir path_to_idir[optional]"
        return
    else:
        odir = sys.argv[1]
        print 'Output Directory = ', odir
        idir = ''
        if len(sys.argv) == 3:
            idir = sys.argv[2]
            print 'Input Directory = ', idir
        
            
        if idir != '' and not os.path.exists(idir):
            raise ValueError("Input Image directory %s does not exist" % idir)
#         print "Idir= %s, Odir=%s" % (idir, odir)
        # cfnames = genConfigFiles(idir, odir)
        cfnames = parseConfigFile(idir, odir)
        for i, fname in enumerate(cfnames):
            print "Running Test case Number %d" % (i + 1)
            sys.argv = ['', '--configfile=' + os.path.join(os.path.pardir, fname)];
            print sys.argv
            # now call the running routine...
            r.main()
    
          
     
if __name__ == "__main__":
#     prepare_list('/home/hussain/datasets/LFW/lfwa', '/tmp')
     print 'ARGV:', sys.argv[1:]
     sys.exit(main())       
