"""
Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
All rights reserved.
Released under BSD License
-------------------------
For license terms please see license.lic
"""
"""This is the main driver class for python wrapper"""
import sys
import argparse
import ComputeFeatures as cf

print 'ARGV:', sys.argv[1:]
parser = argparse.ArgumentParser(description=" Feature Computation ...", version=1.0)
ggroup = parser.add_argument_group('General Options')  # general Group
ggroup.add_argument('-f', '--features', dest="ftype",
                    default="LQP", help="Choose among LBP, LTP, LBP+LTP or LQP feature type")
ggroup.add_argument('-l', '--listfile',
                    dest='listfile',
                    help="list file containing list of cropped images to compute features ")
ggroup.add_argument('-d', '--dir',
                    dest='featdir', default="./",
                    help="directory path for storing feature files")                   
ggroup.add_argument('--width',
                    dest='width',
                    help="Width of Cropped Images ", default=150, type=int)
ggroup.add_argument('--height',
                    dest='height',
                    help="Height of Cropped Images ", default=80, type=int)
ggroup.add_argument('-c', '--cellsize',
                    dest='cellsize',
                    help="cellsize of histogram grid ", default=10, type=int)
ggroup.add_argument('-t', '--tolerance',
                    dest='tol',
                    help="tolerance values used for LTP or LQP features ",
                    default=[5], type=int, action="append")

lqpgroup = parser.add_argument_group('LQP Options')  # general Group


lqpgroup.add_argument('--lqptype',
                    dest='lqptype',
                    help="LQP Type 2. Circular, 9. Hor+Ver+Diag+ADiag",
                    default=2, type=int)
lqpgroup.add_argument('-s', '--size',
                    dest='lqpsize',
                    help="LQP size radius of LQP Disk or other types",
                    default=[5], type=int, action="append")
lqpgroup.add_argument('-w', '--words',
                    dest='cbsize',
                    help="Codebook size (number of visual words) used for LQP computation",
                    default=[150], type=int, action="append")
lqpgroup.add_argument('-e', '--encoding',
                    dest='coding',
                    help="0. Binary, 1. Ternary, 4. Split-Ternary",
                    default=4, type=int)                    

options = parser.parse_args()
# ftype, featdir, listfile
feat = cf.Features(options)

