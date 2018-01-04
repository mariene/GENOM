# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 15:11:38 2018

@author: 3202002
"""

from argparse import ArgumentParser, RawTextHelpFormatter
import os


def launch_simul(n,theta,rep) : 
    
    s = "./SimulTrees/SiteFrequencySpectrum "+n+" "+theta+" "+ rep +" > out.txt"
    print (s)
    os.system(s)

description = \
    "Description:\n\n" + \
    "Permet de lancer le logiciel SimulTrees\n"
    #"via geometric complementarity (sampling step of protein docking).\n"    

parser = ArgumentParser(description=description, formatter_class=RawTextHelpFormatter,
                        add_help=True)

parser.add_argument("-n", "--sample_size", help="total sample size.", 
                    required = True)
parser.add_argument("-t", "--theta",
                    help="Population mutation rate (2pNmu).",required = True)
parser.add_argument("-r", "--rep", 
                    help="number of replicates (independant trees).",required = True)
                    
                    
if __name__ == "__main__":
    args = parser.parse_args()
    #print(args.theta)
    launch_simul(args.sample_size,args.theta,args.rep)