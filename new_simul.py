#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 14 23:04:05 2018

@author: maureen
"""
from new_extract import *

def recup_freq(fichier = "out.txt") : 
    """Recuperation de l'esperance E[X_i] 
    
    Parametres
    ----------
    fichier : str, nom du fichier (par defaut : out.txt)
    
    Retourne
    --------
    col2 : list, [float] liste des E[X_i]
    
    """
    
    
    with open(fichier,"r") as f:
        contenu=f.readlines()

    tmp = 0
    while 'TMRCA' not in contenu[tmp]:
        #print contenu[tmp]
        tmp+=1
    tmp+=2
#    mat = []
#    for l in contenu[tmp:-1]:
#        print l
#        mat.append([float(flo) for flo in l.split()])
    mat = [[float(flo) for flo in l.split()] for l in contenu[tmp:-3]] 
    col2 = [k[1] for k in mat]
    #print col2
    return col2


def compute_error(liste1, liste2):
    """ 
    Pour calculer l'erreur on prend la distance de Manhattan
    """
      
    e=0
    #print len(liste2),len(liste1) 
    if len(liste1)==len(liste2):
        
        for i,j in zip(liste1,liste2):
            #print i,j
            e+=abs(i-j)
            
    return e



p=[]
error=100000000
population='Winters'
n=len(d[population])*2
#dictionary[1/70.]
rep=20

observe=[' ']*n

for freq,occ in dictionary.items():
    observe[int(freq*n)]=occ
observe=observe[1:]

for t in range(41216,412160,50000)   :
    #print t
    s = "./SimulTrees/SiteFrequencySpectrum "+str(n)+" "+str(t)+" "+ str(rep) + 'out.txt'
    # p=2 ?
    # -F pour replier le spectre
    # ./SimulTrees/SiteFrequencySpectrum 70 242268 100 -p 1 > out.txt
    
    os.system(s)
    expected=recup_freq('out.txt')
    #print 'coucou'
    e1= compute_error(expected, observe)
    #print e1
    if e1 <error :
        error=e1
        p.append(t)




recup_freq()

