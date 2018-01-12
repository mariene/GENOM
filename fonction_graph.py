# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 17:40:00 2017

@author: 3600744



Code consu pour lire le fichier hist.txt

Permet de tracer 2 types d'histogramme à partir d'un fichier contenant
les fréquences des SNP et leurs occurences


"""

import os
import numpy as np
import matplotlib.pyplot as plt

def read_list(nom_fichier):
    #freq=cle
    #occ = value

    dico={}
    with open(nom_fichier) as file_:
        for line in file_:
            
            dico[(eval(line.split('\t')[0]))]=(eval(line.split('\t')[1]))
    return dico


def replier_hist(liste_freq, liste_occ):
    dico={}
    
    for x in range(len(liste_freq)) :
        if liste_freq[x] <= 0.5:
            dico[liste_freq[x]]=(liste_occ[x])
        if liste_freq[x] > 0.5:
            if 1-liste_freq[x] in dico.keys():
                dico[1-liste_freq[x]]+=(liste_occ[x])
            else:
                dico[1-liste_freq[x]]=(liste_occ[x])
    return dico




def plot_hist(liste_freq, liste_occ):

    plt.figure(figsize=(10, 8) )
    for lol in range(len(liste_freq)):
        plt.plot( [liste_freq[lol], liste_freq[lol]],[0, liste_occ[lol]]    )
    
    plt.title('Spectre de frequence')
    plt.ylabel('Occurences')
    plt.xlabel('Frequences')
    plt.legend(loc='upper right')
    plt.savefig('figure1.pdf')
    plt.show()
    

def plot_hist_replier(dico):
    liste_freq=dico.keys()
    liste_occ=[dico[k] for k in liste_freq]

    plt.figure(figsize=(10,8) )
    for lol in range(len(liste_freq)):
        plt.plot( [liste_freq[lol], liste_freq[lol]],[0, liste_occ[lol]]    )
    
    plt.title('Spectre de frequence replie')
    plt.ylabel('Occurences')
    plt.xlabel('Frequences')
    plt.legend(loc='upper right')
    plt.savefig('figure1.pdf')
    plt.show()



if __name__ == "__main__":
    
    dico=read_list('hist.txt')

    liste_freq= sorted(dico.keys()) 
    liste_occ=[dico[k] for k in liste_freq]


    plot_hist(liste_freq, liste_occ)
    
    dico_replier=replier_hist(liste_freq, liste_occ)
    plot_hist_replier(dico_replier)
