#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 13 17:39:16 2018

@author: maureen
"""

import json
import vcf
import os
import matplotlib.pyplot as plt
from numpy import around, histogram,sort
import pickle

def pop(vcf_reader):
    """Renvoie pour chaque pays, la liste d'echantillons
    
    Parametres
    ----------
    path_vcf : str, chemin ou nom du fichier vcf
    
    Retourne 
    --------
    dico : dict, 
        dictionnaire contenant pour chaque pays, la liste d'individus
        {str : list(str)} -> {pays : [echantillons]}
    
    """   
    dico = {'France':[],'Cameroon':[],'Winters':[],'Raleigh':[],'Autre':[]}
    
    for p in vcf_reader.samples : 
        if '_' in p :
            tmp = p.split('_')
        elif '-' in p:
            tmp = p.split('-')
        else :
            tmp = p
               
        if tmp[0:2] == 'FR':
            dico['France'].append(p)            
        elif tmp[0] == 'Winters':
            dico['Winters'].append(p)            
        elif tmp[0] == 'Raleigh':
            dico['Raleigh'].append(p)           
        elif tmp[1] == 'HE':
            dico['Cameroon'].append(p)
        else :
            dico['Autre'].append(p)
    return dico



def comptage_frequence(fichier_vcf):
    
    vcf_reader = vcf.Reader(open(path,'r'))
    d = (pop(vcf_reader))
    
    
    dico_freq={}
    for i in d.keys(): ## pour chaque population
        dico_freq[i]=[]
            
    for record in vcf_reader :   ## pour chaque position
        
        for i in d.keys(): ## pour chaque population
            
            alt1=0. ## Allele 1
            alt2=0. ## Allele 2
            alt3=0. ## alele 0
            nb_tot = 0.0
            for j in d[i]: ## pour chaque individu
                gen = record.genotype(j)['GT']
                a = gen.split('|')
                a= [float(ai) for ai in a]
                if a[0]==1:
                    alt1+=1
                if a[1]==1:
                    alt1+=1
                if a[0]==2:
                    alt2+=1
                if a[1]==2:
                    alt2+=1
                if a[0]==0:
                    alt3+=1
                if a[1]==0:
                    alt3+=1
            #print (len(d[i])*2)==alt1+alt2+alt3
            if alt1/(len(d[i])*2)!=1 and alt1/(len(d[i])*2)!=0:
                dico_freq[i].append(alt1/(len(d[i])*2))
            if alt2/(len(d[i])*2)!=1 and alt2/(len(d[i])*2)!=0:
                dico_freq[i].append(alt2/(len(d[i])*2))

    return dico_freq

def estimation_theta( dico_proba, pop_size):
    """
    moyenne des estimations de theta tout tout i in pop_size
    """
    
    theta=0.0
    i=0.0
    for freq,occ in dico_proba.items():
        theta+=freq*occ*pop_size
#    for g in dico_proba.keys() :
#        theta+= g*pop_size * dico_proba[g]
        i+=1
    theta =theta / i
    return theta
    

def applati_hist(dic, pop_size):
        dico_app={}
        for i in dic.keys() :    
            dico_app[i] = dic[i] * (i*pop_size)     
        return dico_app
    
    
def plot_hist(dico,name):
    plt.figure()
    for k,v in dico.iteritems():
        plt.plot([k,k],[0,v], color = 'c')
#    for l in xrange(len(list(dictionary.keys()))):
#        plt.plot( [list(dico.keys())[l], list(dico.keys())[l]],[0, list(dico.values())[l]], color = 'c'  )
#    
    plt.title('Spectre de freq')
    plt.xlabel('freq')
    plt.ylabel('occ')
    plt.savefig(str(name)+'.png')
    plt.show()
    
    
    
#path = os.path.join(os.getcwd(),'Data','test.vcf')
#path = os.path.join(os.getcwd(),'Data','global.pop.GATK.SNP.hard.filters.V3.phased_all.pop.maf.05.recode.vcf')

vcf_reader = vcf.Reader(open(path,'r'))
d = (pop(vcf_reader))
#frequence=comptage_frequence(path)

with open('result.json', 'r') as fp:
    cpt_freq = json.load(fp)

#with open('dico_frequence','w') as filepickle:
#     pickle.dump(frequence,filepickle)
#    
#with open('dico_frequence','r') as savefilepickle:
#  frequence = pickle.load(savefilepickle)


for population in d.keys():
#population='Autre'
    if population !='Cameroon':
        pop_size=len(d[population])*2.
        
        
        freq_1=frequence[population]
        
        dictionary=(dict((x,freq_1.count(x)) for x in set(freq_1)))
        dico_app=applati_hist(dictionary, pop_size)
    else :
        pop_size=len(d[population])
        
        freq_1=frequence[population]
        dictionary=(dict((x,freq_1.count(x)) for x in set(freq_1)))
        dict_tmp = (dict((x,freq_1.count(x)) for x in set(freq_1)))
        for l in dict_tmp.keys():
            if l *pop_size != int(l*pop_size):
                del(dictionary[l])
        
    dico_app= applati_hist(dictionary, pop_size)
    # plot_hist(dictionary, population)
    # plot_hist(dico_app, str(population)+'_app')
    print (population + '\t'+ str(estimation_theta( dictionary, pop_size)))
