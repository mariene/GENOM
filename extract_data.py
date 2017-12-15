# -*- coding: utf-8 -*-
"""
Created on Wed Dec 06 16:52:32 2017

@author: Mariene
"""
from __future__ import print_function

import vcf
import os
import matplotlib.pyplot as plt


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def calcul_prob_mut(fichier_vcf):
    vcf_reader = vcf.Reader(open(fichier_vcf,'r'))
    
    prob = list()
    for record in vcf_reader :
    
        alt = 0.0
        nb_tot = 0.0
        #ref = 0.0
        for i in vcf_reader.samples : 
            
            gen = record.genotype(i).data.GT
            a = gen.split('|')
            #print (a)
            if a[0] != '0':
                alt+=1
            if a[1] != '0' : 
                alt+=1
            """
            if a[0] == '0':
                ref+=1
            if a[1] == '0':
                ref+=1
            """
            nb_tot +=2
        #print (ref,alt,nb_tot)
        prob.append(alt/nb_tot)
    return prob
    
#plt.plot((prob))

def plot_histo(prob):
    dictionary=(dict((x,prob.count(x)) for x in set(prob)))

    plt.bar(list(dictionary.values()),list(dictionary.keys()), color='c')
    
    plt.title('Histogramme')
    plt.xlabel('Occurences')
    plt.ylabel('Frequences')
    plt.legend(loc='upper right')
    plt.show()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TEST%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path = os.path.join(os.getcwd(),'Data','lol.vcf')
prob = calcul_prob_mut(path)
plot_histo(prob)




"""
print (len(vcf_reader.samples))    

print (record)
print (record.CHROM)
print (record.POS)

print (record.REF)
print (record.ALT)

print (record.is_snp)
print (record.is_indel)
print (record.is_transition)
print (record.is_deletion)

print (vcf_reader.samples)

"""
#print (record.QUAL)
#print (record.FILTER)
#print (record.INFO)
#print (record.FORMAT)
#print (record.samples)
#print (record.genotype)
#print (record.ID)


#print (record.QUAL)
#print (record.FILTER)
#print (record.INFO)
#print (record.FORMAT)
#print (record.samples)
#print (record.genotype)
#print (record.ID)
