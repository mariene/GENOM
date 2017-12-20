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


def pop(path_vcf):
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
    vcf_reader = vcf.Reader(open(path_vcf,'r'))
    
    
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




def calc_prob_snp(record,d):
    """Calcule la probabilite de l'allele mutee pour un SNP pour chaque population
    
    Parametres
    ----------
    record : Record
    d : dict
    
    Retourne
    --------
    dico : dict, 
        {str : float} -> {pays : proba}
    
    """
    dico = {'France':None,'Cameroon':None,'Winters':None,'Raleigh':None,'Autre':None}
    
    for i in d.keys():
        alt = 0.0
        nb_tot = 0.0
        for j in d[i]:
            gen = (record.genotype(j).data.GT)
            a = gen.split('|')
            #print (a)
            if a[0] != '0':
                alt+=1
            if a[1] != '0' : 
                alt+=1
            nb_tot +=2
        dico[i] = (alt/nb_tot)
    return dico
    
        
def calc_all_snp(path):
    """Calcule la probabilite de l'allele mutee de tous les SNP pour chaque population
    
    Parametres
    ----------
    path : str, chemin du fichier ou nom vcf 
    
    Retourne 
    --------
    dico : dict
        {str : list(float)} -> {pays : [proba]}
    
    
    """
    vcf_reader = vcf.Reader(open(path,'r'))
    d = (pop(path))
    dico = {'France':[],'Cameroon':[],'Winters':[],'Raleigh':[],'Autre':[]}
    
    for record in vcf_reader :
        dic = calc_prob_snp(record,d)
        
        for dt in dic.keys():
            dico[dt].append(dic[dt])
    
    return (dico)



def plot_histo(prob):
    """Affiche un graph occ en fct de freq
    
    Parametres 
    ----------
    prob :list, liste de probabilites
    
    
    """
    
    dictionary=(dict((x,prob.count(x)) for x in set(prob)))
    
    #print (list(dictionary.keys()))
    #print (list(dictionary.values()))
    #plt.bar(list(dictionary.keys()),list(dictionary.values()), color='c')
    
    plt.figure()
    for l in range(len(list(dictionary.keys()))):
        plt.plot( [list(dictionary.keys())[l], list(dictionary.keys())[l]],[0, list(dictionary.values())[l]], color = 'c'  )
    
    
    plt.title('Spectre de freq')
    plt.xlabel('freq')
    plt.ylabel('occ')
    plt.legend(loc='upper right')
    plt.show()

    
    return dictionary





def plot_hist_replier(prob):
    
    def replier_hist(prob):
        dictionary=(dict((x,prob.count(x)) for x in set(prob)))
        liste_freq = list(dictionary.keys())
        liste_occ = list(dictionary.values())
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
      
    dico = replier_hist(prob)    
    
    liste_freq=dico.keys()
    liste_occ=[dico[k] for k in liste_freq]

    plt.figure()
    for lol in range(len(liste_freq)):
        plt.plot( [liste_freq[lol], liste_freq[lol]],[0, liste_occ[lol]]    )
    
    plt.title('Spectre de frequence replie')
    plt.ylabel('Occurences')
    plt.xlabel('Frequences')
    plt.legend(loc='upper right')
    plt.savefig('figure1.pdf')
    plt.show()
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TEST%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path = os.path.join(os.getcwd(),'Data','test_40000.vcf')


d = (pop(path))
#plot_histo(prob)



d2 = calc_all_snp(path)        
d3 = plot_histo(d2['France'])

plot_hist_replier(d2['France'])


"""
d4 = plot_histo(d2['Raleigh'])
d5 = plot_histo(d2['Cameroon'])
d6 = plot_histo(d2['Autre'])
d7 = plot_histo(d2['Winters'])
"""



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
