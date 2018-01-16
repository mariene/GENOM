# -*- coding: utf-8 -*-
"""
Created on Wed Dec 06 16:52:32 2017

@author: Mariene
"""
from __future__ import print_function

import vcf
import os
import matplotlib.pyplot as plt
import decimal

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
            if a[0] != '0':
                alt+=1
            if a[1] != '0' : 
                alt+=1
            nb_tot +=2
        dico[i] = (alt/nb_tot)
    return dico
 
def comptage_frequence(fichier_vcf):
    """
    
    
    """
    
    def initDico(liste):
        dico = dict()
        for i in range (0,len(liste)+1) : 
            dico[str(i)] = 0.0
        return dico
    
    vcf_reader = vcf.Reader(open(path,'r'))
    d = (pop(path))
    
    
    dico_freq={}
    for i in d.keys(): ## pour chaque population
        dico_freq[i]=[]
            
    for record in vcf_reader :   ## pour chaque position
        liste_alt = list(record.ALT)
        #print (liste_alt)
        for i in d.keys(): ## pour chaque population
            dico = initDico(liste_alt)

            for j in d[i]: ## pour chaque individu
                gen = record.genotype(j)['GT']
                a = gen.split('|')
                
                dico[a[0]] += 1 
                dico[a[1]]+=1
            #print(dico)
            
            #print (len(d[i])*2)==alt1+alt2+alt3
            
            for k in dico.keys():
                if k != '0':
                    if dico[k]/(len(d[i])*2)!=1 and dico[k]/(len(d[i])*2)!=0:
                        dico_freq[i].append(dico[k]/(len(d[i])*2))
            
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


def calc_prob_snp_bis(record,d):
    """Calcule la probabilite de l'allele mutee pour un SNP pour chaque population
    
    
    Parametres
    ----------
    record : Record
    d : dict
    
    Retourne
    --------
    dico : dict, 
        {str : float} -> {pays : proba}
        
    Commentaire
    -----------
    J'ai rajoute une condition pour filtrer les donnees, je recupere juste les 
    genotypes qui n'ont pas un DS 'entier'
    
    """
    dico = {'France':None,'Cameroon':None,'Winters':None,'Raleigh':None,'Autre':None}
    
    for i in d.keys():
        alt = 0.0
        nb_tot = 0.0
        for j in d[i]:
            if record.genotype(j).data.DS == int(record.genotype(j).data.DS):
                gen = (record.genotype(j).data.GT)
                #print (record.genotype(j).data.DS)
                
                #print (record.genotype(j))
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
    """Affiche le spectre de frequence (occ en fct de freq)
    
    Parametres 
    ----------
    prob :list, liste de probabilites
    
    Retourne 
    --------
    dictionary : dict,
        dictionnaire contenant pour chaque proba le nombre d'occurence
    
    """
    
    dictionary=(dict((x,prob.count(x)) for x in set(prob))) 
    plt.figure()
    """
    for l in range(len(list(dictionary.keys()))):
        plt.plot( [list(dictionary.keys())[l], list(dictionary.keys())[l]],[0, list(dictionary.values())[l]], color = 'c'  )
    """
    plt.hist(list(dictionary.values()),bins =100)# len(dictionary.keys()))
    
    plt.title('Spectre de freq')
    plt.xlabel('freq')
    plt.ylabel('occ')
    plt.legend(loc='upper right')
    plt.show()

    
    return dictionary



def replier_hist(prob):
    decimal.getcontext().prec = 1 
    #a = decimal.Decimal(1)/decimal.Decimal(3)
    prob = [round(x,2) for x in prob]
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


def plot_hist_replier(prob):
    """Affiche le spectre replie"""
    dico = replier_hist(prob)    

    liste_freq= list(dico.keys())

    liste_occ=list(dico.values())

    plt.figure()
    for lol in range(len(liste_freq)):

        plt.plot( [liste_freq[lol], liste_freq[lol]],[0, liste_occ[lol]])
    
    plt.title('Spectre de frequence replie')
    plt.ylabel('Occurences')
    plt.xlabel('Frequences')
    plt.legend(loc='upper right')
    plt.savefig('figure1.pdf')
    plt.show()
    
    return dico


def plot_hist_applati(dico1):
    """
    
    Commentaires
    ------------
    c'est pas ça
    """
    
    def applati_hist(dic):
        j = 1
        for i in dic.keys() : 
            dic[i] = dic[i] * (1./j)
            j +=1
            
        return dic
      
    dico = applati_hist(dico1)    
    
    liste_freq=dico.keys()
    liste_occ=[dico[k] for k in liste_freq]

    plt.figure()
    for lol in range(len(liste_freq)):
        plt.plot( [liste_freq[lol], liste_freq[lol]],[0, liste_occ[lol]]    )
    
    plt.title('Spectre de frequence applati')
    plt.ylabel('Occurences')
    plt.xlabel('Frequences')
    plt.legend(loc='upper right')
    plt.savefig('figure1.pdf')
    plt.show()
    
    return dico

def data(prob):
    dictionary=(dict((x,prob.count(x)) for x in set(prob)))
    return dictionary


def premiere_occ(prob):
    dictionary=(dict((x,prob.count(x)) for x in set(prob)))
    return list(dictionary.values())[0]
    

def nb_ech (d,nom):
    return (len(d[nom]))


def applati_hist(dic, pop_size):
    dico_app={}
    for i in dic.keys() :    
        dico_app[i] = dic[i] * (i*pop_size)     
    return dico_app


def conv(m):
    dico = dict()
    for i in np.arange(0.01, 0.51, 0.01).tolist():
        dico[round(i,2)]=0.0
    #print (dico)
    for j in m.keys():
        tmp = m[j]
        dico[round(j,2)]+=tmp
    #print (dico)
    return dico
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TEST%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path = '/media/3202002/PINGOUIN/global.pop.GATK.SNP.hard.filters.V3.phased_all.pop.maf.05.recode.vcf/global.pop.GATK.SNP.hard.filters.V3.phased_all.pop.recode.maf.05.recode.vcf'

#path = os.path.join(os.getcwd(),'Data','test.vcf')
#d = pop(path) 
#plot_histo(prob)
#cpt_freq = comptage_frequence(path)

#print (comptage_frequence(path))
#d2 = calc_all_snp(path)   
     
#d3 = plot_histo(d2['France'])
#d3_replie = replier_hist(d2['France'])

#d4 = plot_hist_replier(d2['France'])
#don = data(d2['France'])



#plot_hist_applati(d4)

"""
d4 = plot_histo(d2['Raleigh'])
d5 = plot_histo(d2['Cameroon'])
d6 = plot_histo(d2['Autre'])
d7 = plot_histo(d2['Winters'])
"""
d3_replie = replier_hist(cpt_freq['Winters'])
d3_bis_replie = conv(d3_replie)
"""
for population in d.keys():
#population='Autre'
    if population !='Cameroon':
        pop_size=len(d[population])*2.
        
        
        freq_1=cpt_freq[population]
        
        dictionary=(dict((x,freq_1.count(x)) for x in set(freq_1)))
        dico_app=applati_hist(dictionary, pop_size)
    else :
        pop_size=len(d[population])
        
        freq_1=cpt_freq[population]
        dictionary=(dict((x,freq_1.count(x)) for x in set(freq_1)))
        dict_tmp = (dict((x,freq_1.count(x)) for x in set(freq_1)))
        for l in dict_tmp.keys():
            if l *pop_size != int(l*pop_size):
                del(dictionary[l])
        
    dico_app= applati_hist(dictionary, pop_size)
    # plot_hist(dictionary, population)
    # plot_hist(dico_app, str(population)+'_app')
    print (population + '\t'+ str(estimation_theta( dictionary, pop_size)))

"""

