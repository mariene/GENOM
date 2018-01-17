# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 15:11:38 2018

@author: 3202002
"""

from argparse import ArgumentParser, RawTextHelpFormatter
import os
import numpy as np
import glob
import math
from multiprocessing import Process,Queue,Pool,cpu_count
import matplotlib.pyplot as plt

def launch_simul(n,theta,rep) : 
    """Permet de lancer le logiciel SimulTrees
    
    Parametres
    ----------
    n : str, nb echantillons
    theta : str, taux de mutation dans la pop
    rep : str, nb repliques
    filename : str, nom du fichier de sortie (par defaut : out.txt)
    
    Output
    ------
    filename : fichier de sortie contenant que le tableau
    
    """
    
    
    path = os.path.join(os.getcwd(),"Sortie")
    
    if not os.path.exists(path):
        os.mkdir(path)
        
    for i in ['-e','-l']:
        for j in np.arange(-1.0, 1.0, 0.1).tolist() : 
            #print (j)
            filename = os.path.join(path,"out"+"_"+str(i[1])+"_"+str(round(j,1))+".txt")
            
            s = "./SimulTrees/SiteFrequencySpectrum "+str(n)+" "+str(theta)+" "+ str(rep) +" "+i+" "+str(j)+" -F " +"> "+filename
            
            os.system(s)
            
            f = open(filename,"r")
            contenu = f.readlines()
            f.close()
            recup = list()
            #print (contenu)
            
            for k in range (len(contenu)-1) :
                

                if contenu[k][:2] != '/*' and contenu[k][:2] != '\n' and contenu[k][:2] != '  ':
                    recup.append(contenu[k])
        
            
            f1 = open(filename,"w")
            tmp = ''.join(recup)
            #print (tmp)
            f1.write(tmp)
            f1.close()
    return path
    


#Multiproc
def lancement(s):
    os.system(s[1])
    filename = s [0]
    f = open(filename,"r")
    contenu = f.readlines()
    f.close()
    recup = list()

    for k in range (len(contenu)-1) :
        if contenu[k][:2] != '/*' and contenu[k][:2] != '\n' and contenu[k][:2] != '  ':
            recup.append(contenu[k])

    
    f1 = open(filename,"w")
    tmp = ''.join(recup)

    f1.write(tmp)
    f1.close()

            
def launch_simul_bis (n,theta,rep,pro=4):
    
    def genere_liste(n,theta,rep):
        path = os.path.join(os.getcwd(),"Sortie")
    
        if not os.path.exists(path):
            os.mkdir(path)
        dico = dict()
        for i in ['-e','-l']:
            for j in np.arange(-1, 1.0, 0.1).tolist() : 
                filename = os.path.join(path,"out"+"_"+str(i[1])+"_"+str(round(j,1))+".txt")
            
                s = "./SimulTrees/SiteFrequencySpectrum "+str(n)+" "+str(theta)+" "+ str(rep) +" "+i+" "+str(j)+" -F" +"> "+filename
                dico[filename]=(filename,s)
        return dico,path
        
    test,path = genere_liste(n,theta,rep)
    pool = Pool(processes = pro)
    results = pool.map(lancement, test.values())
    
    return path
    
    
 

def all_file_freq(path):
    all_file = glob.glob(os.path.join(path,'*'))
    #print (all_file)
    dico = dict()
    for i in all_file :
        res = (recup_freq(i))
        nom = os.path.split(i)[1][:-4]
        dico[nom]=res
    

    return dico



def recup_freq(fichier = "out.txt") : 
    """Recuperation de l'esperance E[Th_i] 
    
    Parametres
    ----------
    fichier : str, nom du fichier (par defaut : out.txt)
    
    Retourne
    --------
    res : list, [float] liste des E[Th_i]
    
    """
    f = open(fichier,"r")
    contenu = f.readlines()
    f.close()
    res = list()
    for i in range (2, len(contenu)):
        tmp = (contenu[i].rstrip().split('\t'))
        res.append(eval(tmp[2]))
    return res


def calc(data,dico):
    def calcul (sfs_obs, sfs_th):
        #print (sfs_obs, sfs_th)
        return(pow((sfs_obs - sfs_th),2) /  sfs_th)
    
    somme = 0.0
    for i in range (len(data)):
        somme += calcul (data[i], dico[i])
    
    return somme
    
    
def meilleur_scena(data,all_freq):
    n = 1000000000000000000000000
    m = 1000000000000000000000000
    
    save_e = None
    save_l = None
    
    for i in all_freq.keys() : 
        
        if i.split('_')[1] == 'l' :
            tmp_n = calc(data,all_freq[i])
            #print (tmp_n)
            if tmp_n < n :
                save_l = i
                n = tmp_n
            
        elif i.split('_')[1] == 'e' :
            tmp_m = calc(data,all_freq[i])
            #print (tmp_m)
            if tmp_m < m :
                save_e = i
                m = tmp_m
    #print (n,m)
    if (n<m) : 
        print ("lineaire")
    if (n>m) : 
        print("expo")
    
    return save_e,save_l
    
    
    
    

                
                
            
        
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Pour faire les tests séparés
#=============================
    
#recup_freq()
#launch_simul('20','2','4')
#d_repli = {0.35: 36, 0.5: 12, 0.44999999999999996: 6, 0.25: 34, 0.175: 56, 0.325: 19, 0.375: 13, 0.025: 106, 0.425: 21, 0.05: 73, 0.19999999999999996: 4, 0.42500000000000004: 3, 0.275: 26, 0.15000000000000002: 2, 0.125: 22, 0.30000000000000004: 4, 0.22499999999999998: 2, 0.0: 237, 0.225: 55, 0.09999999999999998: 1, 0.475: 21, 0.07499999999999996: 3, 0.17500000000000004: 8, 0.15: 31, 0.1: 24, 0.4: 23, 0.32499999999999996: 4, 0.3: 29, 0.075: 32, 0.2: 51, 0.45: 21}


#mini = d_repli [min(d_repli.keys())]

#fichier = launch_simul(100,mini,1000)
"""
fichier = (launch_simul_bis (100,242273.1884057971,1000))
d_all =all_file_freq(fichier)
#donne = d ['out_l_0.9']
plt.plot(  [i for i in range (len(donne))],donne,'b*' )
"""



#print (meilleur_scena(list(d3_bis_replie.values()),d_all))
#plt.plot(list(map(lambda x: x*100,list(d3_bis_replie.keys() ))), list(d3_bis_replie.values()),'c+' )
#plt.plot(list(d3_replie.keys()), list(d3_replie.values()),'c*' )



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Pour lancer le code : python3 simul.py -n 20 -t 2 -r 4
#======================================================
#fichier = launch_simul(20,100,10) 

#donnes_fr = {0.35: 22, 0.5: 12, 0.75: 3, 0.25: 34, 0.425: 21, 0.225: 55, 0.9: 1, 0.125: 20, 0.625: 8, 0.025: 106, 0.3: 29, 0.575: 3, 0.475: 12, 0.075: 32, 0.775: 2, 0.4: 17, 0.65: 14, 0.925: 3, 0.875: 2, 0.8: 4, 0.85: 2, 0.0: 237, 0.05: 73, 0.525: 9, 0.6: 6, 0.675: 4, 0.7: 4, 0.55: 6, 0.15: 31, 0.1: 24, 0.275: 26, 0.825: 8, 0.175: 56, 0.375: 13, 0.2: 51, 0.45: 21, 0.325: 19}

#plt.hist([donne,list(d_repli.values())])
"""
description = \
    "Description:\n\n" + \
    "Permet de lancer le logiciel SimulTrees et d'avoir le fichier de sortie \n"
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
"""
