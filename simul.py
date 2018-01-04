# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 15:11:38 2018

@author: 3202002
"""

from argparse import ArgumentParser, RawTextHelpFormatter
import os


def launch_simul(n,theta,rep, filename = "out.txt") : 
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
    
    s = "./SimulTrees/SiteFrequencySpectrum "+str(n)+" "+str(theta)+" "+ str(rep) +" > out.txt"
    #print (s)
    os.system(s)
    
    f = open(filename,"r")
    contenu = f.readlines()
    f.close()
    recup = list()
    
    for i in range (len(contenu)-1) :
        if contenu[i][:2] != '/*' and contenu[i][:2] != '\n':
            recup.append(contenu[i])
    #print (recup)
            
    f1 = open(filename,"w")
    tmp = ''.join(recup)
    #print (tmp)
    f1.write(tmp)
    f1.close()
    

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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Pour faire les tests séparés
#=============================
    
#recup_freq()
#launch_simul('20','2','4')



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Pour lancer le code : python3 simul.py -n 20 -t 2 -r 4
#======================================================

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
