# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 15:11:38 2018

@author: 3202002
"""

import extract_data
import simul
import json
import vcf
import os, os.path


path = '/media/lahkim/PINGOUIN/global.pop.GATK.SNP.hard.filters.V3.phased_all.pop.maf.05.recode.vcf/global.pop.GATK.SNP.hard.filters.V3.phased_all.pop.recode.maf.05.recode.vcf'
print(path)
if os.path.isfile("result.json"):
  
    with open('result.json', 'r') as fp:
        cpt_freq = json.load(fp)
    print("Loaded result.json")
else:
    print("Computing result.json")
    cpt_freq = extract_data.comptage_frequence(path)
    with open('result.json', 'w') as fp:
        json.dump(cpt_freq, fp)
    print("Saved result.json")
    
occ_autre = extract_data.data(cpt_freq['Autre'])
extract_data.plot_hist_bis(cpt_freq['Autre'])
#plot_hist_applati(d4)
theta_pasSure = (occ_autre[min(occ_autre.keys())])
winter_repli =  (extract_data.replier_hist(cpt_freq['Autre']))
winter_repli_modif = (extract_data.conv(winter_repli))
extract_data.plot_hist_bis(cpt_freq['Autre'],"fig2")
#fichier = (launch_simul_bis (100,242273.1884057971,1000))
#fichier = (simul.launch_simul_bis (100,theta_pasSure,1000))
fichier = os.path.join(os.getcwd(),'Sortie'))
d_all = simul.all_file_freq(fichier)
print (winter_repli_modif)
print (d_all)
print (simul.meilleur_scena(list(winter_repli_modif.values()),d_all))
