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
import matplotlib.pyplot as plt
import numpy as np

path = '/media/3202002/PINGOUIN/global.pop.GATK.SNP.hard.filters.V3.phased_all.pop.maf.05.recode.vcf/global.pop.GATK.SNP.hard.filters.V3.phased_all.pop.recode.maf.05.recode.vcf'
#print(path)
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
    
occ_autre = extract_data.data(cpt_freq['Winters'])
extract_data.plot_hist_bis(cpt_freq['Winters'])
#plot_hist_applati(d4)

winter_repli =  (extract_data.replier_hist(cpt_freq['Winters']))
winter_repli_modif = (extract_data.conv(winter_repli))

theta_pasSure = (winter_repli[min(winter_repli.keys())])

#extract_data.plot_hist_bis(cpt_freq['Winters'],"fig2")
#fichier = (simul.launch_simul_bis (100,242273.1884057971,1000))
#fichier = (simul.launch_simul_bis (100,theta_pasSure,1000))
fichier = os.path.join(os.getcwd(),'Sortie')
d_all = simul.all_file_freq(fichier)
#print (winter_repli_modif)
#print (d_all)
print (theta_pasSure)
expo,lin = (simul.meilleur_scena(list(winter_repli_modif.values()),d_all))
print(expo,lin)
donne_sim = simul.recup_freq( os.path.join(fichier ,lin+".txt"))


dico = winter_repli
#dico=(dict((x,prob.count(x)) for x in set(prob))) 
plt.figure()
for k,v in dico.items():
    #print(k, "->", len(v))
    plt.plot([k,k],[0,v], color = 'c')  
plt.title('Spectre de freq')
plt.xlabel('freq')
plt.ylabel('occ')
#plt.savefig(str(name)+'.png')
#plt.plot( list(dico.keys()),list(dico.values()),'g+' )
#extract_data.plot_hist_bis(donne_sim)
plt.plot( [i for i in np.arange(0.01, 0.51, 0.01).tolist()],donne_sim,'b' )