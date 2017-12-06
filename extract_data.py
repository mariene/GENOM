# -*- coding: utf-8 -*-
"""
Created on Wed Dec 06 16:52:32 2017

@author: Mariene
"""
from __future__ import print_function

import vcf
import os

path = os.path.join(os.getcwd(),'Data','global.pop.GATK.SNP.hard.filters.V3.phased_all.pop.recode.maf.05.recode.vcf')
vcf_reader = vcf.Reader(open(path,'r'))

record = vcf_reader.next() #iterator 
# mais une boucle for ferait l'affaire 
#for record in vcf_reader :
#    print (record)
print (record)
print (record.CHROM)
print (record.POS)
#print (record.ID)
print (record.REF)
print (record.ALT)
#print (record.QUAL)
#print (record.FILTER)
#print (record.INFO)
#print (record.FORMAT)
#print (record.samples)
#print (record.genotype)

print (record.is_snp)
print (record.is_indel)
print (record.is_transition)
print (record.is_deletion)
print (vcf_reader.samples)


