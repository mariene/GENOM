# -*- coding: utf-8 -*-
"""
Created on Wed Dec 06 16:52:32 2017

@author: Mariene
"""
from __future__ import print_function

import vcf
import os


path = os.path.join(os.getcwd(),'Data','test.vcf')
vcf_reader = vcf.Reader(open(path,'r'))

#record = next(vcf_reader) #iterator
# mais une boucle for ferait l'affaire
for record in vcf_reader :
    #print (record.samples)
    #print ((record.samples)[0].sample)
    print (record)
    if len(record.ALT)>1 :
         for i in vcf_reader.samples :
             print (record.genotype(i))

print (len(vcf_reader.samples))
"""
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




