#! /usr/bin/python

class Somatic_mutation:

##### Constructor of a mutation ##############
   def __init__(self, gene, chrom, pos, ref, mut, mut_type):
      self.gene = gene
      self.chrom = chrom
      self.pos = pos
      self.ref = ref
      self.mut = mut
      self.mut_type = mut_type


