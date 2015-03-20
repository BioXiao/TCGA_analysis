#! /usr/bin/python

class Gene:
   'Common base class for Gene'
   gene_count=0
  
   def __init__(self, name, refseq, chrom, start, end):
     self.name = name 
     self.refseq = refseq
     self.chrom = chrom
     self.start = int(start)
     self.end = int(end)

   def returnPos(self):
     pos=self.chrom+":"+str(self.start)+"-"+str(self.end)
     return pos
