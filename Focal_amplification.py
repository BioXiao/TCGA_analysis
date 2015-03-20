#! /usr/bin/python

class Focal_amplification:
  'Class for focal amplifications'
  
  def __init__(self, inputline):
    info = inputline.split("\t")
    self.chrom = info[0]
    self.start = int(info[1])
    self.end = int(info[2])
    self.log2 = float(info[3])
    self.amplified =False
    self.deleted =False
    self.size = self.end - self.start
    if (self.log2 > 0):
      self.amplified =True
    if (self.log2 < 0):
      self.deleted =True
    self.genes=list()
    self.log2_upstream_adj=float(info[10])
    self.log2_downstream_adj=float(info[12])
    self.log2_chrarm_adjust=float(info[8])
    focal_genes = info[5].split(",")
    for gene in focal_genes:
      if gene == "NA":
        continue
      if gene[-9:] == "(partial)":
        self.genes.append(gene[:-9])
      else:
        self.genes.append(gene)
   
##### Print details to stdout ##############
  def printDetails(self):
    print self.chrom+"\t"+str(self.start)+"\t"+str(self.end)+"\t"+str(self.log2)+"\t"+str(self.genes)

##### Check if gene in focal amp ##############
  def checkIfGeneAmp(self,gene):
    if (gene in self.genes) and (self.amplified):
       return True

##### Check if gene in focal del ##############
  def checkIfGeneDel(self,gene):
    if (gene in self.genes) and (self.deleted):
       return True
