#! /usr/bin/python

class CNV_segment:
  'Class for focal amplifications'
  
  def __init__(self, inputline):
    info = inputline.split("\t")
    if info[0] == "23":
      self.chrom = "chrX"
    else:
      self.chrom = "chr"+info[0]
    self.start = int(info[1])
    self.end = int(info[2])
    try:
       self.probes = int(info[3])
    except:
       self.probes = "n/a"
    self.log2 = float(info[4][:-1])
    self.amplified =False
    self.deleted =False
    self.size = self.end - self.start
    if (self.log2 > 0.2):
      self.amplified =True
    if (self.log2 < -0.2):
      self.deleted =True
   
##### Print details to stdout ##############
  def printDetails(self):
    print self.chrom+"\t"+str(self.start)+"\t"+str(self.end)+"\t"+str(self.log2)

##### return Log2 values from position ##############
  def getLog2FromPos(self,chrom,start,end):
    if (chrom == self.chrom) and (self.start < start) and (self.end > end):
      return self.log2
    else:
      return None

##### return size from position ##############
  def getSizeFromPosition(self,chrom,start,end):
    if (chrom == self.chrom) and (self.start < start) and (self.end > end):
      return self.size
    else:
      return None


