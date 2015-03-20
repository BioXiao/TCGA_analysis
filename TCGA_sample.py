#! /usr/bin/python

from Focal_amplification import Focal_amplification
from CNV_segment import CNV_segment
from Somatic_mutation import Somatic_mutation
import os.path

class TCGA_sample:
   'Common base class for all TCGA samples'
   sample_count=0
   sample_ids = list()

##### Constructor of a sample ##############
   def __init__(self, identifier):
      if (identifier in TCGA_sample.sample_ids):
         raise ValueError
      self.id = identifier
      if (self.id[13:14] == "1"):
         self.control = True
      else: 
         self.control = False
      TCGA_sample.sample_count += 1
      TCGA_sample.sample_ids.append(identifier)
      self.focal_amplification_data = False
      self.CNV_data = False
      self.rnaseq_data = False
      self.somatic_mutation_data = False
      self.clinical = False

##### Load Focal amplifications ##############
   def loadFocalOutput(self,focal_file):
      if not self.focal_amplification_data:
            self.focal_amplifications = list()
      self.focal_amplification_count = 0
      try:
         FOCAL=open(focal_file,"r")
      except:
         print "Could not open focal amplification file"
         return
      self.focal_amplification_data = True
      header = FOCAL.readline()
      lines = FOCAL.readlines()
      FOCAL.close()
      for line in lines:
        self.focal_amplification_count += 1
        tmp_focal = Focal_amplification(line)
        self.focal_amplifications.append(tmp_focal)

##### Load CNV data ##############
   def loadCNVData(self,CNV_file):
      if not self.CNV_data:
         self.CNV = list() 
      self.CNV_segments = 0
      try:
         CNV_HANDLE=open(CNV_file,"r")
      except:
         print "Could not open CNV segments file"
         return
      self.CNV_data = True
      header = CNV_HANDLE.readline()
      lines = CNV_HANDLE.readlines()
      CNV_HANDLE.close()
      for line in lines:
        self.CNV_segments += 1
        tmp_cnv = CNV_segment(line)
        self.CNV.append(tmp_cnv)

##### Load RNASeq data ##############
   def loadRNASeq(self,rnaseq_file):
      if not self.rnaseq_data:
         self.rnaseq = dict()
      try:
         RNASEQ=open(rnaseq_file,"r")
      except:
         print "Could not open RNASeq file"
         return
      self.rnaseq_data = True
      header = RNASEQ.readline()
      lines = RNASEQ.readlines()
      RNASEQ.close()
      for line in lines:
        info = line.split("\t")
        gene = info[0].split("|")[0]
        self.rnaseq[gene]=float(info[1][:-1])

##### Load SomaticMutation data ##############
   def loadSomaticMutation(self,maf_file):
      if not self.somatic_mutation_data:
         self.somatic_mutations = list()
         self.genes_affected = list()
      try:
         MAF=open(maf_file,"r")
      except:
         print "Could not open Somatic Mutation (MAF) file"
         return
      self.somatic_mutation_data = True
      header = MAF.readline()
      lines = MAF.readlines()
      MAF.close()
      for line in lines:
        info = line.split("\t")
        mutation = Somatic_mutation(info[0],"chr"+info[4], int(info[5]),info[10], info[12],info[8])
        self.somatic_mutations.append(mutation)
      for mut in self.somatic_mutations:
        if (mut.mut_type != "Silent") and mut.gene not in self.genes_affected:
          self.genes_affected.append(mut.gene)

##### Load Clinical data (BRCA only for now) ##############
   def loadClinicalData(self,biotab_file):
      try:
         BIOTAB=open(biotab_file,"r")
      except:
         print "Could not open Clinical BioTab file"
         return
      header = BIOTAB.readline()
      lines = BIOTAB.readlines()
      BIOTAB.close()
      self.clinical = True
      for line in lines:
        info = line.split("\t")
        if (info[0] == "patient.breast_carcinoma_estrogen_receptor_status"):
           self.er_status = info[1][:-1]
        elif (info[0] == "patient.breast_carcinoma_progesterone_receptor_status"):
           self.pr_status = info[1][:-1]
        elif (info[0] == "patient.lab_proc_her2_neu_immunohistochemistry_receptor_status"):
           self.her2_status = info[1][:-1]
        elif (info[0] == "neoplasm.diseasestage"):
           self.stage = info[1][:-1]
        elif (info[0] == "vitalstatus"):
           self.dead = int(info[1][:-1])
        elif (info[0] == "daystodeath"):
           self.daystodeath = info[1][:-1]
        elif (info[0] == "daystolastfollowup"):
           self.daystolastfollowup = info[1][:-1]
        elif (info[0] == "primarysiteofdesease"):
           self.primarysiteofdesease = info[1][:-1]
        elif (info[0] == "histologicaltype"):
           self.histologicaltype = info[1][:-1]
 
         
##### List genes which are focally amplified ##############
   def listFocalAmpGenes(self):
     gene_list = list()
     for focal_amp in self.focal_amplifications:
        if focal_amp.amplified:
          gene_list.extend(focal_amp.genes)
     return gene_list

##### List genes which are focally deleted ##############
   def listFocalDelGenes(self):
     gene_list = list()
     for focal_amp in self.focal_amplifications:
        if focal_amp.deleted:
          gene_list.extend(focal_amp.genes)
     return gene_list

##### Check if gene is focally amplified ##############
   def checkFocalGeneAmp(self,gene):
     found = False
     for focal_amp in self.focal_amplifications:
        found = focal_amp.checkIfGeneAmp(gene)
        if found:
           return found

##### Get size from focal amp Gene ##############
   def getSizeOfFocal(self,gene):
     found = False
     for focal_amp in self.focal_amplifications:
        found = focal_amp.checkIfGeneAmp(gene)
        if found:
          return focal_amp.size
     return None

##### Print details if gene is focally amplified ##############
   def printDetailsFocalGeneAmp(self,gene):
     found = False
     for focal_amp in self.focal_amplifications:
        found = focal_amp.checkIfGeneAmp(gene)
        if found:
           focal_amp.printDetails()


##### Check if gene is focally deleted ##############
   def checkFocalGeneDel(self,gene):
     found = False
     for focal_amp in self.focal_amplifications:
        found = focal_amp.checkIfGeneDel(gene)
        if found:
           return found

##### Print details if gene is focally deleted ##############
   def printDetailsFocalGeneDel(self,gene):
     found = False
     for focal_amp in self.focal_amplifications:
        found = focal_amp.checkIfGeneDel(gene)
        if found:
           focal_amp.printDetails()


##### Print details if gene is focally deleted ##############
   def getRNASeqFromGene(self,gene):
     found = False
     if hasattr(self,"rnaseq") and (gene in self.rnaseq.keys()):
       return self.rnaseq[gene]
     else:
       return "n/a"

##### Get Log2 ratios of position ##############
   def getLog2FromGene(self,chrom, start, end):
     for segment in self.CNV:
        log2 = segment.getLog2FromPos(chrom, start, end)
        if log2 != None:
           return log2

##### Get Sizes of CNV segment ##############
   def getSizeFromPosition(self,chrom, start, end):
     for segment in self.CNV:
        log2 = segment.getSizeFromPosition(chrom, start, end)
        if log2 != None:
           return log2
