#!/usr/bin/env python
import argparse
import scipy.special
import subprocess as sp
import numpy as np
from Bio import pairwise2
from Bio import SeqIO
from Bio.Seq import Seq
import time
import sys
import os


#######################################################FOR DEBUG###############################################################
Process = "V"
Fastq1 = "/home/zengshuai/fastq_files/A_read1.fastq"
Fastq2 = "/home/zengshuai/fastq_files/A_read2.fastq"
ParamFile = "/home/zengshuai/Splindel_python/Splindel.conf"
OutputPrefix = "/home/zengshuai/Splindel_result/CDTrio/A"
SegLen = 22
MaxSize = 1000
## ------------- get predefined setting from [Splindel Folder]/Splindel.conf file---------------------------------------
if os.path.dirname(sys.argv[0]) != "":
	ParamFile = os.path.dirname(sys.argv[0]) + "/" + "Eindel.conf"
else:
	ParamFile = "Eindel.conf"
if not(os.path.isfile(ParamFile)):
    print "\n##########################################################"
    print "Error: Eindel can't find its configure file (Eindel.conf)"
    print "Please copy Eindel.conf to the same folder with Eindel.py"
    print "##########################################################\n"
    sys.exit()
## ------------------- check if parameter file is exists -----------


##----------------------- get user input -------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("Process", help = "Candidate Indel Detection/ Indel Validation/ Both (C/V/A)")
parser.add_argument("Fq1", help = "Input fastq file, the first paired read (Required)")
parser.add_argument("Fq2", help = "Input fastq file, the second paired read (Required)")
parser.add_argument("OutputPrefix", help = "path with Prefix of output files. (Required)")
parser.add_argument("-MaxSize", type = int, help="Maximum deletion size that Splindel can detect(unit: bp, default value: 1000)", default = 1000)
parser.add_argument("-SegLen", type = int, help="prefix and suffix length (unit: bp, default value: 22)", default = 22)
parser.add_argument("-Core", type = int, help="Max numbers of CPU-cores can be used ([1,2,4], default value: 4)", default = 4)
parser.add_argument("-MinSp", type = int, help="minimal number of reads that identify a candidate indel (default value: 2)", default = 2)

parser.add_argument("-numSRAlt", type = int, help="minimal number of supporting reads for validated indel (default value: 2)", default = 2)
parser.add_argument("-numPSAlt", type = int, help="minimal number of reads that identify a candidate indel (default value: 1)", default = 1)
parser.add_argument("-sprAlt", type = int, help="minimal Spread of supoorting reads for validated indels (default value: 10)", default = 10)
parser.add_argument("-MulRatOri", type = int, help="minimum number of reads for validated indels (default value: 0.3)", default = 0.3)
parser.add_argument("-MulRatAlt", type = int, help="maximal multi ratio for validated indels (default value: 0.3)", default = 0.3)
parser.add_argument("-MulRatAlt2", type = int, help="maximum ratio for xxxxx (default value: 0.9)", default = 0.9)

#parser.add_argument("--IncRepRegion", action="store_false", dest="filCanIndbyRepRegion", default=False, help="Include indels in Repetitative Region. (Default: False)")
#parser.add_argument("-MaxSize", type = int, help="Maximum deletion size that Splindel can detect(unit: bp, default value: 1000)", default = 1000)
parser.add_argument("--IncRepReg", help="Include Indels in repetitive Region", action="store_true")
parser.add_argument("--KeepTmp", help= "Do not delete temp file", action="store_true")
#parser.add_argument("--NoConSeq", help="Not build consensus sequences for remapping", default = 1)
args = parser.parse_args()
Process = args.Process
Fastq1 = args.Fq1
Fastq2 = args.Fq2
OutputPrefix = args.OutputPrefix
SegLen = args.SegLen
MaxSize = args.MaxSize
CPUnum = args.Core
MinSpCDetect = args.MinSp
IncRepRegion = args.IncRepReg
IncRepRegion = args.IncRepReg
###########
CnumSRAlt = args.numSRAlt
CnumPSAlt = args.numPSAlt
CsprAlt = args.sprAlt
CMulRatOri = args.MulRatOri
CMulRatAlt = args.MulRatAlt
CMulRatAlt2 = args.MulRatAlt2
################ Process arguments ############################
KeepTemp = 0
if IncRepRegion == True: KeepTemp = 1
filCanIndbyRepRegion = 1
if IncRepRegion == True: filCanIndbyRepRegion = 0
#######################################################FUNCTIONS###############################################################

def DelfileSeries(OutputPrefix2, filesufix):
    import os.path
    ifile = 1
    #print OutputPrefix2 + str(ifile) + filesufix
    while(os.path.isfile(OutputPrefix2 + str(ifile) + filesufix)):
        RmFile(OutputPrefix2 + str(ifile) + filesufix)
        ifile = ifile + 1
        
def RemoveTempFile(OutputPrefix):
    DelfileSeries(OutputPrefix + ".ori.btRef.", ".bt2")  
    DelfileSeries(OutputPrefix + ".alt2.btRef.", ".bt2")
    DelfileSeries(OutputPrefix + ".ori2.btRef.", ".bt2") 
    DelfileSeries(OutputPrefix + ".ori2.btRef.rev.", ".bt2") 
    DelfileSeries(OutputPrefix + ".ori.btRef.rev.", ".bt2") 
    DelfileSeries(OutputPrefix + ".alt.btRef.rev.", ".bt2") 
    DelfileSeries(OutputPrefix + ".alt.btRef.", ".bt2") 
    DelfileSeries(OutputPrefix + ".alt2.btRef.rev.", ".bt2")
    DelfileSeries(OutputPrefix + "_Con_alt.", ".sam") 
    DelfileSeries(OutputPrefix + "_Con_ori.", ".sam")
    DelfileSeries(OutputPrefix + "_alt.", ".bam")
    DelfileSeries(OutputPrefix + "_ori.", ".bam")
    DelfileSeries(OutputPrefix + "_alt.", ".sam")
    DelfileSeries(OutputPrefix + "_ori.", ".sam")
    DelfileSeries(OutputPrefix + ".", ".sam")
    DelfileSeries(OutputPrefix + ".", ".seg.alin")
    DelfileSeries(OutputPrefix + ".", ".seg.fastq")
    DelfileSeries(OutputPrefix + ".", ".unmapped")
    RmFile(OutputPrefix + ".alt.reference.fasta")
    RmFile(OutputPrefix + ".ori.reference.fasta")
    RmFile(OutputPrefix + ".alt.detail.form")
    RmFile(OutputPrefix + ".ori.detail.form")
    RmFile(OutputPrefix + "_ori.sort.bam")
    RmFile(OutputPrefix + "_alt.sort.bam")
    RmFile(OutputPrefix + "_ori.sort.bam.bai")
    RmFile(OutputPrefix + "_alt.sort.bam.bai")
    RmFile(OutputPrefix + "_ori.bam")
    RmFile(OutputPrefix + "_alt.bam")
    RmFile(OutputPrefix + "_ori.sort.pileup")
    RmFile(OutputPrefix + "_alt.sort.pileup")
    RmFile(OutputPrefix + "_ori.cons.fasta")
    RmFile(OutputPrefix + "_alt.cons.fasta")
    RmFile(OutputPrefix + "_PSpair.csv")
    
def CheckFastqValid(fastq1, fastq2):
	##~~~~~~~ find out read length and check if the fastq file are valid
	FastqValid = 1
	p = sp.Popen("head " + Fastq1, stdout=sp.PIPE, shell=True)
	(output1, err) = p.communicate()
	readLength = len(output1.split("\n")[1])
	return readLength, FastqValid

def CheckInputValid(Process, Fastq1, Fastq2, OutputPrefix, SegLen, MaxSize, CPUnum, IncRepRegion, MinSpCDetect):
	##~~~~~~~~ check completeness and validation of user input
	ErrorInfo = ""
	if (Process == "A"): 
		ProcessStep = "All(Candidate Indel Detect + Candidate Indel Validate)"
	elif (Process == "C"): 
		ProcessStep = "Candidate Indel Detect"
	elif (Process == "V"):
		ProcessStep = "Candidate Indel Validate"
	else:
		ErrorInfo = ErrorInfo + "Error: PROCESS session is missing\n"
	##~~~~~~~~ Check if input fastq files are missing
	if (not os.path.isfile(Fastq1)):
		ErrorInfo = ErrorInfo + "Error: fastq1 file can't be found\n"
	if (not os.path.isfile(Fastq2)):
		ErrorInfo = ErrorInfo + "Error: fastq2 file can't be found\n"
	if (not os.path.isdir("/".join(OutputPrefix.split("/")[:-1]))):
		ErrorInfo = ErrorInfo + "Error: output file is not exist\n"
	##~~~~~~~~ check if fasta file of allele sequence are exists, for validation step
	if (Process == "V") and ((not os.path.isfile(OutputPrefix + ".ori.reference.fasta")) or (not os.path.isfile(OutputPrefix + ".alt.reference.fasta"))):
		ErrorInfo = ErrorInfo + "Candidate Indel Detect Step probably not fnished\n"
	if SegLen <15:
		ErrorInfo = ErrorInfo + "Segment length can't be less than 15 bps\n"
	##~~~~~~~~ Output general information to screen
	print "Process Step: " + ProcessStep
	print "Fastq file1: " + Fastq1
	print "Fastq file2: " + Fastq2
	print "Result file prefix: " + OutputPrefix
	print "Split segment length: " + str(SegLen)
	print "Include Indels in repetitive Region: " + str(IncRepRegion)
	if (Process == "A" or Process == "C"):
		print "minimal number of reads that identify a candidate indel: " + str(MinSpCDetect)
	print "Max number of CPU cores used:" + str(CPUnum)
	return ErrorInfo
	
#################################################################################################
#################################################################################################
def LoadParameters(ParamFile):
	## load general setting from configure file, which is located in [Endel folder]/Endel.conf
	samtoolsPath, bowtie2Path, BowRefIndex, Bow2RefIndex, GenomeReffasta = ("", "", "", "", "")
	with open(ParamFile) as f:
	    for line in f:
		if ("=" in line):
			Prefix = line.split("=")[0]
			if Prefix == "bowtie_path": bowtiePath = line.split("=")[1].rstrip()				## path of bowtie program 
			if Prefix == "bowtie2_path": bowtie2Path = line.split("=")[1].rstrip()				## path of bowtie2 program
			if Prefix == "BowRef_Index": BowRefIndex = line.split("=")[1].rstrip()				## reference index for bowtie
			if Prefix == "Bow2Ref_Index": Bow2RefIndex = line.split("=")[1].rstrip()			## reference index of bowtie2 
			if Prefix == "GenomeRef_fasta": GenomeReffasta = line.split("=")[1].rstrip()		## genome reference in fasta format
			if Prefix == "samtools_path": samtoolsPath = line.split("=")[1].rstrip()			## path of samtools program
			if Prefix == "repetitiveRegion_path": repRegionPath = line.split("=")[1].rstrip()	## path to files for repetitive regions
	f.close()
	if bowtiePath=="" or bowtie2Path =="" or BowRefIndex =="" or Bow2RefIndex =="" or GenomeReffasta =="":
		return (-1)					 ## if path are empty, return error
	return(bowtiePath, bowtie2Path, BowRefIndex, Bow2RefIndex, GenomeReffasta, samtoolsPath, repRegionPath)

#################################################################################################
#################################################################################################
def cut_Unmapped_Reads(unmapped, outputFileName, segLen):
	# cut unmapped reads into prefix and suffix
	# the length of prefix and suffix are defined by SegLen
	fileWrite = open(outputFileName, "w")
	fileRead = open(unmapped, 'r')
	Numline =0
	for line in fileRead:			## read fastq file line by line
		line = line.rstrip()
		Numline =1 + Numline
		if Numline % 4 ==1:			## build prefix/suffix ID
			SeqID1 = line + "/1"; SeqID2 = line + "/2"
		if Numline % 4 ==2:			## get prefix and suffix sequences
			Seq1 = line[:segLen]; Seq2 = line[-segLen:]
		if Numline % 4 ==0:			## get prefix and suffix quality sequences
			SeqQ1 = line[:segLen]; SeqQ2 = line[-segLen:]
			fileWrite.write(SeqID1 + "\n" + Seq1 + "\n" + "+" + "\n" + SeqQ1 + "\n")
			fileWrite.write(SeqID2 + "\n" + Seq2 + "\n" + "+" + "\n" + SeqQ2 + "\n")

#################################################################################################
#################################################################################################
def SegAlignment(Seqseg1, Seqseg2, Bow1, Bow2, Max_Muti_loc, CPUnum):
	##~~~~~~~~ map prefix and suffix segment to genome reference
	##~~~~~~~~ one_ex is 1 + (max number of locations allowed for mapping)
	##~~~~~~~~ -v 2, max 2 mismatches are allowed
	if CPUnum >1:							## if user allow more than one CPU core to process mapping
		one_ex = Max_Muti_loc+1;
		p1 = sp.Popen(bowtiePath + "/bowtie -k " + str(one_ex) + " -v 2 " + BowRefIndex + " " + Seqseg1 + " --suppress 5,6,7 > " + Bow1, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
		p2 = sp.Popen(bowtiePath + "/bowtie -k " + str(one_ex) + " -v 2 " + BowRefIndex + " " + Seqseg2 + " --suppress 5,6,7 > " + Bow2, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
		(output1, err1) = p1.communicate()
		(output2, err2) = p2.communicate()
	else:
		one_ex = Max_Muti_loc+1;
		p1 = sp.Popen(bowtiePath + "/bowtie -k " + str(one_ex) + " -v 2 " + BowRefIndex + " " + Seqseg1 + " --suppress 5,6,7 > " + Bow1, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
		(output1, err1) = p1.communicate()
		p2 = sp.Popen(bowtiePath + "/bowtie -k " + str(one_ex) + " -v 2 " + BowRefIndex + " " + Seqseg2 + " --suppress 5,6,7 > " + Bow2, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
		(output2, err2) = p2.communicate()

#################################################################################################
#################################################################################################
def ParseBowTie1Res(Bow1line):
	line = Bow1line.rstrip()
	words = line.split("\t")
	words2 = words[0].split("/")
	header = words2[0]
	ReadsID = words2[1] + words2[2]
	MM = 0
	if len(words)>4:
		words3 = words[4].split(",")
		MM = len(words3)
	Chr = words[2].replace("Chr", "").replace("chr", "").replace("X", "23").replace("Y", "24")
	Location = words[3]
	Chain = -1
	if words[1] =="+": Chain = 1
	return header, ReadsID, Chr, Location, Chain, MM


#################################################################################################
#~~~~~~ read and Parse sam file
def ParseSam(SamLine):
	line = SamLine.rstrip()
	words = line.split("\t")
	words2 = words[0].split("/")
	words3 = words[11].split(":")
	header = words2[0]
	ReadsID = words2[1]
	AS = words3[2]
	Chr = words[2].replace("Chr", "").replace("chr", "").replace("X", "23").replace("Y", "24")
	Location = words[3]
	Chain = 0
	#if words[1] =="+": Chain = 1
	return header, ReadsID, Chr, Location, Chain, AS


#################################################################################################
def form_Presuf_Matrix(SegAlin1):
	# read Segment alignment file, and store it as a dictionary
	AlinInfo = {}
	fileRead = open(SegAlin1, 'r')
	for line in fileRead:
		line = line.rstrip()
		header, ReadsID, Chr, Location, Chain, MM = ParseBowTie1Res(line)
		StIf = ReadsID + "," + Chr + "," + Location + "," + str(Chain) +", "+  str(MM)
		if Chr.isdigit():
			AlinInfo[header] = StIf + "\n" + str(AlinInfo.get(header, ""))
	return AlinInfo

#################################################################################################
##~~~~~ match genome reference sequence and candidate reads sequences, and file best split point
##~~~~~ pairwise2 function from BioPython is used.
def FindBestSplit(genomeSeq, ReadSeq, infoCollect, ReadID):
	#global n
	#print ReadID
	#n = n +1
	IndelLst =""
	gapOpen = -10
	gapExt = -.1
	GenomeRef, Reads, S, B, E = pairwise2.align.globalxs(genomeSeq.upper(), ReadSeq.upper(), gapOpen, gapExt)[0]
	if len(genomeSeq) > len(ReadSeq):
		indelType = "D"
		relatedindelLoc = Reads.find('-')
	else:
		indelType = "I"
		relatedindelLoc = GenomeRef.find('-')
	#if S<30:
	#	print "==============================="
	#	print infoCollect
	#	print genomeSeq
	#	print ReadSeq
	#	print "+++++++++++++++++++++++++++++++"
	#	print GenomeRef
	#	print Reads
	#	print S
	EndLocsIsGap1 = GenomeRef[0] =="-" or GenomeRef[-1] =="-" or Reads[0] =="-" or  Reads[-1] =="-"
	EndLocsIsGap2 = GenomeRef[1] =="-" or GenomeRef[-2] =="-" or Reads[1] =="-" or  Reads[-2] =="-"
	EndLocsIsGap3 = GenomeRef[2] =="-" or GenomeRef[-3] =="-" or Reads[2] =="-" or  Reads[-3] =="-"
	EndLocsIsGap4 = GenomeRef[3] =="-" or GenomeRef[-4] =="-" or Reads[3] =="-" or  Reads[-4] =="-"
	EndLocsIsGap5 = GenomeRef[4] =="-" or GenomeRef[-5] =="-" or Reads[4] =="-" or  Reads[-5] =="-"
	if (EndLocsIsGap1 or EndLocsIsGap2 or EndLocsIsGap3 or EndLocsIsGap4 or EndLocsIsGap5):
		#print "==============================="
		#print infoCollect
		#print genomeSeq
		#print ReadSeq
		#print "+++++++++++++++++++++++++++++++"
		#print GenomeRef
		#print Reads
		#print S
		return -1, "", ""
	else:
		Allele1 = ""; Allele2 = ""
		for i in range(1,(len(Reads)-1)):
			if (Reads[i] == '-') or (GenomeRef[i] == '-'): Allele1 = Allele1 + Reads[i]; Allele2 = Allele2 + GenomeRef[i]
			if (Reads[i] == '-') and (Reads[i-1] != '-'): startIndex = i
			if (Reads[i] == '-') and (Reads[i+1] != '-'): 
				IndelLst = ",".join([IndelLst, str(startIndex),Allele1,Allele2])
				Allele1 = ""; Allele2 = ""
			if (GenomeRef[i] == '-') and (GenomeRef[i-1] != '-'):startIndex = i
			if (GenomeRef[i] == '-') and (GenomeRef[i+1] != '-'): 
				IndelLst = ",".join([IndelLst, str(startIndex),Allele1,Allele2])
				Allele1 = ""; Allele2 = ""
		return relatedindelLoc, indelType, IndelLst

#################################################################################################
##~~~~~ retrieve genome reference, my_data contains chromosome/locations information
def RetreveGenomeReference(GenomeReffasta, my_data):
	#############################
	def retrieveSeq(locPair): return str(seqRecord.seq[locPair[0]:locPair[1]])
	startLoc = np.array([0 for x in range(my_data.shape[0])])
	EndLoc = np.array([0 for x in range(my_data.shape[0])])
	startLoc[my_data[:,4]=='-1'] = my_data[my_data[:,4]=='-1',3].astype(np.int)
	startLoc[my_data[:,4]=='1'] = my_data[my_data[:,4]=='1',8].astype(np.int)
	EndLoc[my_data[:,4]=='-1'] = my_data[my_data[:,4]=='-1',8].astype(np.int) + SegLen
	EndLoc[my_data[:,4]=='1'] = my_data[my_data[:,4]=='1',3].astype(np.int) + SegLen
	#############################
	preSeq = np.array(["" for x in range(my_data.shape[0])], dtype="|S500")
	posSeq = np.array(["" for x in range(my_data.shape[0])], dtype="|S500")
	midSeq = np.array(["" for x in range(my_data.shape[0])], dtype="|S1000")
	ValidES = EndLoc -startLoc >0
	for seqRecord in SeqIO.parse(GenomeReffasta, "fasta"):
		#print "Retrieve Genome Sequence from " + seqRecord.id
		ChrID = seqRecord.id.replace("chr", "").replace("X", "23").replace("Y", "24")
		Locs = zip(startLoc[my_data[:,2] ==ChrID]-readLength, startLoc[my_data[:,2] ==ChrID])
		preSeq[my_data[:,2] ==ChrID] = map(retrieveSeq, Locs)
		Locs = zip(EndLoc[my_data[:,2] ==ChrID], EndLoc[my_data[:,2] ==ChrID]+readLength)
		posSeq[my_data[:,2] ==ChrID] = map(retrieveSeq, Locs)
		Selected = np.logical_and(my_data[:,2] ==ChrID, ValidES)
		Locs = zip(startLoc[Selected], EndLoc[Selected])
		midSeq[Selected] = map(retrieveSeq, Locs)
	return preSeq, posSeq, midSeq, startLoc, EndLoc, ValidES

#################################################################################################
def RetreveFastqSeq(fastqfile1, fastqfile2, readsIDsub):
	uniqueReadsID = np.unique(readsIDsub)
	hashtest = {}
	for x in uniqueReadsID:
		hashtest[x] = ""
	output = np.array(["" for x in range(readsIDsub.shape[0])], dtype="|S250")
	Noline = 3
	print " Loading sequence from 1st unmapped reads..."
	with open(fastqfile1) as f:
		for line in f:
			Noline +=1
			if Noline==4:
				Noline =0
				RecordSignal = 0
				ReadIDs = line
				if hashtest.has_key(ReadIDs[1:-1]):
					thisReadID = ReadIDs[1:-1]
					RecordSignal = 1
			if Noline ==1 and RecordSignal==1: output[np.where(np.array(thisReadID == readsIDsub))]=line[:-1]
	print " Loading sequences from 2nd unmapped reads..."
	Noline = 3
	with open(fastqfile2) as f:
		for line in f:
			Noline +=1
			if Noline==4:
				Noline =0
				RecordSignal = 0
				ReadIDs = line
				if hashtest.has_key(ReadIDs[1:-1]):
					thisReadID = ReadIDs[1:-1]
					RecordSignal = 1
			if Noline ==1 and RecordSignal==1: output[np.where(np.array(thisReadID == readsIDsub))]=line[:-1]
	return output

#########################################################################################################################################################################
###~~~~~ test if two indels can be consolidated
def TwoRefAreIdentical(unIn, indelID, SeqPre, Seq2, SeqPos, RelatedLoc):
	# intialize
	IndelInfoArray = []
	NumIndel =0
	NumTotal = len(unIn)
	# for each candidate indel, build matrix to record if we want to keep this indel
	for x in unIn:
		Inddetail = indelID[x].split("-")
		IndelInfoArray.append([x, Inddetail[0], Inddetail[1], Inddetail[2], Inddetail[3], "P"])			# P = Pass, keep
	IndelInfoArray = np.array(IndelInfoArray)
	IndelInfoArray = IndelInfoArray[IndelInfoArray[:,1].argsort(),:]									# sort by chr
	for i1 in range(IndelInfoArray.shape[0]):
		hitChr = 0
		NumIndel = NumIndel + 1
		sys.stdout.write("\r%f%%" % (float(NumIndel)/float(NumTotal)*100))
		for i2 in range(i1+1, IndelInfoArray.shape[0]):
			Chr1, Loc1, Typ1, Len1 = IndelInfoArray[i1, 1:5]
			Chr2, Loc2, Typ2, Len2 = IndelInfoArray[i2, 1:5]
			unInIndex1 = int(IndelInfoArray[i1, 0]); unInIndex2 = int(IndelInfoArray[i2, 0]);
			if Chr1 == Chr2:
				hitChr = 1
			if (Chr1 != Chr2 and hitChr == 1):# or IndelInfoArray[i2, 5] == "F":
				#print "break in advanced"
				break
			if Chr1 == Chr2 and Len1 == Len2 and abs(int(Loc1)-int(Loc2))<30 and IndelInfoArray[i2, 5] == "P":
				gapOpen = -20
				gapExt = -.1
				RefSeq1 = "".join((SeqPre[RelatedLoc >-1][unInIndex1], Seq2[RelatedLoc >-1][unInIndex1], SeqPos[RelatedLoc >-1][unInIndex1]))	# build ref1
				RefSeq2 = "".join((SeqPre[RelatedLoc >-1][unInIndex2], Seq2[RelatedLoc >-1][unInIndex2], SeqPos[RelatedLoc >-1][unInIndex2]))	# build ref2
				#break
				AlignRef1, AlignRef2, S, B, E = pairwise2.align.globalxs(RefSeq1.upper(), RefSeq2.upper(), gapOpen, gapExt)[0]
				AlignRef1=AlignRef1[10:-10]; AlignRef2=AlignRef2[10:-10]
				SharedStartPos = max(min(AlignRef1.find("G"), AlignRef1.find("T"), AlignRef1.find("C"), AlignRef1.find("A")), min(AlignRef2.find("G"), AlignRef2.find("T"), AlignRef2.find("C"), AlignRef2.find("A")))
				AlignRev1 = AlignRef1[::-1]; AlignRev2 = AlignRef2[::-1]
				SharedEndPos = max(min(AlignRev1.find("G"), AlignRev1.find("T"), AlignRev1.find("C"), AlignRev1.find("A")), min(AlignRev2.find("G"), AlignRev2.find("T"), AlignRev2.find("C"), AlignRev2.find("A")))
				if not ("-" in AlignRef1[SharedStartPos:(len(AlignRef1)-SharedEndPos)] or "-" in AlignRef2[SharedStartPos:(len(AlignRef1)-SharedEndPos)]):
					IndelInfoArray[i2, 5] = "F"
	return IndelInfoArray[IndelInfoArray[:,5]=="P",0]

#########################################################################################################################################################################


def formatAltRefSeq(OutputPrefix, fastq1, fastq2, SPRthreshold = 1):
	#######
	global my_data
	def GetChoosedList(indelIDtest):
		IndexSet = np.where(np.array(indelID) == indelIDtest)#; np.random.shuffle(IndexSet)
		return IndexSet[0][0]
	#######
	readsIDsub = np.core.defchararray.add(my_data[:,0], map(lambda x: "/" + x[0], my_data[:,1]))
	ReadsSeqList = RetreveFastqSeq(OutputPrefix + ".1.unmapped", OutputPrefix + ".2.unmapped", readsIDsub)
	MinValid = abs(my_data[:,3].astype(int)- my_data[:,8].astype(int)) >SegLen
	my_data = my_data[MinValid,:]
	preSeq, posSeq, midSeq, StartLoc, EndLoc, ValidES = RetreveGenomeReference(GenomeReffasta, my_data)		##ValidES end-start>0
	##my_data = my_data[ValidES]
	Seq1 = midSeq[ValidES]; Seq2 = ReadsSeqList[MinValid][ValidES]
	SeqPre = preSeq[ValidES]; SeqPos = posSeq[ValidES]
	Seq2[my_data[ValidES,4] == '-1']=map(lambda x: str(Seq(x).reverse_complement()),Seq2[my_data[ValidES,4] == '-1'])
	#n = 0
	TempResults = map(FindBestSplit, Seq1, Seq2, my_data[ValidES], range(len(Seq1)))
	RelatedLoc = np.array(TempResults)[:,0].astype(np.int)
	TypeList = np.array(TempResults)[:,1][RelatedLoc >-1]
	ChrList = my_data[ValidES,2][RelatedLoc >-1]
	Loclist = (StartLoc[ValidES] + RelatedLoc)[RelatedLoc >-1]
	LenList = (abs(abs(my_data[ValidES,3].astype(np.int) - my_data[ValidES,8].astype(np.int))- (readLength - SegLen)))[RelatedLoc >-1]
	indelID = np.array(map(lambda a,b,c,d:a+"-"+b+"-"+c+"-"+d, ChrList, Loclist.astype(np.str), TypeList, LenList.astype(np.str)))
	DetailedIndel = np.array(TempResults)[RelatedLoc >-1,2]
	SRread = NumbersSR4Candidate(indelID)
	unIn = map(GetChoosedList, np.unique(indelID[SRread>SPRthreshold]))
	##################consolidate###################
	unInCon = TwoRefAreIdentical(unIn, indelID, SeqPre, Seq2, SeqPos, RelatedLoc)
	unIn = np.array(unInCon, dtype ='int')
	##################consolidate###################
	#StartLoc[ValidES][unIn]
	OriginalSeq = map(lambda a,b,c : a.upper()+b.upper()+c.upper(), SeqPre[RelatedLoc >-1][unIn], Seq1[RelatedLoc >-1][unIn], SeqPos[RelatedLoc >-1][unIn])
	AlternativeSeq = map(lambda a,b,c : a.upper()+b.upper()+c.upper(), SeqPre[RelatedLoc >-1][unIn], Seq2[RelatedLoc >-1][unIn], SeqPos[RelatedLoc >-1][unIn])
	f=open(OutputPrefix + ".alt.reference.fasta",'w')
	test = map(lambda a,b,c : f.write(">"+a+"-alt_"+str(b-readLength)+"\n"+c+"\n"), indelID[unIn], StartLoc[ValidES][RelatedLoc >-1][unIn], AlternativeSeq)
	f=open(OutputPrefix + ".ori.reference.fasta",'w')
	test = map(lambda a,b,c : f.write(">"+a+"-ori_"+str(b-readLength)+"\n"+c+"\n"), indelID[unIn], StartLoc[ValidES][RelatedLoc >-1][unIn], OriginalSeq)
	f=open(OutputPrefix + ".alt.detail.form",'w')
	test = map(lambda a,b,c: f.write(a+":" + str(c) + ","+b+"\n"), indelID[unIn], DetailedIndel[unIn], StartLoc[ValidES][RelatedLoc >-1][unIn])
	f.close()


#################################################################################################
##~~~~~ calculate number of supporting reads for indel allele
def NumbersSR4Candidate(indelID):
	NumbersSR = {}
	for indelid in indelID:
		if NumbersSR.has_key(indelid):
			NumbersSR[indelid] = NumbersSR[indelid] + 1
		else:
			NumbersSR[indelid] = 1
	NumSR = np.zeros(len(indelID))
	idata = -1
	for indelid in indelID:
		idata = idata + 1
		NumSR[idata] = NumbersSR[indelid]
	return NumSR

#################################################################################################
##~~~~~ align all fastq files to newly-built references
##~~~~~ raw alternative reference and original reference, when maptime ==1
##~~~~~ consensus alternative reference and original reference, when maptime ==2
def WholeAlignOriAlt(OutputPrefix, fastq1, fastq2, maptime):
    #### asign different file name according to maptime
	if maptime ==1:
		alt = OutputPrefix + ".alt.reference.fasta"
		dbAlt = OutputPrefix + ".alt.btRef"
		ori = OutputPrefix + ".ori.reference.fasta"
		dbOri = OutputPrefix + ".ori.btRef"
		OutputPrefix1 = OutputPrefix
	elif maptime ==2:
		alt = OutputPrefix + "_alt.cons.fasta"
		dbAlt = OutputPrefix + ".alt2.btRef"
		ori = OutputPrefix + "_ori.cons.fasta"
		dbOri = OutputPrefix + ".ori2.btRef"
		OutputPrefix1 = OutputPrefix + "_Con"
	p1 = sp.Popen(bowtie2Path +"/bowtie2-build --quiet " + alt + " " + dbAlt, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
	p2 = sp.Popen(bowtie2Path +"/bowtie2-build --quiet " + ori + " " + dbOri, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
	(output1, err1) = p1.communicate(); (output2, err2) = p2.communicate()
	#cutoff = READ_QUALITY_cutoff1/read_length			#####cut off 
	p1 = sp.Popen(bowtie2Path +"/bowtie2 -k 10 --gbar " + str(readLength) + " --score-min L,-3,-0.12 --no-unal --quiet -x "+ dbAlt +" -U " + fastq1 + " -S "+OutputPrefix1+"_alt.1.sam", stdout=sp.PIPE, shell=True)
	p2 = sp.Popen(bowtie2Path +"/bowtie2 -k 10 --gbar " + str(readLength) + " --score-min L,-3,-0.12 --no-unal --quiet -x "+ dbOri +" -U " + fastq1 + " -S "+OutputPrefix1+"_ori.1.sam", stdout=sp.PIPE, shell=True)
	#cutoff = READ_QUALITY_cutoff2/read_length;
	p3 = sp.Popen(bowtie2Path +"/bowtie2 -k 10 --gbar " + str(readLength) + " --score-min L,-3,-0.12 --no-unal --quiet -x "+ dbAlt +" -U " + fastq2 + " -S "+OutputPrefix1+"_alt.2.sam", stdout=sp.PIPE, shell=True)
	p4 = sp.Popen(bowtie2Path +"/bowtie2 -k 10 --gbar " + str(readLength) + " --score-min L,-3,-0.12 --no-unal --quiet -x "+ dbOri +" -U " + fastq2 + " -S "+OutputPrefix1+"_ori.2.sam", stdout=sp.PIPE, shell=True)
	#print "################"
	#print output1
	#print err1
	#print "################"
	(output1, err1) = p1.communicate(); (output2, err2) = p2.communicate()
	(output3, err3) = p3.communicate(); (output4, err4) = p4.communicate()
	return OutputPrefix1+"_alt.1.sam", OutputPrefix1+"_ori.1.sam", OutputPrefix1+"_alt.2.sam", OutputPrefix1+"_ori.2.sam"

##################################################################################################
##################################################################################################
def buildPileup(samtoolsPath, OutputPrefix):
	# Pileup file is used to build consensus sequences. 
	# pileup file is built using samtools
	print " building Pileup file..."
	p1 = sp.Popen(samtoolsPath +"/samtools view -Sb " + OutputPrefix+ "_ori.1.sam > "+ OutputPrefix + "_ori.1.bam", stdout=sp.PIPE, shell=True)
	p2 = sp.Popen(samtoolsPath +"/samtools view -Sb " + OutputPrefix+ "_ori.2.sam > "+ OutputPrefix + "_ori.2.bam", stdout=sp.PIPE, shell=True)
	p3 = sp.Popen(samtoolsPath +"/samtools view -Sb " + OutputPrefix+ "_alt.1.sam > "+ OutputPrefix + "_alt.1.bam", stdout=sp.PIPE, shell=True)
	p4 = sp.Popen(samtoolsPath +"/samtools view -Sb " + OutputPrefix+ "_alt.2.sam > "+ OutputPrefix + "_alt.2.bam", stdout=sp.PIPE, shell=True)
	(output1, err1) = p1.communicate(); (output2, err2) = p2.communicate()
	(output3, err3) = p3.communicate(); (output4, err4) = p4.communicate()
	p1 = sp.Popen(samtoolsPath +"/samtools merge -f "+OutputPrefix+"_ori.bam " + OutputPrefix+ "_ori.1.bam "+ OutputPrefix + "_ori.2.bam", stdout=sp.PIPE, shell=True)
	p2 = sp.Popen(samtoolsPath +"/samtools merge -f "+OutputPrefix+"_alt.bam " + OutputPrefix+ "_alt.1.bam "+ OutputPrefix + "_alt.2.bam", stdout=sp.PIPE, shell=True)
	(output1, err1) = p1.communicate(); (output2, err2) = p2.communicate()
	p1 = sp.Popen(samtoolsPath +"/samtools sort "+OutputPrefix+"_ori.bam " + OutputPrefix+ "_ori.sort", stdout=sp.PIPE, shell=True)
	p2 = sp.Popen(samtoolsPath +"/samtools sort "+OutputPrefix+"_alt.bam " + OutputPrefix+ "_alt.sort", stdout=sp.PIPE, shell=True)
	(output1, err1) = p1.communicate(); (output2, err2) = p2.communicate()
	p1 = sp.Popen(samtoolsPath +"/samtools index "+OutputPrefix+ "_ori.sort.bam", stdout=sp.PIPE, shell=True)
	p2 = sp.Popen(samtoolsPath +"/samtools index "+OutputPrefix+ "_alt.sort.bam", stdout=sp.PIPE, shell=True)
	(output1, err1) = p1.communicate(); (output2, err2) = p2.communicate()
	#p1 = sp.Popen(samtoolsPath +"/samtools mpileup " +OutputPrefix+ "_ori.sort.bam >" +OutputPrefix+ "_ori.sort.pileup -f " +OutputPrefix+ ".ori.reference.fasta", stdout=sp.PIPE, shell=True)
	#p2 = sp.Popen(samtoolsPath +"/samtools mpileup " +OutputPrefix+ "_alt.sort.bam >" +OutputPrefix+ "_alt.sort.pileup -f " +OutputPrefix+ ".alt.reference.fasta", stdout=sp.PIPE, shell=True)
	p1 = sp.Popen(samtoolsPath +"/samtools mpileup " +OutputPrefix+ "_ori.sort.bam >" +OutputPrefix+ "_ori.sort.pileup", stdout=sp.PIPE, shell=True)
	p2 = sp.Popen(samtoolsPath +"/samtools mpileup " +OutputPrefix+ "_alt.sort.bam >" +OutputPrefix+ "_alt.sort.pileup", stdout=sp.PIPE, shell=True)
	(output1, err1) = p1.communicate(); (output2, err2) = p2.communicate()

##################################################################################################
###~~~~~build consensus reference sequences, reading from raw sequences and pileup file
def buildConsensusSeq(OutputPrefix, AltOri):
	fh = open(OutputPrefix+ "_"+AltOri+".sort.pileup")
	fo = open(OutputPrefix+ "_"+AltOri+".cons.fasta", "w")
	RefSeqRecord= {}
	# read original/alternative reference sequences  
	for seqRecord in SeqIO.parse(OutputPrefix+ "."+AltOri+".reference.fasta", "fasta"):
		RefSeqRecord[seqRecord.id] = str(seqRecord.seq)
	lastreadID = ''
	while True:
		# read pileup file line by line
		lineContent = fh.readline()
		if not lineContent: break
		readID = lineContent.split("\t")[0]
		LocOnRef = lineContent.split("\t")[1]
		OriNC = lineContent.split("\t")[2].upper()
		pileUp = lineContent.split("\t")[4].upper().replace(",",".").replace(".", OriNC)
		#print pileUp
		if readID != lastreadID:
			if lastreadID != '': fo.write(">" + lastreadID + "\n" + "".join(SeqChaList) + "\n")
			lastreadID = readID
			SeqChaList = [RefSeqRecord[lastreadID][i] for i in range(len(RefSeqRecord[lastreadID]))]
		indexLar = np.array([pileUp.count('A'), pileUp.count('T'), pileUp.count('C'), pileUp.count('G')])
		maxIndex = np.where(indexLar == max(indexLar))[0][0]
		if maxIndex ==0: SeqChaList[int(LocOnRef)-1]='A'
		if maxIndex ==1: SeqChaList[int(LocOnRef)-1]='T'
		if maxIndex ==2: SeqChaList[int(LocOnRef)-1]='C'
		if maxIndex ==3: SeqChaList[int(LocOnRef)-1]='G'
	fo.write(">" + lastreadID + "\n")
	fo.write("".join(SeqChaList) + "\n")
	fo.close()

####################################################################################################
def LoadSam(samFile):
	# read sam file, and output a list contains [readid, chr, loc, mapping Quality]
	fh = open(samFile)
	ReadsID = []; RefID = []; LocOnRef = []; MappingQ = []
	while True:
		lineContent = fh.readline()
		if not lineContent: break
		if not lineContent[0] =='@':
			ReadsID.append(lineContent.split("\t")[0])
			MappingQ.append(lineContent.split("\t")[11].split(":")[2])
			RefID.append(lineContent.split("\t")[2])
			LocOnRef.append(lineContent.split("\t")[3])
	return(zip(ReadsID, RefID, LocOnRef, MappingQ))

###################################################################################################
def ParseFirstRunSam(SamFile1, mapReadID):
	# read sam file, and output a list contains [readid, chr, loc, mapping Quality]
	fh = open(SamFile1)
	ReadsID = []; RefID = []; LocOnRef = []; MappingQ = []
	while True:
		lineContent = fh.readline()
		if not lineContent: break
		if not lineContent[0] =='@':
			if lineContent.split("/")[0] in mapReadID:
				ReadsID.append(lineContent.split("\t")[0])
				MappingQ.append(lineContent.split("\t")[11].split(":")[2])
				RefID.append(lineContent.split("\t")[2])
				LocOnRef.append(lineContent.split("\t")[3])
	return(zip(ReadsID, RefID, LocOnRef, MappingQ))

#######################################################################################################################################
####~~~~~ get alternative/original/Genome alignment from previous aligments
def MakeAlignArray(RefAlignFile, GenAlignFile):
	RefAlign = np.array(LoadSam(RefAlignFile))
	################################
	AltRefID = RefAlign[:,1]
	ALTtag = np.array(map(lambda x: x.find("alt")>0, AltRefID))
	AltAlignment =RefAlign[ALTtag]
	OriAlignment =RefAlign[np.logical_not(ALTtag)]
	################################
	ReadsIDInc = np.concatenate((np.array(map(lambda x: x.split("/")[0], AltAlignment[:,0])), np.array(map(lambda x: x.split("/")[0], OriAlignment[:,0]))))
	ReadsIDInc = np.unique(ReadsIDInc)
	mapReadID =  dict(zip(ReadsIDInc, np.array([1 for x in range(ReadsIDInc.shape[0])])))
	RefAlignment =np.array(ParseFirstRunSam(GenAlignFile, mapReadID))
	return AltAlignment, OriAlignment, RefAlignment


###############################################################################################
### ~~~~~ 
def dbbinom(x, N, a, b):
	n = N
	k = x
	temp = scipy.special.beta(x+a, n-x+b)/scipy.special.beta(a, b)
	if 0 <= k <= n:
		ntok = 1; ktok = 1
		for t in xrange(1, min(k, n - k) + 1):
			ntok *= n; ktok *= t; n -= 1
		return  temp *(ntok // ktok) 
	else:
		return 0

###################################################################################################
def dbbinomAr(x, n, a, b):
	output = [0 for y in range(len(x))]
	for i in range(len(x)):
		output[i] = dbbinom(x[i], n[i], a, b)
	return np.array(output)

####################################################################################################
def ParameterEstimate_bBinomial(x, n):
	######
	t_a=0.0001; t_b=0.0001; t_a2=0.0001; t_b2=0.0001; t_a3=0.0001; t_b3=0.0001;
	deri_a=0.01; deri_b=0.01; deri_a2=0.01; deri_b2=0.01; deri_a3=0.01; deri_b3=0.01;
	w1=0.5; w2=0.3; w3=1-w1-w2; w0 = 1
	a=4; b=4; a_2=10; b_2=40; a_3=95; b_3=5;
	def Llikelyhood(x, n, a, b, a_2, b_2, a_3, b_3):
		return sum(np.log(w1* dbbinomAr(x,n,a,b)+w2* dbbinomAr(x,n,a_2,b_2)+(1-w1-w2)* dbbinomAr(x,n,a_3,b_3)))
	######
	undateIteration = 0
	while (abs(w1 -w0)>0.0001) and undateIteration <50:
		undateIteration +=1 
		w0=w1; No_Mid=0; No_left=0; No_right=0;
		TotalProb = w1*dbbinomAr(x,n,a,b) + w2*dbbinomAr(x,n,a_2,b_2) + (1-w1-w2)*dbbinomAr(x,n,a_3,b_3)
		P_mid=w1*dbbinomAr(x,n,a,b)/TotalProb
		P_left=w2*dbbinomAr(x,n,a_2,b_2)/TotalProb
		P_right=1-P_mid-P_left
		genotype = np.array(map(lambda x: np.where(x ==max(x))[0][0],np.column_stack((P_left,P_mid,P_right))))
		w1_0 = str(w1); w2_0 = str(w2); w3_0 = str(w3)
		w1 = float(sum(genotype ==0)) / float(genotype.shape[0])
		w2 = float(sum(genotype ==1)) / float(genotype.shape[0])
		w3 = float(sum(genotype ==2)) / float(genotype.shape[0])
		start=0;
		temp0=Llikelyhood(x, n, a, b, a_2, b_2, a_3, b_3)
		temp00 = 0; no_iteration_M = 0
		print " update weight from " + w1_0[:5] + ", " + w2_0[:5] + ", " + w3_0[:5] + " To " + str(w1)[:5] + ", " + str(w2)[:5] + ", " + str(w3)[:5]
		while((abs(temp0-temp00)>float(len(n))/400)) and no_iteration_M <100:
			no_iteration_M = no_iteration_M + 1
			temp00=temp0;
			start=1;
			deri_a=(Llikelyhood(x, n, a*1.01, b, a_2, b_2, a_3, b_3)-Llikelyhood(x, n, a, b, a_2, b_2, a_3, b_3))/(a*1.01)
			a=a+t_a*deri_a
			if (a<=0):a=0.000001
			deri_b=(Llikelyhood(x, n, a, b*1.01, a_2, b_2, a_3, b_3)-Llikelyhood(x, n, a, b, a_2, b_2, a_3, b_3))/(b*0.01)
			b=b+t_b*deri_b;
			if (b<=0):b=0.00001
			deri_a2=(Llikelyhood(x, n, a, b, a_2*1.01, b_2, a_3, b_3)-Llikelyhood(x, n, a, b, a_2, b_2, a_3, b_3))/(a_2*0.01);
			a_2=a_2+t_a2*deri_a2; 
			if (a_2<=0):a_2=0.00001
			deri_b2=(Llikelyhood(x, n, a, b, a_2, b_2*1.01, a_3, b_3)-Llikelyhood(x, n, a, b, a_2, b_2, a_3, b_3))/(b_2*0.01);
			b_2=b_2+t_b2*deri_b2;
			if (b_2<=0):b_2=0.00001
			deri_a3=(Llikelyhood(x, n, a, b, a_2, b_2, a_3*1.01, b_3)-Llikelyhood(x, n, a, b, a_2, b_2, a_3, b_3))/(a_3*0.01);
			a_3=a_3+t_a3*deri_a3;
			if (a_3<=0):a_3=0.00001
			deri_b3=(Llikelyhood(x, n, a, b, a_2, b_2, a_3, b_3*1.01)-Llikelyhood(x, n, a, b, a_2, b_2, a_3, b_3))/(b_3*0.01);
			b_3=b_3+t_b3*deri_b3;
			if (b_3<=0):b_3=0.00001
			w3=1-w1-w2;
			temp0 = Llikelyhood(x, n, a, b, a_2, b_2, a_3, b_3)
			#print temp0, abs(temp0-temp00), float(len(n))/400
	return a, b, a_2, b_2, a_3, b_3, w1, w2, w3

####################################################################################################
def prediction_genotype(x, n, a, b, a_2, b_2, a_3, b_3, w1, w2, w3):
	# predict genotypes, given:
	# 	x: vector of total number of supporting reads for alternative allele
	#	n: bector of total number of supporting reads for both alternative allele and original allele
	#	a, a_2, a_3: alpha for first, second and third beta-binomial component 
	#	b, b_2, b_3: beta for first, second and third beta-binomial component
	#	w1, w2, w3 weights for first, second and third components
	TotalProb = w1*dbbinomAr(x,n,a,b) + w2*dbbinomAr(x,n,a_2,b_2) + (1-w1-w2)*dbbinomAr(x,n,a_3,b_3)
	# calculate probability for each genotype
	P_mid=w1*dbbinomAr(x,n,a,b)/TotalProb		# first component(homozygous reference)
	P_left=w2*dbbinomAr(x,n,a_2,b_2)/TotalProb	# second component (hetezygous indel)
	P_right=1-P_mid-P_left						# third component (homozygous indel)
	# predict genotype, which is the one with largest probability
	genotype = np.array(map(lambda x: np.where(x ==max(x))[0][0],np.column_stack((P_left,P_mid,P_right))))	
	return genotype, P_left, P_mid, P_right

###################################################################################################

def judgeSegPairValid(OutputPrefix, AlinInfo, fileID):
	#
	print " file ID: " + str(fileID) + ", processed"
	ProcessID = 0
	if fileID == 1:
		f=open(OutputPrefix + "_PSpair.csv",'w')
	else:
		f=open(OutputPrefix + "_PSpair.csv",'a')
	TotalKeyNum = len(AlinInfo)
	#f.write("########################" + str(fileID))		#####for debug only
	finishedTag = 0
	for key, value in AlinInfo.iteritems():
		ProcessID +=1
		if (ProcessID*100/TotalKeyNum > finishedTag):
			#sys.stdout.write("\r%f%%" % (ProcessID*100/TotalKeyNum))
			finishedTag = ProcessID*100/TotalKeyNum
		allCombine = []
		readsalines = np.array(value.split("\n")[:-1])
		repValid = np.array([False for x in range(len(readsalines))])
		ReadsID = np.array(map(lambda x: x[:2], readsalines))
		if sum(ReadsID=="21") * sum(ReadsID=="22")>0 and sum(ReadsID=="21") <21 and sum(ReadsID=="22")<21:	##filter those segment with more than
			repValid[np.logical_or(ReadsID=="21", ReadsID=="22")] = True				##21 locs on genome
		
		if sum(ReadsID=="11") * sum(ReadsID=="12")>0 and sum(ReadsID=="11") <21 and sum(ReadsID=="12")<21:
			repValid[np.logical_or(ReadsID=="11", ReadsID=="12")] = True
		
		readsalines = readsalines[repValid]
		for i in range(0,len(readsalines)-1):
				for r in range(i+1,len(readsalines)):
					if readsalines[i][:1] == readsalines[r][:1] and readsalines[i][1:2] != readsalines[r][1:2]:	##readsID and SegID valid
						allCombine.append(readsalines[i]+ "," + readsalines[r])
		if len(readsalines)>0:
			allCombineAr = np.array([x.split(",") for x in allCombine])
			Chr1 = np.array(allCombineAr[:,1]); Chr2 = np.array(allCombineAr[:,6])
			Chain1 = np.array(allCombineAr[:,3]); Chain2 = np.array(allCombineAr[:,8])
			loc1 = np.array(allCombineAr[:,2],dtype='int32'); loc2 = np.array(allCombineAr[:,7],dtype='int32')
			Criterion = (Chain1 == Chain2) & (abs(loc1 - loc2) <1200) & (Chr1 == Chr2)# & (abs(loc2 - loc1)-(readLength-SegLen))!=0
			if sum(Criterion) >0 and min(abs(abs(loc2 - loc1)[Criterion]-(readLength-SegLen))) != 0:
				CriterionMin = abs(abs(loc2 - loc1)[Criterion]-(readLength-SegLen))==min(abs(abs(loc2 - loc1)[Criterion]-(readLength-SegLen)))
				ReadID = np.repeat(key, sum(CriterionMin), axis=0)
				tempRes = np.c_[ReadID, allCombineAr[Criterion,][CriterionMin]]
				for Seq in tempRes:
					f.write(",".join(Seq) + "\n")
	f.close()


#######################################################################################################################################
###	output data to vcf format
#######################################################################################################################################
def OutputVCFfile(OutputPrefix, IndelSupport, Genotype, P0, P1, P2, ShowGT):
	PassedIndelTag = map(lambda x: x[0]+"-"+x[1]+"-"+x[2]+"-"+x[3], IndelSupport)
	if ShowGT == 0:
		f = open(OutputPrefix + ".vcf",'w')
		detailForm = np.genfromtxt(OutputPrefix+ '.alt.detail.form', delimiter=':', dtype="|S10000")
		detailFormDict = dict(detailForm)
		modiIndelForm = map(lambda x:detailFormDict[x.replace("X", "23").replace("Y", "24")], PassedIndelTag)
		headerShow = "##fileformat=VCFv4.0\n\
##source=Eindel V1.0\n\
##INFO=<ID=SPO,Number=1,Type=Integer,Description=\"number of supporting reads on original reference\">\n\
##INFO=<ID=PSO,Number=1,Type=Integer,Description=\"number of pair-supporting reads in original reference\">\n\
##INFO=<ID=MGO,Number=1,Type=Integer,Description=\"number of supporting reads on original reference which can map to other locations of genomic reference\">\n\
##INFO=<ID=SprO,Number=1,Type=Integer,Description=\"spread of supporting reads on original reference\">\n\
##INFO=<ID=SPA,Number=1,Type=Integer,Description=\"number of pair-supporting reads on alternative reference\">\n\
##INFO=<ID=PSA,Number=1,Type=Integer,Description=\"number of paired supporting reads on original reference\">\n\
##INFO=<ID=MGA,Number=1,Type=Integer,Description=\"number of supporting reads on alternative reference which can map to other locations of genomic reference\">\n\
##INFO=<ID=MAA,Number=1,Type=Integer,Description=\"number of supporting reads on alternative reference which can map to other alternative reference(s)\">\n\
##INFO=<ID=SprA,Number=1,Type=Integer,Description=\"spread of supporting reads on original reference\">\n\
##FILTER=<ID=PASS,Number=1,Type=String,Description=\"Pass\">\n\
##FILTER=<ID=FAIL,Number=1,Type=String,Description=\"Non-indel\">\n\
##FORMAT=<ID=GT,Number=1,Type=Integer,Description=\"Genotype(0:non-indel, 1:heterozygous indel, 2:homozygous indel)\">\n\
##FORMAT=<ID=P1,Number=1,Type=Float,Description=\"probability of non-indel genotype\">\n\
##FORMAT=<ID=P2,Number=1,Type=Float,Description=\"probability of heterozygous indel\">\n\
##FORMAT=<ID=P3,Number=1,Type=Float,Description=\"probability of homozygous indel\">\n\
"
		sample_tag = OutputPrefix.split("/")[-1]
		IndelHeader = "#CHROM	POS	ID	REF	ALT	FILTER	INFO	FORMAT	" + sample_tag + "\n"
		f.write(headerShow)
		f.write(IndelHeader)
		IndelSupport[IndelSupport[:,0] =="23",0] = "X"
		IndelSupport[IndelSupport[:,0] =="24",0] = "Y"
		for iPhe in range(len(IndelSupport)):
			StartLoc = modiIndelForm[iPhe].split(",")[0]
			GT = Genotype[iPhe]; Prob0 = P0[iPhe]; Prob1 = P1[iPhe]; Prob2 = P2[iPhe]; 
			for iIndel in range((len((modiIndelForm[iPhe]).split(","))-2)/3):
				relatedLoc = str(int(modiIndelForm[iPhe].split(",")[iIndel*3 + 2]) + int(StartLoc))
		 		Allele2 = modiIndelForm[iPhe].split(",")[iIndel*3+3].replace("-", "")
		 		Allele1 = modiIndelForm[iPhe].split(",")[iIndel*3+4].replace("-", "")
				Filter = "FAIL"
				if GT>0:
					Filter = "PASS"
				if (Allele1==""): Allele1="-"
				if (Allele2==""): Allele2="-"
				f.write(IndelSupport[iPhe,0] + "\t" + relatedLoc + "\t.\t" + Allele1 + "\t" + Allele2 + "\t")
				f.write(Filter + "\t")
				f.write("SPO=" + IndelSupport[iPhe,4]+";PSO="+IndelSupport[iPhe,5]+";MGO="+IndelSupport[iPhe,6]+";SprO="+IndelSupport[iPhe,7]+";SPA="+IndelSupport[iPhe,8]+";PSA=")
				f.write(IndelSupport[iPhe,9]+";MGA="+IndelSupport[iPhe,10]+";MAA="+IndelSupport[iPhe,11]+";SprA="+IndelSupport[iPhe,12]+"\t")
				f.write("GT:P1:P2:P3\t")
				f.write(str(GT)+":"+str(Prob0)+":"+str(Prob1)+":"+str(Prob2))
				f.write("\n")
		f.close()
	if ShowGT == 0:
		f = open(OutputPrefix + ".txt",'w')
		f.write("CHROM\tPOS\tCOMT\tTYPE\tSIZE\tSPO\tPSO\tMGO")
		f.write("\tSprO\tSPA\tPSA\tMGA\tMAA\tSprA\tGT\tP1\tP2\tP3\n")
		IndelSupport[IndelSupport[:,0] =="23",0] = "X"
		IndelSupport[IndelSupport[:,0] =="24",0] = "Y"
		for iPhe in range(len(IndelSupport)):
			relatedLoc = IndelSupport[iPhe, 1]
			GT = Genotype[iPhe]; Prob0 = P0[iPhe]; Prob1 = P1[iPhe]; Prob2 = P2[iPhe];
			f.write(IndelSupport[iPhe,0] + "\t" + relatedLoc + "\t.\t" + IndelSupport[iPhe,2]+"\t"+IndelSupport[iPhe,3] + "\t")
			f.write(IndelSupport[iPhe,4]+"\t"+IndelSupport[iPhe,5]+"\t"+IndelSupport[iPhe,6]+"\t"+IndelSupport[iPhe,7]+"\t"+IndelSupport[iPhe,8]+"\t")
			f.write(IndelSupport[iPhe,9]+"\t"+IndelSupport[iPhe,10]+"\t"+IndelSupport[iPhe,11]+"\t"+IndelSupport[iPhe,12]+"\t")
			f.write(str(GT)+"\t"+str(Prob0)+"\t"+str(Prob1)+"\t"+str(Prob2))
			f.write("\n")
		f.close()
	print " Output vcf file: " + OutputPrefix + ".vcf"
	print " Output txt file: " + OutputPrefix + ".txt"


################# for rep regions#################################################################
def RetrieveRepRegionDict(ChrID, ConFile):
	data0 = np.genfromtxt(ConFile, delimiter="\t", dtype='str')
	DictRegion = {}
	ChrTag = data0[:,0]==ChrID
	for regionTag in range(sum(ChrTag)):
		for r in range(data0[ChrTag,1][regionTag].astype(int)/10,data0[ChrTag,2][regionTag].astype(int)/10):
			DictRegion[r]=''
	return DictRegion

def KeepNonRepIndel(ConFile):
	RepRegion = np.array([False for x in range(len(my_data))])
	for ChrID in np.unique(my_data[:,2]):
		ChrQuery = ChrID.replace("23", "Y").replace("24", "X")
		#print "	Removing Chr" + ChrQuery
		RepDict = RetrieveRepRegionDict(ChrQuery, ConFile)
		RepRegionTag = map(lambda x,y: RepDict.has_key(int(x)/10) or RepDict.has_key(int(y)/10), my_data[my_data[:,2]== ChrID,3], my_data[my_data[:,2]== ChrID,8])
		#RepRegionTag = np.logical_or(RepDict[], RepDict[my_data[my_data[:,1]== ChrID,2]])
		RepRegion[my_data[:,2]== ChrID] = np.array(RepRegionTag)
	return RepRegion

#######################################################################################################################################
###	Convert array to Dictionary
#######################################################################################################################################
def ConvertArray2ArrayDict(NumpyMatrix, Column):
	ResultDict = {}
	for x in NumpyMatrix:
		if ResultDict.has_key(x[Column]):
			ResultDict[x[Column]] = np.vstack((ResultDict[x[Column]], [x]))
		else:
			ResultDict[x[Column]] = np.array([x])
	return ResultDict

#######################################################################################################################################
###	Convert array to Dictionary
#######################################################################################################################################
def GetPairedSubReadsID(SubReadsID):
	if SubReadsID.find("/2") >0: countpart = SubReadsID.replace("/2", "/1")
	if SubReadsID.find("/1") >0: countpart = SubReadsID.replace("/1", "/2")
	return countpart

#######################################################################################################################################
###	Get Alignment information dictionary, with Reads ID as the keys
#######################################################################################################################################
def Get_AlignInfoDict_ByReadsID(Alignment):
	AlignmentHReads = []
	if len(Alignment)>0: AlignmentHReads = np.hstack((np.array(map(lambda x: [x[0].split("/")[0]], Alignment)), Alignment))
	AlignDict = ConvertArray2ArrayDict(AlignmentHReads, 0); AlignmentHReads = []
	return AlignDict

#######################################################################################################################################
###	Convert array to Dictionary
#######################################################################################################################################
def filterSupportReads(AltAlignment, OriAlignment, RefAlignment):
	global AllDataMatrix
	OriAlignDict = Get_AlignInfoDict_ByReadsID(OriAlignment)
	AltAlignDict = Get_AlignInfoDict_ByReadsID(AltAlignment)
	RefAlignDict = Get_AlignInfoDict_ByReadsID(RefAlignment)
	uniqueReads = np.unique(np.concatenate((AltAlignDict.keys(), OriAlignDict.keys())))
	NumReadsProcessed = 0
	ProcessTag = 0
	for read in uniqueReads:
		NumReadsProcessed = NumReadsProcessed+1
		if (NumReadsProcessed * 100 /len(uniqueReads)) > ProcessTag:
			sys.stdout.write("\r%f%%" % (NumReadsProcessed * 100 /len(uniqueReads)))
			ProcessTag = NumReadsProcessed * 100 /len(uniqueReads)
		if AltAlignDict.has_key(read):
			ReadsAlign = np.array(AltAlignDict[read])
			if OriAlignDict.has_key(read): ReadsAlign = np.vstack((ReadsAlign, OriAlignDict[read]))
		elif OriAlignDict.has_key(read): ReadsAlign = np.array(OriAlignDict[read])
		ReadAlignKeep = np.ones(ReadsAlign.shape[0])			###	this column record which alignment record should be kept (0 remove)
		################judge if the reads can mapped to both refs
		if ReadsAlign.shape[0] >1:					#########only judge those reads can be mapped to both reference 
			for subReads in np.unique(ReadsAlign[:,1]):		#########
				for RefSeqID in ReadsAlign[ReadsAlign[:,1]==subReads,2]:
					if RefSeqID.find("ori") > 0: countpart = RefSeqID.replace("ori", "alt")
					if RefSeqID.find("alt") > 0: countpart = RefSeqID.replace("alt", "ori")
					if any(ReadsAlign[:,2]==countpart):
						ReadAlignKeep[np.logical_or(ReadsAlign[:,2]==countpart, ReadsAlign[:,2]==RefSeqID)] = 0
		if sum(ReadAlignKeep==1)>0 :
			################judge if the support reads of genome reference has paired support
			GenRefsPairs = 0
			GenRefs = np.array([])
			if RefAlignDict.has_key(read): GenRefs = RefAlignDict[read]
			GePaired = np.zeros(len(GenRefs))
			for i in xrange(len(GenRefs)):
				GenRef = GenRefs[i]
				if GenRef[1].find("/2") >0: countpart = GenRef[1].replace("/2", "/1")
				if GenRef[1].find("/1") >0: countpart = GenRef[1].replace("/1", "/2")
				PairIndex = np.logical_and(np.logical_and(GenRefs[:,1] ==countpart, GenRefs[:,2] == GenRef[2]), np.abs(GenRefs[:,3].astype(int) - int(GenRef[3]))<1500)
				if any(PairIndex): GePaired[i] = min(np.abs(GenRefs[:,3].astype(int) - int(GenRef[3]))[PairIndex])
			###====================================================
			SPPaired = np.zeros(sum(ReadAlignKeep==1))
			IndRefs = ReadsAlign[ReadAlignKeep==1]
			if len(GenRefs) >0:
				for i in xrange(len(IndRefs)):
					IndRef = IndRefs[i]
					##---- find out id of paired reads
					if IndRef[1].find("/2") >0: countpart = IndRef[1].replace("/2", "/1")
					if IndRef[1].find("/1") >0: countpart = IndRef[1].replace("/1", "/2")
					##----
					IndLoc = int(IndRef[2].split("_")[1]) + int(IndRef[3])
					GenRefs_Chr = np.array(map(lambda x: x.replace("chr", "").replace("Chr", ""), GenRefs[:,2]))
					IndChr = IndRef[2].split("-")[0].replace("chr", "").replace("Chr", "")
					#################
					PairIndex = np.logical_and(np.logical_and(GenRefs[:,1] == countpart, GenRefs_Chr == IndChr), np.abs(GenRefs[:,3].astype(int) - IndLoc)<1500)
					if any(PairIndex): SPPaired[i] = min(np.abs(GenRefs[:,3].astype(int) - IndLoc)[PairIndex])
				###########find the most sutable pairs
				PairedDist = np.reshape(np.concatenate([GePaired, SPPaired]), (-1,1))
				ValidByMostPossible = np.hstack((np.vstack((GenRefs, IndRefs)), PairedDist))
				Distnum = ValidByMostPossible[:,5].astype("float").astype("int")
				IndGenalignBothRfPr = ValidByMostPossible[max(Distnum) ==Distnum]
			else:
				IndGenalignBothRfPr = np.hstack((IndRefs, np.reshape(SPPaired, (-1,1))))
			################find out reads status###################
			if any(map(lambda x: ("alt" in x) or ("ori" in x), IndGenalignBothRfPr[:,2])):
				#break
				ReadsChain = np.array(map(lambda x: x[1][-1], IndGenalignBothRfPr))
				AltIndex = np.array(map(lambda x: ("alt" in x[2]), IndGenalignBothRfPr))
				OriIndex = np.array(map(lambda x: ("ori" in x[2]), IndGenalignBothRfPr))
				GenIndex = np.array(map(lambda x: (not ("ori" in x[2])) and  (not ("alt" in x[2])), IndGenalignBothRfPr))
				###########judge for alternative reference
				if sum(AltIndex) >0:
					for alignRead in IndGenalignBothRfPr[AltIndex]:
						GenolocType = "UniGenLocs"
						if sum(np.logical_and(alignRead[1].split("/")[1] == ReadsChain, GenIndex)) > 0:
							GenolocType = "MulGenLocs"
						MulIndelType = "UniIndelLocs"
						if sum(np.logical_and(alignRead[1].split("/")[1] == ReadsChain, AltIndex)) > 1:
							MulIndelType = "MulIndelLocs"
						AllDataMatrix.append(alignRead.tolist() + [GenolocType, MulIndelType])
				###########judge for original reference
				MulIndelType = "NA"
				for alignRead in IndGenalignBothRfPr[OriIndex]:
					GenolocType = "UniGenLocs"
					if sum(np.logical_and(alignRead[1].split("/")[1] == ReadsChain, GenIndex)) > 1:
						GenolocType = "MulGenLocs"
					AllDataMatrix.append(alignRead.tolist() + [GenolocType, MulIndelType])

##########################################################################

##########################################################################
def OutputIndelSupportMatrix(AllDataMatrix):
	###
	tempDICT = {}
	for data in AllDataMatrix[:,2:]: tempDICT[":".join(data)] = ""
	###
	UniqDataMatrix = []
	for key in tempDICT.keys(): UniqDataMatrix.append(key.split(":"))
	###
	UniqDataMatrix = np.array(UniqDataMatrix)
	filterData = 1- float(UniqDataMatrix.shape[0])/AllDataMatrix.shape[0]
	print "\n	Removed PCR duplications.... " + str(filterData*100) + "% reads are filtered"
	###
	tempDICT = {}
	for data in UniqDataMatrix:
		IndelID = data[0].split("_")[0].replace("-alt", "").replace("-ori", "")
		if tempDICT.has_key(IndelID):
			tempDICT[IndelID] = np.vstack((tempDICT[IndelID], data))
		else:
			tempDICT[IndelID] = np.array([data.tolist()])
	###
	IndelSupport = []
	for IndelID in tempDICT.keys():
		#print IndelID
		Chr, Loc, Typ, Len = IndelID.split("-")
		IndelAlign = tempDICT[IndelID]
		#StartLocs = indelAlt.split("_")[1]
		AltSupReadTag = np.core.defchararray.rfind(IndelAlign[:,0], "alt") >0
		OriSupReadTag = np.core.defchararray.rfind(IndelAlign[:,0], "ori") >0
		if sum(OriSupReadTag) ==0:
			SpOri = 0; PsOri = 0; MulOri = 0; SprOri = 0
		else:
			SpOri = sum(OriSupReadTag)
			PsOri = sum(IndelAlign[OriSupReadTag,3] != "0.0")
			MulOri = sum(IndelAlign[OriSupReadTag,4] == "MulGenLocs")
			AlignLocs = IndelAlign[OriSupReadTag,1].astype(np.float32)
			SprOri = max(AlignLocs) - min(AlignLocs)
		if sum(AltSupReadTag) > 0:
			#SpOri = 0; PsOri = 0; MulOri = 0; SprOri = 0
			SpAlt = sum(AltSupReadTag)
			PsAlt = sum(IndelAlign[AltSupReadTag,3] != "0.0")
			MulAlt = sum(IndelAlign[AltSupReadTag,4] == "MulGenLocs")
			MMAlt = sum(IndelAlign[AltSupReadTag,5] == "MulIndelLocs")
			AlignLocs = IndelAlign[AltSupReadTag,1].astype(np.float32)
			SprAlt = max(AlignLocs) - min(AlignLocs)
			IndelSupport.append([Chr, Loc, Typ, Len, SpOri, PsOri, MulOri, SprOri, SpAlt, PsAlt, MulAlt, MMAlt, SprAlt])
	return np.array(IndelSupport)

#######################################################################################################################################
def removeDuplicateAlign(my_data):
	PlusChainPass = np.logical_and(my_data[:,3].astype(int) -my_data[:,8].astype(int) >0, my_data[:,4] =="1")
	MinusChainPass = np.logical_and(my_data[:,3].astype(int) -my_data[:,8].astype(int) <0, my_data[:,4] =="-1")
	ChainPass = np.logical_or(PlusChainPass, MinusChainPass)
	my_data = my_data[ChainPass]
	tempLocRec = {}
	removingItem = np.zeros(len(my_data))
	i = -1
	for Align in my_data:
		i = i + 1
		if (tempLocRec.has_key(".".join((Align[2], Align[3], Align[8]))) or tempLocRec.has_key(".".join((Align[2], Align[8], Align[3])))):
			removingItem[i] = 1
		else:
			tempLocRec[".".join((Align[2], Align[3], Align[8]))] = 1
			removingItem[i] = 0
	return my_data[removingItem==0]

#######################################################################################################################################
###	Remove Files
#######################################################################################################################################
def RmFile(filePath):
	try:
		os.remove(filePath); return 1
	except OSError:
		pass; return 0

#######################################################################################################################################
###	this function divide support reads alignment files into sub files
#######################################################################################################################################
def DivSupReadAlignFile(samFileAlt1, samFileAlt2, samFileOri1, samFileOri2, samFilePair1, samFilePair2, RefAlignPre, GenAlignPre):
	RefAlign = [samFileAlt1, samFileAlt2, samFileOri1, samFileOri2]
	GenAlign = [samFilePair1, samFilePair2]
	################
	AlinFNIndex = {}
	CurrentFileID = 1
	MaxID = CurrentFileID
	ReadsInCurrentFile = 0
	LastFileID = 0
	###	load genome reference aligment file	###
	removeTag = RmFile(RefAlignPre+"_"+str(CurrentFileID))							#Remove previous files
	for AlinFile in RefAlign:
		fileRead = open(AlinFile, 'r')
		for line in fileRead:
			line = line.rstrip()
			ReadsIDMain = line.split("\t")[0].split("/")[0]
			if AlinFNIndex.has_key(ReadsIDMain):
				FileID = AlinFNIndex[ReadsIDMain]
			else:
				if ReadsInCurrentFile >500000:
					ReadsInCurrentFile = 0
					CurrentFileID = CurrentFileID + 1
					MaxID = CurrentFileID
					removeTag = RmFile(RefAlignPre+"_"+str(CurrentFileID))
				ReadsInCurrentFile = ReadsInCurrentFile + 1
				FileID = CurrentFileID
				AlinFNIndex[ReadsIDMain] = FileID
			if FileID != LastFileID: fileWrite = open(RefAlignPre+"_"+str(FileID), "a")
			fileWrite.write(line + "\n")
			LastFileID = FileID
	###	load genome reference aligment file	###
	LastFileID = 0
	for iFile in range(1,(MaxID+1)): removeTag = RmFile(GenAlignPre+"_"+str(iFile))				#Remove previous files
	for AlinFile in GenAlign:
		fileRead = open(AlinFile, 'r')
		for line in fileRead:
			line = line.rstrip()
			ReadsIDMain = line.split("\t")[0].split("/")[0]
			if AlinFNIndex.has_key(ReadsIDMain):
				FileID = AlinFNIndex[ReadsIDMain]
				if FileID != LastFileID: fileWrite = open(GenAlignPre+"_"+str(FileID), "a")
				fileWrite.write(line + "\n")
				LastFileID = FileID
	del AlinFNIndex
	return MaxID

#######################################################################################################################################
### split and sort aligment file in candidate detection step
#######################################################################################################################################

def SplitAndSortSeqAlinFile(SegAlin1, SegAlin2, outputFileName):
	AlinFNIndex = {}
	CurrentFileID = 1
	ReadsInCurrentFile = 0
	LastFileID = 0
	fileRead = open(SegAlin1, 'r')
	for line in fileRead:
		line = line.rstrip()
		ReadsIDMain = line.split("\t")[0].split("/")[0]
		if AlinFNIndex.has_key(ReadsIDMain):
			FileID = AlinFNIndex[ReadsIDMain]
		else:
			if ReadsInCurrentFile >1000000:
				ReadsInCurrentFile = 0
				CurrentFileID = CurrentFileID + 1
				MaxID = CurrentFileID
				rmTag = RmFile(outputFileName+"_"+str(CurrentFileID))		##remove files
			ReadsInCurrentFile = ReadsInCurrentFile + 1
			FileID = CurrentFileID
			AlinFNIndex[ReadsIDMain] = FileID
		if FileID != LastFileID:
			fileWrite = open(outputFileName+"_"+str(FileID), "a")
		fileWrite.write(line + "\n")
		LastFileID = FileID
	fileRead = open(SegAlin2, 'r')
	for line in fileRead:
		line = line.rstrip()
		ReadsIDMain = line.split("\t")[0].split("/")[0]
		if AlinFNIndex.has_key(ReadsIDMain):
			FileID = AlinFNIndex[ReadsIDMain]
		else:
			if ReadsInCurrentFile >1000000:
				ReadsInCurrentFile = 0
				CurrentFileID = CurrentFileID + 1
				MaxID = CurrentFileID
				rmTag = RmFile(outputFileName+"_"+str(CurrentFileID))		##remove files
			ReadsInCurrentFile = ReadsInCurrentFile + 1
			FileID = CurrentFileID
			AlinFNIndex[ReadsIDMain] = FileID
		if FileID != LastFileID:
			fileWrite = open(outputFileName+"_"+str(FileID), "a")
		fileWrite.write(line + "\n")
		LastFileID = FileID
	return MaxID

def GetGenotype(IndelSupport, CnumSRAlt = 1, CnumPSAlt = 0, CsprAlt = 10, CMulRatOri = 1, CMulRatAlt = 1, CMulRatAlt2 = 1):
	##~~~~~~~ calculate genotype for each indel
	numSRAlt = IndelSupport[:,8].astype(float); numSROri = IndelSupport[:,4].astype(float)
	numPSAlt = IndelSupport[:,9].astype(float); numPSOri = IndelSupport[:,5].astype(float)
	sprAlt = IndelSupport[:,12].astype(float); sprOri = IndelSupport[:,7].astype(float)
	MulRatOri = np.zeros(IndelSupport.shape[0])
	MulRatAlt = np.zeros(IndelSupport.shape[0])
	MulRatAlt2 = np.zeros(IndelSupport.shape[0])
	MulRatOri[numSROri>0] = IndelSupport[numSROri>0,6].astype(float)/numSROri[numSROri>0]
	MulRatAlt[numSRAlt>0] = IndelSupport[numSRAlt>0,10].astype(float)/numSRAlt[numSRAlt>0]
	MulRatAlt2[numSRAlt>0] = IndelSupport[numSRAlt>0,11].astype(float)/numSRAlt[numSRAlt>0]
	##-------- filter candidate indels by a series critieria
	ValidHardFilter = np.logical_and.reduce((numSRAlt>CnumSRAlt, numPSAlt>CnumPSAlt, sprAlt>CsprAlt, \
						MulRatOri<CMulRatOri, MulRatAlt<CMulRatAlt, MulRatAlt2<CMulRatAlt2))
	##-------- for those indel passed criterion, get number of supoorting read
	IndelSupportValid = IndelSupport[ValidHardFilter,]
	Genotype = np.array([-1 for x in range(IndelSupport.shape[0])])
	P0 = np.array([0.0 for x in range(IndelSupport.shape[0])]); P1 = P0.copy(); P2 = P0.copy()
	StatusList = np.array(["Failed" for x in range(IndelSupport.shape[0])])
	##-------- only use indels with size less than 6 bps as training data
	x = map(lambda x: int(x[8]), IndelSupport[ValidHardFilter][IndelSupport[ValidHardFilter,3].astype(int)<=6,])
	n = map(lambda x: int(x[4])+int(x[8]), IndelSupport[ValidHardFilter][IndelSupport[ValidHardFilter,3].astype(int)<=6,])
	##-------- random sample 4000 indels from last step
	xTrain, nTrain = (np.array(x), np.array(n))
	xTrain = xTrain[nTrain<500]
	nTrain = nTrain[nTrain<500]
	xTrain = xTrain[:4000]
	nTrain = nTrain[:4000]
	xTrain = xTrain.tolist()
	nTrain = nTrain.tolist()
	##-------- estimate parameter for beta-binomial model
	a, b, a_2, b_2, a_3, b_3, w1, w2, w3 = ParameterEstimate_bBinomial(xTrain, nTrain)
	QGTag = np.logical_and(IndelSupport[:,3].astype(int)<=6, ValidHardFilter)
	################
	x = np.array(x)
	n = np.array(n)
	x[n>1000] = (x[n>1000] / ((n[n>1000]).astype(float)/1000)).astype(int)
	n[n>1000] = 1000
	x = x.tolist()
	n = n.tolist()
	##-------- predict genotype based on beta-binomial model
	Genotype[QGTag], P0[QGTag], P1[QGTag], P2[QGTag] = prediction_genotype(x, n, a, b, a_2, b_2, a_3, b_3, w1, w2, w3)
	Genotype[ValidHardFilter == False] =  0
	##-------- for those indels with size > 6 bps, predict genotype
	QGCtr = np.logical_and(ValidHardFilter, IndelSupport[:,3].astype(int)>6)
	Genotype[QGCtr] = 2
	##-------- if supporting reads of both original references and alternative reference are qualified
	##-------- the indel will be asigned as a heterozygous
	Genotype[np.logical_and(QGCtr, np.logical_and.reduce((numSROri>2, numPSOri>1, sprOri>10)))] = 1
	return Genotype, P0, P1, P2

############################### Functions ############################################################################################################################################################
######################################################################################################################################################################################################
#Process = "C"
#Fastq1 = "/media/stephen/Elements/Eindel/fasta/dw_generated_40x_100bps.20141112.read1.fastq"
#Fastq2 = "/media/stephen/Elements/Eindel/fasta/dw_generated_40x_100bps.20141112.read2.fastq"
#ParamFile = "/home/zengshuai/dist/Eindel.conf"
#OutputPrefix = "/media/stephen/Elements/Eindel/simDataEindelResult/SimulatedData/T2014111240x100bp/T2014111240x100bp"
#SegLen = 22
#MaxSize = 1000
#filCanIndbyRepRegion = 1
#################### for debug use only #####################
debug = 0
if debug ==1:
	Process = "V"
	Fastq1 = "/media/stephen/Elements/Eindel/fasta/CDtrio/A_read1.fastq"
	Fastq2 = "/media/stephen/Elements/Eindel/fasta/CDtrio/A_read2.fastq"
	OutputPrefix = "/media/stephen/Elements/Eindel/simDataEindelResult/CDTrio/A"
	SegLen = 22
	MaxSize = 1000
	ParamFile = "/media/stephen/Elements/Eindel/Eindel.conf"

#################################################
print "================ Parameter & Option ================"
print "Spindel version: 0.1.0"
InputValid = 0
ErrorInfo = CheckInputValid(Process, Fastq1, Fastq2, OutputPrefix, SegLen, MaxSize, CPUnum, IncRepRegion, MinSpCDetect)
if ErrorInfo !="":
	print "\n############## Error Info ####################\n" + ErrorInfo
else:
	InputValid = 1
	readLength, FastqValid = CheckFastqValid(Fastq1, Fastq2)
	bowtiePath, bowtie2Path, BowRefIndex, Bow2RefIndex, GenomeReffasta, samtoolsPath, RepRgionPath = LoadParameters(ParamFile)
	print "Read Length is " + str(readLength);

print "================ Parameter & Option ================\n"

if (Process=="C" or Process=="A") and (InputValid == 1):
	##############intial mapping, user can specify number of core(CPU) to process initial mapping
	print "############intial mapping############################"
	if CPUnum>1:	#mapping with parallel process
		p = sp.Popen(bowtie2Path +"/bowtie2 -k 10 --gbar " + str(readLength) + " --no-unal --score-min L,-3,-0.12 --no-hd --quiet --omit-sec-seq -x " + Bow2RefIndex + " -U "+ Fastq1 + " -S " + OutputPrefix + ".1.sam  --un " + OutputPrefix + ".1.unmapped", stdout=sp.PIPE, shell=True)
		p2 = sp.Popen(bowtie2Path +"/bowtie2 -k 10 --gbar " + str(readLength) + " --no-unal --score-min L,-3,-0.12 --no-hd --quiet --omit-sec-seq -x " + Bow2RefIndex + " -U "+ Fastq2 + " -S " + OutputPrefix + ".2.sam  --un " + OutputPrefix + ".2.unmapped", stdout=sp.PIPE, shell=True)
		(output1, err1) = p.communicate()
		(output2, err2) = p2.communicate()
	else:					#mapping reads to genome one by one
		p = sp.Popen(bowtie2Path +"/bowtie2 -k 10 --gbar " + str(readLength) + " --no-unal --score-min L,-3,-0.12 --no-hd --quiet --omit-sec-seq -x " + Bow2RefIndex + " -U "+ Fastq1 + " -S " + OutputPrefix + ".1.sam  --un " + OutputPrefix + ".1.unmapped", stdout=sp.PIPE, shell=True)
		(output1, err1) = p.communicate()
		p2 = sp.Popen(bowtie2Path +"/bowtie2 -k 10 --gbar " + str(readLength) + " --no-unal --score-min L,-3,-0.12 --no-hd --quiet --omit-sec-seq -x " + Bow2RefIndex + " -U "+ Fastq2 + " -S " + OutputPrefix + ".2.sam  --un " + OutputPrefix + ".2.unmapped", stdout=sp.PIPE, shell=True)
		(output2, err2) = p2.communicate()
		
	##############cut unmapped reads into three segments, and length of first segments and last segments is SegLen. 
	print "############unmapped reads cutting####################"
	cut_Unmapped_Reads(OutputPrefix + ".1.unmapped", OutputPrefix + ".1.seg.fastq", SegLen)
	cut_Unmapped_Reads(OutputPrefix + ".2.unmapped", OutputPrefix + ".2.seg.fastq", SegLen)
	
	##############Second mapping, map segments(prefix and suffix) to genome references
	print "############segment alignment#########################"
	SegAlignment(OutputPrefix + ".1.seg.fastq", OutputPrefix + ".2.seg.fastq", OutputPrefix + ".1.seg.alin", OutputPrefix + ".2.seg.alin", 20, CPUnum)
	
	##############Due to limitation of memory, split files and sort aligments
	MaxID = SplitAndSortSeqAlinFile(OutputPrefix + ".1.seg.alin", OutputPrefix + ".2.seg.alin", OutputPrefix + ".SortSegAlign")
	
	##############
	print "############read-seg alignment validating#############"
	for Fileid in range(1,(MaxID+1)):
		AlinInfo = form_Presuf_Matrix(OutputPrefix + ".SortSegAlign_" + str(Fileid))
		rmfileTag = RmFile(OutputPrefix + ".SortSegAlign_" + str(Fileid))
		judgeSegPairValid(OutputPrefix, AlinInfo, Fileid)
	del AlinInfo
	#print "\n"
	my_data = np.genfromtxt( OutputPrefix+ '_PSpair.csv', delimiter=',', dtype="|S200")
	############## by default, indels within repetitive regions are removed
	if filCanIndbyRepRegion == 1:
		print "############Remove Indels in Repetitive Region########"
		IndelInRepRegion = KeepNonRepIndel(RepRgionPath)						## find out which indels are within repetitive region
		my_data = my_data[IndelInRepRegion == False]								## filter those indels within repetitive region
	#################################
	print "############Alt/Ori reference sequencing building#####"
	my_data = removeDuplicateAlign(my_data)												## remove duplications and consolidate indels
	#################################
	formatAltRefSeq(OutputPrefix, Fastq1, Fastq2, MinSpCDetect-1)	## build reference sequences for alternative and original alleles
	del my_data
	
	########################### this is the end of candidate-detection step ################
	
	
if (Process == "V" or Process == "A") and (InputValid == 1):
	###
	print "\n############Second run mapping########################"
	AlignAlt1, AlignOri1, AlignAlt2, AlignOri2 = WholeAlignOriAlt(OutputPrefix, Fastq1, Fastq2, 1)
    
	print "############Build consensus sequences#################"
	buildPileup(samtoolsPath, OutputPrefix)
	buildConsensusSeq(OutputPrefix, "ori")
	buildConsensusSeq(OutputPrefix, "alt")
    
	print "############Third run mapping#########################"
	WholeAlignOriAlt(OutputPrefix, Fastq1, Fastq2, 2)
	###	load and divide files into sub-files	###
	samFileAlt1 = OutputPrefix + "_Con_alt.1.sam"; samFileAlt2 = OutputPrefix + "_Con_alt.2.sam"
	samFileOri1 = OutputPrefix + "_Con_ori.1.sam"; samFileOri2 = OutputPrefix + "_Con_ori.2.sam"
	samFilePair1 = OutputPrefix + ".1.sam"; samFilePair2 = OutputPrefix + ".2.sam"
	RefAlignPre = OutputPrefix + "_Con_Align.sam"; GenAlignPre = OutputPrefix + "_SPRgenome.sam"
	MaxID = DivSupReadAlignFile(samFileAlt1, samFileAlt2, samFileOri1, samFileOri2, samFilePair1, samFilePair2, RefAlignPre, GenAlignPre)
	AllDataMatrix = []
	for iFile in range(1,(MaxID+1)):
		print " Processing the " + str(iFile) + " file. Totally " + str(MaxID) + " files: "
		AltAlignment, OriAlignment, RefAlignment = MakeAlignArray(RefAlignPre + "_" + str(iFile), GenAlignPre + "_" + str(iFile))
		filterSupportReads(AltAlignment, OriAlignment, RefAlignment)
		os.remove(RefAlignPre + "_" + str(iFile))
		os.remove(GenAlignPre + "_" + str(iFile))
	IndelSupport = OutputIndelSupportMatrix(np.array(AllDataMatrix))
	del AllDataMatrix
	#############
	print "############Predict genotype##########################"
	Genotype, P0, P1, P2 = GetGenotype(IndelSupport, CnumSRAlt-1, CnumPSAlt-1, CsprAlt, CMulRatOri, CMulRatAlt, CMulRatAlt2)
	print "############Output result files#######################"
	OutputVCFfile(OutputPrefix, IndelSupport, Genotype, P0, P1, P2, 0)
	if KeepTemp == 0:
		RemoveTempFile(OutputPrefix)
