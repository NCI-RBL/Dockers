# This file holds the programs used in the part of the workflow invovled in creating ribosome density data structures.
import argparse
import pysam 
import twobitreader
import GFF
from Bio.Seq import Seq
#from rpgetreads import SeqReadGetter

### For hg19, hg38, and mm10
# Build density file based on
# 1. transcripts (e.g. canonical transcripts etc. )
# 2. GTF file
# 3. bam file input
# 4. .2bit file for genome sequence. This can be commented out. 

####### Note GFFparser parse 1-based GFF file coordinates start position into 0-based. BUT!!!!, end position still 1-based !!!!
# BAM files give 0-based coordinates. 
# totalreads is the total mapped reads from alignment.
# mapped reads for normalization could be either reads that mapped to known transcripts or the total mapped reads from alignment. 

class densebuilder(object):
	def __init__(self,bamfile,GTFgen,genome,riboshiftdict,threshold,totreads,outputdata,assignment,bamfileoutput):
		self.bamfile= bamfile
		self.GTFgen= GTFgen
		self.genome= genome
		self.assignment= assignment
		self.riboshiftdict= riboshiftdict
		self.threshold= threshold
		self.totreads= totreads
		self.outputdata= outputdata
		self.bamfileout= bamfileout

	def builddense(self):
		transcriptdict= {}
		mappedlocalreads= 0
		dumpedreads= 0
		illegalreads= 0
		tooshortlongreads= 0
		wrongstrandreads= 0
		totalreads= 0	# not totreads
		GFFlist= self.makeGFFlist(self.GTFgen)
		
		for chrom in GFFlist:
			if len(chrom)> 6:	continue	# Quick and dirty 
			transcriptnum= -1
			for transcript in GFFlist[chrom].features:	
				transcriptnum+= 1
				trsp_id= transcript.id # it is a number 
				trsp_strand= transcript.strand
				trsp_chromstart= int(transcript.location.start.position)  # 0-based
				trsp_chromend= int(transcript.location.end.position) 
				transcriptlist= [0.0 for x in range(abs(trsp_chromend- trsp_chromstart))] # a list for transcript (pre-mRNA), not CDS

				gb= self.getbam_5or3counts(self.bamfile,transcriptlist,chrom,transcriptnum,trsp_chromstart,trsp_chromend,trsp_id,trsp_strand,self.riboshiftdict,self.assignment,self.bamfileout) # return a riboshifted list (0-based) of unspliced counts 

				mappedlocalreads+= gb[1]
				dumpedreads+= gb[2]
				illegalreads+= gb[3]
				tooshortlongreads+= gb[4]
				wrongstrandreads+= gb[5]
				totalreads+= gb[6]
				exonsplicedseq= Seq('')
				exonsplicedcounts= []
				transcriptseq= Seq(genome[chrom][trsp_chromstart: trsp_chromend])

				startcodonspos= stopcodonpos= 'absence'
				for item in GFFlist[chrom].features[transcriptnum].sub_features:
					if trsp_strand== 1:
						if item.type== 'exon': # or item.type== 'CDS':	# For yeast, use 'CDS'
							exonstart= int(item.location.start.position)  # 0-based position
							exonend= int(item.location.end.position) # not 0-based
							exonstart_feat= exonstart- trsp_chromstart
							exonend_feat= exonend- trsp_chromstart # Not 0-based, it is fine for length....next line. 
							exonsplicedcounts+= gb[0][exonstart_feat:exonend_feat]	# takes from exonstart to exonend-1
							exonsplicedseq+= transcriptseq[exonstart_feat:exonend_feat] # takes from exonstart to exonend-1
						if item.type== 'start_codon':
							startcodonpos= item.location.start.position # 0-based position
							startcodonmrnapos=  self.chrpostomrnapos(startcodonpos,chrom,transcriptnum,GFFlist)	# spliced mRNA position
						if item.type== 'stop_codon':
							stopcodonpos= item.location.end.position- 1 # 0-based position
							stopcodonmrnapos= self.chrpostomrnapos(stopcodonpos,chrom,transcriptnum,GFFlist)

					if trsp_strand== -1:
						transcriptseq_rev= transcriptseq.reverse_complement()
						if item.type== 'exon':
							exonstart= int(item.location.start.position)  # 0-based position
							exonend= int(item.location.end.position)	# not 0-based 
							exonstart_feat= (trsp_chromend-1)- (exonend- 1) 		# 0-based
							exonend_feat= (trsp_chromend-1)- exonstart 		# 0-based
							exoncounts= gb[0][exonstart_feat:exonend_feat+ 1]	# both 0-based, need to +1 for length
							exonsplicedcounts= exoncounts+ exonsplicedcounts  # exoncounts added to the upstream of existing counts, so don't flip again.
							exonseq= transcriptseq_rev[exonstart_feat:exonend_feat+ 1] 
							exonsplicedseq= exonseq+ exonsplicedseq
						if item.type== 'start_codon':
							startcodonpos= item.location.end.position- 1 # Need to -1 to be 0-based.
							startcodonmrnapos= self.chrpostomrnapos(startcodonpos,chrom,transcriptnum,GFFlist)
						if item.type== 'stop_codon':
							stopcodonpos= item.location.start.position	# start.position is 0-based already. 
							stopcodonmrnapos= self.chrpostomrnapos(stopcodonpos,chrom,transcriptnum,GFFlist) 

				#if startcodonpos!= 'absence' or stopcodonpos!= 'absence':
				cdsseq= exonsplicedseq[startcodonmrnapos: stopcodonmrnapos+ 1] # take from startcodonmrnapos to stopcodonmrnapos
				cdscounts= exonsplicedcounts[startcodonmrnapos: stopcodonmrnapos+ 1] # take from startcodonmrnapos to stopcodonmrnapos
					### Sanity check 
				if str(cdsseq[:3].upper())!= "ATG":	continue	# ignore non-AUG start codons
				stopcodon= str(cdsseq[-3:].upper())
				if stopcodon!= "TGA" and stopcodon!= "TAG" and stopcodon!= "TAA":	continue	# ignore weird stop codons
					#utr5len= startcodonmrnapos
					#utr3len= len(exonsplicedseq)- stopcodonmrnapos- 1

				if sum(cdscounts)>= float(self.threshold): # thresholding minimal reads per CDS. 
					transcriptdict[trsp_id]= exonsplicedcounts
		if self.totreads== '-1':
			print(str(totalreads)+" total mapped reads used for normalization.")
			self.norm_m(transcriptdict,totalreads) # Normalzied by total reads mapped to transcriptdict only... but not total mapped reads.
		else:
			print(str(self.totreads)+" total mapped reads from alignment used for normalization.")
			self.norm_m(transcriptdict,self.totreads)

		self.writecountsf(transcriptdict, self.outputdata)

		# Write output file of comments.
		fc= open(self.outputdata+"output.txt","w")
		fc.write("Density was built with parameters:\n")
		fc.write("riboshiftdict="+str(self.riboshiftdict)+"\n")
		fc.write("threshold="+str(self.threshold)+"\n")
		fc.write("assignment="+str(self.assignment)+"\n")
		fc.write("totreads="+str(self.totreads)+"\n")
		fc.write("reads mapped to known canonical coding transcripts: " +str(mappedlocalreads)+ "\n")
		fc.write("reads are dumped, due to weird cigar codes: " +str(dumpedreads)+"\n")
		fc.write("reads are illegal, mapped outside of annotated transcripts: " +str(illegalreads)+ "\n")	
		fc.write("reads are too short/long: " +str(tooshortlongreads)+ "\n")
		fc.write("reads are on the wrong strand: " +str(wrongstrandreads)+ "\n")
		fc.write("total mapped reads from aligner: "+ str(totalreads))
		fc.close()

		print(str(mappedlocalreads)+" reads within length limitation mapped to known canonical coding transcripts. ")
		print(str(dumpedreads)+" reads are dumped, due to weird cigar codes.")
		print(str(illegalreads)+" reads are illegal, mapped outside of annotated transcripts.")
		print(str(tooshortlongreads)+" reads are too short/long.")
		print(str(wrongstrandreads)+" reads are on the wrong strand.")
		print(str(totalreads)+" total mapped reads from aligner. ")


	def getbam_5or3counts(self,bamfile,transcriptlist,chrom,transcriptnum,trsp_chromstart,trsp_chromend,trsp_id,trsp_strand,riboshiftdict,assignment,bamfileout):
		mappedlocalreads= 0
		dumpedreads= 0
		illegalreads= 0
		tooshortlongreads= 0
		wrongstrandreads= 0
		totalreads= 0

		for read in bamfile.fetch(chrom,trsp_chromstart,trsp_chromend):
			totaljunctlen= 0
			readlen= len(read.seq)
			quitsignal= 0
			totalreads+= 1
			for element in read.cigar:
				if element[0]== 0: # (M) alignemnt match
					pass
				elif element[0]== 1: # (I) insertion to reference 
					readlen-= element[1] 
				elif element[0]== 2: # (D) deletion from the reference
					readlen+= element[1]
				elif element[0]== 3: # (N) skipped region from the reference, --> junction(s) in mapped reads
					totaljunctlen+= element[1]
				elif element[0]== 4:    # soft clipping from ends, discard softclipped nts 
					pass
				else:
					quitsignal= 1

			if quitsignal== 1:
				dumpedreads+= 1
				continue

			if readlen in riboshiftdict:  # filter read length of interest
				riboshift= riboshiftdict[readlen][0]
			else:
				tooshortlongreads+= 1
				continue

			if read.is_reverse== False: read_strand= 1
			else: read_strand= -1

			if read_strand!= trsp_strand:   
				wrongstrandreads+= 1
				continue

			if self.assignment== '5':
				if not read.is_reverse:  # mapped to forward strand
					junctlen= self.junctlen_for_riboshift(read.cigar, riboshift) # Determine the length of junctions required for riboshift
					read5end= int(read.pos)- trsp_chromstart+ junctlen+ riboshift   #riboshifted 5'end of read on pre-mRNA 
					
					if read5end >= len(transcriptlist) or read5end< 0:  # possible due to annotation of different isoforms, a failsafe just for now....
						illegalreads+= 1
						continue 

				else: # mapped to reverse strand
					junctlen= self.junctlen_for_riboshift(read.cigar[::-1], riboshift) # Determine the length of junctions required for riboshift. Reverse read.cigar, couting form 3'end of reads
					read5end_chrompos= read.pos+ totaljunctlen+ (readlen- 1)  # the chrom pos of read5end 
					read5end= (trsp_chromend-1)- read5end_chrompos+ junctlen+ riboshift  # trsp_chrmend -1 refers to the start position (0-based) where trsp_chromstart:trsp_chromend takes from trsp_chromstart to trsp_chromend-1 position.

					if read5end >= len(transcriptlist) or read5end< 0:  # possible due to annotation of different isoforms 
						illegalreads+= 1
						continue 

				transcriptlist[read5end]+= 1

			if self.assignment== '3':
				if not read.is_reverse:  # mapped to forward strand
					junctlen= self.junctlen_for_riboshift(read.cigar[::-1], abs(riboshift)) # Determine the length of junctions required for riboshift 
					read3end_chrompos= int(read.pos)+ totaljunctlen+ (readlen- 1)   # the chrom position of read3end 
					read3end= read3end_chrompos- trsp_chromstart- junctlen+ riboshift  # riboshift is a negative value in 3'end alignment.

					if read3end >= len(transcriptlist) or read3end< 0:  # possible due to annotation of different isoforms, a failsafe just for now....
						illegalreads+= 1
						continue 

				else: # mapped to reverse strand
					junctlen= self.junctlen_for_riboshift(read.cigar, abs(riboshift)) # Determine the length of junctions required for riboshift. Reverse read.cigar, couting form 3'end of reads
					read3end= (trsp_chromend-1)- int(read.pos)+ junctlen+ riboshift  # riboshift is a negative value in 3'end alignment. 

					if read3end >= len(transcriptlist) or read3end< 0:  # possible due to annotation of different isoforms 
						illegalreads+= 1	
						continue 

				transcriptlist[read3end]+= 1

			mappedlocalreads+= 1
			bamfileout.write(read)

		return [transcriptlist,mappedlocalreads,dumpedreads,illegalreads,tooshortlongreads,wrongstrandreads,totalreads]


    # Function to determine the combined junction length for riboshift. 
    # For plus strand only.
	def junctlen_for_riboshift(self,readcigar,riboshift):
		junctnum= 0
		mappedexonlen= 0
		insert= 0
		for x in range(len(readcigar)):
			if readcigar[x][0]== 0:	mappedexonlen+= readcigar[x][1]
			if mappedexonlen-1 < riboshift:	junctnum+= 1

		for y in range(junctnum):
			if readcigar[y][0]== 3:	insert+= readcigar[y][1]

		return insert

	# Simple program, give it chrom pos and it tells you the mrna position.
	# 0 is the start of spliced transcripts.
	# input position should be 0-based!!!
	# output is also 0-based. 
	def chrpostomrnapos(self,chrpos,chrom,featnum,GFFlist):
	    trsp_id= GFFlist[chrom].features[featnum].id
	    trsp_strand= GFFlist[chrom].features[featnum].strand
	    trsp_chromstart= int(GFFlist[chrom].features[featnum].location.start.position)  # 0-based
	    trsp_chromend= int(GFFlist[chrom].features[featnum].location.end.position)
	    sublist=[]

	    for subfeature in GFFlist[chrom].features[featnum].sub_features:     # Make list of features
	        if subfeature.type== 'exon':
	            start= subfeature.location.start.position
	            end= subfeature.location.end.position
	            sublist.append([start,end])

	    if trsp_strand== -1:    
	        sublist.reverse()
	    assert len(sublist)!= 0
	    
	    prevexonlen= 0 
	    for item in sublist:
	        exonstart= item[0]
	        exonend= item[1]
	        exonlen= exonend- exonstart
	        
	        if trsp_strand== 1:
	            if chrpos>= exonstart and chrpos< exonend:      
	                mrnapos= prevexonlen+ chrpos- exonstart
	                return mrnapos
	            else:	prevexonlen+= exonlen
	        else:
	            if chrpos< exonend and chrpos>= exonstart:
	                mrnapos= prevexonlen+ (exonend-1)- chrpos       # Need -1 because end is not actual end, it is 1 beyond end.
	                return mrnapos
	            else:	prevexonlen+= exonlen 


	# Divide every read value by number of mapped reads. Values then in rpm.
	# No need to return readcounts since it is a mutable type being a list.
	def norm_m(self,readcounts,reads):
	    for chrom in readcounts.keys():
	        for position in range(len(readcounts[chrom])):
	            readcounts[chrom][position]/=float(reads)
	            readcounts[chrom][position]*=1E6  


	# Function to write counts files to disk - Use for float.
	def writecountsf(self,resultlist,filestring):  #Resultlist it is the list of counts, filestring is the file prefix for each chr to be written.
	    import struct
	    f2= open(filestring+"keys","w")
	    for trsp in resultlist.keys():
	        f= open(filestring+ trsp, "wb")
	        for position in resultlist[trsp]:
	            f.write(struct.pack("f",position))
	        f.close()   
	        f2.write(trsp+ "\n")
	    f2.close() 


	# Tool for loading the entire genome into memory.
	def makeGFFlist(self, GTFgen):
		GTFlist={}
		for chr in self.GTFgen:
			GTFlist[chr.id]=chr
		return GTFlist

if __name__== '__main__':
	parser= argparse.ArgumentParser()
	parser.add_argument('--bamfileinput', help= 'sorted and indexed bam file from alignment')
	parser.add_argument('--GTFfile', help= 'coding GTF file')
	parser.add_argument('--twobitfile', help= 'twobitfile of genome')
	parser.add_argument('--assignment', help= '5 or 3 end', required= True)
	parser.add_argument('--riboshiftdict', help= 'dictionary of riboshifts')
	parser.add_argument('--threshold', default= 1, help= 'thresholding read counts')
	parser.add_argument('--totreads', default= -1, help= 'reads for normalization')
	parser.add_argument('--outputdata', help= 'output data filepath')
	parser.add_argument('--bamfileoutput', help= 'output bam file')
	parser.add_argument
	args = parser.parse_args()

	import ast
	riboshiftdict= ast.literal_eval(args.riboshiftdict) #convert string into dictionary
	#print "parsing gff..."
	GTFgen= GFF.parse(args.GTFfile)
	#print "loading bam file..."
	bamfile= pysam.AlignmentFile(args.bamfileinput, "rb")
	#print "loading genome..."
	genome= twobitreader.TwoBitFile(args.twobitfile)
	print("writing bam out file...")
	bamfileout= pysam.AlignmentFile(args.bamfileoutput, "wb", template= bamfile)


	rfpdense= densebuilder(bamfile,GTFgen,genome,riboshiftdict,int(args.threshold),args.totreads,args.outputdata,args.assignment,bamfileout)
	rfpdense.builddense()
