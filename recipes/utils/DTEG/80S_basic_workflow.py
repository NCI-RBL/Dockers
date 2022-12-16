#!/usr/bin/env python

__author__ = 'Colin Wu'
__modified_by__ = 'Wilfried Guiblet'

# A basic workflow for 80S profiling
# Only contains 4 basic functions

#sys.path.append("/mnt/gridftp/wuc11/code/80S") # add path containing other scripts that might need to be imported

import os, shutil, sys, argparse, subprocess, csv
from pathos.multiprocessing import ProcessingPool as Pool  # A workaround that makes multiprocessing possible in a class:  https://www.py4u.net/discuss/12126   ### pip install --user pathos


class SubunitProfiling(object):

	def __init__(self,datapath,exp,ncrnaDir,genomeDir,GTFfile,twobitfile,mismatch,fiveprimetrim,linker,threshold,filtermodule,exclusionmodule,equalweight,UTRfilestring,ftsize,regionlength5,regionlength3,normalization,Rpath,sizedistr_outpath,outpath):
		self.datapath= datapath
		self.exp= exp
		self.ncrnaDir= ncrnaDir
		self.genomeDir= genomeDir
		self.GTFfile= GTFfile
		self.twobitfile= twobitfile
		self.mismatch= mismatch
		self.fiveprimetrim= fiveprimetrim
		self.linker= linker
		self.threshold= threshold
		self.filtermodule= filtermodule
		self.exclusionmodule= exclusionmodule
		self.equalweight= equalweight
		self.UTRfilestring= UTRfilestring
		self.ftsize= ftsize
		self.regionlength5= regionlength5
		self.regionlength3= regionlength3
		self.normalization= normalization
		self.Rpath= Rpath
		self.outpath= outpath
		self.sizedistr_outpath= sizedistr_outpath
		print('Init Success')

	def StarAlignment(self):
		rootpath= '%s/%s' %(self.datapath, self.exp)
		#commandstring= 'python StarAlign_main.py --rootpath %s --file %s --ncrnaDir %s --genomeDir %s --mismatch %s --fiveprimetrim %s --linker %s' %(rootpath, self.exp, self.ncrnaDir, self.genomeDir, self.mismatch, self.fiveprimetrim, self.linker)
		commandstring= 'StarAlign_main.py --rootpath %s --file %s --ncrnaDir %s --genomeDir %s --mismatch %s --fiveprimetrim %s --linker %s' %(rootpath, self.exp, self.ncrnaDir, self.genomeDir, self.mismatch, self.fiveprimetrim, self.linker)
		print(commandstring)
		os.system(commandstring)
		print('STAR Align Success')

	def SizeDistr(self):
		import collections
		rootpath= '%s/%s' %(self.datapath, self.exp)
		mismatch= str(self.mismatch).replace('.','')
		bamfileinput= '%s/%s_starM%s/%s_match.bam' %(rootpath,self.exp,mismatch,self.exp)
		minlen,maxlen= self.ftsize.strip("'").split(",")[0],self.ftsize.strip("'").split(",")[-1]
		distribution= collections.OrderedDict()
		for i in range(int(minlen),int(maxlen)+1): distribution[i]= 0

		self.count_bam2(distribution,bamfileinput)   # dont count soft clipped nt

		totcounts= 0
		for size,counts in distribution.items():	totcounts+= int(counts)
		
		outlist=[]
		headers= ["size","counts","fraction"]
		outlist.append(headers)

		for key,value in distribution.items():
			if totcounts!= 0:	fraction= float(value)/float(totcounts)
			######## ADDED BY WIL ########
			else: fraction = 0 
			##############################
			line= [key,value,fraction]
			outlist.append(line)

		if not os.path.exists(self.sizedistr_outpath): os.makedirs(self.sizedistr_outpath)
		csvoutfile= '%s/%s_%sto%s_sizedistr_M%s.csv' %(self.sizedistr_outpath,self.exp,minlen,maxlen,mismatch)
		self.writerows(outlist,csvoutfile)

		# Plot size distribution using R
		Rcode= '%s/size_distribution.R' %(self.Rpath)
		if not os.path.exists(self.outpath): os.makedirs(self.outpath)
		pdfoutfile= '%s/%s_%sto%s_sizedistr_M%s.pdf' %(self.outpath,self.exp,minlen,maxlen,mismatch)

		input= 'Rscript %s %s %s %s %s %s' %(Rcode,csvoutfile,minlen,maxlen,self.exp,pdfoutfile)
		subprocess.call(input, shell= True)
		print('Size Distrib Success')

	# input file is bam file
	# Count soft clipped nt
	def count_bam(self,distribution,bamfile):
		import pysam
		totalreads= 0
		bam= pysam.AlignmentFile(bamfile, "rb")
		for read in bam.fetch():  # iterate through every read
			if read.is_unmapped:	continue # Ideally, tophat outputs only mapped reads 
			#if str(read.flag)!= "0" and str(read.flag)!= "16":	continue # Filter out secondary alignments here if needed. 
			length= len(read.seq)
			if length in distribution:
				distribution[length]+= 1
			totalreads+=1
		return distribution


	# Don't count soft clipped nt 
	def count_bam2(self,distribution,bamfile):
		import pysam
		totalreads= 0
		bam= pysam.AlignmentFile(bamfile, "rb")
		for read in bam.fetch():
			if read.is_unmapped:	continue 
			length= 0
			for element in read.cigar:
				if element[0]!= 4:	length+= element[1]	
			if length in distribution:	distribution[length]+= 1
			totalreads+=1
		return distribution


	def writerows(self,intable,outfilestring):
		import csv
		fc= open(outfilestring, "w")
		writer= csv.writer(fc)
		writer.writerows(intable)
		fc.close()


	def GetMappedReads(self):
		#countsfile= '/mnt/gridftp/wuc11/data/count_dict.csv'
		f= open(self.countsfile,'U')
		table= self.readOrderedDict(f)
		rootpath= '%s/%s' %(self.datapath, self.exp)
		mismatch= str(self.mismatch).replace('.','')
		logfile= '%s/%s_starM%s/Log.final.out' %(rootpath,self.exp,mismatch)
		infile= open(logfile,'U')
		reader= csv.reader(infile, delimiter= '\t')
		for line in reader: 
			if len(line)<= 1:	continue
			if line[0]!= '                   Uniquely mapped reads number |':	continue

			mappedreads= line[1]
			if sample not in table:	table[sample]= [mappedreads]
		infile.close()
		f.close()
		self.writedicttocsv2(table,self.countsfile)


	# Make a csv back into a dictionary. f is open file handle. Handles multiple hits for same gene by putting in a _# term for them.
	def readOrderedDict(self,f):
		from collections import OrderedDict
		previousgene= ""
		counter= 1
		filegen= csv.reader(f,delimiter=',')
		output= OrderedDict()
		for gene in filegen:
			if gene[0]== previousgene:
				modgenename= gene[0]+"_"+str(counter)
				counter+= 1
			else:
				modgenename= gene[0]
				counter= 1
			output[modgenename]= []
			for element in gene[1:]:
				output[modgenename].append(element)
			previousgene= gene[0]
		return output

	# Write dictionary to csv file, output the first item of the list
	def writedicttocsv2(self,indict,outfilestring):
		import csv
		fout= open(outfilestring,"w")
		writer= csv.writer(fout)
		for key,value in indict.items():
			writer.writerow([key,value[0]])


	def Densebuilder(self,assignment):
		rootpath= '%s/%s' %(self.datapath, self.exp)
		mismatch= str(self.mismatch).replace('.','')
		ftsize= self.ftsize.strip("'").split(",")
		bamfileinput= '%s/%s_starM%s/%s_match.bam' %(rootpath,self.exp,mismatch,self.exp)
		if self.normalization== '1':	totreads= '-1'
		else:	totreads= '1E6'      # This gives raw reads mapped to each nucleotide.

		def DenseMultitask(parameters):
			size, assignment= parameters[0],parameters[1]
			bamfileoutput= '%s/%s_starM%s/%s_%s_match.bam' %(rootpath,self.exp,mismatch,self.exp,size)
			densityfilepath= '%s/DensityUnnormalized/density%sM%s_%s' %(rootpath,assignment,mismatch,size)
			densityfileout= '%s/%s_%sf_' %(densityfilepath,self.exp,size)
			if not os.path.exists(densityfilepath): os.makedirs(densityfilepath)
			riboshiftdict= "{"+str(size)+':[0]'+"}"
			#commandstring= 'python densebuilder_main.py --bamfileinput %s --GTFfile %s --twobitfile %s --assignment %s --riboshiftdict %s --threshold %s --totreads %s --outputdata %s --bamfileoutput %s' %(bamfileinput,self.GTFfile,self.twobitfile,assignment,riboshiftdict,self.threshold,totreads,densityfileout,bamfileoutput)
			commandstring= 'densebuilder_main.py --bamfileinput %s --GTFfile %s --twobitfile %s --assignment %s --riboshiftdict %s --threshold %s --totreads %s --outputdata %s --bamfileoutput %s' %(bamfileinput,self.GTFfile,self.twobitfile,assignment,riboshiftdict,self.threshold,totreads,densityfileout,bamfileoutput)
			print(commandstring)
			os.system(commandstring)

		def NormDenseMultitask(parameters):
			size, assignment= parameters[0],parameters[1]
			bamfileoutput= '%s/%s_starM%s/%s_%s_match.bam' %(rootpath,self.exp,mismatch,self.exp,size)
			densityfilepath= '%s/DensityNormalized/density%sM%s_%s' %(rootpath,assignment,mismatch,size)
			densityfileout= '%s/%s_%sf_' %(densityfilepath,self.exp,size)
			if not os.path.exists(densityfilepath): os.makedirs(densityfilepath)
			riboshiftdict= "{"+str(size)+':[0]'+"}"
			#commandstring= 'python densebuilder_main.py --bamfileinput %s --GTFfile %s --twobitfile %s --assignment %s --riboshiftdict %s --threshold %s --totreads %s --outputdata %s --bamfileoutput %s' %(bamfileinput,self.GTFfile,self.twobitfile,assignment,riboshiftdict,self.threshold,totreads,densityfileout,bamfileoutput)
			commandstring= 'densebuilder_main.py --bamfileinput %s --GTFfile %s --twobitfile %s --assignment %s --riboshiftdict %s --threshold %s --totreads %s --outputdata %s --bamfileoutput %s' %(bamfileinput,self.GTFfile,self.twobitfile,assignment,riboshiftdict,self.threshold,totreads,densityfileout,bamfileoutput)
			print(commandstring)
			os.system(commandstring)


		### TRY TO FIX THIS - MULTITHREADING SEEMS NOT TO BE WORKING ###
		params_in= [[x,str(assignment)] for x in ftsize]
		pool= Pool(len(ftsize))
		if self.normalization== '1':	pool.map(NormDenseMultitask,params_in)
		else:	pool.map(DenseMultitask,params_in)
		################################################################

	print('DenseBuilder Success')

	def Makeavggene(self,assignment,alignpos):
		ftsize= self.ftsize.strip("'").split(",")
		rootpath= '%s/%s' %(self.datapath, self.exp)
		mismatch= str(self.mismatch).replace('.','')
		thresh= str(0)  # rpkm threshold   ### currently still hard-coded. 

		def Avggene(parameters):
			size,assignment,alignpos= parameters[0],parameters[1],parameters[2]
			trspdictfilestring= '%s/DensityUnnormalized/density%sM%s_%s/%s_%sf_' %(rootpath,assignment,mismatch,size,self.exp,size)
			outfilebase= '%s/%s_%sf'%(outfolder,self.exp,size)
			#commandstring= 'python makeavggene_main.py --regionlength5 %s --regionlength3 %s --trspdictfilestring %s --UTRfilestring %s --filtermodule %s --exclusionmodule %s --threshold %s --alignpos %s --equalweight %s --outfilebase %s' %(self.regionlength5,self.regionlength3,trspdictfilestring,self.UTRfilestring,self.filtermodule,self.exclusionmodule,thresh,alignpos,self.equalweight,outfilebase)
			commandstring= 'makeavggene_main.py --regionlength5 %s --regionlength3 %s --trspdictfilestring %s --UTRfilestring %s --filtermodule %s --exclusionmodule %s --threshold %s --alignpos %s --equalweight %s --outfilebase %s' %(self.regionlength5,self.regionlength3,trspdictfilestring,self.UTRfilestring,self.filtermodule,self.exclusionmodule,thresh,alignpos,self.equalweight,outfilebase)
			print(commandstring)
			os.system(commandstring)

		def NormAvggene(parameters):
			size,assignment,alignpos= parameters[0],parameters[1],parameters[2]
			trspdictfilestring= '%s/DensityNormalized/density%sM%s_%s/%s_%sf_' %(rootpath,assignment,mismatch,size,self.exp,size)
			outfilebase= '%s/%s_%sf'%(outfolder,self.exp,size)
			#commandstring= 'python makeavggene_main.py --regionlength5 %s --regionlength3 %s --trspdictfilestring %s --UTRfilestring %s --filtermodule %s --exclusionmodule %s --threshold %s --alignpos %s --equalweight %s --outfilebase %s' %(self.regionlength5,self.regionlength3,trspdictfilestring,self.UTRfilestring,self.filtermodule,self.exclusionmodule,thresh,alignpos,self.equalweight,outfilebase)
			commandstring= 'makeavggene_main.py --regionlength5 %s --regionlength3 %s --trspdictfilestring %s --UTRfilestring %s --filtermodule %s --exclusionmodule %s --threshold %s --alignpos %s --equalweight %s --outfilebase %s' %(self.regionlength5,self.regionlength3,trspdictfilestring,self.UTRfilestring,self.filtermodule,self.exclusionmodule,thresh,alignpos,self.equalweight,outfilebase)
			print(commandstring)
			os.system(commandstring)

		params_in= [[x,str(assignment),str(alignpos)] for x in ftsize]
		pool= Pool(len(ftsize))
		if self.normalization== '1':	
			outfolder= '%s/DensityNormalized/avggene%s_ORF_0shift_%s' %(rootpath,alignpos,assignment)
			if not os.path.exists(outfolder):   os.makedirs(outfolder)
			pool.map(NormAvggene,params_in)
		else:	
			outfolder= '%s/DensityUnnormalized/avggene%s_ORFraw_0shift_%s' %(rootpath,alignpos,assignment)
			if not os.path.exists(outfolder):   os.makedirs(outfolder)
			pool.map(Avggene,params_in)


	# Plot heatmap using R
	def plot_output_START(self,assignment,alignpos,normalization):
		Rcode= '%s/START_heatmap%s_40S_local.R' %(self.Rpath,str(assignment))
		rootpath= '%s/%s' %(self.datapath, self.exp)
		if not os.path.exists(self.outpath):   os.makedirs(self.outpath)
		if self.normalization== '1':
			wd= '%s/DensityNormalized/avggene%s_ORF_0shift_%s' %(rootpath,alignpos,assignment)
			outfile= '%s/%s_START_heatmap_normalized%s.pdf' %(self.outpath, self.exp,assignment)
		else:
			wd= '%s/DensityUnnormalized/avggene%s_ORFraw_0shift_%s' %(rootpath,alignpos,assignment)
			outfile= '%s/%s_START_heatmap_rawreads%s.pdf' %(self.outpath, self.exp,assignment)
		ftmin,ftmax= self.ftsize.strip("'").split(",")[0],self.ftsize.strip("'").split(",")[-1]

		input= 'Rscript %s %s %s %s %s %s %s %s' %(Rcode,self.exp,wd,outfile,self.regionlength5,self.regionlength3,ftmin,ftmax)
		print(str(Rcode), str(self.exp), str(wd))
		print(str(input))
		subprocess.call(input, shell= True)


def parse_setting(setting_file):
	filein= open(setting_file)
	params= {}
	for line in filein:
		#line= line.strip()
		key_value= line.split('=')
		if len(key_value)== 2:	params[key_value[0].strip()]= key_value[1].strip()

	datapath= params['datapath']
	fastq_filein= params['fastq_filein']
	exp_name= params['exp_name']
	fastq_rearrangement= params['fastq_rearrangement']
	fastq_filein_list= list(fastq_filein.split(','))
	exp_name_list= list(exp_name.split(','))

	if fastq_rearrangement== '1':
		for i in range(len(exp_name_list)):
			folderoutpath= '%s/%s' %(datapath,exp_name_list[i].strip("'"))
			if fastq_filein_list[i]!= exp_name_list[i] and not os.path.exists(folderoutpath): 
				fastq_gz_in= '%s/%s.fastq.gz' %(datapath,fastq_filein_list[i].strip("'"))
				fastq_gz_exp= '%s/%s.fastq.gz' %(datapath,exp_name_list[i].strip("'"))
				try:	os.rename(fastq_gz_in, fastq_gz_exp)
				except	IOError:	sys.exit("No such a file... ")

			folderoutpath= '%s/%s' %(datapath,exp_name_list[i].strip("'"))
			if not os.path.exists(folderoutpath):
				os.makedirs(folderoutpath)
				shutil.move(datapath+'/%s.fastq.gz' %(exp_name_list[i].strip("'")),folderoutpath+'/%s.fastq.gz' %(exp_name_list[i].strip("'")))
	return params


def main():
	params= parse_setting(sys.argv[1])

	datapath= params['datapath']
	fastq_filein= params['fastq_filein']
	exp_name= params['exp_name']
	ncrnaDir= params['ncrnaDir']
	genomeDir= params['genomeDir']
	GTFfile= params['GTFfile']
	twobitfile= params['twobitfile']
	UTRfilestring= params['UTRfilestring']
	mismatch= params['mismatch']
	fiveprimetrim= params['fiveprimetrim']
	linker= params['linker']
	threshold= params['threshold']
	filtermodule= params['filtermodule']
	exclustionmodule= params['exclustionmodule']
	equalweight= params['equalweight']
	ftsize= params['ftsize']
	regionlength5= params['regionlength5']
	regionlength3= params['regionlength3']
	fastq_filein_list= list(fastq_filein.split(','))
	exp_name_list= list(exp_name.split(','))
	normalization= params['normalization']
	Rpath= params['Rpath']
	sizedistr_outpath= params['sizedistr_outpath']
	outpath= params['outpath']


	for exp in exp_name_list:
		Run= SubunitProfiling(datapath, exp.strip("'"), ncrnaDir, genomeDir, GTFfile, twobitfile, mismatch, fiveprimetrim, linker, threshold, filtermodule, exclustionmodule, equalweight, UTRfilestring, ftsize, regionlength5, regionlength3, normalization, Rpath, sizedistr_outpath, outpath)
		
		Run.StarAlignment()
		Run.SizeDistr()
		Run.Densebuilder(5)
		Run.Makeavggene(5,1)
		Run.plot_output_START(5,1,normalization)


if __name__ == '__main__':
	main()
