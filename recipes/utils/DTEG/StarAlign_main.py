#!/usr/bin/env python

import os, subprocess, sys
import argparse

class Alignment(object):
    def __init__(self, rootpath, file, ncrnaDir, genomeDir, mismatch, fiveprimetrim, linker):
        self.rootpath= rootpath
        self.file= file
        self.ncrnaDir= ncrnaDir
        self.genomeDir= genomeDir
        self.mismatch= mismatch
        self.fiveprimetrim= fiveprimetrim
        self.linker= linker


    def StarAlign(self):
        inputfile= '%s/%s.fastq.gz' %(self.rootpath,self.file)
        tempfolder= '%s/temp' %(self.rootpath)
        if not os.path.exists(tempfolder):    os.makedirs(tempfolder)
        ncrna_outpath= '%s/ncrna_out' %(self.rootpath)
        if not os.path.exists(ncrna_outpath):    os.makedirs(ncrna_outpath)
        ncrna_outfile= '%s/%s_Unmapped.fastq.gz' %(ncrna_outpath,self.file)

        #if not os.path.exists(ncrna_outfile):
            # Trim 5'end of reads
        trimmedfastq= '%s/%s_5trimmed.fastq' %(tempfolder,self.file)
        starInput= '%s/%s_5trimmed-trimmed.fastq' %(tempfolder,self.file)
        CMMD1= 'seqtk trimfq -b %s %s > %s' %(self.fiveprimetrim,inputfile,trimmedfastq)
        #if not os.path.exists(trimmedfastq):
        subprocess.Popen(CMMD1, shell=True).wait()
        
        # Remove adaptor
        #if not os.path.exists(starInput):
        CMMD2= 'skewer -t 26 -x %s -l 15 %s' %(self.linker,trimmedfastq)
        subprocess.Popen(CMMD2, shell=True).wait()

        # Remove ncrna
        #if os.path.exists(starInput):
        CMMD3= 'STAR --runThreadN 26 --readFilesIn %s --genomeDir %s --outReadsUnmapped Fastx --outSAMtype None --outFilterMismatchNoverLmax 0.3 --outFileNamePrefix %s/' %(starInput,self.ncrnaDir,ncrna_outpath)
        subprocess.Popen(CMMD3, shell=True).wait()
            
        compression= 'pigz -p 26 -c %s/Unmapped.out.mate1 > %s && rm %s/Unmapped.out.mate1' %(ncrna_outpath,ncrna_outfile,ncrna_outpath)
        subprocess.Popen(compression, shell=True).wait()

        ### --outFilterMismatchNmax 2
        ### --outFilterMismatchNoverLmax

        outfolder= '%s/%s_starM%s/' %(self.rootpath,self.file,str(self.mismatch).replace('.',''))
        if not os.path.exists(outfolder):  os.makedirs(outfolder)
        outfile= '%sAligned.sortedByCoord.out.bam' %(outfolder)
        newoutfile= '%s%s_match.bam' %(outfolder,self.file)
        outfileindex= '%s%s_match.bai' %(outfolder,self.file)

        # Genomic alignment
        CMMD4='STAR --runThreadN 26 --readFilesIn %s --readFilesCommand gunzip -c --limitBAMsortRAM 600000000000 --genomeDir %s --outFilterType BySJout --outFilterIntronMotifs RemoveNoncanonicalUnannotated --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --outFilterMultimapNmax 1 --outFilterMismatchNoverLmax %s --outWigType wiggle read1_5p --outWigNorm RPM --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --outFileNamePrefix %s' %(ncrna_outfile,self.genomeDir,self.mismatch,outfolder)
        if os.path.exists(newoutfile)== False or os.path.exists(outfileindex)== False:
            subprocess.Popen(CMMD4, shell=True).wait()
            if os.path.exists(outfile):    os.rename(outfile,newoutfile)
            CMMD= 'samtools index %s %s' %(newoutfile,outfileindex)
            subprocess.Popen(CMMD, shell=True).wait()
            os.rename('%sSignal.UniqueMultiple.str1.out.wig' %(outfolder), '%s%s_Signal.UniqueMultiple.str1.out.wig' %(outfolder,self.file))
            os.rename('%sSignal.UniqueMultiple.str2.out.wig'%(outfolder), '%s%s_Signal.UniqueMultiple.str2.out.wig' %(outfolder,self.file))
            os.rename('%sSignal.Unique.str1.out.wig' %(outfolder), '%s%s_Signal.Unique.str1.out.wig' %(outfolder,self.file))
            os.rename('%sSignal.Unique.str2.out.wig' %(outfolder), '%s%s_Signal.Unique.str2.out.wig' %(outfolder,self.file))

        if os.path.exists('%sAligned.toTranscriptome.out.bam' %(outfolder))== True: 
            outtrxfile= '%s%s_Transcriptome_match.bam' %(outfolder,self.file)
            outtrxfileindex= '%s%s_Transcriptome_match.bai' %(outfolder,self.file)
            sort_transcriptome= 'samtools sort -@ 10 %sAligned.toTranscriptome.out.bam -o %s' %(outfolder,outtrxfile)
            subprocess.Popen(sort_transcriptome, shell=True).wait()

            index_transcriptome= 'samtools index %s %s' %(outtrxfile,outtrxfileindex)
            subprocess.Popen(index_transcriptome, shell=True).wait()

            if os.path.exists(outtrxfile):  subprocess.Popen('rm %sAligned.toTranscriptome.out.bam' %(outfolder), shell = True).wait()

        f= open('%sStarAlignment_log.txt' %(outfolder),'w')
        f.write("\nSTAR alignment was done with following parameters:\n")
        try:    f.write(CMMD1+"\n")
        except NameError:   pass
        try:    f.write(CMMD2+"\n")
        except NameError:   pass
        try:    f.write(CMMD3+"\n")
        except NameError:   pass
        try:    f.write(CMMD4+"\n")
        except NameError:   pass
        f.close()



if __name__== '__main__':
    parser= argparse.ArgumentParser()
    parser.add_argument('--rootpath', help= 'working directory')
    parser.add_argument('--file', help= 'input dataset')
    parser.add_argument('--ncrnaDir', default= '/mnt/gridftp/wuc11/index/rna/star_hg19_ncrna', help= 'index for ncRNA depletion')
    parser.add_argument('--genomeDir', default= '/mnt/gridftp/wuc11/index/hg19/star_hg19', help= 'index for genomic alignment')
    parser.add_argument('--mismatch', default= 0.1, help= 'percent of mismatch allowed in STAR')
    parser.add_argument('--fiveprimetrim', default= 4, help= 'bases to trim off from the 5-end')
    parser.add_argument('--linker', default= 'NNNNNNCACTCGGGCACCAAGGA', help= 'linker sequence to remove from 3-end')
    parser.add_argument
    args = parser.parse_args()

    RunAlignment= Alignment(args.rootpath,args.file,args.ncrnaDir,args.genomeDir,args.mismatch,args.fiveprimetrim,args.linker)
    RunAlignment.StarAlign()

