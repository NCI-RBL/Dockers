#!/usr/bin/env python

__author__ = 'Colin Wu'
__modified_by__ = 'Wilfried Guiblet'

import os, sys

folder, file, genomeDir, mismatch= sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

rootpath= '/%s/%s' % (folder,file)
#genomeDir= '/index/hg19/star_hg19'
# mismatch        #outFilterMismatchNoverLmax

############################################################################################################################################

starInput1= '%s/%s_R1_001.fastq.gz' %(rootpath,file)
starInput2= '%s/%s_R2_001.fastq.gz' %(rootpath,file)


# Align to genome
outfolder= '%s/%s_starM%s/' %(rootpath,file,str(mismatch).replace('.',''))
if not os.path.exists(outfolder):  os.makedirs(outfolder)
outfile= '%sAligned.sortedByCoord.out.bam' %(outfolder)
newoutfile= '%s%s_match.bam' %(outfolder,file)
outfileindex= '%s%s_match.bai' %(outfolder,file)
unmapped1= '%s/%s_unmapped1.fastq.gz' %(outfolder,file)
unmapped2= '%s/%s_unmapped2.fastq.gz' %(outfolder,file)

if not os.path.exists(outfileindex):
    CMMD='STAR --runThreadN 45 --readFilesIn %s %s --readFilesCommand gunzip -c --limitBAMsortRAM 6000000000 --alignIntronMax 1000 --outSJfilterIntronMaxVsReadN 1000 1000 1000 --genomeDir %s --outFilterType BySJout --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outFilterMultimapNmax 1 --outFilterMismatchNoverLmax %s --outWigType wiggle read1_5p --outWigNorm RPM --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --alignEndsType Local --outFileNamePrefix %s' %(starInput1,starInput2,genomeDir,mismatch,outfolder)
    os.system(CMMD)

    os.rename('%sSignal.UniqueMultiple.str1.out.wig' %(outfolder), '%s%s_Signal.UniqueMultiple.str1.out.wig' %(outfolder,file))
    os.rename('%sSignal.UniqueMultiple.str2.out.wig'%(outfolder), '%s%s_Signal.UniqueMultiple.str2.out.wig' %(outfolder,file))
    os.rename('%sSignal.Unique.str1.out.wig' %(outfolder), '%s%s_Signal.Unique.str1.out.wig' %(outfolder,file))
    os.rename('%sSignal.Unique.str2.out.wig' %(outfolder), '%s%s_Signal.Unique.str2.out.wig' %(outfolder,file))
    
    if os.path.exists(outfile):
        os.rename(outfile,newoutfile)
        CMMD= 'samtools index %s %s' %(newoutfile,outfileindex)
        os.system(CMMD)

    if os.path.exists(unmapped1):
        command= 'pigz -c %s/Unmapped.out.mate1 > %s' %(outfolder,unmapped1)
        print(command)
        os.system(command)

        CMMD= 'rm %s/Unmapped.out.mate1' % outfolder
        os.system(CMMD)
        
    if os.path.exists(unmapped2):
        command= 'pigz -c %s/Unmapped.out.mate2 > %s' %(outfolder,unmapped2)
        print(command)
        os.system(command)

        CMMD= 'rm %s/Unmapped.out.mate2' % outfolder
        os.system(CMMD)

