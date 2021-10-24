import sys
import gzip
import pysam
import re
from contextlib import closing
import argparse
from numpy import inf ##for the start dist 

#arguments
argparser = argparse.ArgumentParser(description = 'VCF annotation of gene proximity information')
argparser.add_argument('-i', '--in-vcf', metavar = 'file', dest = 'in_VCF', required = True, help = 'Input VCF/BCF file. File may be compressed with gzip/bgzip.')
argparser.add_argument('-c', '--in-GTF', metavar = 'file', dest = 'in_GTF_file', required = True, help = 'Input GTF sorted, compressed bgzip and indexed with tabix')
argparser.add_argument('-o', '--out-vcf', metavar = 'file', dest = 'out_VCF', required = True, help = 'Output VCF file. File will be compressed with bgzip.')


#function
def addGTF_info(in_VCF, in_GTF_file, out_VCF):
  tabixfile= pysam.TabixFile(in_GTF_file) 
  with pysam.VariantFile(in_VCF, 'r') as ifile, pysam.BGZFile(out_VCF, 'w') as ofile:
##Header addition
    for x in ['GENES_IN', 'GENES_200KB','GENE_NEAREST']:                 #exception if already present
      if x in ifile.header.info:                                         #
        raise Exception('{} already exists in input VCF/BCF.'.format(x)) #
    ifile.header.add_line('##INFO=<ID=GENES_IN,Number=.,Type=String,Description=" comma-separated list of identifiers (gene_id key from GENCODE GTF) of overlapping genes. If there are no overlapping genes, the empty value must be denoted with .">')
    ifile.header.add_line('##INFO=<ID=GENES_200KB,Number=.,Type=String,Desctiption="comma-separated list of identifiers (gene_id key from GENCODE GTF) of genes which are within +/-200000 base pairs from the variant. If there are no such genes, the empty value must be denoted with .">')
    ifile.header.add_line('##INFO=<ID=GENE_NEAREST,Number=.,Type=String,Desctiption="identifier (gene_id key from GENCODE GTF) of the nearest gene">')
    ofile.write('{}'.format(ifile.header).encode()) #writing
    last=[] #kept in memory for nearest determination optimization

    for record in ifile :
      genes_in = []
      genes_200kb = []
      genes_nearest = []
      if record.pos-200000<0 : #start exceptions -200000kb 
        for gene in tabixfile.fetch(record.chrom,0,parser=pysam.asGTF()):
          genes_200kb.append(gene.gene_id)
          if not last:
            last=gene
          if min(abs(gene.start-record.pos),abs(gene.end-record.pos)) < min(abs(last.end-record.pos),abs(last.start-record.pos)):
            genes_nearest=gene.gene_id
            last=gene
          if min(abs(gene.start-record.pos),abs(gene.end-record.pos))>200000: #200kb cap if 200kb
            break 
      else :                   #regular +/-200kb skip beguining
        for gene in tabixfile.fetch(record.chrom,record.pos-200000, parser=pysam.asGTF()):
          if min(abs(gene.start-record.pos),abs(gene.end-record.pos))>200000: #200kb cap
            break 
          else:
            genes_200kb.append(gene.gene_id)
            if not last:
              last=gene
            if min(abs(gene.start-record.pos),abs(gene.end-record.pos)) < min(abs(last.end-record.pos),abs(last.start-record.pos)): # nearest gene
              genes_nearest=gene.gene_id                                                                                            #
              last=gene                                                                                                             #
      for gene in tabixfile.fetch(record.chrom,record.pos, record.pos+len(record.ref),parser=pysam.asGTF()):  #gene IN
        genes_in.append(gene.gene_id)                                                                         #
        genes_nearest="."                                                                                     #
        last=gene                                                                                             #
      if not last:                                                                                                                 #nearest gene no 200kb no gene_in
        dist=inf   
      if not genes_200kb and not genes_in:
        if not last: 
          for gene in tabixfile.fetch(record.chrom,0,parser=pysam.asGTF()):                                                          #
            if min(abs(gene.start-record.pos),abs(gene.end-record.pos)) < dist:                                                      #
              dist=min(abs(gene.start-record.pos),abs(gene.end-record.pos))                                                          #
              genes_nearest=gene.gene_id                                                                                             #
            else:                                                                                                                    #
              last=gene                                                                                                              #
              break 
         if last :
            for gene in tabixfile.fetch(record.chrom, max(last.end,last.start),parser=pysam.asGTF()):                                                          #
             if min(abs(gene.start-record.pos),abs(gene.end-record.pos)) < min(abs(last.end-record.pos),abs(last.start-record.pos)):                                                      #
               dist=min(abs(gene.start-record.pos),abs(gene.end-record.pos))                                                          #
               genes_nearest=gene.gene_id                                                                                             #
             else:                                                                                                                    #
               last=gene                                                                                                              #
               break 
      record.info['GENES_IN'] = ",".join(set(genes_in))
      if not genes_200kb:
        genes_200kb.append(".")
      record.info['GENES_200KB'] = ",".join(set(genes_200kb))
      record.info['GENE_NEAREST'] = genes_nearest
      ofile.write('{}'.format(record).encode())
if __name__ == "__main__":
  args = argparser.parse_args()
  addGTF_info(args.in_VCF, args.in_GTF_file, args.out_VCF)
