# Test_2
Test_2 for job application


## Requiered software
python >3

tabix

  python packages
  
    sys, gzip, pysam, re, contextlib ,argparse, numpy 

## Requiered files
    sorted bgzip tabix indexed GTF file**
    indexed VCF file
**regular GTF file can be prepared using :

    (grep ^"#" in.gtf; grep -v ^"#" in.gtf | sort -k1,1 -k4,4n) | bgzip  > in.sorted.gtf.gz ; tabix -gff in.sorted.gtf.gz
   
## Usage

python script
  -i/--in-vcf VCF input 
  -c/--in-GTF GTF input 
  -o/--out-vcf VCF output
