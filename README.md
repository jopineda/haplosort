# haplosort

A script to filter out haplotype discordant variants given a candidate VCF file and a phased BAM file. The primary purpose is to remove candidate somatic variants that do not segregate by haplotype.

# quickstart

```
python haplosort.py input.vcf phased.bam highconf.bed tumour_id --filters HD
```

# Dependencies

* python3
* pyvcf
* pysam

# Installation

```
virtualenv venv -p [/path/to/python3]
source venv activate
git clone https://github.com/jopineda/haplosort.git
pip install -r requirements.txt
```

