'''
Reads a candidate somatic vcf file
and applies high confidence, min AF, and haplotype discodant filter

python haplosort.py input.vcf phased.bam highconf.bed tumour_id --filters HD
'''
import argparse
import math
import subprocess
import re
import pysam
import vcf
import os

def main():
    parser = argparse.ArgumentParser(description='Rewrite vcf files with haplodiscordant variants and haploconcordant variants')
    parser.add_argument('vcffile', help='Name of vcf file with candidate somatic mutations')
    parser.add_argument('bamfile', help='Name of phased bam file')
    parser.add_argument('highconffile', help='Name of high confidence bed file')
    parser.add_argument('tumour_id', help='Tumour ID from BAM')
    parser.add_argument('--filters', nargs='+', required=True, help="Choose either high confidence (HC), haplotype discordant (HD) or allele frequency (AF)")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s beta')
    args = parser.parse_args()

    # get filters
    avail_filters = ["HC", "HD", "AF"]
    filters = []
    for f in args.filters:
        if f not in avail_filters:
            print("Choose HC HD and/or AF")
            exit(0)
        else:
            filters.append(f)

    # get old vcf file name
    print("[haplosort] reading candidate vcf file " + args.vcffile)
    basename = os.path.basename(args.vcffile).strip(".gz").strip(".vcf")
    newname = basename +  "." + "_".join(filters) + ".clean.vcf"
    removed_filename =  basename +  "." + "_".join(filters) + ".removed.vcf"
    print("[haplosort] storing all passing variants in vcf: " + newname)
    print("[haplosort] storing all removed variants in vcf: " + removed_filename)

    # read vcffiles for reading and writing
    vcf_reader = vcf.Reader(open(args.vcffile, 'r'))
    vcf_writer = vcf.Writer(open(newname, 'w'), vcf_reader)
    vcf_removed_writer = vcf.Writer(open(removed_filename, 'w'), vcf_reader)

    # get high confidence file
    high_conf_intervals = {}
    if "HC" in filters:
        file = open(args.highconffile, 'r')
        for line in file:
            info = line.strip('\n').split()
            if info[0] not in high_conf_intervals.keys():
                high_conf_intervals[info[0]] = list()
            high_conf_intervals[info[0]].append((int(info[1]), int(info[2])))

    # open phased BAM file
    print("[haplosort] reading phased BAM file " + args.bamfile)
    open_bamfile = pysam.AlignmentFile(args.bamfile, "rb")

    # read vcffile with candidate somatic mutations
    print("[haplosort] filtering ...")
    num_HC_filtered = 0
    num_HD_filtered = 0
    num_minAF_filtered = 0
    for rec in vcf_reader:
        # filters
        # only filter SNVs
        if len(rec.ALT[0]) > 1 or len(rec.REF) > 1:
            continue
        try:
            call = rec.genotype(args.tumour_id)
            gt = call.data.GT
            af = -1.0
            try:
                af = float(call.data.AD[1])/(float(call.data.AD[1]) + float(call.data.AD[0]))
            except:
                af = float(call.data.AF)
        except:
            try:
                af = rec.INFO['AF'][0]
            except:
                try:
                    af = float(rec.INFO['AC'][1])/rec.INFO['DP']
                except:
                    af = rec.INFO['AD'][1]/(rec.INFO['AD'][1] + rec.INFO['AD'][0])
        assert af >= 0.0
        if "AF" in filters and af < 0.06:
            num_minAF_filtered += 1
            vcf_removed_writer.write_record(rec)
            continue
        if "HC" in filters and not is_in_highconfregion(rec.CHROM, rec.POS, high_conf_intervals):
            num_HC_filtered += 1
            vcf_removed_writer.write_record(rec)
            continue
        if "HD" in filters and (gt == "0/1" or gt == "0|1") and is_haplo_discordant(open_bamfile, rec.CHROM, rec.POS, rec.ALT[0], rec.REF):
            num_HD_filtered += 1
            vcf_removed_writer.write_record(rec)
            continue
        vcf_writer.write_record(rec)
    print("[haplosort] number filtered = {}".format(num_minAF_filtered + num_HD_filtered + num_HC_filtered)) 
    print("[haplosort] done")

def is_haplo_discordant(open_bamfile, chr, pos, alt, ref):
    '''Return True if the variant at chr:pos does not
    segregate properly into haplotypes
    '''
    # stores counts
    for pileupcolumn in open_bamfile.pileup(chr, pos-1, pos+1, stepper='all'):
        h1a = 0
        h1r = 0
        h2a = 0
        h2r = 0
        ma = 0
        mr = 0
        nrnab = 0
        if (pileupcolumn.pos + 1) == pos :
            bases = ''
            # go through all reads at specific position
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    query_len = pileupread.alignment.query_length
                    base = pileupread.alignment.query_sequence[pileupread.query_position]
                    hp = -9 if not pileupread.alignment.has_tag("HP") else int(pileupread.alignment.get_tag("HP"))
                    bases += base
                    if hp == 1 and (base == alt):
                        h1a += 1
                    elif hp == 1 and (base == ref):
                        h1r += 1
                    elif hp == 2 and (base == alt):
                        h2a += 1
                    elif hp == 2 and (base == ref):
                        h2r += 1
                    elif hp == -9 and (base == alt):
                        ma += 1
                    elif hp == -9 and (base == ref):
                        mr += 1
                    else:
                        nrnab += 1
                    if (h1a > 0 and h2a > 0):
                        return True
    return False


def is_in_highconfregion(chr, pos, high_conf_intervals):
    '''Return True if chr:pos is in a high conf region
    defined by high_conf_intervals list
    '''
    if chr not in high_conf_intervals:
        return True
    for (lower, upper) in high_conf_intervals[chr]:
        if lower <= pos <= upper:
            return True
    return False

if __name__ == '__main__':
	main()
        

