#!/usr/bin/env python
import re
import sys
import pysam
import regex
import parasail
import statistics
from itertools import islice

def revComp(s):
        revCDict = {'G':'C','C':'G','T':'A','A':'T','N':'N','-':'-'}
        return ''.join([revCDict[x] for x in s])[::-1]

def removeAdapterEnd(adapter, seq, qual):
        matches = regex.finditer(adapter, seq, overlapped = False)
        for it in matches:
                end = it.end()
        try:
                seq = seq[end:]
                qual = qual[end:]	
        except:
                seq = seq
                qual = qual
        return seq, qual

def removeAdapterStart(adapter, seq, qual):
        seq = revComp(seq)
        qual = qual[::-1]
        seq, qual = removeAdapterEnd(adapter, seq, qual)
        seq = revComp(seq)
        qual = qual[::-1]
        return seq, qual

def findSpikein(ccs_name, ccs_pass, barcode_dict, strand, seq, qual):
        spikein = []
        all_bc_loci = [(0, len(seq) + 10)]

        for barcode_name, barcode_seq in barcode_dict.items():
                bc_loci = []
                pattern = "("+ barcode_seq + "){e<=2}(ATGGTGAGCAAGGGCG){e<=2}"
                matches = regex.finditer(pattern, seq, overlapped = False)
                for it in matches:
                        bc_loci.append((barcode_name,it.end()))
                all_bc_loci.extend(bc_loci)
        sort_loci = sorted(all_bc_loci, key=lambda all_bc_loci: all_bc_loci[1])

        for i in range(len(sort_loci)-1):
                barcode_name = sort_loci[i][0]
                barcode_seq = barcode_dict[barcode_name]
                read_name_serial = ccs_name + '_' + str(i+1)
                read_id = read_name_serial + ':' + barcode_name + ':' + ccs_pass
                target_seq = seq[sort_loci[i][1] - 16:sort_loci[i+1][1]]
                target_qual = qual[sort_loci[i][1] - 16:sort_loci[i+1][1]]
                read = (ccs_name, read_id, barcode_name, strand, target_seq, target_qual)
                spikein.append(read)
        return spikein

infile = sys.argv[1]
passfile = sys.argv[2]
barcodefile = sys.argv[3]

if len(sys.argv) != 4:
        print('\tneed: infile, passfile, barcodefile, pass_thres, qual_thres')
        print('\tPlease contact WJQ if you have any questions')
        quit()

barcode_dict = {}
barcodefile = pysam.FastxFile(barcodefile)
for read in islice(barcodefile, None):
        barcode_dict[read.name] = read.sequence

pass_dict = {}
passfile = open(passfile,'r')
for line in passfile:
        (ccs_name, ccs_pass) = line.strip('\n').split('\t')
        pass_dict[ccs_name] = ccs_pass

spikeins = []
ccs_count = 0
infile = pysam.FastxFile(infile)
for ccs in islice(infile, None):
        ccs_count = ccs_count + 1
        ccs_name = ccs.name
        ccs_pass = pass_dict[ccs_name]

        seq = ccs.sequence
        qual = ccs.quality
        strand = '+'
        spikein = findSpikein(ccs_name, ccs_pass, barcode_dict, strand, seq, qual)

        seq = revComp(ccs.sequence)
        qual = ccs.quality[::-1]
        strand = '-'
        spikein.extend(findSpikein(ccs_name, ccs_pass, barcode_dict, strand, seq, qual))
        spikeins.extend(spikein)

for read in spikeins:
        end = ''
        (read_name, read_id, spikein_barcode, strand, seq, qual) = read
        (read_name_serial, barcode_name, ccs_pass) = read_id.split(':')
        ccs_pass = int(ccs_pass)        

        TGA = '(GACGAGCTGTACAAGTGA){e<=3}'
        matches = regex.finditer(TGA, seq, overlapped = False)
        for it in matches:
                end = it.end()
                break
        if  end < 1000:
                tso5 = '(AAGCAGTGGTATCAACGCAGAGTAC){e<=3}'
                seq1, qual = removeAdapterStart(tso5, seq, qual)
                if len(seq1) != len(seq):
                        seq = seq1
                        tag = 'expected'
                else:
                        tag = 'NO_suff3'
        else:
                tag = 'NO_suff5'

        if tag == 'expected':
                out = [read_name, read_id, spikein_barcode, tag, strand, seq, qual]
                out = '\t'.join([str(x) for x in out])
                print(out)

        else:
                err = [read_name, read_id, spikein_barcode, tag, strand, seq, qual]
                err = '\t'.join([str(x) for x in err])
                sys.stderr.write(err + '\n')

logfile = sys.argv[1] + '.log.txt'
log_file = open(logfile, 'w')
log_file.write(str(ccs_count) + '\n')
log_file.close()
