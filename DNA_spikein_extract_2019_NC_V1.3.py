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

def toDigit(ascii_list):
        digit_list = [(ord(i)-33) for i in ascii_list]
        return digit_list

def findSpikein(ccs_name, ccs_pass, barcode_dict, seq, qual):
        spikein = []
        all_bc_loci = [(0, len(seq) + 10)]

        for barcode_name, barcode_seq in barcode_dict.items():
                bc_loci = []
                pattern = "("+ barcode_seq + ")(ATGGTGAGCAAGGGCG){e<=2}"
                matches = regex.finditer(pattern, seq, overlapped = False)
                for it in matches:
                        bc_loci.append((barcode_name,it.end()))
                all_bc_loci.extend(bc_loci)
        sort_loci = sorted(all_bc_loci, key=lambda all_bc_loci: all_bc_loci[1])

        for i in range(len(sort_loci)-1):
                barcode_name = sort_loci[i][0]
                barcode_seq = barcode_dict[barcode_name]
                read_name_serial = ccs_name + '_' + str(i+1)
                read_id = read_name_serial + ':' + ccs_pass
                target_seq = seq[sort_loci[i][1] - 16:sort_loci[i+1][1]]
                target_qual = qual[sort_loci[i][1] - 16:sort_loci[i+1][1]]
                read = (ccs_name, read_id, barcode_name, target_seq, target_qual)
                spikein.append(read)
        return spikein

infile = sys.argv[1]
passfile = sys.argv[2]
barcodefile = sys.argv[3]
pass_thres = int(sys.argv[4])
qual_thres = float(sys.argv[5])

if len(sys.argv) != 6:
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
        spikein = findSpikein(ccs_name, ccs_pass, barcode_dict, seq, qual)

        seq = revComp(ccs.sequence)
        qual = ccs.quality[::-1]
        spikein.extend(findSpikein(ccs_name, ccs_pass, barcode_dict, seq, qual))
        spikeins.extend(spikein)

for read in spikeins:
        end = start1 = start2 = start3 = ''
        (read_name, read_id, spikein_barcode, seq, qual) = read
        (read_name_serial, ccs_pass) = read_id.split(':')
        ccs_pass = int(ccs_pass)

        suff5 = '(GACGAGCTGTACAAGTG){e<=2}([CGT]{0,3})(A)'
        matches = regex.finditer(suff5, seq, overlapped = False)
        for it in matches:
                end = it.end()
                break
        if  end < 1000:
                cds = seq[:end]
                cds_qual = qual[:end]
                seq = seq[end:]
                qual = qual[end:]

                suff3 = '(GTACTCTGCGTTGATAC){e<=2}'
                matches = regex.finditer(suff3, seq, overlapped = False)
                for it in matches:
                        start1 = it.start()
                        break
                if start1 != '':
                        suff3 = '(G[GT]{0,4})(TACTCTGCGTTGATAC){e<=2}'
                        matches = regex.finditer(suff3, seq, overlapped = False)
                        for it in matches:
                                start2 = it.start()
                                break
                        if start2 != '':
                                seq = seq[:start2]
                                qual = qual[:start2]
                        else:
                                seq = seq[:start1]
                                qual = qual[:start1]
                        tag = 'expected'
                else:
                        tag = 'NO_suff3'

        elif spikein_barcode == 'DNA_spikein_0A':
                suff3 = '(GTACTCTGCGTTGATAC){e<=2}'
                matches = regex.finditer(suff3, seq, overlapped = False)
                for it in matches:
                        start3 = it.start()
                        break
                if start3 != '':
                        seq = seq[:start3]
                        qual = qual[:start3]
                        if seq == '':
                                tag = 'expected'
                        else:
                                tag = 'NO_suff3'
                else:
                        tag = 'NO_suff3'
        else:
                tag = 'NO_suff5'

        qual2 = toDigit(qual)
        if qual2 != []:
                qual_avg2 = qual_avg = statistics.mean(qual2)
        else:
                seq = qual = qual_avg = 'NA'
                qual_avg2 = 93

        if tag == 'expected' and ccs_pass >= pass_thres and (qual_avg >= qual_thres or qual_avg2 >= qual_thres):
                out = [read_name, read_id, spikein_barcode, tag, seq, qual, qual_avg, cds, cds_qual]
                out = '\t'.join([str(x) for x in out])
                print(out)

        elif tag == 'NO_suff3' or tag == 'NO_suff5':
                err = [read_name, read_id, spikein_barcode, tag, seq, qual, qual_avg, cds, cds_qual]
                err = '\t'.join([str(x) for x in err])
                sys.stderr.write(err + '\n')

logfile = sys.argv[1] + '.log.txt'
log_file = open(logfile, 'w')
log_file.write(str(ccs_count) + '\n')
log_file.close()

