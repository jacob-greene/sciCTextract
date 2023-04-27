#!/usr/bin/env python
import argparse
import csv
import gzip
import logging
import os
import sys
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from collections import defaultdict, namedtuple

"""
Simple demultiplexing for sciCUT&Tag runs. Assumes Illumina sequencers,
the process has been tested on HiSeq 2500 and NextSeq 2000 output.
"""
logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')

BARCODE_LEN = 8 # Hardwired for now
UNDETERMINED = "Undetermined"
Prefix = namedtuple("Prefix", "tag j7_name j7_seq j5_name j5_seq")
Suffix = namedtuple("Suffix", "tag i7_name i7_seq i5_name i5_seq")

def revcomp(seq, t=str.maketrans('ACGTacgt', 'TGCAtgca')):
    """Reverse-complement of a constrained sequence string. Marginally faster
    than Bio.Seq but not general purpose (no IUPAC, works on strings)"""
    return seq.translate(t)[::-1]

def fastq_reader(filename):
    """Wrap Biopython streaming fastq reader. Assumes gzip'ed input"""
    with gzip.open(filename, "rt") as infile:
        for (hdr, seq, qual) in FastqGeneralIterator(infile):
            yield(hdr, seq, qual)

def emit_fastq(outfile, hdr, i7, i5, j7, j5, seq, qual):
    acc, flag = hdr.split(" ")
    _, _, flowcell, lane, tile, x, y = acc.split(":")
    print(f"@{flowcell}:{lane}:{tile}:{x}:{y}_{i7}_{i5}_{j7}_{j5} {flag}",
          seq, "+", qual, sep="\n", file=outfile)

class SampleBarcodes:
    """Track individual barcodes with optional, but highly recommended,
    single subsitiution error tolerance. Includes utilities to write 
    process reads and write to sample-specific files."""
    # barcode lookup tables optionally supporting one mismatch
    def __init__(self, exact_mode=False, forward_mode=False):
        self.exact_mode = exact_mode
        self.forward_mode = forward_mode
        self.i7_lookup = {}
        self.i5_lookup = {}
        self.j7_lookup = {}
        self.j5_lookup = {}
        self.sample_names = set()
        self.sample_lookup = {}
        self.outfiles = {}

    def enumerate_substitutions(self, lookup):
        result = {}
        for seq in lookup:
            if self.exact_mode:
                result[seq] = seq # Exact matching; not recommended
                continue
            u = set()
            for i in range(len(seq)):
                for k in "ACGTN":
                    s = seq[:i] + k + seq[i+1:]
                    u.add(s)
            for s in u:
                if s in result:
                    print("COLLISION", seq + " " + s, file=sys.stderr)
                result[s] = seq
        return result

    def enumerate_sample_tables(self, prefixes, suffixes):

        # separate corrections for each of the four types of barcode
        i7_lookup = set()
        i5_lookup = set()
        j7_lookup = set()
        j5_lookup = set()

        for suffix in suffixes:
            i7_lookup.add(suffix.i7_seq)
            i5_lookup.add(suffix.i5_seq)
            for prefix in prefixes:
                j7_lookup.add(prefix.j7_seq)
                j5_lookup.add(prefix.j5_seq)
                key = (suffix.i7_seq, suffix.i5_seq, prefix.j7_seq, prefix.j5_seq)
                name = prefix.tag + suffix.tag
                assert key not in self.sample_lookup, f"Collision on sample {name} for barcodes {key}"
                self.sample_names.add(name)
                self.sample_lookup[key] = name

        self.i7_lookup = self.enumerate_substitutions(i7_lookup)
        self.i5_lookup = self.enumerate_substitutions(i5_lookup)
        self.j7_lookup = self.enumerate_substitutions(j7_lookup)
        self.j5_lookup = self.enumerate_substitutions(j5_lookup)

    def open_outfiles(self, outdir):
        assert os.path.exists(outdir), f"ERROR: Output direcory {outdir} does not exist"
        for sample in self.sample_names:
            self.outfiles[sample] = (
                gzip.open(os.path.join(outdir, f"{sample}_R1.fq.gz"), "wt"),
                gzip.open(os.path.join(outdir, f"{sample}_R2.fq.gz"), "wt")
            )
        self.outfiles[UNDETERMINED] = (
            gzip.open(os.path.join(outdir, f"Undetermined_R1.fq.gz"), "wt"),
            gzip.open(os.path.join(outdir, f"Undetermined_R2.fq.gz"), "wt")
        )

    def process_fastq(self, outdir, r1_filename, r2_filename, i1_filename, i2_filename):

        self.open_outfiles(outdir)

        total, matched = 0, 0

        r1_reader = fastq_reader(r1_filename)
        r2_reader = fastq_reader(r2_filename)
        i1_reader = fastq_reader(i1_filename)
        i2_reader = fastq_reader(i2_filename)

        while True:
            try:
                h1, s1, q1 = next(r1_reader)
                h2, s2, q2 = next(r2_reader)
                _h1, i1, _ = next(i1_reader)
                _h2, i2, _ = next(i2_reader)
            except StopIteration:
                break

            assert _h1 == h1, h1
            assert _h2 == h2, h2

            total += 1
            j7 = i1[:BARCODE_LEN]
            i7 = i1[-BARCODE_LEN:]

            # Need to flip for NextSeq & NovaSeq
            if not self.forward_mode:
                i2 = revcomp(i2)

            i5 = i2[:BARCODE_LEN]
            j5 = i2[-BARCODE_LEN:]

            if i7 in self.i7_lookup \
               and i5 in self.i5_lookup \
               and j7 in self.j7_lookup \
               and j5 in self.j5_lookup:
                matched += 1

                i7 = self.i7_lookup[i7]
                i5 = self.i5_lookup[i5]
                j7 = self.j7_lookup[j7]
                j5 = self.j5_lookup[j5]

                key = (i7, i5, j7, j5)
                if key in self.sample_lookup:
                    sample = self.sample_lookup[key]
                else:
                    # Should not reach this
                    log.warn(f"Valid barcodes were not mapped to a sample {i7}_{i5}_{j7}_{j5}")
                    sample = UNDETERMINED
            else:
                sample = UNDETERMINED
            o1, o2 = self.outfiles[sample]
            emit_fastq(o1, h1, i7, i5, j7, j5, s1, q1)
            emit_fastq(o2, h2, i7, i5, j7, j5, s2, q2)

            if total % 100000 == 0:
                logging.info(f"Matched {matched} reads of {total}")

        for sample in self.outfiles:
            o1, o2 = self.outfiles[sample]
            o1.close()
            o2.close()

        logging.info(f"Matched {matched} reads of {total}")

def read_tn5_barcodes(filename):
    prefixes = []
    with open(filename) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            prefix = row["Sample Name"]
            assert prefix != ""
            j7_name = row["Tn5_s7"]
            j7_seq = row["Tn5_s7_seq"]
            j5_name = row["Tn5_s5"]
            j5_seq = row["Tn5_s5_seq"]

            if "sci_i7_Ad_seq" in row and row["sci_i7_Ad_seq"] != "":
                assert False
            if "sci_i5_Ad_seq" in row and row["sci_i5_Ad_seq"] != "":
                assert False

            prefixes.append(Prefix(prefix, j7_name, j7_seq, j5_name, j5_seq))
    return prefixes

def read_primer_barcodes(filename):
    suffixes = []
    with open(filename) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            suffix = row["ID"]
            assert suffix != ""
            i7_name = row["i7_index_id"]
            i7_seq = row["i7_index_seq"]
            i5_name = row["i5_index_id"]
            i5_seq = row["i5_index_seq"]
            assert i7_name.startswith("P7_") and i5_name.startswith("P5_"), i7_name
            i7_name = "Ad2." + i7_name[3:] 
            i5_name = "Ad1." + i5_name[3:] 
            suffixes.append(Suffix(suffix, i7_name, i7_seq, i5_name, i5_seq))
    return suffixes

def parse_args():
    parser = argparse.ArgumentParser(
                    prog='sciCTextract',
                    description='Simple sciCUT&Tag demultiplexer')

    parser.add_argument('--Tn5_Barcode',
                        type=str,
                        required=True,
                        help='CSV format table of Tn5 barcodes and sample name prefixes')
    parser.add_argument('--Primer_Barcode',
                        type=str,
                        required=True,
                        help='CSV format table of Primer barcodes and sample name suffixes')
    parser.add_argument('--outdir',
                        type=str,
                        required=True,
                        help='Name of pre-existing output directory for Fastq files')
    parser.add_argument('--exact-mode',
                        action='store_true',
                        help='Require exact matching of barcode sequences rather than default one-mismatch tolerance. Not generally recommended')
    parser.add_argument('--forward-mode',
                        action='store_true',
                        help='Use forward-strand (Workflow A) processing for MiSeq, HiSeq, etc.')
    parser.add_argument('r1_filename',
                        type=str,
                        help='First read filename')
    parser.add_argument('r2_filename',
                        type=str,
                        help='Second read filename')
    parser.add_argument('i1_filename',
                        type=str,
                        help='First index read filename (i7)')
    parser.add_argument('i2_filename',
                        type=str,
                        help='Second index read filename (i5)')
    args = parser.parse_args()
    return args

def main():
    args = parse_args()

    tn5_barcode_filename = args.Tn5_Barcode
    primer_barcode_filename = args.Primer_Barcode
    outdir = args.outdir

    r1_filename = args.r1_filename
    r2_filename = args.r2_filename
    i1_filename = args.i1_filename
    i2_filename = args.i2_filename

    prefixes = read_tn5_barcodes(tn5_barcode_filename)
    suffixes = read_primer_barcodes(primer_barcode_filename)

    sample_barcodes = SampleBarcodes(exact_mode=args.exact_mode, forward_mode=args.forward_mode)
    sample_barcodes.enumerate_sample_tables(prefixes, suffixes)
    sample_barcodes.process_fastq(outdir, r1_filename, r2_filename, i1_filename, i2_filename)

if __name__ == "__main__":
    sys.exit(main())
