#!/usr/bin/env python
import argparse
import csv
import gzip
import logging
import os
import sys
import queue
import threading
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from collections import defaultdict, namedtuple
from concurrent.futures import ThreadPoolExecutor, as_completed

# Initialize logging
logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,  # Set level to INFO for the log output we need
    datefmt='%Y-%m-%d %H:%M:%S')

# A shared log queue for storing log messages in a FIFO order
log_queue = queue.Queue()

# Define constants
BARCODE_LEN = 8  # Hardwired for now
UNDETERMINED = "Undetermined"
Prefix = namedtuple("Prefix", "tag j7_name j7_seq j5_name j5_seq")
Suffix = namedtuple("Suffix", "tag i7_name i7_seq i5_name i5_seq")

def revcomp(seq, t=str.maketrans('ACGTacgt', 'TGCAtgca')):
    """Reverse-complement of a constrained sequence string."""
    return seq.translate(t)[::-1]

def fastq_reader(filename):
    """Wrap Biopython streaming fastq reader. Assumes gzip'ed input"""
    with gzip.open(filename, "rt") as infile:
        for (hdr, seq, qual) in FastqGeneralIterator(infile):
            yield (hdr, seq, qual)

def emit_fastq(outfile, hdr, i7, i5, j7, j5, seq, qual):
    acc, flag = hdr.split(" ")
    _, _, flowcell, lane, tile, x, y = acc.split(":")
    print(f"@{flowcell}:{lane}:{tile}:{x}:{y}_{i7}_{i5}_{j7}_{j5} {flag}",
          seq, "+", qual, sep="\n", file=outfile)

def log_worker():
    """Thread for handling log messages from the queue."""
    while True:
        msg = log_queue.get()
        if msg is None:  # None is a signal to stop the logging thread
            break
        logging.info(msg)

class SampleBarcodes:
    """Track individual barcodes and process reads."""
    
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
                result[seq] = seq  # Exact matching
                continue
            u = set()
            for i in range(len(seq)):
                for k in "ACGTN":
                    s = seq[:i] + k + seq[i + 1:]
                    u.add(s)
            for s in u:
                if s in result:
                    logging.warning(f"COLLISION {seq} {s}")
                result[s] = seq
        return result

    def enumerate_sample_tables(self, prefixes, suffixes):
        """Enumerate sample tables for barcodes."""
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
        """Open files for writing output."""
        assert os.path.exists(outdir), f"ERROR: Output directory {outdir} does not exist"
        for sample in self.sample_names:
            self.outfiles[sample] = (
                gzip.open(os.path.join(outdir, f"{sample}_R1.fq.gz"), "wt"),
                gzip.open(os.path.join(outdir, f"{sample}_R2.fq.gz"), "wt")
            )
        self.outfiles[UNDETERMINED] = (
            gzip.open(os.path.join(outdir, f"Undetermined_R1.fq.gz"), "wt"),
            gzip.open(os.path.join(outdir, f"Undetermined_R2.fq.gz"), "wt")
        )

    def process_fastq(self, outdir, r1_filename, r2_filename, i1_filename, i2_filename, chunk_size=10000):
        """Process the fastq files using multithreading."""
        self.open_outfiles(outdir)

        total, matched = 0, 0
        r1_reader = fastq_reader(r1_filename)
        r2_reader = fastq_reader(r2_filename)
        i1_reader = fastq_reader(i1_filename)
        i2_reader = fastq_reader(i2_filename)

        chunk = []
        num_threads = os.cpu_count()  # Detect available cores
        logging.info(f"Detected {num_threads} CPU threads for processing")

        with ThreadPoolExecutor(max_workers=num_threads) as executor:
            futures = []

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
                chunk.append((h1, s1, q1, h2, s2, q2, i1, i2))

                if len(chunk) >= chunk_size:
                    futures.append(executor.submit(process_chunk, chunk, self))
                    chunk = []

                # Log percentage every 1 million reads processed
                if total % 1000000 == 0 and total > 0:
                    percentage_matched = (matched / total) * 100
                    log_queue.put(f"Processed {total} reads, matched {matched} ({percentage_matched:.2f}%)")

            if chunk:
                futures.append(executor.submit(process_chunk, chunk, self))

            # Aggregate matched counts from all completed threads
            for future in as_completed(futures):
                matched += future.result()

        percentage_matched = (matched / total) * 100 if total > 0 else 0
        log_queue.put(f"Total matched {matched} of {total} reads ({percentage_matched:.2f}%)")

def process_chunk(chunk, sample_barcodes):
    """Process a chunk of reads."""
    matched_chunk = 0
    for h1, s1, q1, h2, s2, q2, i1, i2 in chunk:
        j7 = i1[:BARCODE_LEN]
        i7 = i1[-BARCODE_LEN:]

        if not sample_barcodes.forward_mode:
            i2 = revcomp(i2)

        i5 = i2[:BARCODE_LEN]
        j5 = i2[-BARCODE_LEN:]

        if i7 in sample_barcodes.i7_lookup and i5 in sample_barcodes.i5_lookup and j7 in sample_barcodes.j7_lookup and j5 in sample_barcodes.j5_lookup:
            i7 = sample_barcodes.i7_lookup[i7]
            i5 = sample_barcodes.i5_lookup[i5]
            j7 = sample_barcodes.j7_lookup[j7]
            j5 = sample_barcodes.j5_lookup[j5]

            key = (i7, i5, j7, j5)
            if key in sample_barcodes.sample_lookup:
                sample = sample_barcodes.sample_lookup[key]
                matched_chunk += 1  # Increment matched reads within this chunk
            else:
                sample = UNDETERMINED
        else:
            sample = UNDETERMINED

        o1, o2 = sample_barcodes.outfiles[sample]
        emit_fastq(o1, h1, i7, i5, j7, j5, s1, q1)
        emit_fastq(o2, h2, i7, i5, j7, j5, s2, q2)

    return matched_chunk  # Return the number of matched reads for this chunk

def read_tn5_barcodes(filename):
    """Read Tn5 barcodes from a CSV file."""
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

            prefixes.append(Prefix(prefix, j7_name, j7_seq, j5_name, j5_seq))
    return prefixes

def read_primer_barcodes(filename):
    """Read primer barcodes from a CSV file."""
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
            suffixes.append(Suffix(suffix, i7_name, i7_seq, i5_name, i5_seq))
    return suffixes

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
                    prog='sciCTextract',
                    description='Simple sciCUT&Tag demultiplexer')

    parser.add_argument('--Tn5_Barcode', type=str, required=True,
                        help='CSV format table of Tn5 barcodes and sample name prefixes')
    parser.add_argument('--Primer_Barcode', type=str, required=True,
                        help='CSV format table of Primer barcodes and sample name suffixes')
    parser.add_argument('--outdir', type=str, required=True,
                        help='Name of pre-existing output directory for Fastq files')
    parser.add_argument('--exact-mode', action='store_true',
                        help='Require exact matching of barcode sequences rather than default one-mismatch tolerance. Not generally recommended')
    parser.add_argument('--forward-mode', action='store_true',
                        help='Use forward-strand (Workflow A) processing for MiSeq, HiSeq, etc.')
    parser.add_argument('r1_filename', type=str, help='First read filename')
    parser.add_argument('r2_filename', type=str, help='Second read filename')
    parser.add_argument('i1_filename', type=str, help='First index read filename (i7)')
    parser.add_argument('i2_filename', type=str, help='Second index read filename (i5)')
    args = parser.parse_args()
    return args

def start_logging_thread():
    """Start the logging thread."""
    log_thread = threading.Thread(target=log_worker)
    log_thread.daemon = True  # Set as a daemon thread so it exits when the main thread does
    log_thread.start()
    return log_thread

def stop_logging_thread(log_thread):
    """Stop the logging thread."""
    log_queue.put(None)  # Sending None to signal the log worker to stop
    log_thread.join()

def main():
    """Main execution function."""
    log_thread = start_logging_thread()

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

    stop_logging_thread(log_thread)

if __name__ == "__main__":
    sys.exit(main())
