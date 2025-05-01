import argparse
import gzip
import sys
import numpy as np

# ASCII values of DNA characters.
bA = ord('A')
bC = ord('C')
bG = ord('G')
bT = ord('T')
bN = ord('N')


# ASCII value equivalent to a quality score of 0.
MIN_QUAL = 33


def reverse_complement(s):
    t = np.zeros_like(s)
    v = t[::-1]
    v[s == bA] = bT
    v[s == bC] = bG
    v[s == bG] = bC
    v[s == bT] = bA
    v[s == bN] = bN
    return t


def read_line(f):
    """Read a line from f with trailing whitespace stripped.

    Unlike f.readline(), this function raises an EOFError instead of returning
    an empty string at the end of the file.
    """
    line = f.readline()
    if len(line) == 0:
        raise EOFError
    return line.rstrip()


def merge_reads(s1, s2, q1, q2, amplen):
    """Merge paired end reads of an amplicon and return sequence,
    quality, and number of mismatches.

    s1 -- Sequence of first read (ndarray of int).
    s2 -- Sequence of second read (ndarray of int).
    q1 -- Quality of first read (ndarray of int).
    q2 -- Quality of second read (ndarray of int).
    amplen -- Length of amplicon (int).

    """
    # If the amplicon is of length L and the reads are lengths l1, l2 then:
    # - read 1 from 0 to L-l2-1 inclusive doesn't overlap
    # - read 1 from L-l2 to l1-1 inclusive overlaps with read 2
    # - read 2 from 0 to l1+l2-L-1 inclusive overlaps with read 1
    # - read 2 from l1+l2-L to its end doesn't overlap

    # A picture for clarity:
    # s1 coords: 0                                      l1-1
    #            |                                      |
    #            ----------------------------------------
    #                                          ------------------------------
    #                                          |        |                   |
    # s1 coords:                               L-l2     |                   L-1
    # s2 coords:                               0        l1+l2-L-1

    # Reverse complement read 2 and reverse its quality scores.
    s2 = reverse_complement(s2)
    q2 = q2[::-1]

    # This is where we'll put the merged sequence and quality score.
    s = np.zeros(amplen, dtype=np.int8)
    q = np.zeros(amplen, dtype=np.int8)

    # If the reads overlap correctly, then s1[offset+i] == s2[i], assuming s2 is
    # the reverse complement of the reverse read.
    offset = amplen - len(s2)

    # Fill in the parts of the merged sequence where the reads don't overlap.
    s[:offset] = s1[:offset]
    q[:offset] = q1[:offset]
    s[len(s1):] = s2[len(s1)+len(s2)-amplen:]
    q[len(s1):] = q2[len(s1)+len(s2)-amplen:]

    # Create a set of views into the overlapping region. We can directly compare
    # vs1[i] to vs2[i] and use that to fill in vs[i] with all indexing taken
    # care of.
    vs1 = s1[offset:]
    vq1 = q1[offset:]
    vs2 = s2[:len(vs1)]
    vq2 = q2[:len(vs1)]
    vs = s[offset:len(s1)]
    vq = q[offset:len(s1)]

    # Quality score of matching bases is the larger of the two quality
    # scores. Quality score of mismatched bases is the difference of the two
    # quality scores. If the mismatched bases have equal quality scores, the
    # base is written as an N with the minimum possible quality.

    # Positions where the reads agree.
    ieq = vs1 == vs2
    vs[ieq] = vs1[ieq]
    vq[ieq] = np.maximum(vq1[ieq], vq2[ieq])

    # Positions where the reads disagree.
    ineq = vs1 != vs2
    mismatches = ineq.sum()

    # Positions where the reads disagree and read 1 has the higher quality.
    ir1 = np.logical_and(ineq, vq1 > vq2)
    vs[ir1] = vs1[ir1]
    vq[ir1] = MIN_QUAL + vq1[ir1] - vq2[ir1]

    # Positions where the reads disagree and read 2 has the higher quality.
    ir2 = np.logical_and(ineq, vq2 > vq1)
    vs[ir2] = vs2[ir2]
    vq[ir2] = MIN_QUAL + vq2[ir2] - vq1[ir2]

    # Positions where the reads disagree and they have equal qualities.
    irn = np.logical_and(ineq, vq1 == vq2)
    vs[irn] = bN
    vq[irn] = MIN_QUAL

    return s, q, mismatches


def bytes_to_str(v):
    return bytearray(v).decode('ascii')


def read_seqs(f):
    while True:
        # Read the sequence ID. If there's nothing to read, then we're done.
        try:
            seq_id = read_line(f)
        except EOFError:
            return

        # If we successfully read a sequence ID, then running out of stuff to
        # read means a truncated record.
        try:
            seq = np.array(bytearray(read_line(f), 'ascii'), dtype=np.int8)
            qual_id = read_line(f)
            qual = np.array(bytearray(read_line(f), 'ascii'), dtype=np.int8)
        except EOFError:
            raise EOFError('EOF while reading sequence.')

        # Some simple checks of the data.
        if seq_id[0] != '@':
            raise ValueError("Sequence ID doesn't begin with '@'.")
        if qual_id[0] != '+':
            raise ValueError("Quality ID doesn't begin with '+'.")
        if len(seq) != len(qual):
            raise ValueError("Sequence and quality are different lengths.")

        yield (seq_id, seq, qual_id, qual)


def write_seq(seq_id, seq, qual, f):
    """Write a FASTQ entry to an open file handle."""
    print(seq_id, file=f)
    print(bytes_to_str(seq), file=f)
    print('+', file=f)
    print(bytes_to_str(qual), file=f)
    
def compare_seq_ids(id1, id2):
	id1_parts = id1.split()
	id2_parts = id2.split()
	l1 = id1_parts[0]
	r1 = " ".join(id1_parts[1:])
	l2 = id2_parts[0]
	r2 = " ".join(id2_parts[1:])
	match = l1 == l2
	prefix = l1 if match else None
	return match, prefix



def run(f1, f2, outf, amplen, max_mm, check_ids):

    written = 0
    discarded = 0

    for i, (r1, r2) in enumerate(zip(read_seqs(f1), read_seqs(f2)), 1):
        seq_id1, seq1, qual_id1, qual1 = r1
        seq_id2, seq2, qual_id2, qual2 = r2

        if check_ids:
            match, prefix = compare_seq_ids(seq_id1, seq_id2)
            if not match:
                raise ValueError('Reads do not appear to match.')
            seq_id = prefix + ' merged'
        else:
            seq_id = 'merged-%08i' % i

        s, q, n_mm = merge_reads(seq1, seq2, qual1, qual2, amplen)

        if max_mm is None or n_mm <= max_mm:
            write_seq(seq_id, s, q, outf)
            written += 1
        else:
            discarded += 1

        total = written + discarded
        if total > 0 and total % 10000 == 0:
            print('written:   %8i' % written, file=sys.stderr)
            print('discarded: %8i' % discarded, file=sys.stderr)

    return written, discarded


def bounded_int(low=None, high=None, msg=None):
    def f(s, msg=msg):
        try:
            i = int(s)
        except:
            raise argparse.ArgumentTypeError("Not an integer: '%s'" % s)

        if msg is None:
            msg = "Invalid value: '%i'"
        msg = msg % i

        if (low is not None and i < low) or (high is not None and i > high):
            raise argparse.ArgumentTypeError(msg)

        return i
    return f

positive_int = bounded_int(low=1, msg="Not a positive integer: '%i'.")
non_negative_int = bounded_int(low=1, msg="Not a non-negative integer: '%i'.")

def parse_args(argv):
    parser = argparse.ArgumentParser(
        description='Merge paired end reads of amplicons from FASTQ files.')

    parser.add_argument(
        'file1',
        help='FASTQ file containing the forward reads (possibly gzipped).')

    parser.add_argument(
        'file2',
        help='FASTQ file containing the reverse reads (possibly gzipped).')

    parser.add_argument(
        'output',
        nargs='?',
        help='Output file. If omitted, stdout will be used.')

    #input the required value for length of amplicon.
    parser.add_argument(
        '-l', '--amplicon-length',
        help='Length of the amplicon.',
        required=True,
        type=positive_int)

    parser.add_argument(
        '-n', '--no-id-check',
        help=('Don\'t check that corresponding forward and reverse reads have'
              'matching sequence IDs.'),
        action='store_true')
    
    #input the required value of max mismatches, usually 4.
    parser.add_argument(
        '-m', '--max-mismatches',
        help=('Maximum number of mismatches a read can have '
              'without being discarded.'),
        type=non_negative_int,
        default=None)

    parser.add_argument(
        '-i', '--input-format',
        help='Is input gzipped or not? Guessed by extension if not specified.',
        choices=['gzip', 'uncompressed'],
        default=None)

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args(sys.argv)

    if args.input_format is None:
        args.input_format = 'gzip' if args.file1.endswith('.gz') else 'uncompressed'

    check_ids = not args.no_id_check

    open_in = gzip.open if args.input_format == 'gzip' else open
    f1 = open_in(args.file1, 'rt')
    f2 = open_in(args.file2, 'rt')

    if args.output is not None:
        outf = open(args.output, 'wt')
    else:
        outf = sys.stdout

    run(f1, f2, outf, args.amplicon_length, args.max_mismatches, check_ids)

    f1.close()
    f2.close()
    if outf is not sys.stdout: outf.close()







