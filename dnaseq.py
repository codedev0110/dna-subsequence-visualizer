#!/usr/bin/env python3

import sys
import unittest
from dnaseqlib import *

### Utility classes ###

# Maps integer keys to a set of arbitrary values.
class Multidict:
    def __init__(self, pairs=None):
        if pairs is None:
            pairs = []
        self.d = {}
        for (k,v) in pairs:
            self.put(k, v)
    # Associates the value v with the key k.
    def put(self, k, v):
        if k in self.d:
            self.d[k].append(v)
        else:
            self.d[k] = [v]
    # Gets any values that have been associated with the key k; or, if
    # none have been, returns an empty sequence.
    def get(self, k):
        if k in self.d:
            return self.d[k]
        return []

# Given a sequence of nucleotides, return all k-length subsequences
# and their hashes.  (What else do you need to know about each
# subsequence?)
def subsequenceHashes(seq, k):
    try:
        subseq = ''
        pos = 0
        while len(subseq) < k:
            subseq += next(seq)
        rh = RollingHash(subseq)
        yield (rh.current_hash(), pos, subseq)
        while True:
            prev = subseq[0]
            nxt = next(seq)
            subseq = subseq[1:] + nxt
            pos += 1
            rh.slide(prev, nxt)
            yield (rh.current_hash(), pos, subseq)
    except StopIteration:
        return

# Similar to subsequenceHashes(), but returns one k-length subsequence
# every m nucleotides.  (This will be useful when you try to use two
# whole data files.)
def intervalSubsequenceHashes(seq, k, m):
    try:
        subseq = ''
        pos = 0
        while len(subseq) < k:
            subseq += next(seq)
        rh = RollingHash(subseq)
        if pos % m == 0:
            yield (rh.current_hash(), pos, subseq)
        while True:
            prev = subseq[0]
            nxt = next(seq)
            subseq = subseq[1:] + nxt
            pos += 1
            rh.slide(prev, nxt)
            if pos % m == 0:
                yield (rh.current_hash(), pos, subseq)
    except StopIteration:
        return

# Searches for commonalities between sequences a and b by comparing
# subsequences of length k.  The sequences a and b should be iterators
# that return nucleotides.  The table is built by computing one hash
# every m nucleotides (for m >= k).
def getExactSubmatches(a, b, k, m):
    table = Multidict()
    for (h, pos, subseq) in intervalSubsequenceHashes(a, k, m):
        table.put(h, (pos, subseq))
    for (h, posb, subseqb) in subsequenceHashes(b, k):
        for (posa, subseqa) in table.get(h):
            if subseqa == subseqb:
                yield (posa, posb)

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print('Usage: {0} [file_a.fa] [file_b.fa] [output.png]'.format(sys.argv[0]))
        sys.exit(1)

    # The arguments are, in order: 1) Your getExactSubmatches
    # function, 2) the filename to which the image should be written,
    # 3) a tuple giving the width and height of the image, 4) the
    # filename of sequence A, 5) the filename of sequence B, 6) k, the
    # subsequence size, and 7) m, the sampling interval for sequence
    # A.
    compareSequences(getExactSubmatches, sys.argv[3], (500,500), sys.argv[1], sys.argv[2], 8, 100)
