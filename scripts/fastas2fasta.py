'''
fastas2fasta.py - concatenate sequences from multiple fasta files
==========================================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Genomics Sequences MultipleAlignments FASTA Manipulation

Purpose
-------

This script reads sequences from two or more :term:`fasta` formatted
files and outputs a new file with the sequences concatenated per
entry.

All files must have the same number of sequences and the id of
the first file is output.

Usage
-----

Example::

   python fastas2fasta.py a.fasta b.fasta > c.fasta

If a.fasta is::

  >1
  AAACC
  >2
  CCCAA

and b.fasta is::

  >a
  GGGGTTT
  >b
  TTTTGGG

then the output will be::

  >1
  AAACCGGGGTTT
  >2
  CCCAATTTTGGG


Type::

   python fastas2fasta.py --help

for command line help.

Command line options
--------------------

'''
import sys
import re

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.FastaIterator as FastaIterator

# ------------------------------------------------------------


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(version="%prog version: $Id:",
                            usage=globals()["__doc__"])

    parser.add_option("-s", "--stack", action="store_true", dest="stack",
                      help=("instead of concatenating reads, write out read 1"
                            " and read 2 alternately"))

    (options, args) = E.Start(parser)

    if len(args) < 2:
        raise ValueError(
            "please supply at least two filenames to concatenate.")

    iterators = []
    for a in args:
        iterators.append(FastaIterator.FastaIterator(IOTools.openFile(a, "r")))

    ninput, noutput, nerrors = 0, 0, 0

    while 1:

        sequences = []
        ids = []

        for iterator in iterators:
            try:
                cur_record = iterator.next()
            except StopIteration:
                break

            sequences.append(re.sub(" ", "", cur_record.sequence))
            ids.append(cur_record.title)

        if not sequences:
            break
        ninput += 1

        if len(sequences) != len(iterators):
            raise "unequal number of sequences in files."

        noutput += 1

        if options.stack:
            options.stdout.write("\n".join(['>' + x[0] + '\n' + x[1] for x in
                                            zip(ids, sequences)]) + "\n")
        else:
            options.stdout.write(">%s\n%s\n" % (ids[0],
                                                "".join(sequences)))

            E.info("ninput=%i, noutput=%i, nerrors=%i" % (ninput, noutput,
                                                          nerrors))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
