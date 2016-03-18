'''
tax2tax.py - parse metagenomic taxonomy files
=============================================

:Author: Jethro Johnson
:Release: $Id$
:Date: |today|
:Tags: Metagenomics

Purpose
-------

Several  databases exist to provide taxonomic information for archived
metagenomic sequences. Information downloaded from these databases
commonly consists of a fasta file and a taxonomy file. Taxonomy files
have a basic structure of sequence ID followed by a ';' separated list
of taxonomic levels. However, beyond this, there is limited consistency
with regard to formatting.

This script exists to parse taxonomy files downloaded from databases
such as Greengenes, RDP & Sliva, or output by programmes such as
metaphlabn, into a common format. Various different options exist for
defining the format of the output files.

'''

import CGAT.Experiment as E
import sys
import re


def main(argv=None):

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id$",
        usage=globals()["__doc__"])

    parser.add_option("-f", "--field-split", dest="field_split",
                      help="regular expression for defining the separation"
                      " between sequence ID and phylogeny information")
    parser.add_option("-t", "--tax-split", dest="tax_split",
                      help="regular expression for defining the separation"
                      " between the tax fields in phylogeny")
    parser.add_option("-e", "--empty-value", dest="empty_value",
                      help="value to include in outfile when tax data"
                      " are missing")
    parser.add_option("-a", "--add-prefix", dest="add_prefix",
                      action='store_true', help="add single character prefix,"
                      " with double underscore-separation, to taxonomic field")
    parser.add_option("-g", "--clean-greengenes", dest="clean_greengenes",
                      action="store_true", help="remove additional brackets"
                      "([]) added in greengenes databases")
    parser.add_option("-s", "--remove-whitespace", dest="remove_white_space",
                      action="store_true", help="remove spaces from tax field")
    parser.add_option("-p", "--prefixes", dest="prefixes",
                      help="comma-separated list ofprefixes to retain in the"
                      " output file")
    parser.add_option("-o", "--out-format", type='choice', dest="out_format",
                      choices=['table', 'taxonomy'],
                      help="the chosen output format")
    parser.add_option("-c", "--headers", dest="headers", help="comma-separated"
                      " list of headers to include in output tables")
    parser.add_option("-r", "--replace-species", dest="replace_species",
                      action="store_true", help="replace species with seq ID")
    parser.add_option("-m", "--min-tax-levels", dest="min_tax_levels",
                      help="minimum number of tax levels to include in output")
    parser.add_option("-d", "--drop-empty", dest="drop_empty",
                      action="store_true", help="keep empty fields at"
                      " end of line in taxonomic output")
    parser.add_option("-x", "--add-seq-preifx", dest="add_seq_prefix",
                      help="add a prefix to the output sequence ID")
    parser.add_option("-b", "--remove-confidence-thresholds",
                      dest="remove_confidence", action="store_true",
                      help="remove the confidence thresholds added by"
                      " RDP classifier")
    parser.add_option("-z", "--final-column", dest="final_column", help="Name"
                      " of a final column containing counts (if present)")
    parser.add_option("-y", "--remove-final-column",
                      dest="remove_final_column", action="store_true",
                      help="drop final column from output, if present")

    parser.set_defaults(field_split='\t|\s{4,}',
                        tax_split=';',
                        empty_value='unknown',
                        add_prefix=False,
                        clean_greengenes=False,
                        remove_white_space=False,
                        prefixes='k__,p__,c__,o__,f__,g__,s__',
                        out_format='taxonomy',
                        headers='SeqID,Kingdom,Phylum,Class,Order,Family,Genus,Species',
                        replace_species=False,
                        min_tax_levels=0,
                        drop_empty=False,
                        add_seq_prefix='',
                        remove_confidence=False,
                        final_column=False,
                        remove_final_column=False
                        )

    (options, args) = E.Start(parser)

    prefixes = options.prefixes.split(',')

    if options.headers:
        if options.out_format == 'table':
            hdrs = options.headers.split(',')
            if options.final_column and not options.remove_final_column:
                hdrs.append(options.final_column)
            options.stdout.write('\t'.join(options.headers.split(',')) + '\n')
        else:
            E.warn('Outputting taxonomy file, so headers are ignored')

    fc = None
    for record in options.stdin:
        record = record.strip()
        seqID, phylogeny = re.split(options.field_split, record)

        # there may be a terminal ';' in file
        if phylogeny.endswith(';'):
            phylogeny = phylogeny.split(options.tax_split)[:-1]
        else:
            phylogeny = phylogeny.split(options.tax_split)

        # pop final column, if present
        if options.final_column:
            fc = phylogeny.pop()

        # fetch output taxonomic fields
        tax_out = []
        # allow ambiguous number of tax fields if prefixes are present
        if None not in [re.match('[a-z]+__', x) for x in phylogeny]:
            nf = 0
            # retrieve only the main tax levels
            for p in prefixes:
                loc = [i for i, elem in enumerate(phylogeny)
                       if elem.startswith(p)]
                assert len(loc) <= 1, \
                    'Multiple instances of {} in {}'.format(p,
                                                            "\t".join(
                                                                phylogeny))
                if len(loc) == 1:
                    p_out = phylogeny[loc[0]].split('__')[1]
                else:
                    p_out = None
                if p_out:
                    tax_out.append(p_out)
                    nf += 1
                else:
                    tax_out.append(options.empty_value)

            # only output fields with a minimum number of taxonomic levels
            if nf < int(options.min_tax_levels):
                continue

        else:
            # if tax level ids are absent expect at most 7 fields
            assert len(phylogeny) <= 7, \
                'Ambiguous tax fields: {}'.format("\t".join(phylogeny))
            # only output fields with a minimum number of taxonomic levels
            if len(phylogeny) < int(options.min_tax_levels):
                continue
            for i in range(7):
                try:
                    tax_out.append(phylogeny[i])
                except IndexError:
                    tax_out.append(options.empty_value)

        # option to add a prefix to the sequence ID
        if options.add_seq_prefix:
            seqID = options.add_seq_prefix + seqID

        # clean up taxonomic fields
        if options.replace_species:
            # guess empty fields based on prefixes
            tax_out = [prefixes[x[0]] + seqID if x[1] == options.empty_value
                       else x[1] for x in enumerate(tax_out)]
            # for simplicity, the last entry will just be seqID
            tax_out[-1] = seqID
        if options.clean_greengenes:
            tax_out = [re.sub('\[|\]', '', x) for x in tax_out]
        if options.remove_white_space:
            tax_out = [re.sub('\s+', '_', x) for x in tax_out]
        if options.add_prefix:
            tax_out = ["".join(x) for x in zip(prefixes, tax_out)]
        if options.remove_confidence:
            tax_out = [re.sub('\(\d+\.?\d*\)', '', x) for x in tax_out]

        # add final column, if present and requested
        if fc and not options.remove_final_column:
            tax_out.append(fc)

        if options.out_format == 'taxonomy':
            line_out = seqID + '\t' + ';'.join(tax_out) + ';\n'
            if options.drop_empty:
                regex = re.compile('(' + options.empty_value + ';)+$')
                line_out = re.sub(regex, '', line_out)
        else:
            line_out = seqID + '\t' + '\t'.join(tax_out) + '\n'

        options.stdout.write(line_out)

if __name__ == "__main__":
    sys.exit(main(sys.argv))
