import sys
import CGAT.Experiment as E
import ast as ast
import CGAT.IOTools as IOTools
import os


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option(
        "-m", "--method", dest="methods", type="choice", action="append",
        choices=("standardise",),
        help="actions to perform on annotations")

    parser.add_option("-p", "--prefix", dest="prefix", type="str",
                      help="prefix for output files")

    parser.add_option("-l", "--location", dest="loc", type="str",
                      help="path to annotation file")

    parser.add_option("-k", "--keycol", dest="kcol", type="str",
                      help="column name or index of ID")

    parser.add_option("-o", "--othercol", dest="ocol", type="str",
                      help="column name or index of description")

    parser.add_option("-d", "--delim1", dest="d1", type="str",
                      help="column delimiter in input")

    parser.add_option("--delim2", dest="d2", type="str",
                      help="within column delimiter in input")

    parser.add_option("--ignore_delim_key", dest="idk", type="int",
                      help="ignore delimiters within the key column")

    parser.add_option("--ignore_delim_other", dest="ido", type="int",
                      help="ignore delimiters within the description column")

    parser.add_option("--translate_ids", dest="translate", type="int",
                      help="translate ids to another type? 1=yes, 0=no")

    parser.add_option("--translate_file", dest="tfile", type="str",
                      help="path to translation file")

    parser.add_option("--translate_from", dest="tfrom", type="str",
                      help="translate IDs from this type")

    parser.add_option("--translate_to", dest="tto", type="str",
                      help="translate IDs to this type")

    parser.set_defaults(
        methods=[],
        prefix="",
        loc="",
        kcol="0",
        ocol="1",
        d1="\t",
        d2="",
        idk=0,
        ido=0,
        translate=0,
        tto="",
        tfrom="",
        tfile="",
    )

    (options, args) = E.Start(parser, add_pipe_options=True)
    options.d1 = ast.literal_eval("'%s'" % options.d1)
    options.d2 = ast.literal_eval("'%s'" % options.d2)
    for method in options.methods:
        if method == "standardise":
            assert len(options.loc) != 0, """Location of annotation file is
            required"""
            pairs = readFile(options)
            if options.translate == 1:
                tDict = makeTranslationDict(options)
            else:
                tDict = dict()
            GenesToTerms = mapGenesAndTerms(options, tDict, pairs)
            outf = "%s_mapped.tsv" % options.prefix
            saveDict(GenesToTerms, outf)


def readFile(options):

    if options.kcol.isdigit():
        if int(options.kcol) < int(options.ocol):
            colnames = ["c1", "c2"]
        else:
            colnames = ["c2", "c1"]

        usecols = [int(s) for s in [options.kcol, options.ocol]]
        options.kcol = "c1"
        options.ocol = "c2"

    else:
        colnames = options.kcol, options.ocol

        allcols = IOTools.openFile(options.loc).readline().split(options.d1)
        k = allcols.index(options.kcol)
        o = allcols.index(options.ocol)

        usecols = [k, o]

        if k > o:
            colnames = colnames[1], colnames[0]

    i = 0
    res = []
    with IOTools.openFile(options.loc) as inf:
        for line in inf:
            line = line.strip().split(options.d1)
            if i != 0:
                res.append((line[usecols[0]], line[usecols[1]]))
            i += 1
    return res


def makeTranslationDict(options):
    '''
    translates gene IDs into ID types used by different databases
    '''
    idtype = options.tfrom
    transtype = options.tto
    idtsv = options.tfile
    tID = dict()

    assert os.path.exists(idtsv), "translation file does not exist"
    assert len(idtype) != 0, "starting ID type not specified"
    assert len(transtype) != 0, "final ID type not specified"

    i = 0
    with IOTools.openFile(idtsv) as inf:
        for line in inf:
            line = line.strip().split("\t")
            if i == 0:
                assert idtype in line, """%s column not in
                translation file""" % idtype
                ind1 = line.index(idtype)
                assert transtype in line, """%s column not in
                translation file""" % transtype
                ind2 = line.index(transtype)

            else:
                if len(line) > ind1:
                    if line[ind1] not in tID:
                        tID[line[ind1]] = set()
                    if len(line) > ind2:
                        tID[line[ind1]].add(line[ind2])

            i += 1

    return tID


def mapGenesAndTerms(options, tDict, pairlist):
    '''
    Builds a GenesToTerms dictionary based on the MappingDF from the
    mapping file.
    '''
    GenesToTerms = dict()

    for pair in pairlist:
        #  split columns by delim2 if needed
        if options.d2 == "" or options.idk == 1:
            genes = [str(pair[0])]
        else:
            genes = str(pair[0]).split(options.d2)

        if options.d2 == "" or options.ido == 1:
            terms = [str(pair[1])]
        else:
            terms = str(pair[1]).split(options.d2)

        #  build a dictionary - one key for each item in keycol with a set
        #  of every associated othercol term
        for gene in genes:
            gene = gene.replace(" ", "")
            if gene in tDict or len(tDict.keys()) == 0:
                if gene in tDict:
                    t_gene = tDict[gene]
                else:
                    t_gene = set([gene])

                for onegene in t_gene:
                    if onegene not in GenesToTerms:
                        GenesToTerms[onegene] = set()

                    for term in terms:
                        if term != '':
                            term = term.replace(" ", "_").replace(",", "")
                            GenesToTerms[onegene].add(term)

    # remove empty sets from dicts
    for item in GenesToTerms.items():
        if len(item[1]) == 0:
            del GenesToTerms[item[0]]
    return GenesToTerms


def saveDict(adict, outfile):
    '''
    Stores a parsed mapping file in a flat file.
    '''
    out = IOTools.openFile(outfile, "w")
    out.write("id1\tid2\n")

    for id1, id2 in adict.items():
        if len(id1) != 0:
            out.write("%s\t%s\n" % (id1,
                                    ",".join(id2)
                                    if len(id2) != 0 else "-"))
    out.close()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
