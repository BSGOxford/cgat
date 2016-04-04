'''
_annot - direct annotation of ensemblg to annotation term - something to test
for enrichment

_details - further details of an annotation e.g definitions of go terms

_geneid - alternative gene id which may be used in other databases

_other - another annotation (e.g. gene name, position etc.) which are of
interest but are not relevent in enrichment testing)

_ont - hierarchical ontology

'''

import CGATPipelines.Pipeline as P
from Bio import Entrez
import numpy as np
import httplib2
import json as json
import sqlite3
from intermine.webservice import Service
import string
import re
import os
import datetime
import xml.etree.ElementTree as ET
import pandas as pd

# load options from the config file
PARAMS = P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])

dt = datetime.datetime.now()
dtstring = "%s%s%s_%sh%sm%ss" % (
    dt.year, dt.strftime("%B")[:3], dt.day, dt.hour, dt.minute, dt.second)
os.mkdir("%s_%s" % (PARAMS['db_hostname'], dtstring))

dbname = "%s_%s/%s" % (PARAMS['db_hostname'], dtstring, PARAMS['db_name'])
Entrez.email = PARAMS['entrez_email']
ohost = PARAMS['entrez_host']


def testAPIs():
    '''
    Tests the APIs used by the script with simple queries to check that
    they are running as expected.
    If adding a new annotation from a different API add a new test here
    and return 0 if the test fails.
    '''
    #  Test Entrez API - Look up symbol APOBEC3G, Entrez ID 60489 should be
    #  amongst the results
    E = Entrez.esearch(db="gene", term="(APOBEC3G[Preferred+Symbol])",
                       retmode="text", retmax=1000000)
    res = Entrez.read(E)
    E.close()
    if "60489" in res['IdList']:
        print "Entrez API test OK"
    else:
        print "Entrez API test failed"
        return 0

    # Test MyGeneInfo API - Look up symbol for Entrez ID 60489, result
    # should be APOBEC3G.
    h = httplib2.Http()
    res, con = h.request('http://mygene.info/v2/gene', 'POST',
                         'ids=60489&fields=symbol',
                         headers={'content-type':
                                  'application/x-www-form-urlencoded'})
    # convert json format into Python
    jcon = json.loads(con)
    if jcon[0]['symbol'] == "APOBEC3G":
        print "MyGeneInfo test OK"
    else:
        print "MyGeneInfo test failed"
        return 0

    # Test MouseMine API - Look up symbol for Apoe gene, should return Apoe.
    if int(PARAMS['homologues_mousephenotype']) == 1 or int(
            PARAMS['homologues_mousepathway']) == 1:

        service = Service("http://www.mousemine.org/mousemine/service")
        query = service.new_query("Gene")
        query.add_view("symbol")
        query.add_constraint("Gene", "LOOKUP", "Apoe", "M. musculus", code="A")
        for row in query.rows():
            symbol = row['symbol']
        if symbol == "Apoe":
            print "MouseMine test OK"
        else:
            print "MouseMine test failed"

    # Test HumanMine API - Look up symbol for APOBEC3G, should return APOBEC3G.
    if int(PARAMS['homologues_hpo']) == 1:
        service = Service("http://www.humanmine.org/humanmine/service")
        query = service.new_query("Gene")
        query.add_view("symbol")
        query.add_constraint("Gene", "LOOKUP", "APOBEC3G", code="A")
        for row in query.rows():
            symbol = row['symbol']
        if symbol == "APOBEC3G":
            print "HumanMine test OK"
        else:
            print "HumanMine test failed"
            return 0

    # If the code reaches here, everything is OK.
    return 1


def add_db_table(dbname, tablename, cnames, ctypes, zipped, PS=False):
    '''
    Adds a database table, named tablename, to the database dbname.
    Cnames - column names as a Python list
    Ctypes - column types as strings in a Python list
    e.g. ['varchar (250)', 'varchar(25)', 'int']
    zipped - a list of tuples to fill the table, each item in the list should
    correspond to a row and should be a tuple of the value for each column
    in that row

    E.g.
    add_db_table(
    "csvdb", "entrez2symbol",
    ["entrez_id", "gene_symbol"],
    ['int', 'varchar (25)'],
    [(1111, 'aaa'), (2222, 'bbb'), (3333, 'ccc)])

    would create the table "entrez2symbol" in the database "csvdb" with the
    layout
    entrez_id    gene_symbol
    1111         aaa
    2222         bbb
    3333         ccc

    If PS is True SQL statements are printed before they are run - used for
    debugging.
    '''
    #  check there are the right number of column names for the dataset
    assert len(cnames) == len(zipped[0]), "length of column names doesn't \
    match number of fields"

    # connect to the database
    dbh = sqlite3.connect(dbname)
    c = dbh.cursor()

    # build parts of the SQL statement
    cnamestring = ",".join(cnames)
    cnamestypes = ",".join([" ".join(z) for z in zip(cnames, ctypes)])
    qs = ",".join("?" * len(cnames))

    # if the table already exists in the database, delete it
    statement = '''DROP TABLE IF EXISTS '%s' ;''' % (tablename)
    c.execute(statement)

    # create the table with the appropriate columns and column types
    statement = 'CREATE TABLE %s (%s);' % (tablename, cnamestypes)
    if PS is True:
        print statement
    c.execute(statement)

    # populate the table with the data
    statement = '''INSERT INTO %s (%s) VALUES (%s);''' % (tablename,
                                                          cnamestring, qs)
    if PS is True:
        print statement
    c.executemany(statement, zipped)
    dbh.commit()

    # disconnect
    c.close()
    dbh.close()


def pd2sql(df, tabname):
    dbh = sqlite3.connect(dbname)
    df.to_sql(tabname, dbh, index=False, if_exists='replace')
    dbh.close()


def host2symbol(hostids):
    '''
    Translates entrez ids to gene symbols for different hosts, using the
    MyGeneInfo API.

    hostids should be a Python list of entrez IDs to translate.

    A gene found in multiple hosts will have a different entrez ID in each host
    but may have a common symbol across hosts.
    '''
    symboldict = dict()

    # If the list of entrez ids is long, cut it into chunks of 500 ids to send
    # to the API in each query (this size seems to run the fastest overall).
    allgenes = np.unique(np.array(hostids))
    if len(allgenes) > 500:
        a = len(allgenes) / 500
        segs = np.array_split(allgenes, a)
    else:
        segs = [allgenes]
    for seg in segs:
        # Retrieve symbols
        h = httplib2.Http()
        headers = {'content-type': 'application/x-www-form-urlencoded'}
        params = 'ids=%s&fields=symbol' % (",".join(seg))
        res, con = h.request('http://mygene.info/v2/gene', 'POST', params,
                             headers=headers)
        # Parse Json output into Python format and put it into a dictionary.
        jcon = json.loads(con)
        for line in jcon:
            if 'symbol' in line:
                symboldict[line['query']] = line['symbol']
    return symboldict


def mineOntology(host_abb, genes, views, constraints=[], host="",
                 name=None, hostpd=None, ind=None):
    '''
    Simplifies mining from mousemine, humanmine and other databases with
    APIs which use intermine.

    host_abb - host abbreviation in the database name - mouse for mousemine,
    human for humanmine etc.  These seem to be standardised for all hosts.
    genes - list of genes to pull from the database
    views - which fields to pull from the database about these genes
    constraints - how to filter the data
    host - some databases require a host name in the query - this is formatted
    as e.g. "M. musculus" or "H. sapiens" (the space is needed).
    table_name - name of database table to put the data into

    E.g. if m is a list of mouse gene symbols
    views = ["ontologyAnnotations.ontologyTerm.description",
             "ontologyAnnotations.ontologyTerm.name"]
    constraints = [
        ("Gene.ontologyAnnotations.ontologyTerm.namespace=MPheno.ontology")]

    mineOntology('mouse',
                 m,
                 views,
                 constraints,
                 "M. musculus",
                 "mouse_phenotype")

    Would add a table to the database called "mouse_phenotype" containing
    the ontologyAnnotations.ontologyTerm.name and
    ontologyAnnotations.ontologyTerm.description fields from MPheno.ontology
    set of annotations for genes in list m.

    The "QueryBuilder" on the database websites make it easier to put
    together these queries, e.g.
    http://www.mousemine.org/mousemine/customQuery.do

    '''
    # If the list of gene ids is long, cut it into chunks of 1000 ids to send
    # to the API in each query (this size seems to run the fastest overall).

    glist = np.array(genes)
    if len(glist) > 1000:
        a = len(glist) / 1000
        segs = np.array_split(glist, a)
    else:
        segs = [glist]

    # store the data in here
    z = []

    # API uses letters to distinguish between constraints
    alpha = list(string.ascii_uppercase)

    for seg in segs:
        # Connect to the API
        service = Service("http://www.%smine.org/%smine/service" % (host_abb,
                                                                    host_abb))
        query = service.new_query("Gene")
        query.add_view(",".join(views))
        # Some databases require a host name
        if host != "":
            query.add_constraint("Gene", "LOOKUP", ",".join(seg), host,
                                 code="A")
        else:
            query.add_constraint("Gene", "LOOKUP", ",".join(seg), code="A")

        # Apply the constraints
        if len(constraints) != 0:
            i = 1
            for constraint in constraints:
                letter = alpha[i]
                if len(constraint.split("=")) == 2:
                    L = constraint.split("=")
                    query.add_constraint(L[0], "=", L[1], code=letter)
                elif re.search("IS NOT NULL", constraint):
                    p1 = constraint.replace(" IS NOT NULL", "")
                    query.add_constraint(p1, "IS NOT NULL", code=letter)
                i = i + 1

        # Parse the output into a list of tuples
        j = 0
        for row in query.rows():
            t = [row['symbol']]
            for v in views:
                t.append(row[v])
            z.append(tuple(t))
            j += 1

    symbols = []
    terms = []
    other = []
    for k in z:
        symbols.append(k[0])
        terms.append(k[ind+1])
        L = list(k)
        L.remove(k[0])
        L.remove(k[ind+1])
        k2 = [k[ind+1]]
        k2 += L
        other.append(tuple(k2))

    cols = [v.replace(".", "_") for v in views[0:ind] + views[(ind+1):]]

    cname = hostpd.columns[1]
    dfs = pd.DataFrame(zip(symbols, terms),
                       columns=[cname, name])

    results = hostpd.merge(dfs, left_on=cname, right_on=cname)
    results = results.drop(cname, 1)
    pd2sql(results, "ensemblg2%s$annot" % name)

    # add the new table to the database
    add_db_table(dbname,
                 "%s$details" % name,
                 ([name] + cols),
                 ["varchar(250)"] * (len(views)),
                 other, PS=True)


def parseOwlOntology(ontologypath, tablename):
    '''
    Parses an OWL format ontology to find the hierarchy of terms.
    Generates a database table, tablename with a list of terms in column
    1 and all the terms of which they are a direct subclass in column 2, with
    one term per row.
    For example, from the human phenotype ontology:
    Glomerulonephritis (HP:0000099) is a subclass of nephritis (HP:0000123)
    and of "Abnormality of the Glomerulus" (HP:0000095).
    This would be recorded as:
    id    subclass_of
    HP:0000099    HP:0000095
    HP:0000099    HP:0000123
    '''
    # download an up to date ontology file, parse the xml data into a
    # Python "ElementTree" and delete the ontology file.
    ontologyfile = P.getTempFilename(".")
    os.system("wget -O %s %s" % (ontologyfile, ontologypath))
    tree = ET.parse(ontologyfile)
    os.remove(ontologyfile)

    # Traverse the tree and identify the classes and subclasses
    # The ns dictionary links "prefixes" of xml tags to the full tag
    # names, see here for more info -
    # https://docs.python.org/2/library/xml.etree.elementtree.html
    ns = {'owl': 'http://www.w3.org/2002/07/owl#',
          'rdfs': 'http://www.w3.org/2000/01/rdf-schema#'}
    root = tree.getroot()
    r = dict()
    for aclass in root.findall('owl:Class', ns):
        nam = aclass.attrib.values()[0].split("/")[-1]
        for subclassof in aclass.findall('rdfs:subClassOf', ns):
            if nam not in r:
                r[nam] = set()
            if len(subclassof.attrib.values()) == 1:
                r[nam].add(subclassof.attrib.values()[0].split("/")[-1])

    klist = []
    vlist = []
    for item in r.items():
        key = item[0]
        vals = item[1]
        for v in vals:
            klist.append(key.replace("_", ":"))
            vlist.append(v.replace("_", ":"))

    add_db_table(dbname,
                 tablename,
                 ["id", "subclass_of"],
                 ["varchar(250)", "varchar(250)"],
                 zip(klist, vlist))


def getEntrezGeneIds():
    '''
    Gets all the Gene IDs for a particular host, specified in PARAMS, from
    Entrez Gene and returns them as a list.
    '''

    print 'retrieving entrez gene ids'
    host = PARAMS['entrez_host']

    # Limited to IDs which are current and not obsolete (i.e. "alive")
    term = '("alive"[Properties]) AND %s[Taxonomy ID]' % host

    E = Entrez.esearch(db="gene", term=term, retmode="text", retmax=1000000)
    res = Entrez.read(E)
    E.close()

    return res['IdList']


def getMyGeneInfo(allids):
    '''
    Gets information from MyGeneInfo.org for the list of genes pulled from
    the Entrez Gene database.  This information is combined into a big
    python dictionary.
    The fields to get from MyGeneInfo are specified in pipeline.ini, the list
    of all possible fields is at
    http://docs.mygene.info/en/latest/doc/data.html#available-fields
    Each MyGeneInfo 'object' is formatted slightly differently, depending on
    the field and host, so adding another field means writing a new function
    to parse the information retrieved (similar to the ParseAndLoad functions
    below).
    '''

    print 'retrieving data from mygeneinfo'
    fields = PARAMS['my_gene_info_annotations'].split(",")
    if 'ensembl' not in fields:
        fields.append('ensembl')
    DB = dict()
    # If the list of entrez ids is long, cut it into chunks of 500 ids to send
    # to the API in each query (this size seems to run the fastest overall).
    if len(allids) >= 500:
        ids = np.array(allids)
        a = len(ids) / 500
        idsets = np.array_split(ids, a)
    else:
        idsets = [allids]
    for idset in idsets:
        h = httplib2.Http()
        headers = {'content-type': 'application/x-www-form-urlencoded'}

        # get all the ids and fields specified
        params = 'ids=%s&fields=%s' % (",".join(idset), ",".join(fields))
        res, con = h.request('http://mygene.info/v2/gene', 'POST',
                             params, headers=headers)

        # convert JSON format to nested Python dictionaries.
        jcon = json.loads(con)
        for line in jcon:
            for f in line.keys():
                if f not in DB:
                    DB[f] = dict()
                DB[f][line['query']] = line[f]
    return DB


def ParseAndLoadEnsembl(DB):
    '''
    Parses and loads the Ensembl dataset retrieved from MyGeneInfo (DB)
    Always creates a database table, entrez2ensemblg, converting entrez gene
    ids to ensembl gene ids.
    If 'transcript' or 'protein' (or 'all') is specified in the pipeline.ini,
    generates ensemblg2ensemblt and ensemblg2ensemblp tables respectively,
    translating between different types of ensembl ID.
    Where there is not a 1-1 mapping, a row is created for each combination,
    e.g. if entrezgene '111' corresponds to both ensembl gene 'aaa' and
    ensembl gene 'bbb' the layout would be:
    entrez    ensemblg
    111       aaa
    111       bbb
    '''
    print 'parsing ensembl data'
    ensembldict = DB['ensembl']

    # options from pipeline.ini for which ensembl id types to retrieve
    options = PARAMS['my_gene_info_ensembl'].split(",")

    # store the relationships between various gene ids as lists of tuples
    entrez2ensembl = []
    ensembl2entrez = []
    g2t = []
    t2g = []
    g2p = []
    p2g = []

    # for each entrez ID in the dictionary of ensembl results from mygeneinfo
    for q in ensembldict:
        # each ID leads to another dictionary with keys 'gene', 'transcript'
        # and 'protein' containing ensembl ids of these types.
        res = ensembldict[q]
        # in the "ensembldict" created from the mygeneinfo output
        # 1:1 mappings have dictionary values stored as strings
        # and 1 to many mappings have dictionary values stored as
        # lists.  This makes a list with one element if the value is a
        # string, so that subsequent steps can be applied in both cases.
        if type(res) is not list:
            res = [res]
        for L in res:
            # L['gene'], L['transcript'] and L['protein'] are also a
            # mixture of lists and strings depending if there is one or more
            # than one value, most of the following is to deal with this
            if type(L['gene']) is list:
                entrez2ensembl += [q] * len(L['gene'])
                ensembl2entrez += L['gene']
                glist = L['gene']
            else:
                entrez2ensembl.append(q)
                ensembl2entrez.append(L['gene'])
                glist = [L['gene']]
            for gene in glist:
                if 'transcript' in L:
                    if type(L['transcript']) is list:
                        g2t += [gene] * len(L['transcript'])
                        t2g += L['transcript']
                    else:
                        g2t.append(gene)
                        t2g.append(L['transcript'])
                if 'protein' in L:
                    if type(L['protein']) is list:
                        g2p += [gene] * len(L['protein'])
                        p2g += L['protein']
                    else:
                        g2p.append(gene)
                        p2g.append(L['protein'])

    # load the entrez2ensemblg table
    add_db_table(dbname,
                 'ensemblg2entrez$geneid',
                 ['ensemblg', 'entrez'],
                 ['int', 'varchar(25)'],
                 zip(ensembl2entrez, entrez2ensembl))

    # load transcript IDs if requested
    if 'all' in options or 'transcript' in options:
        add_db_table(dbname,
                     'ensemblg2ensemblt$other',
                     ['ensemblg', 'ensemblt'],
                     ['varchar(25)', 'varchar(25)'],
                     zip(g2t, t2g))

    # load protein IDs if requested
    if 'all' in options or 'protein' in options:
        add_db_table(dbname,
                     'ensemblg2ensemblp$other',
                     ['ensemblg', 'ensemblp'],
                     ['varchar(25)', 'varchar(25)'],
                     zip(g2p, p2g))
    ensdf = pd.DataFrame(zip(ensembl2entrez, entrez2ensembl),
                         columns=['ensemblg', 'entrez'])
    return ensdf


def ParseAndLoadGO(DB, ensdf):
    '''
    Parses and loads the GO dataset retrived from MyGeneInfo (in DB).
    Which GO namespace to store (any combination of BP (biological process),
    CC (cellular component), or MF (molecular function) or 'all' can be
    specified in the pipeline.ini.
    Generates two database tables:
    entrez2go maps entrez IDs to GO IDs
    goid2goinfo contains further information on this GO id - term name,
    description, evidence and namespace.
    '''

    print 'parsing go data'
    if 'go' not in DB:
        print "GO annotations not available for this host"
        return 0
    gores = DB['go']
    godict = dict()
    # options from pipeline.ini for which go namespaces are of interest
    options = PARAMS['my_gene_info_go'].split(",")
    goids = []
    goset = set()
    allgoids = []
    entrez2go = []
    gotypes = []
    # for each entrez ID in the go dictionary
    for q in gores:
        # res is a dictionary with keys BP, CC and MF
        res = gores[q]
        for cat in ['BP', 'MF', 'CC']:
            if cat in res and (cat in options or 'all' in options):
                # As for ensembl above, this deals with the mixture of
                # lists and strings as dictionary values returned by
                # mygeneinfo
                if type(res[cat]) is not list:
                    res[cat] = [res[cat]]
                for L in res[cat]:
                    allgoids.append(L['id'])
                    entrez2go.append(q)
                    if L['id'] not in goset:
                        goset = goset | set([L['id']])
                        goids.append(L['id'])
                        gotypes.append(cat)
                        for subterm in ['term', 'evidence']:
                            if subterm not in godict:
                                godict[subterm] = []
                            if subterm in L:
                                godict[subterm].append(L[subterm])
                            else:
                                godict[subterm].append("")
    dfs = pd.DataFrame(zip(entrez2go, allgoids), columns=['entrez', 'go'])
    results = ensdf.merge(dfs, left_on='entrez', right_on='entrez')
    results = results.drop("entrez", 1)
    pd2sql(results, "ensemblg2go$annot")

    # load the go id to go details table
    add_db_table(dbname,
                 'go$details',
                 ['go', 'gotype', 'goterm', 'goevidence'],
                 ['varchar(25)', 'varchar(25)', 'varchar(255)', 'varchar(25)'],
                 zip(goids, gotypes, godict['term'], godict['evidence']))


def ParseAndLoadPathway(DB, ensdf):
    '''
    Parses the "pathway" dictionary extracted from MyGeneInfo above.
    The pathway annotation sets of interest are specified in pipeline.ini, or
    if 'all' is specified all pathway annotations available in the database
    are used.
    Each annotation set is used to generate 2 database tables:
    entrez2xxx - translates entrez ids to annotation ids
    xxx - annotation ids and the pathways they correspond to
    '''
    print 'parsing pathway data'
    if 'pathway' not in DB:
        print "pathway annotations not available for this species"
        return 0
    pathwaydict = DB['pathway']
    D = dict()

    # If specific annotation sets to use are provided in pipeline.ini use those
    # otherwise get all the keys from the dictionaries in pathwaydict.
    typelist = PARAMS['my_gene_info_pathway'].split(",")

    if 'all' in typelist:
        allkeys = set()
        for k in DB['pathway'].keys():
            allkeys = allkeys | set(DB['pathway'][k].keys())
        typelist = list(allkeys)

    # Make empty lists to store the lists of tuples to load into the database
    # The D dictionary contains four lists per host, named host_ids2entrez,
    # host_entrez2ids, host_id and host_name.
    for db in typelist:
        D["%s_ids2entrez" % db] = []
        D["%s_entrez2ids" % db] = []
        D["%s_id" % db] = []
        D["%s_name" % db] = []

    # for each entrez id in the pathway dictionary
    for q in pathwaydict.keys():
        # for each of the requested annotation types
        for db in typelist:
            idset = set()
            if db in pathwaydict[q]:
                results = pathwaydict[q][db]
                if type(results) is not list:
                    results = [results]
                for b in results:
                    D['%s_ids2entrez' % db].append(b['id'])
                    D["%s_entrez2ids" % db] .append(q)
                    if b['id'] not in idset:
                        D["%s_id" % db].append(b['id'])
                        D["%s_name" % db].append(b['name'])
                        idset.add(b['id'])

    # put these into the database
    for db in typelist:
        ids2entrez = D["%s_ids2entrez" % db]
        entrez2ids = D["%s_entrez2ids" % db]
        id = D["%s_id" % db]
        name = D["%s_name" % db]
        identrez = zip(entrez2ids, ids2entrez)
        idname = zip(id, name)
        dfs = pd.DataFrame(identrez, columns=['entrez', db])
        results = ensdf.merge(dfs, left_on='entrez', right_on='entrez')
        results = results.drop("entrez", 1)
        pd2sql(results, "ensemblg2%s$annot" % db)
        add_db_table(dbname,
                     "%s$details" % db,
                     [db, 'term'],
                     ['varchar(25)', 'varchar(25)'],
                     idname)


def ParseAndLoadHomologene(DB, ensdf):
    '''
    Uses information from the homologene dataset retrieved from MyGeneInfo
    to link entrez IDs from the host organism to those of homologous genes
    in other organisms of interest.
    The organisms of interest can be specified in pipeline.ini or, if 'all',
    is specified, every available organism is used (~25 organisms).
    A summary table "species_ids" is generated showing which
    host taxonomy id corresponds to which host scientific name.
    For each taxon a table is generated, taxonname_entrez, linking host entrez
    ids to those of the taxon.
    Entrez ids for each taxon are translated into gene symbols, these are
    stored in tables as entrez_taxonname_symbol.
    '''
    print 'parsing homologene data'
    if 'homologene' not in DB:
        print "homologene not available for this host"
        return 0
    homodb = DB['homologene']
    taxa = str(PARAMS['my_gene_info_homologene']).split(",")
    if '9606' not in taxa:
        taxa.append('9606')

    # retrieve IDs for all the taxa listed in PARAMS["my_gene_info_homologene"]
    # or if this is "all" check which taxa are in the input
    if 'all' in taxa:
        allspp = set()
        for q in homodb.keys():
            geneset = homodb[q]
            if type(geneset) is not list:
                geneset = [geneset]
            for L in geneset:
                spp = [str(line[0]) for line in L['genes']]
                allspp = allspp | set(spp)
        allspp = list(allspp)
    else:
        allspp = taxa

    # retrieve species names of the taxa from entrez taxonomy
    asp = ",".join(allspp)
    E = Entrez.efetch(db="taxonomy", id=asp, retmode="xml", retmax=1000000)
    res = Entrez.read(E)
    names = ["_".join(r['ScientificName'].split(" ")) for r in res]
    E.close()

    if str(ohost) not in taxa:
        taxa.append(ohost)
        names.append(PARAMS['entrez_sciname'])

    # load species names into the species_ids database table
    # add_db_table(dbname,
    #              "species_ids",
    #              ['id', 'species'],
    #              ['varchar(25)', 'varchar(250)'],
    #              zip(allspp, names))

    # parse the list of entrez ids in each host

    hostdict = dict()
    for host in names:
        hostdict[host] = []
        hostdict["%s_entrez" % host] = []
    for q in homodb.keys():
        for item in homodb[q]['genes']:
            hostid = str(item[0])
            if hostid in allspp:
                hostgeneid = str(item[1])
                host = names[allspp.index(hostid)]
                hostdict[host].append(hostgeneid)
                hostdict["%s_entrez" % host].append(q)
    pdd = dict()
    # load database tables
    for host in names:
        hh = host.split("_")
        hid = "%s%s" % (hh[0][0], hh[1][0:3])
        ids = host2symbol(hostdict[host])
        df1 = pd.DataFrame(zip(hostdict['%s_entrez' % host], hostdict[host]),
                           columns=['entrez', '%s_entrez' % hid])
        df2 = pd.DataFrame(zip(ids.keys(), ids.values()),
                           columns=['%s_entrez' % hid, '%s_symbol' % hid])
        dfs = df1.merge(df2, left_on='%s_entrez' % hid,
                        right_on='%s_entrez' % hid)
        print dfs.columns
        dfs = dfs.drop("%s_entrez" % hid, 1)

        results = ensdf.merge(dfs, left_on='entrez', right_on='entrez')
        results = results.drop("entrez", 1)
        pd2sql(results, "ensemblg2%s_symbol$other" % hid)
        pdd[host] = results

    return hostdict, pdd


def addMousePhenotypes(hostdict, pdd):
    '''
    Using the mineOntology function above, adds data from the mouse
    phenotype ontology to the database.
    Creates a database table, mouse_phenotype, showing gene symbols in
    mouse and the corresponding phenotypes.
    '''
    print 'adding mouse phenotypes'
    assert "Mus_musculus" in hostdict, "Can't add mouse phenotypes without \
    mining homologene for mouse"
    constraints = [
        ("Gene.ontologyAnnotations.ontologyTerm.namespace=MPheno.ontology")]
    views = ["ontologyAnnotations.ontologyTerm.description",
             "ontologyAnnotations.ontologyTerm.name",
             "ontologyAnnotations.ontologyTerm.namespace",
             "ontologyAnnotations.ontologyTerm.identifier"]
    m = hostdict["Mus_musculus"]
    mousepd = pdd['Mus_musculus']
    mineOntology('mouse',
                 m,
                 views,
                 constraints,
                 "M. musculus",
                 "MPI",
                 mousepd, 3)


def addMousePathways(hostdict, pdd):
    '''
    Using the mineOntology function above, adds data about mouse
    pathways to the database.
    Creates a database table, mouse_pathway, showing gene symbols in
    mouse and the corresponding pathways.
    '''
    print 'adding mouse pathways'
    assert "Mus_musculus" in hostdict, "Can't add mouse pathways without \
    mining homologene for mouse"
    constraints = []
    m = hostdict["Mus_musculus"]
    views = ["pathways.identifier", "pathways.name"]
    mousepd = pdd['Mus_musculus']
    mineOntology('mouse',
                 m,
                 views,
                 constraints,
                 "M. musculus",
                 "mousepathway",
                 mousepd, 0)


def addHumanPhenotypes(hostdict, pdd):
    '''
    Using the mineOntology function above, adds data from the human
    phenotype ontology to the database.
    Creates a database table, hpo, showing gene symbols in
    human and the corresponding phenotypes.
    '''

    print 'adding human phenotypes'
    assert "Homo_sapiens" in hostdict, "Can't add human phenotypes without \
    mining homologene for human"
    constraints = []
    h = hostdict["Homo_sapiens"]
    humanpd = pdd['Homo_sapiens']
    views = ["diseases.hpoAnnotations.hpoTerm.description",
             "diseases.hpoAnnotations.hpoTerm.name",
             "diseases.hpoAnnotations.hpoTerm.namespace",
             "diseases.hpoAnnotations.hpoTerm.identifier"]
    mineOntology('human',
                 h,
                 views,
                 constraints,
                 "",
                 "hpo",
                 humanpd, 3)


def main():
    t = testAPIs()
    if t == 0:
        print "Failed - Check APIs"
        return 0
    ids = getEntrezGeneIds()
    print "done"

    mgi = getMyGeneInfo(ids)
    print "Loading ensembl ids"
    ensdf = ParseAndLoadEnsembl(mgi)
    print "done"
    symbols = host2symbol(ids)
    dfs = pd.DataFrame(zip(symbols.keys(), symbols.values()),
                       columns=['entrez', 'symbol'])
    results = ensdf.merge(dfs, left_on='entrez', right_on='entrez')
    results = results.drop("entrez", 1)
    pd2sql(results, "ensemblg2symbol$geneid")
    annots = PARAMS['my_gene_info_annotations'].split(",")
    if 'go' in annots and (len(PARAMS['my_gene_info_go']) != 0):
        ParseAndLoadGO(mgi, ensdf)
        parseOwlOntology(PARAMS['my_gene_info_goont'], "go$ont")
        print "done"
    if 'pathway' in annots and len(PARAMS['my_gene_info_pathway']) != 0:
        ParseAndLoadPathway(mgi, ensdf)
        print "done"
    if 'homologene' in annots and PARAMS['my_gene_info_homologene'] != '':
        D, pdd = ParseAndLoadHomologene(mgi, ensdf)
        print "done"
    try:
        D.values()
    except:
        return 0

    if int(PARAMS['homologues_mousephenotype']) == 1:
        addMousePhenotypes(D, pdd)
        print "done"
    if int(PARAMS['homologues_mousepathway']) == 1:
        addMousePathways(D, pdd)
        print "done"
    if int(PARAMS['homologues_hpo']) == 1:
        addHumanPhenotypes(D, pdd)
        parseOwlOntology(PARAMS['homologues_hpoont'], "hpo$ont")
        print "done"

if __name__ == "__main__":
    main()
