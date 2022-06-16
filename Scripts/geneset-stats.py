###
### Read GMT file and list of genes in
### the gene sets.
### as well as statitics for each genesets
###
### Mestivier Denis - June 2022
###

import sys

###
### Command line 
###

if len( sys.argv )!= 2:
    print( "Syntax: %s file.gmt" % sys.argv[0] )
    sys.exit(1)

###
### Import GMT 
###

d = {}
for lig in open( sys.argv[1], "rt" ):
    cols = lig[:-1].split("\t")
    gs   = cols[0]
    d[ gs ] = {}
    d[ gs ]["link"] = cols[1]
    d[ gs ]["GS"]   = cols[2:]

###
### List of genes
###

listOfGenes = []

for gs in d.keys():
    for gene in d[gs]["GS"]:
        if gene not in listOfGenes:
            listOfGenes.append( gene )

###
### report stats at gene level
###

sep  = "\t"

sys.stdout.write("genename\tNbGenesets\tlistOfGenesets\n" )

for gene in listOfGenes:
    # list of all geneset to which 'gene' belong to
    listOfM = []
  
    # inspect each geneset
    for gs in d.keys():
        if gene in d[gs]["GS"]:
            listOfM.append( gs )

    # output
    oo  = gene
    oo += sep + str( len(listOfM) )
    oo += sep + ";".join(listOfM)
    sys.stdout.write(oo + "\n" )

###
### report stats at geneset level
###

sys.stderr.write("geneset\tsize\n" )
for gs in d.keys():
    oo  = gs + sep
    oo += str( len(d[gs]["GS"]) )
    sys.stderr.write(oo + "\n" )
