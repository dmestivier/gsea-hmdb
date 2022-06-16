# parse un fichier XML de la HMDB
#
# Je veux extraire les elements + les proteines associees
#

import sys
import gzip

################################################################
###
### compte le nombre d'espace en debut de ligne
###
################################################################

def nlevel( s ):
    n1 = len( s )
    n2 = len( s.lstrip() )
    return n1-n2

#################################################
###
### Recupere l'info entre deux balises
### On peut avoir des info contenant des ';' !!!
### que l'on remplace par "_"
###
#################################################

def getval( m, tag ):
    m = m.replace("<"+tag+">", "" )
    m = m.replace("</"+tag+">", "" )
    m = m.replace(";", "_" )
    return m.strip()

#################################################
###
### Print metabolite
###
#################################################

def print_metabolite_stdout( m, lt ):
    tmp = m["accession"] 
    # je ne prend qu'une accession
    tmp = tmp.split(";")
    tmp = tmp[0]
    
    oo  = tmp 
    sep = "\t"

    for tag in lt:
        if tag not in m.keys():
            oo += sep + "NONE"
        else:
            if len(m[tag])>0:
                if m[tag][0] == ";":
                    oo += sep + m[tag][1:]
                else:
                    oo += sep + m[tag]

    sys.stdout.write( oo + "\n" )

##########################################################
###
### Handle tag
###
##########################################################

def handleTag( m, nlvl, lig, tg, nsp, adding ):
    if lig.find ("<" + tg + ">" )>=0:
        if tg not in m.keys():
            m[tg] = ""

        if nlvl == nsp:
            if adding==True:
                m[ tg ] += ";" + getval( lig[:-1], tg )
            else:
                m[ tg ] = getval( lig[:-1], tg )
    else:
        # on peut alors le cas vide : <TAG/>
        if lig.find("<" + tg + "/>")>=0:
            if tg not in m.keys():
                m[tg] = ""

            if nlvl == nsp:
                if adding==True:
                    m[tg] += ";" + "NOT-PROVIDED"
                else:
                    m[tg] += "NONE:" + m["accession"]
    return m
        
##################################################
###
### Affiche sous forme d'interaction 
###   metabolite - gene
###
##################################################

def print_interaction_stderr( met, lt ):
    # test si les deux infos sont là

    accession = met["accession"]
    metaboName = met["name"]

    if "gene_name" in met.keys():
        # commence toujours par un ";" ( corriger plus tard)
        # et donc une case vide quand on splitte si on ne 
        # commence par en `1:`
        listOfGenes = met["gene_name"][1:].split(";")
    else:
        listOfGenes = []

    if "uniprot_id" in met.keys():
        # commence toujours par un ";" ( corriger plus tard)
        # et donc une case vide quand on splitte si on ne 
        # commence par en `1:`
        listOfProt = met["uniprot_id"][1:].split(";")
    else:
        listOfProt = []

    if len( listOfProt )==len( listOfGenes ):
        for i in range(0,len(listOfGenes) ):
            sys.stderr.write( accession + "\t" + metaboName + "\t" + listOfGenes[i] + "\t" + listOfProt[i] + "\n" )
    else:
        sys.stderr.write("Error: pas le meme nombre de uniprot que de genes\n" )
        sys.stderr.write(met)
        sys.stderr.write("\n")
        sys.exit(1)

#################################################
###
### MAIN
###
#################################################

###
### ligne de commande
###

if len( sys.argv )!=2:
    print( "Syntax: %s file.xml.gz" % sys.argv[0] )
    sys.exit(1)

#fn = "uu.xml"
fn = sys.argv[1]

###
### liste des tags a rechercher
###

listOfTags = [ "accession", "name", "synonym", "kegg_id", "chebi_id", "protein_accession", "gene_name", "uniprot_id"]

###
### Output/header
###

hh = "metabolite"
sep = "\t"
for tag in listOfTags:
    hh += sep + tag
sys.stdout.write( hh + "\n" )
sys.stderr.write( "HMDBID\tHMDBname\tgenename\tuniprot\n" )

met = {}
for lig in gzip.open( fn, "rt" ):

    # quel niveau (=nb espace en debut)
    nident = nlevel(lig)

    # chaque item
    met = handleTag( met, nident, lig, "accession", 2, False )
    #met = handleTag( met, nident, lig, "accession", 4, True ) # les secondary_accessions
    met = handleTag( met, nident, lig, "name", 2, False )
    met = handleTag( met, nident, lig, "synonym", 4, True )
    met = handleTag( met, nident, lig, "kegg_id", 2, False )
    met = handleTag( met, nident, lig, "chebi_id", 2, False )
    met = handleTag( met, nident, lig, "protein_accession", 6, True )
    met = handleTag( met, nident, lig, "uniprot_id", 6, True )
    met = handleTag( met, nident, lig, "gene_name", 6, True )
   
    # fin de description d'un metabo. Affiche et remet a 0
    if lig.find("</metabolite>")>=0:
        print_metabolite_stdout( met, listOfTags) 
        print_interaction_stderr( met, listOfTags )
        met = {}

