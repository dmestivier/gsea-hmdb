###
### Lit un fichier GMT et calcule un overlap entre deux GS
### On agrege les geneset identiques (ie meme liste de genes)
### en une classe.
### On produit egalement un fichier annexe qui associe le
###  nom de la classe avec la liste de ses elements (hmdb)
###
### Mestivier Denis - Avril 2022
###

import sys
        
#####################################################
###
### Construit une table de contingence
###
# > getContbl(go.obj)
#        not1  in1
#   not2 1647 10865
#   in2  6101  2583
###
#####################################################

def get_crosstab( L1, L2 ):
    # list des items
    LL = set(L1).union( L2)
    
    c11 = sum( [ g not in L1 and g not in L2 for g in LL ] )
    c12 = sum( [ g     in L1 and g not in L2 for g in LL ] )
    c21 = sum( [ g not in L1 and g     in L2 for g in LL ] )
    c22 = sum( [ g     in L1 and g     in L2 for g in LL ] )

    return [ [  c11, c12 ],  [ c21, c22 ] ]

###

##################################################
###
### Compute Jaccard similarity index
### URL: https://datascienceparichay.com/article/jaccard-similarity-python
###
##################################################
        
def get_jaccard_sim( a, b ):
    # convert to set
    a = set(a)
    b = set(b)
    # calcucate jaccard similarity
    j = float(len(a.intersection(b))) / len(a.union(b))
    return j
    

#################################################################
###
### Teste si une liste de gene est identique a une liste
### deja connue.
###
##################################################################

def is_strictly_overlapping( gs2name, gs2list, d, THR_POVER=99.9 ):
    noverlap = 0
    poverlap1 = 0
    poverlap2 = 0

    for gs1 in d.keys():
        noverlap = len( list( set( d[gs1]["GS"] ) & set( gs2list ) ) )
        poverlap1 = noverlap*100.0/len( d[gs1]["GS"] )
        poverlap2 = noverlap*100.0/len( gs2list  )
        
        #if poverlap2 > THR_POVER: 
        #return gs1, noverlap, poverlap1, poverlap2

        # Si les deux listes sont identiques, agrege
        if noverlap == len( d[gs1]["GS"] ) and noverlap == len( gs2list ):
            return gs1, noverlap, poverlap1, poverlap2

    return None, noverlap, poverlap1, poverlap2

##########################################################
###
### Cree une nouvelle entree
###
##########################################################

def geneset_new( d, gid, gsname, gslink, gslist ):
    gsnew = "GS" + str( gid )
    d[ gsnew ]          = {}
    d[ gsnew ]["link"]  = [ gslink ]
    d[ gsnew ]["GS"]    = gslist
    d[ gsnew ]["names"] = [ gsname ] 
    
    return d

##########################################################
###
### Clusterise un geneset dans un deja present
###
##########################################################

def geneset_merge( d, gsid, gsname, gslink, gslist ):
    d[gsid]["link" ].append( gslink )
    d[gsid]["names"].append(gsname)
    d[gsid]["GS"] = list( set( d[gsid]["GS"] + gslist )  )
    return d

###
### Print geneset
###

def geneset_print( d, gs, sep="\t" ):
    oo  = gs
    oo += sep + str( len( d[gs]["names"] ))
    oo += "\n"
    for n in d[gs]["names"]:
        oo += sep + n + "\n" 
    #oo += sep + ";".join(d[gs]["names"]) 
    oo += sep + "Geneset:" + sep + ";".join(d[gs]["GS"])
    oo += sep + "sizeGS:" + sep + str(len(d[gs]["GS"]))
    print(oo)

##########################################################
###
### Ligne de commande
###
##########################################################

if len( sys.argv )!= 2:
    print( "Syntax: %s file.gmt" % sys.argv[0] )
    sys.exit(1)

###
### Import GMT 
###

d = {}
L = []

for lig in open( sys.argv[1], "rt" ):
    cols = lig[:-1].split("\t")
    gsname   = cols[0]
    gslink   = cols[1]
    gslist   = cols[2:]

    d[ gsname ]          = {}
    d[ gsname ]["link"]  = gslink
    d[ gsname ]["GS"]    = gslist

    L.append( [ len(gslist), gsname ] )

L.sort(reverse=False)

###
### re-order. Start with smallest genesets
###

gmt = {}
for s, gsname  in L:
    gmt[ gsname ]          = {}
    gmt[ gsname ]["link"]  = d[ gsname ]["link"]
    gmt[ gsname ]["GS"]    = d[ gsname ]["GS"]

###
### 'clustering'
###

# Step 1 = merge les genesets strictement equivalent
#          i.e. meme liste de genes

gsid = 1
d = {}

for gsname in gmt.keys():
    gslink   = gmt[gsname]["link"] 
    gslist   = gmt[gsname]["GS"] 

    # est-ce un nouveau geneset qui a un equivalent
    # (ie un autre geneset avec le meme ensemble de genes)
    st, no, po1, po2 = is_strictly_overlapping( gsname, gslist, d )
    if st != None:
        sys.stderr.write("Merging:\t" + gsname + "\t" + str(len(gslist)) + "\tinto\t" + str(st) + "\t" + str(len(d[st]["GS"] )) + "\t:\t" + str(no) + "\t" + str(po1) + "\t" + str(po2) + "\n" )
        d = geneset_merge( d, st, gsname, gslink, gslist )
    else:
        sys.stderr.write("Create:\t" + gsname + "\t" + str(len(gslist)) + "\tinto\tGS" + str( gsid ) + "\n" )
        d    = geneset_new( d, gsid, gsname, gslink, gslist )
        gsid = gsid + 1

# Sort by decreasing size (ie. occurences)
L = [ [ len(d[gs]["names"]), gs ] for gs in d.keys() ]
L.sort(reverse=True)

###
### output.
###

fnoutGMT = sys.argv[1].replace(".gmt", "") + "-merged.gmt"
fnoutDIC = sys.argv[1].replace(".gmt", "") + "-dico.csv" 

outGMT = open( fnoutGMT, "wt" )
outDIC = open( fnoutDIC, "wt" )
outDIC.write("#ClassID\tHMDBName\n" )

sep = "\t"

for s, gs in L:
    nelts = len( d[gs]["names"] )
    newgs = gs + "_M" + str(nelts) + "_G" + str( len(d[gs]["GS"]) )
    # Nouveau GMT
    if nelts>1:
        oo  = newgs 
        oo += sep + "SEE_DICO_FILES_ASSOCIATED"
    else:
        oo  = d[gs]["names"][0]
        oo += sep + d[gs]["link"][0]
    
    oo += sep + "\t".join(d[gs]["GS"])
    outGMT.write(oo + "\n" )

    # Dico : equivalence nomClass <-> nomHMDB
    #        si necessaire
    if nelts>1:
        for n in d[gs]["names"]:
            oo = newgs + sep + n
            outDIC.write( oo + "\n" )

outDIC.close()
outGMT.close()
