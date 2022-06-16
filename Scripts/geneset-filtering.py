###
### Read GMT file.
### Implement filtering
###
### Rappel: 
###
###  "1-Methylhistidine" https://hmdb.ca/metabolites/HMDB0000001 CNDP1   PRMT3
###
### Mestivier Denis - Mars 2022
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
### Filtering
###

THR_SIZE=10

# remove "NOT-PROVIDED" genenames
tagOFF = "NOT-PROVIDED"

for gs in d.keys():
    if tagOFF in d[gs]["GS"]:
        # remove
        d[gs]["GS"] = [ gene for gene in d[gs]["GS"] if gene != tagOFF ]

# remove duplicated genenames
for gs in d.keys():
    d[gs]["NEWGS"] = list( set( d[gs]["GS"] ) )

# filter two small GS (THR=10)
for gs in d.keys():
    if len(d[gs]["NEWGS"]) < THR_SIZE:
        d[ gs ]["STATUS"] = "remove"
        d[ gs ]["NEWGS" ] = []
    else:
        d[ gs ]["STATUS"] = "keep"

# Save new GMT
fnout = sys.argv[1]
fnout = fnout.replace( ".gmt", "" )
fnout = fnout + "-filtered.gmt"

fid = open( fnout, "wt" )
sep = "\t"

for gs in d.keys():
    if d[gs]["NEWGS"] != []:
        oo  = gs
        oo += sep + d[gs]["link"]
        oo += sep + "\t".join( d[gs]["NEWGS"] )
        fid.write( oo + "\n" )

fid.close()
