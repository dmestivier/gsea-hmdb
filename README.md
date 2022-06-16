# gsea-hmdb

We develop a Metabolomic Signature for GSEA analysis using data from the HMDB
(Human Metabolome Database).

For the impatient, our Metabolomic Signature can be downloaded in the 
`./GMT` directory. Be carefull to download also the class description
file (see below) :

- hmdb_metabolites-filtered-merged.gmt
- hmdb_metabolites-filtered-dico.csv

-----------------------------------------------------------

## Reminder

Source: https://www.gsea-msigdb.org/gsea/index.jsp:

"Gene Set Enrichment Analysis (GSEA) is a computational method that determines whether 
an a priori defined set of genes shows statistically significant, concordant differences 
between two biological states (e.g. phenotypes)."

In case of a metabolomic signature, each metabolite is associated with a set
of genes (Gene Set, GS).

This Metabolomic Signature is provided as a GMT file in order to be
used with classical pipelines, in this case the GSEA software.

------------------------------------------------------------

- gmt: gene sets file
- gene set name => metabolite
- for each gene set :
	- genese name = metabolite
	- list of genes in that gene set
	- using HGNC gene symbol.
- the GSEA team recommends that the file name include the 
		  gene identifier format you used to list the genes

-----------------------------------------------------------

## Metabolomics Signature creation

If you want to build from scratch our Metabolomic Signature, 
please follow these steps :

### Clone this repository :

```bash
$ git@github.com:dmestivier/gsea-hmdb.git
$ cd gsea-hmdb
```

### Download the latest XML file from HMDB :

Be aware that the XML file is 6Go uncompressed.

```bash
$ wget https://hmdb.ca/system/downloads/current/hmdb_metabolites.zip
$ unzip hmdb_metabolites.zip
$ mv hmdb_metabolites.xml  ./Data/.
$ mv hmdb_metabolites.xml Data/.
$ gzip Data/hmdb_metabolites.xml
$ rm hmdb_metabolites.zip
```

### Extract Metabolite-gene associations :

We parse the XML file (compressed) in order to extract every interactions
recorded between a metabolite and genes. This script is very slow... be patient.

```bash
$ python ./Scripts/parsor.py ./Data/hmdb_metabolites.xml.gz > ./HMDB/hmdb_metabolites_metadata.csv 2> ./HMDB/hmdb_metabolites_interactions.csv
$ gzip ./HMDB/hmdb_metabolites_interactions.csv
$ gzip ./HMDB/hmdb_metabolites_metadata.csv 
```

Note :
- Metabolites with no gene associations are excluded.
- files created are (zip-compressed) TEXT format files with tab-separated columns 
  and could be imported in every software (after decompression).
- `hmdb_metabolites_metadata.csv` list metabolites for which
  metabolite-gene associations have been found. It reports some common identifiers
  (ex: keggid) for further processing/validation if needed.
- `hmdb_metabolites_interactions.csv` list every metabolite-gene interactions.


### Create Genesets (GS) from Metabolite-Genes associations

We can now merge all genes associated to the same metabolite into a gene set (GS) :

```bash
$ python ./Scripts/interactions2gmt.py ./HMDB/hmdb_metabolites_interactions.csv.gz > ./GMT/hmdb_metabolite_raw.symbols.gmt

```

### Filtering steps

Two steps :

1. remove duplicated gene name in the same (GS) (due to different uniprot identifiers not taken into account here)
2. remove GS with less than 10 genes.

```bash
$ python ./Scripts/geneset-filtering.py ./GMT/hmdb_metabolite_raw.symbols.gmt
$ wc -l ./GMT/*.gmt
19178 ./GMT/hmdb_metabolite_raw.symbols-filtered.gmt
22849 ./GMT/hmdb_metabolite_raw.symbols.gmt
```

A new file `hmdb_metabolite_raw.symbols-filtered.gmt` is automatically generated.

### Merge redundant GS

Some metabolites are associated with exaclty the *same* list of genes.
These GS are merged into a class. Because it is very difficule to infer a class name
form the names of the metabolites, we decide to use a generic class name
such as : `GSxxx_number-of-metabolites_size-of-the-geneset` 

For example `GS445_M13621_G34` represent a class of 13621 metabolites 
are described by the same geneset of 34 genes. These is a familly 
of triradyglycerols. One member is `TG(14:0/14:0/14:0)`.

A dictionnary file (TEXT tab-sep) is provided with the members of all
the class generated.

```bash
$ python ./Scripts/geneset-merge-identical.py ./GMT/hmdb_metabolite_raw.symbols-filtered.gmt
```

Two files are automatically generated :

- `./GMT/hmdb_metabolite_raw.symbols-filtered-dico.csv` : description of the classes when redundant GS
- `./GMT/hmdb_metabolite_raw.symbols-filtered.gmt` : the GMT file.

Rename :

```bash
$ mv hmdb_metabolite_raw.symbols-filtered-merged.gmt hmdb_metabolites.symbols.gmt
$ mv hmdb_metabolite_raw.symbols-filtered-dico.csv hmdb_metabolites-dico.csv
```

----------------------------------------------------------------------

## Some statistics

How many metabolites

```bash
$ zcat ./HMDB/hmdb_metabolites_metadata.csv.gz | awk '{print $1}' | sort | wc -l
217921
```

How many genes without uniprot identifiers:

	- combien de genes sans *uniprot* ?
```bash
$ zcat ./HMDB/hmdb_metabolites_metadata.csv.gz | grep "NOT-PROVIDED" | wc -l
3344
```

Number of metabolites-gene associations:
	- 
```
$ zcat ./HMDB/hmdb_metabolites_interactions.csv.gz | wc -l
863760
```

```bash
$ wc -l ./GMT/hmdb_metabolite*
  19178 ./GMT/hmdb_metabolite_raw.symbols-filtered.gmt
  22849 ./GMT/hmdb_metabolite_raw.symbols.gmt
  18635 ./GMT/hmdb_metabolites-dico.csv
    604 ./GMT/hmdb_metabolites.symbols.gmt

dmestivier@labo ~/Travail/gsea-metabolites/GMT (master) $ wc -l hmdb_metabolites*
 18646 hmdb_metabolites-filtered-dico.csv
 19178 hmdb_metabolites-filtered.gmt
   571 hmdb_metabolites-filtered-merged.gmt
 22849 hmdb_metabolites.gmt
```

So starting with 863.760 metabolite-gene interactions for 217.921 metabolites,
we create 22.849 GS: 3671 GS had a size less than 10 genes, and we keep
19.178 GS.

From these 19.178 GS, 18.574 GS ($96.85\%$) are "redundant" meaning that the same GS is
associated to a different identifier (metabolite). 604 genesets are retained.

Let's draw an histogram of the size...

```bash
$ python ./Scripts/geneset-stats.py ./GMT/hmdb_metabolites.symbols.gmt > genes.csv 2> GS.csv
```

... using R :

```R
> gs = read.csv("GS.csv", h=T, sep="\t" )
> genes = read.csv("genes.csv", h=T, sep="\t" )
> dim(gs)
[1] 604   2
> dim(genes)
[1] 5769    3
> hist( gs$size, xlab="GS size", ylab="#Occ", col="blue" )
```

The metabolites associated with the largest GS :

```R
> gs[ order( gs$size, decreasing=T )[1:10], ]
                    geneset size
101                   Water 1137
102                 Calcium 1080
103  Adenosine triphosphate 1027
104                     ADP  781
105               Phosphate  777
108  Guanosine triphosphate  557
109               Magnesium  463
110           Pyrophosphate  335
111 Adenosine monophosphate  265
112                     NAD  258
```

Be aware that GSEA sofware use a default upper limit of 500 for the genesets size.

Genes associated in many GS :

```R
> genes[ order( genes$NbGenesets, decreasing=T )[1:10], 1:2]
    genename NbGenesets
474   CYP3A4        168
466   CYP2D6        120
470   CYP2C9        118
472   CYP1A2        109
461  CYP2C19        103
467   CYP2C8         94
682    ABCB1         90
460   CYP2B6         86
482   CYP3A5         85
463   CYP2E1         79
```



