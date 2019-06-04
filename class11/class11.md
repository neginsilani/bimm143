Class11: structural Bioinformatics
================
Negin Silani
5/7/2019

\#\#The PDBB database

The [PDB](http://www.rcsb.org/) is the main repository for biomolecular
structural data.

Here we examine the contents of the PDB:

Q1: Download a CSV file from the PDB site (accessible from “Analyze” -\>
“PDB Statistics” \> “by Experimental Method and Molecular Type”. Move
this CSV file into your RStudio project and determine the percentage of
structures solved by X-Ray and Electron Microscopy. From the website
what proportion of structures are protein? Aim to have a rendered GitHub
document with working code that yields your answers

``` r
db<- read.csv("Data Export Summary.csv", row.names=1)
head(db)
```

    ##                     Proteins Nucleic.Acids Protein.NA.Complex Other  Total
    ## X-Ray                 126880          2012               6547     8 135447
    ## NMR                    11062          1279                259     8  12608
    ## Electron Microscopy     2277            31                800     0   3108
    ## Other                    256             4                  6    13    279
    ## Multi Method             129             5                  2     1    137

How many are X-Ray, etc..

``` r
(db$Total/sum(db$Total)) * 100
```

    ## [1] 89.35736481  8.31777489  2.05041595  0.18406244  0.09038191

What percent are proteins..

``` r
(sum(db$Proteins)/sum(db$Total)) * 100
```

    ## [1] 92.75955

we could also try the datapasta package and copy from website and
“Addins” \> "Paste as data.frame..

``` r
library(datapasta)
        
tmp<- data.frame(stringsAsFactors=FALSE,
                     Experimental.Method = c("X-Ray", "NMR", "Electron Microscopy", "Other",
                                             "Multi Method", "Total"),
                                Proteins = c(126880, 11062, 2277, 256, 129, 140604),
                           Nucleic.Acids = c(2012, 1279, 31, 4, 5, 3331),
                      ProteinComplex = c(6547, 259, 800, 6, 2, 7614),
                                   Other = c(8, 8, 0, 13, 1, 30),
                                   Total = c(135447, 12608, 3108, 279, 137, 151579)
                  )
```

Q2: Type HIV in the PDB website search box on the home page and
determine how many HIV-1 protease structures are in the current PDB?

There are 1157 as of 2019-05-07 see:
<http://www.rcsb.org/pdb/results/results.do?tabtoshow=Current&qrid=880B2268>

\#\#Section 3 Using Bio3D

``` r
library(bio3d)

pdb<- read.pdb("1hsg.pdb")
pdb
```

    ## 
    ##  Call:  read.pdb(file = "1hsg.pdb")
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1686,  XYZs#: 5058  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 172  (residues: 128)
    ##      Non-protein/nucleic resid values: [ HOH (127), MK1 (1) ]
    ## 
    ##    Protein sequence:
    ##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
    ##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
    ##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
    ##       VNIIGRNLLTQIGCTLNF
    ## 
    ## + attr: atom, xyz, seqres, helix, sheet,
    ##         calpha, remark, call

``` r
atom.select(pdb, "protein", value= TRUE)
```

    ## 
    ##  Call:  trim.pdb(pdb = pdb, sele)
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1514,  XYZs#: 4542  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 0  (residues: 0)
    ##      Non-protein/nucleic resid values: [ none ]
    ## 
    ##    Protein sequence:
    ##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
    ##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
    ##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
    ##       VNIIGRNLLTQIGCTLNF
    ## 
    ## + attr: atom, helix, sheet, seqres, xyz,
    ##         calpha, call

Q6. How many amino acid residues are there in this pdb object and what
are the two nonprotein residues?

amino acid residues: 1686 two non protein residues: HOH (127), MK1 (1)

``` r
attributes(pdb)
```

    ## $names
    ## [1] "atom"   "xyz"    "seqres" "helix"  "sheet"  "calpha" "remark" "call"  
    ## 
    ## $class
    ## [1] "pdb" "sse"

``` r
head(pdb$atom)
```

    ##   type eleno elety  alt resid chain resno insert      x      y     z o
    ## 1 ATOM     1     N <NA>   PRO     A     1   <NA> 29.361 39.686 5.862 1
    ## 2 ATOM     2    CA <NA>   PRO     A     1   <NA> 30.307 38.663 5.319 1
    ## 3 ATOM     3     C <NA>   PRO     A     1   <NA> 29.760 38.071 4.022 1
    ## 4 ATOM     4     O <NA>   PRO     A     1   <NA> 28.600 38.302 3.676 1
    ## 5 ATOM     5    CB <NA>   PRO     A     1   <NA> 30.508 37.541 6.342 1
    ## 6 ATOM     6    CG <NA>   PRO     A     1   <NA> 29.296 37.591 7.162 1
    ##       b segid elesy charge
    ## 1 38.10  <NA>     N   <NA>
    ## 2 40.62  <NA>     C   <NA>
    ## 3 42.64  <NA>     C   <NA>
    ## 4 43.40  <NA>     O   <NA>
    ## 5 37.87  <NA>     C   <NA>
    ## 6 38.40  <NA>     C   <NA>

# Print a subset of $atom data for the first two atoms

``` r
aa321(pdb$atom$resid)
```

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ##    [1] "P" "P" "P" "P" "P" "P" "P" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "I"
    ##   [18] "I" "I" "I" "I" "I" "I" "I" "T" "T" "T" "T" "T" "T" "T" "L" "L" "L"
    ##   [35] "L" "L" "L" "L" "L" "W" "W" "W" "W" "W" "W" "W" "W" "W" "W" "W" "W"
    ##   [52] "W" "W" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "R" "R" "R" "R" "R" "R"
    ##   [69] "R" "R" "R" "R" "R" "P" "P" "P" "P" "P" "P" "P" "L" "L" "L" "L" "L"
    ##   [86] "L" "L" "L" "V" "V" "V" "V" "V" "V" "V" "T" "T" "T" "T" "T" "T" "T"
    ##  [103] "I" "I" "I" "I" "I" "I" "I" "I" "K" "K" "K" "K" "K" "K" "K" "K" "K"
    ##  [120] "I" "I" "I" "I" "I" "I" "I" "I" "G" "G" "G" "G" "G" "G" "G" "G" "Q"
    ##  [137] "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "L" "L" "L" "L" "L" "L" "L" "L" "K"
    ##  [154] "K" "K" "K" "K" "K" "K" "K" "K" "E" "E" "E" "E" "E" "E" "E" "E" "E"
    ##  [171] "A" "A" "A" "A" "A" "L" "L" "L" "L" "L" "L" "L" "L" "L" "L" "L" "L"
    ##  [188] "L" "L" "L" "L" "D" "D" "D" "D" "D" "D" "D" "D" "T" "T" "T" "T" "T"
    ##  [205] "T" "T" "G" "G" "G" "G" "A" "A" "A" "A" "A" "D" "D" "D" "D" "D" "D"
    ##  [222] "D" "D" "D" "D" "D" "D" "D" "D" "D" "D" "T" "T" "T" "T" "T" "T" "T"
    ##  [239] "V" "V" "V" "V" "V" "V" "V" "L" "L" "L" "L" "L" "L" "L" "L" "E" "E"
    ##  [256] "E" "E" "E" "E" "E" "E" "E" "E" "E" "E" "E" "E" "E" "E" "E" "E" "M"
    ##  [273] "M" "M" "M" "M" "M" "M" "M" "S" "S" "S" "S" "S" "S" "L" "L" "L" "L"
    ##  [290] "L" "L" "L" "L" "P" "P" "P" "P" "P" "P" "P" "G" "G" "G" "G" "R" "R"
    ##  [307] "R" "R" "R" "R" "R" "R" "R" "R" "R" "W" "W" "W" "W" "W" "W" "W" "W"
    ##  [324] "W" "W" "W" "W" "W" "W" "K" "K" "K" "K" "K" "K" "K" "K" "K" "P" "P"
    ##  [341] "P" "P" "P" "P" "P" "K" "K" "K" "K" "K" "K" "K" "K" "K" "M" "M" "M"
    ##  [358] "M" "M" "M" "M" "M" "I" "I" "I" "I" "I" "I" "I" "I" "G" "G" "G" "G"
    ##  [375] "G" "G" "G" "G" "I" "I" "I" "I" "I" "I" "I" "I" "G" "G" "G" "G" "G"
    ##  [392] "G" "G" "G" "F" "F" "F" "F" "F" "F" "F" "F" "F" "F" "F" "I" "I" "I"
    ##  [409] "I" "I" "I" "I" "I" "K" "K" "K" "K" "K" "K" "K" "K" "K" "V" "V" "V"
    ##  [426] "V" "V" "V" "V" "R" "R" "R" "R" "R" "R" "R" "R" "R" "R" "R" "Q" "Q"
    ##  [443] "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Y" "Y" "Y" "Y" "Y" "Y" "Y" "Y" "Y" "Y"
    ##  [460] "Y" "Y" "D" "D" "D" "D" "D" "D" "D" "D" "Q" "Q" "Q" "Q" "Q" "Q" "Q"
    ##  [477] "Q" "Q" "I" "I" "I" "I" "I" "I" "I" "I" "L" "L" "L" "L" "L" "L" "L"
    ##  [494] "L" "I" "I" "I" "I" "I" "I" "I" "I" "E" "E" "E" "E" "E" "E" "E" "E"
    ##  [511] "E" "I" "I" "I" "I" "I" "I" "I" "I" "C" "C" "C" "C" "C" "C" "G" "G"
    ##  [528] "G" "G" "H" "H" "H" "H" "H" "H" "H" "H" "H" "H" "K" "K" "K" "K" "K"
    ##  [545] "K" "K" "K" "K" "A" "A" "A" "A" "A" "I" "I" "I" "I" "I" "I" "I" "I"
    ##  [562] "G" "G" "G" "G" "T" "T" "T" "T" "T" "T" "T" "V" "V" "V" "V" "V" "V"
    ##  [579] "V" "L" "L" "L" "L" "L" "L" "L" "L" "V" "V" "V" "V" "V" "V" "V" "G"
    ##  [596] "G" "G" "G" "P" "P" "P" "P" "P" "P" "P" "T" "T" "T" "T" "T" "T" "T"
    ##  [613] "P" "P" "P" "P" "P" "P" "P" "V" "V" "V" "V" "V" "V" "V" "N" "N" "N"
    ##  [630] "N" "N" "N" "N" "N" "I" "I" "I" "I" "I" "I" "I" "I" "I" "I" "I" "I"
    ##  [647] "I" "I" "I" "I" "G" "G" "G" "G" "R" "R" "R" "R" "R" "R" "R" "R" "R"
    ##  [664] "R" "R" "N" "N" "N" "N" "N" "N" "N" "N" "L" "L" "L" "L" "L" "L" "L"
    ##  [681] "L" "L" "L" "L" "L" "L" "L" "L" "L" "T" "T" "T" "T" "T" "T" "T" "Q"
    ##  [698] "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "I" "I" "I" "I" "I" "I" "I" "I" "G"
    ##  [715] "G" "G" "G" "C" "C" "C" "C" "C" "C" "T" "T" "T" "T" "T" "T" "T" "L"
    ##  [732] "L" "L" "L" "L" "L" "L" "L" "N" "N" "N" "N" "N" "N" "N" "N" "F" "F"
    ##  [749] "F" "F" "F" "F" "F" "F" "F" "F" "F" "P" "P" "P" "P" "P" "P" "P" "Q"
    ##  [766] "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "I" "I" "I" "I" "I" "I" "I" "I" "T"
    ##  [783] "T" "T" "T" "T" "T" "T" "L" "L" "L" "L" "L" "L" "L" "L" "W" "W" "W"
    ##  [800] "W" "W" "W" "W" "W" "W" "W" "W" "W" "W" "W" "Q" "Q" "Q" "Q" "Q" "Q"
    ##  [817] "Q" "Q" "Q" "R" "R" "R" "R" "R" "R" "R" "R" "R" "R" "R" "P" "P" "P"
    ##  [834] "P" "P" "P" "P" "L" "L" "L" "L" "L" "L" "L" "L" "V" "V" "V" "V" "V"
    ##  [851] "V" "V" "T" "T" "T" "T" "T" "T" "T" "I" "I" "I" "I" "I" "I" "I" "I"
    ##  [868] "K" "K" "K" "K" "K" "K" "K" "K" "K" "I" "I" "I" "I" "I" "I" "I" "I"
    ##  [885] "G" "G" "G" "G" "G" "G" "G" "G" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q"
    ##  [902] "L" "L" "L" "L" "L" "L" "L" "L" "K" "K" "K" "K" "K" "K" "K" "K" "K"
    ##  [919] "E" "E" "E" "E" "E" "E" "E" "E" "E" "A" "A" "A" "A" "A" "L" "L" "L"
    ##  [936] "L" "L" "L" "L" "L" "L" "L" "L" "L" "L" "L" "L" "L" "D" "D" "D" "D"
    ##  [953] "D" "D" "D" "D" "T" "T" "T" "T" "T" "T" "T" "G" "G" "G" "G" "A" "A"
    ##  [970] "A" "A" "A" "D" "D" "D" "D" "D" "D" "D" "D" "D" "D" "D" "D" "D" "D"
    ##  [987] "D" "D" "T" "T" "T" "T" "T" "T" "T" "V" "V" "V" "V" "V" "V" "V" "L"
    ## [1004] "L" "L" "L" "L" "L" "L" "L" "E" "E" "E" "E" "E" "E" "E" "E" "E" "E"
    ## [1021] "E" "E" "E" "E" "E" "E" "E" "E" "M" "M" "M" "M" "M" "M" "M" "M" "S"
    ## [1038] "S" "S" "S" "S" "S" "L" "L" "L" "L" "L" "L" "L" "L" "P" "P" "P" "P"
    ## [1055] "P" "P" "P" "G" "G" "G" "G" "R" "R" "R" "R" "R" "R" "R" "R" "R" "R"
    ## [1072] "R" "W" "W" "W" "W" "W" "W" "W" "W" "W" "W" "W" "W" "W" "W" "K" "K"
    ## [1089] "K" "K" "K" "K" "K" "K" "K" "P" "P" "P" "P" "P" "P" "P" "K" "K" "K"
    ## [1106] "K" "K" "K" "K" "K" "K" "M" "M" "M" "M" "M" "M" "M" "M" "I" "I" "I"
    ## [1123] "I" "I" "I" "I" "I" "G" "G" "G" "G" "G" "G" "G" "G" "I" "I" "I" "I"
    ## [1140] "I" "I" "I" "I" "G" "G" "G" "G" "G" "G" "G" "G" "F" "F" "F" "F" "F"
    ## [1157] "F" "F" "F" "F" "F" "F" "I" "I" "I" "I" "I" "I" "I" "I" "K" "K" "K"
    ## [1174] "K" "K" "K" "K" "K" "K" "V" "V" "V" "V" "V" "V" "V" "R" "R" "R" "R"
    ## [1191] "R" "R" "R" "R" "R" "R" "R" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Y"
    ## [1208] "Y" "Y" "Y" "Y" "Y" "Y" "Y" "Y" "Y" "Y" "Y" "D" "D" "D" "D" "D" "D"
    ## [1225] "D" "D" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "I" "I" "I" "I" "I" "I"
    ## [1242] "I" "I" "L" "L" "L" "L" "L" "L" "L" "L" "I" "I" "I" "I" "I" "I" "I"
    ## [1259] "I" "E" "E" "E" "E" "E" "E" "E" "E" "E" "I" "I" "I" "I" "I" "I" "I"
    ## [1276] "I" "C" "C" "C" "C" "C" "C" "G" "G" "G" "G" "H" "H" "H" "H" "H" "H"
    ## [1293] "H" "H" "H" "H" "K" "K" "K" "K" "K" "K" "K" "K" "K" "A" "A" "A" "A"
    ## [1310] "A" "I" "I" "I" "I" "I" "I" "I" "I" "G" "G" "G" "G" "T" "T" "T" "T"
    ## [1327] "T" "T" "T" "V" "V" "V" "V" "V" "V" "V" "L" "L" "L" "L" "L" "L" "L"
    ## [1344] "L" "V" "V" "V" "V" "V" "V" "V" "G" "G" "G" "G" "P" "P" "P" "P" "P"
    ## [1361] "P" "P" "T" "T" "T" "T" "T" "T" "T" "P" "P" "P" "P" "P" "P" "P" "V"
    ## [1378] "V" "V" "V" "V" "V" "V" "N" "N" "N" "N" "N" "N" "N" "N" "I" "I" "I"
    ## [1395] "I" "I" "I" "I" "I" "I" "I" "I" "I" "I" "I" "I" "I" "G" "G" "G" "G"
    ## [1412] "R" "R" "R" "R" "R" "R" "R" "R" "R" "R" "R" "N" "N" "N" "N" "N" "N"
    ## [1429] "N" "N" "L" "L" "L" "L" "L" "L" "L" "L" "L" "L" "L" "L" "L" "L" "L"
    ## [1446] "L" "T" "T" "T" "T" "T" "T" "T" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q"
    ## [1463] "I" "I" "I" "I" "I" "I" "I" "I" "G" "G" "G" "G" "C" "C" "C" "C" "C"
    ## [1480] "C" "T" "T" "T" "T" "T" "T" "T" "L" "L" "L" "L" "L" "L" "L" "L" "N"
    ## [1497] "N" "N" "N" "N" "N" "N" "N" "F" "F" "F" "F" "F" "F" "F" "F" "F" "F"
    ## [1514] "F" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X"
    ## [1531] "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X"
    ## [1548] "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X"
    ## [1565] "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X"
    ## [1582] "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X"
    ## [1599] "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X"
    ## [1616] "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X"
    ## [1633] "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X"
    ## [1650] "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X"
    ## [1667] "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X"
    ## [1684] "X" "X" "X"

``` r
atom.select(pdb, "protein")
```

    ## 
    ##  Call:  atom.select.pdb(pdb = pdb, string = "protein")
    ## 
    ##    Atom Indices#: 1514  ($atom)
    ##    XYZ  Indices#: 4542  ($xyz)
    ## 
    ## + attr: atom, xyz, call

``` r
pdb$atom[1:2, c("eleno", "elety", "x","y","z")]
```

    ##   eleno elety      x      y     z
    ## 1     1     N 29.361 39.686 5.862
    ## 2     2    CA 30.307 38.663 5.319

# Note that individual $atom records can also be accessed like this

``` r
pdb$atom$elety[1:2]
```

    ## [1] "N"  "CA"

# Which allows us to do the following

``` r
#plot.bio3d(pdb$atom$b[pdb$calpha], sse=pdb, typ="l", ylab=“B-factor”)
```

Q7. What type of R object is pdb$atom? HINT: You can always use the
str() function to get a useful summery of any R object.

# Print a summary of the coordinate data in $xyz

``` r
pdb$xyz
```

    ## 
    ##    Total Frames#: 1
    ##    Total XYZs#:   5058,  (Atoms#:  1686)
    ## 
    ##     [1]  29.361  39.686  5.862  <...>  30.112  17.912  -4.791  [5058] 
    ## 
    ## + attr: Matrix DIM = 1 x 5058

Atom selection is done via the function **atom.select()**

``` r
prot.pdb<- atom.select(pdb, "ligand", value= TRUE)
write.pdb(prot.pdb, file="1hsg_protein.pdb")
```

``` r
#pdb$atom[inds$atom, ]
```

``` r
lig.pdb<- atom.select(pdb, "ligand", value= TRUE)
write.pdb(lig.pdb, file="1hsg_protein.pdb")
```

``` r
aa <- get.seq("1ake_A")
```

    ## Warning in get.seq("1ake_A"): Removing existing file: seqs.fasta

``` r
aa
```

    ##              1        .         .         .         .         .         60 
    ## pdb|1AKE|A   MRIILLGAPGAGKGTQAQFIMEKYGIPQISTGDMLRAAVKSGSELGKQAKDIMDAGKLVT
    ##              1        .         .         .         .         .         60 
    ## 
    ##             61        .         .         .         .         .         120 
    ## pdb|1AKE|A   DELVIALVKERIAQEDCRNGFLLDGFPRTIPQADAMKEAGINVDYVLEFDVPDELIVDRI
    ##             61        .         .         .         .         .         120 
    ## 
    ##            121        .         .         .         .         .         180 
    ## pdb|1AKE|A   VGRRVHAPSGRVYHVKFNPPKVEGKDDVTGEELTTRKDDQEETVRKRLVEYHQMTAPLIG
    ##            121        .         .         .         .         .         180 
    ## 
    ##            181        .         .         .   214 
    ## pdb|1AKE|A   YYSKEAEAGNTKYAKVDGTKPVAEVRADLEKILG
    ##            181        .         .         .   214 
    ## 
    ## Call:
    ##   read.fasta(file = outfile)
    ## 
    ## Class:
    ##   fasta
    ## 
    ## Alignment dimensions:
    ##   1 sequence rows; 214 position columns (214 non-gap, 0 gap) 
    ## 
    ## + attr: id, ali, call

# Blast or hmmer search

``` r
b <- blast.pdb(aa)
```

    ##  Searching ... please wait (updates every 5 seconds) RID = FDDHSFFM015 
    ##  ..
    ##  Reporting 97 hits
