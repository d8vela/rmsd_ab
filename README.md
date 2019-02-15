
# Description:

This script will calculate the mutliple RMSD values for a given antibody model vs reference (X-ray) structure.

Install JSON module from CPAN (for decoding JSON output from antibody.cc)
 > cpan install JSON

Install ProFit (for superpositioning and RMSD calculation)
Download and install: http://www.bioinf.org.uk/software/swreg.html

Define path to ProFit executable (edit script appropriately if necessary):
Modify line 13 of the rmsd_ab.ama2.pl define the path of the local ProFit program.

# Notes

All PDB files used as input must be consistent in numbering using 1 to the last residue for chains L and H.

# Example Execution:

## Example Input:

Model PDB file (Model Ouput Stucture):
> model/model-0.relaxed.renum.pdb

Crystal Structure (Correct Structure)
> xtal/4kmt.pdb

JSON file that defines all CDR Ranges
> json/4kmt.json


## Example Command:
> ./rmsd_ab.ama2.pl model/model-0.relaxed.renum.pdb xtal/4kmt.pdb json/4kmt.json


## STDOUT (using input and command above):
# Header description: Superposition Range, RMSD Calculation Range, "RMS:", RMSD Value
```
All All RMS: 0.821
FR FR RMS: 0.587
FR Loops RMS: 1.494
H H RMS: 0.873
FRH FRH RMS: 0.329
FRH LoopH RMS: 1.862
FRH H1 RMS: 1.133
FRH H2 RMS: 0.850
FRH H3 RMS: 2.472
L L RMS: 0.589
FRL FRL RMS: 0.554
FRL LoopL RMS: 0.813
FRL L1 RMS: 0.773
FRL L2 RMS: 0.918
FRL L3 RMS: 0.804
H L RMS: 0.950
L H RMS: 1.213
Can Can RMS: 0.777
FR ALL RMS: 0.827
```


# Definitions

## Terms taken from AMA-1 Article (Almagro, Proteins 2011), unless stated otherwise:

```
FR: Framework Regions (The b-sheets and nonhypervariable loops serve as structural support to the antigen-binding site)
HVL: Hypervariable loops (L1,L2,L3,H1,H2,H3)


All: whole FV; 
LOOPS: all hyperveriable loops;
Can: Loops with canonical structures (NOTE: Not exactly like AMA defintion.  This considers all CDR loops, except CDR H3)
L: V_L (Light variable domain)
H: V_H (Heavy variable domain)
```

## Terms taken from the AMA-2 Article (Almagro, Proteins 2014), unless stated otherwise:

```
All:  The Entire Fv domains (V_L + L_H)
FR:   Framework only for Fv [V_L + V_H] (do not include HVL)


L:   Entire V_L domain (entire light variable domain)
FRL:  Framework for light variable domain (do not include HVL)
LoopL:  All HVL Loops of the light varible domain
L1:   HVL 1 of the light variable domain
L2:  HVL 2 of the light variable domain
L3:  HVL 3 of the light variable domain


H:   V_H domain (entire heavy variable domain)
FRH:  Framework for heavy variable domain (do not include HVL)
LoopH: All HVL Loops of the heavy variable domain
H1:  HVL 1 of the heavy variable domain
H2:  HVL 2 of the heavy variable domain
H3:  HVL 3 of the heavy varaible domain

Can:  HVL with canonical structures (Note: Not exactly like AMA Defintion.  This considers all CDR loops, except CDR H3)
FW:   Should this be FR (Framework)?  Possible typo in the AMA2 article for Table II (Almagro, Proteins 2014)
```
