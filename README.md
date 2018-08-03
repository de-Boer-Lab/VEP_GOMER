# VEP_GOMER
VEP Plug-in for calculating differential TF binding to SNPs


## Installation
```
cd /path/to/where/you/want/to/store/files
git clone https://github.com/Carldeboer/VEP_GOMER/
cd VEP_GOMER
./SETUP.sh
#this should generate  Homo_sapiens.VEP.CISBP.txt for use with VEP and Homo_sapiens.TF_Information.txt which contains the Motif_ID <-> TF mapping.
#link the GOMER VEP plugin to VEP. Usually VEP modules are found here: ~/.vep/Plugins/
ln -s ~/.vep/Plugins/ `pwd`/GOMER.pm
```

Now GOMER VEP should be usable.  Here is an example:
```
./vep -i EXAMPLE.vcf  --plugin GOMER,/path/to/where/you/want/to/store/files/VEP_GOMER/Homo_sapiens.VEP.CISBP.txt -o EXAMPLE_GOMER_VEP.out --cache --force_overwrite --pick
```
Notes: 
* Since this is not gene-specific, I include `--pick` so that there is only one output per variant
* Binding scores are output as follows: Motif_ID/scoreAlternate/scoreReference/(scoreAlternate - scoreReference),
* Scores are log occupancy, as in the GOMER, but instead of log(P(bound)) we report log(E(bound)) [expected binding, summing the predicted binding events across the locus]
* `Homo_sapiens.VEP.CISBP.txt` is of the form: 1 motif per line, `Motif_ID  motifFile scoreScale  scoreCenter scoreCutoff deltaCutoff`
* Default scoreScale and scoreCenter are 1 and 0, respectively (i.e. occupancy scores are not scaled or centered)
* scoreCutoff determines the minimal binding score for a locus to get before we consider it "potentially occupied" by the TF. The included defaults in `Homo_sapiens.VEP.CISBP.txt` are set such that only the top 5% of binding sites will be included (using a +/-50 bp window). 
* deltaCutoff is set to 0.1 in `Homo_sapiens.VEP.CISBP.txt` so that there has to be at least a 0.1 difference in binding score between alleles - this, combined with the SNP being within the top 5% of binding sites has the effect of most motifs having "hits" at ~1% of SNPs.
