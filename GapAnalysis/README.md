# Gap Analysis
---

Scripts and command line code to replicate the Gap and Repeat intersection analysis in the manuscript.

## Requirements
* Bedtools
* Samtools
* RepeatMasker
* BWA
* Perl version 5+
* Python version 3.6+

## Gap Fill Analysis

```bash
sbatch --nodes=1 --mem=12000 --ntasks-per-node=1 -p msn --wrap="bwa index umd3_reference_genome.fasta"

# Identifying gaps
sbatch --nodes=1 --mem=25000 --ntasks-per-node=3 -p medium --wrap="perl identifyFilledGaps.pl -o umd3_reference_genome.fasta -s ../ARS-UCD1.2_Btau5.0.1Y/ARS-UCD1.2_Btau5.0.1Y.fa -g ~/rumen_longread_metagenome_assembly/binaries/GetMaskBedFasta/store/GetMaskBedFasta.jar -j /software/7/apps/java/1.8.0_121/bin/java -d umd3_gaps_onarsucd.newlogic.tab"

python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f umd3_gaps_on_arsucd.tab -c 0 -m -d '\t'
|Entry    | Value|
|:--------|-----:|
|Closed   | 58110|
|Trans    |  8508|
|Unmapped |  5570|

# Now to print out a table for plotting in R
perl -lane 'if($F[0] =~ /Trans/){next;} print "$F[0]\t$F[4]\t$F[10]";' < umd3_gaps_onarsucd.newlogic.tab > umd3_gaps_onarsucd.newlogic.lens
```

## Repeat association analysis

```bash
sbatch --nodes=1 --mem=45000 --ntasks-per-node=50 -p msn --wrap="RepeatMasker -pa 50 -q -species cow -no_is -gff ../ARS-UCD1.2_Btau5.0.1Y/ARS-UCD1.2_Btau5.0.1Y.fa"

# Here is my one-liner for greping out the repeats in bed format
perl -e '<>; <>; <>; while(<>){ $_ =~ s/^\s+//; @s = split(/\s+/); $orient = ($s[8] eq "+")? "+" : "-"; $qlen = $s[12] - $s[11]; print "$s[4]\t$s[5]\t$s[6]\t$orient\t$s[9]\t$s[10]\t$qlen\n";}' < ARS-UCD1.2_Btau5.0.1Y.fa.out > ARS-UCD1.2_Btau5.0.1Y.fa.out.bed

python3 intersectGapFlanksWithRepeats.py -g umd3_gaps_onarsucd.newlogic.tab -r ARS-UCD1.2_Btau5.0.1Y.fa.out.bed -o umd3_gaps_repeat_intersections
```