# Cattle Centromere Analysis
---

These are simplified Unix commandline commands used to generate summaries of centromeric repeats in the [ARS-UCD1.2](https://www.ncbi.nlm.nih.gov/assembly/GCF_002263795.1/) and the [UMD3.1](https://ccb.jhu.edu/bos_taurus_assembly.shtml) cattle reference genome assemblies.

Actual filenames were changed in the creation of this note so as to be more easy to recognize. 

## Software Prerequisites
* [RepeatMasker](http://www.repeatmasker.org/)
* [Bedtools](https://bedtools.readthedocs.io/en/latest/)
* [R](https://www.r-project.org/)

## RepeatMasking

First, we need to generate a list of repetitive sites in the references.

```bash
RepeatMasker -pa 20 -q -species cow -no_is -gff ARS-UCD1.2.fa
RepeatMasker -pa 20 -q -species cow -no_is -gff UMD3.1.fa

# Next, we need to convert the output into bed format for easier parsing
perl -e '<>; <>; <>; while(<>){ $_ =~ s/^\s+//; @s = split(/\s+/); $orient = ($s[8] eq "+")? "+" : "-"; $qlen = $s[12] - $s[11]; print "$s[4]\t$s[5]\t$s[6]\t$orient\t$s[9]\t$s[10]\t$qlen\n";}' < ARS-UCD1.2.fa.out > ARS-UCD1.2.fa.out.bed
perl -e '<>; <>; <>; while(<>){ $_ =~ s/^\s+//; @s = split(/\s+/); $orient = ($s[8] eq "+")? "+" : "-"; $qlen = $s[12] - $s[11]; print "$s[4]\t$s[5]\t$s[6]\t$orient\t$s[9]\t$s[10]\t$qlen\n";}' < UMD3.1.fa.out > UMD3.1.fa.out.bed
```

## Centromere analysis

Now that we have a list of repetitive sites in the reference, we need to extract the centromeric satellite sequence for further analysis. 

```bash
bedtools sort -i ARS-UCD1.2.fa.out.bed | grep 'Satellite/centr' | bedtools merge -i stdin -c 5 -o distinct -d 5000 > ARS-UCD1.2.fa.out.centromeres.bed
bedtools sort -i UMD3.1.fa.out.bed | grep 'Satellite/centr' | bedtools merge -i stdin -c 5 -o distinct -d 5000 > UMD3.1.fa.out.centromeres.bed

# Next, we're going to convert them into something that is useable in an R dataframe
perl -lane '$a = "ARSUCD"; $l = $F[2] - $F[1]; $p = ($F[0] =~ m/^NKLS/)? "UNPLACED" : "CHRSCAFF";if($p ne "UNPLACED" && $F[1] <= 500000){ $p = "CHRSTART";} print join("\t", @F) . "\t$l\t$a\t$p";' < ARS-UCD1.2.fa.out.centromeres.bed > ARS-UCD1.2.fa.out.centromeres.bed.rtab
perl -lane '$a = "UMD3"; $l = $F[2] - $F[1]; $p = ($F[0] =~ m/^718/)? "UNPLACED" : "CHRSCAFF"; if($p ne "UNPLACED" && $F[1] <= 500000){ $p = "CHRSTART";} print join("\t", @F) . "\t$l\t$a\t$p";' < UMD3.1.fa.out.centromeres.bed > UMD3.1.fa.out.centromeres.bed.rtab

cat ARS-UCD1.2.fa.out.centromeres.bed.rtab UMD3.1.fa.out.centromeres.bed.rtab > combined_centromere.rtab
```

From here, it's just a matter of reading the data into R and using the built-in functions to calculate the statistics needed.

```R
library(ggplot2)
library(dplyr)
data <- read.delim("combined_centromere.rtab", header=FALSE)

colnames(data) <- c("chr", "start", "end", "cats", "Length", "Assembly", "Scaff")
data.sum <- group_by(data, Assembly, Scaff) %>% summarize(Num = n(), Min= min(Length), Max= max(Length), Mean = mean(Length), SD = sd(Length), SE = (sd(Length)/ sqrt(length(Length))))

# Reordering levels to make it neater
data.sum$Scaff = factor(data.sum$Scaff, levels(data.sum$Scaff)[c(2,1,3)])

pdf("centromere_size_histogram.pdf", useDingbats=FALSE)
ggplot(data.sum, aes(x=Assembly, y=Mean, fill=Assembly)) + scale_fill_brewer(palette="Dark2") + geom_bar(stat="identity", colour="black", size=0.2) + geom_errorbar(aes(ymin = Mean - (2 * SE), ymax = Mean + (2 * SE)), width=0.2) + xlab("Assembly") + ylab("Centromeric repeat length (bp)") + theme_bw() + facet_wrap(~Scaff)
dev.off()

# The following will show a table of the summary statistics for each assembly, grouped by centromere category
data.sum
```