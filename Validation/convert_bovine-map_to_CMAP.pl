#!/usr/bin/perl
# This script was written by the Computomix team
# The intended use of this script is to convert optical map formats into a file format that works with BioNano software


my $usage="$0 <mapfile> [<factor, def: 1000>]\n";
my $file = shift or die $usage;
my $factor;
if ($factor eq "") { $factor = 1000; }
open OUT, ">$file.cmap";
open KEY, ">$file.key";
open F,"<$file";
my %markers=(); my %scaffsize=(); my $chr=""; my $chr_nr=1; my $chr_transl=0;
while (<F>) {
        chomp;
        if ($_ =~ m/^chr/) {
                $chr_transl++;
                $chr = substr($_, 3, index($_, "_", 3));
                $chr_nr = 1;
                $chr .= "_".$chr_nr if ($_ =~ m/gap/);
                print KEY "$chr\t$chr_transl\n";
        }
        elsif ($_ !~ m/^$/) {
                my @a=split" ";
                my @tmp=(); $markers{$chr_transl} = \@tmp;
                for (my $i=2; $i<@a; ++$i) {
                        if ($a[$i] == "500.00") { # gap -> start new map:
                                my $tmp = $chr_nr+1;
                                $chr =~ s/${chr_nr}$/$tmp/;
                                $chr_nr++;
                                $chr_transl++;
                                my @tmp=(); $markers{$chr_transl} = \@tmp;
                                print KEY "$chr\t$chr_transl\n";
                                next;
                        }
                        $scaffsize{$chr_transl} += $a[$i] * $factor;
                        push @{$markers{$chr_transl}}, $a[$i] * $factor;
                }
        }
}
close F;
close KEY;
print OUT "# CMAP File Version: 0.1\n";
print OUT "# Label Channels:    1\n";
print OUT "# Nickase Recognition Site 1:        GGATCC\n";
print OUT "#h CMapId\tContigLength\tNumSites\tSiteID\tLabelChannel\tPosition\tStdDev\tCoverage\tOccurrence\n";
print OUT "#f int\tfloat\tint\tint\tint\tfloat\tfloat\tint\tint\n";
foreach my $chr (sort{$a<=>$b} keys %scaffsize) {
        my $c=1;
        my $cum_sum = 0;
        foreach my $value (@{$markers{$chr}}) {
                $cum_sum += $value;
                my $channel = ( $c == scalar(@{$markers{$chr}}) ? "0" : "1" );
                print OUT $chr . "\t" . $scaffsize{$chr} . "\t" . (scalar(@{$markers{$chr}})-1) . "\t" .
                          $c++ . "\t" . $channel . "\t" . $cum_sum . "\t".sprintf("%.1f", ($value*.1))."\t1\t1\n";#\t0\t0\t0\n";
        }
}
close OUT;