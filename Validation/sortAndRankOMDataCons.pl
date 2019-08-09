#!/usr/bin/perl
#SBATCH --nodes=1
#SBATCH --mem=15000
#SBATCH --ntasks-per-node=1
#SBATCH -p short
# This script is designed to take a BNG optical map and generate rank order estimates for it based on coordinate
# order and map order
# BtOMcontig_36-chr10_OM-1_035    NC_037337.1     F       0       49      566     619     0       442050  4958992 5406759
# BtOMcontig_36-chr10_OM-1_036    NC_037337.1     F       0       49      584     638     0       396960  5097917 5495426
# BtOMcontig_36-chr10_OM-1_037    NC_037337.1     F       0       49      599     653     0       291100  5279543 5566475
# BtOMcontig_36-chr10_OM-1_038    NC_037337.1     F       0       49      615     668     0       276080  5389288 5662516
# BtOMcontig_36-chr10_OM-1_039    NC_037337.1     F       0       49      634     683     0       370690  5477536 5849865

use strict;
my $usage = "perl $0 <input om aln file> <output rank file> <output concensus file>\n";
chomp @ARGV;

unless(scalar(@ARGV) ==3){
	print $usage;
	exit;
}

my %data; # {refchr} = [omcontig, omcontigintnum, refstart, refend]
my %omctglow; # {omcontigname} = lowest number observed
my %omctgmax; # {omcontigname} = total number observed
my %OmRefCtgCts; # {omcontig}->{refchr}
open(my $IN, "< $ARGV[0]") || die "Could not find aln file!\n$usage";
my $head = <$IN>;
while(my $line = <$IN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	my @omsegs = split(/[_-]/, $segs[0]);
	my $omnum = $omsegs[-1];
	my $omctgname = join("_", @omsegs[0..1]);
	push(@{$data{$segs[1]}}, [$omctgname, $omnum, $segs[5], $segs[6]]);

	if(!exists($omctglow{$omctgname})){
		$omctglow{$omctgname} = $omnum;
	}
	# Reset OM window number down to the lowest one observed
	$omctglow{$omctgname} = ($omnum < $omctglow{$omctgname})? $omnum : $omctglow{$omctgname};
	
	# Total number of OM windows observed
	$omctgmax{$omctgname} += 1;

	# For later stats file
	$OmRefCtgCts{$omctgname}->{$segs[1]} += 1;
}
close $IN;

# define the list of ordered reference chromosomes by lexicographic sorting
my @orderedChrs = sort{$a cmp $b} keys(%data);

# define the list of ordered optical maps by their position in the list
# Addition: make sure that the OM consensus places the OM squarely in the chromosome!
my $cons = determineCons(\%OmRefCtgCts);
my @orderedOM;
my %seen;
foreach my $chr (@orderedChrs){
	foreach my $row (@{$data{$chr}}){
		if(!exists($seen{$row->[0]}) && $cons->{$row->[0]} eq $chr){
			push(@orderedOM, $row->[0]);
			$seen{$row->[0]} = 1;
		}
	}
}

my %chrRank = map{$orderedChrs[$_] => $_} 0..$#orderedChrs;

# Adjust OMRank order because of total number of entries
my %omAdjust; my $prevOM = 1;
foreach my $om (@orderedOM){
	$omAdjust{$om} = $prevOM;
	$prevOM += $omctgmax{$om};
}

open(my $OUT, "> $ARGV[1]");
print {$OUT} "RefName\tOMName\tRefRank\tOMRank\n";
my $chrRank = 1;
# now run through and generate rank order
foreach my $chr (@orderedChrs){
	foreach my $row (sort{$a->[2] <=> $b->[2]} @{$data{$chr}}){
		my $adjOMRank = $omAdjust{$row->[0]} + ($row->[1] - $omctglow{$row->[0]});
		print {$OUT} "$chr\t$row->[0]\t$chrRank\t$adjOMRank\n";
		$chrRank += 1;
	}
}
close $OUT;

# Now print the number of OM conflicts
open(my $OUT, "> $ARGV[2]");
foreach my $om (sort{scalar(keys(%{$OmRefCtgCts{$b}})) <=> scalar(keys(%{$OmRefCtgCts{$a}}))} keys(%OmRefCtgCts)){
	my $count = scalar(keys(%{$OmRefCtgCts{$om}}));
	print {$OUT} "$om\t$count";
	foreach my $k (sort{$OmRefCtgCts{$om}->{$b} <=> $OmRefCtgCts{$om}->{$a}} keys(%{$OmRefCtgCts{$om}})){
		print {$OUT} "\t$k\:" . $OmRefCtgCts{$om}->{$k};
	}
	print {$OUT} "\n";
}
close $OUT;

exit;

sub determineCons{
	my ($href) = @_;
	# {OM}->{Ref} = num
	my %cons;
	foreach my $om (keys(%{$href})){
		my $winner = "NA";
		my $val = -1;
		foreach my $ref (keys(%{$href->{$om}})){
			if($href->{$om}->{$ref} > $val){
				$val = $href->{$om}->{$ref};
				$winner = $ref;
			}
		}
		$cons{$om} = $winner;
	}
	return \%cons;
}
