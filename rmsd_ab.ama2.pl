#!/usr/bin/perl

# Written by David La
# Updated: Mon Jul 16 09:55:19 PDT 2018

# Description:
# This script will calculate the mutliple RMSD values for a given antibody model vs reference (X-ray) structure.

# Install JSON module from CPAN (for decoding JSON output from antibody.cc)
# Install ProFit (for superpositioning and RMSD calculation)

# Define path to ProFit
my $profit_exe = "$ENV{'HOME'}/bin/profit";

# File requirements
# Reference PDB (X-ray) have chain H and L
# Mobile PDB (Model) have chain H and L

my $usage = "Usage: rmsd_ab.pl <model_pdb_file> <reference_pdb_file> <report.json> [chain_range_pdb_file]\n";

use strict;
use JSON qw( decode_json );
use Data::Dumper; 

my $model_pdb_file = $ARGV[0] or die $usage;
my $ref_pdb_file = $ARGV[1] or die $usage;
my $report_file = $ARGV[2] or die $usage;

# PDB file used to determine chain ranges to assess for the model PDB and reference PDB 
# This is useful for full length structure outputs from rebuilding loops in xtal structures, 
# but only assessing partial domains from predicted model outputs.
my $chain_pdb_file = $ARGV[3];

# ------------- Check Data -------------

# Get PDB Chain Ranges
my %model_chain_range;
if ($chain_pdb_file) {
	%model_chain_range = getChainRange($chain_pdb_file);
}
else {
	%model_chain_range = getChainRange($model_pdb_file);
}

my %ref_chain_range = getChainRange($ref_pdb_file);

my $chain_H_begin = ($model_chain_range{'H'}{'begin'} > $ref_chain_range{'H'}{'begin'}) ? $model_chain_range{'H'}{'begin'} : $ref_chain_range{'H'}{'begin'};
my $chain_H_end = ($model_chain_range{'H'}{'end'} < $ref_chain_range{'H'}{'end'}) ? $model_chain_range{'H'}{'end'} : $ref_chain_range{'H'}{'end'};
my $chain_L_begin = ($model_chain_range{'L'}{'begin'} > $ref_chain_range{'L'}{'begin'}) ? $model_chain_range{'L'}{'begin'} : $ref_chain_range{'L'}{'begin'};
my $chain_L_end = ($model_chain_range{'L'}{'end'} < $ref_chain_range{'L'}{'end'}) ? $model_chain_range{'L'}{'end'} : $ref_chain_range{'L'}{'end'};

my %chain_range;
$chain_range{'H'}{'begin'} = $chain_H_begin;
$chain_range{'H'}{'end'} = $chain_H_end;
$chain_range{'L'}{'begin'} = $chain_L_begin;
$chain_range{'L'}{'end'} = $chain_L_end;

# Read entire JSON file into string
local $/ =  undef; # Unset $/, the Input Record Separator, to make <> give you the whole file at once.
open my $json_file, "<", $report_file or die "Cannot open $report_file\n";
my $json = <$json_file>;

# Decode JSON report file
my $report_ref = decode_json($json);
my %report =%{$report_ref};
#print Dumper(\%report);

# Get Loop Ranges
my $L1_begin = $report{l1}{'begin'};
my $L1_end = $report{l1}{'end'};

my $L2_begin = $report{l2}{'begin'};
my $L2_end = $report{l2}{'end'};

my $L3_begin = $report{l3}{'begin'};
my $L3_end = $report{l3}{'end'};

my $H1_begin = $report{h1}{'begin'};
my $H1_end = $report{h1}{'end'};

my $H2_begin = $report{h2}{'begin'};
my $H2_end = $report{h2}{'end'};

my $H3_begin = $report{h3}{'begin'};
my $H3_end = $report{h3}{'end'};

# Output Data
#print "Chain H: $chain_H_begin $chain_H_end\n";
#print "Chain L: $chain_L_begin $chain_L_end\n";

#print "L1 $L1_begin $L1_end\n";
#print "L2 $L2_begin $L2_end\n";
#print "L3 $L3_begin $L3_end\n";

#print "H1 $H1_begin $H1_end\n";
#print "H2 $H2_begin $H2_end\n";
#print "H3 $H3_begin $H3_end\n";


# ------------- Main -------------

# Generate ProFit Code
my $profit_code;
my $rms;

# ---- All All ----
$profit_code = genProFitCode(\%chain_range,\%chain_range);
#print "\n\n";
#print $profit_code;
$rms = runProFit($ref_pdb_file,$model_pdb_file,$profit_code,$profit_exe);
print "All All RMS: $rms\n";

# ---- FR FR ----
# Select all loop ranges
my %all_loops = selectRanges(\%report,'h1','h2','h3','l1','l2','l3');

#print Dumper(\%all_loops);
#print "\n\n";

# Select framework ranges (inverse loop selection)
my %FR_range = inverseRanges(\%chain_range,\%all_loops);
#print "\n\n";
#print Dumper(\%FR_range);

# Generate ProFit code
$profit_code = genProFitCode(\%FR_range,\%FR_range);
#print $profit_code;
$rms = runProFit($ref_pdb_file,$model_pdb_file,$profit_code);
print "FR FR RMS: $rms\n";

# ---- FR Loops ----
# Generate ProFit code
$profit_code = genProFitCode(\%FR_range,\%all_loops);
#print $profit_code;
$rms = runProFit($ref_pdb_file,$model_pdb_file,$profit_code);
print "FR Loops RMS: $rms\n";

# ---- H H ----
# Generate ProFit code
my %H_range;
$H_range{'H'}{'begin'} = $chain_range{'H'}{'begin'};
$H_range{'H'}{'end'} = $chain_range{'H'}{'end'};

$profit_code = genProFitCode(\%H_range,\%H_range);
#print $profit_code;
$rms = runProFit($ref_pdb_file,$model_pdb_file,$profit_code);
print "H H RMS: $rms\n";

# ---- FRH FRH ----
# Select all loop ranges
my %loopH_range = selectRanges(\%report,'h1','h2','h3');

# Heavy Framework ranges
my %FRH_range = inverseRanges(\%H_range,\%loopH_range);

# Generate ProFit code
$profit_code = genProFitCode(\%FRH_range,\%FRH_range);
#print $profit_code;
$rms = runProFit($ref_pdb_file,$model_pdb_file,$profit_code);
print "FRH FRH RMS: $rms\n";

# ---- FRH LoopH ----

# Generate ProFit code
$profit_code = genProFitCode(\%FRH_range,\%loopH_range);
#print $profit_code;
$rms = runProFit($ref_pdb_file,$model_pdb_file,$profit_code);
print "FRH LoopH RMS: $rms\n";

# ---- FRH H1 ----
my %H1_range = selectRanges(\%report,'h1');
$profit_code = genProFitCode(\%FRH_range,\%H1_range);
#print $profit_code;
$rms = runProFit($ref_pdb_file,$model_pdb_file,$profit_code);
print "FRH H1 RMS: $rms\n";

# ---- FRH H2 ----
my %H2_range = selectRanges(\%report,'h2');
$profit_code = genProFitCode(\%FRH_range,\%H2_range);
#print $profit_code;
$rms = runProFit($ref_pdb_file,$model_pdb_file,$profit_code);
print "FRH H2 RMS: $rms\n";

# ---- FRH H3 ----
my %H3_range = selectRanges(\%report,'h3');
$profit_code = genProFitCode(\%FRH_range,\%H3_range);
#print $profit_code;
$rms = runProFit($ref_pdb_file,$model_pdb_file,$profit_code);
print "FRH H3 RMS: $rms\n";

# ---- L L ----
my %L_range;
$L_range{'L'}{'begin'} = $chain_range{'L'}{'begin'};
$L_range{'L'}{'end'} = $chain_range{'L'}{'end'};

$profit_code = genProFitCode(\%L_range,\%L_range);
#print $profit_code;
$rms = runProFit($ref_pdb_file,$model_pdb_file,$profit_code);
print "L L RMS: $rms\n";

# ---- FRL FRL ----
# Select all loop ranges
my %loopL_range = selectRanges(\%report,'l1','l2','l3');

# Heavy Framework ranges
my %FRL_range = inverseRanges(\%L_range,\%loopL_range);

# Generate ProFit code
$profit_code = genProFitCode(\%FRL_range,\%FRL_range);
#print $profit_code;
$rms = runProFit($ref_pdb_file,$model_pdb_file,$profit_code);
print "FRL FRL RMS: $rms\n";


# ---- FRL LoopL ----
$profit_code = genProFitCode(\%FRL_range,\%loopL_range);
#print $profit_code;
$rms = runProFit($ref_pdb_file,$model_pdb_file,$profit_code);
print "FRL LoopL RMS: $rms\n";


# ---- FRL L1 ----
my %L1_range = selectRanges(\%report,'l1');
$profit_code = genProFitCode(\%FRL_range,\%L1_range);
#print $profit_code;
$rms = runProFit($ref_pdb_file,$model_pdb_file,$profit_code);
print "FRL L1 RMS: $rms\n";

# ---- FRL L2 ----
my %L2_range = selectRanges(\%report,'l2');
$profit_code = genProFitCode(\%FRL_range,\%L2_range);
#print $profit_code;
$rms = runProFit($ref_pdb_file,$model_pdb_file,$profit_code);
print "FRL L2 RMS: $rms\n";


# ---- FRL L3 ----
my %L3_range = selectRanges(\%report,'l3');
$profit_code = genProFitCode(\%FRL_range,\%L3_range);
#print $profit_code;
$rms = runProFit($ref_pdb_file,$model_pdb_file,$profit_code);
print "FRL L3 RMS: $rms\n";


# ---- H L ----
$profit_code = genProFitCode(\%H_range,\%L_range);
#print $profit_code;
$rms = runProFit($ref_pdb_file,$model_pdb_file,$profit_code);
print "H L RMS: $rms\n";

# ---- L H ----
$profit_code = genProFitCode(\%L_range,\%H_range);
#print $profit_code;
$rms = runProFit($ref_pdb_file,$model_pdb_file,$profit_code);
print "L H RMS: $rms\n";


# ---- Can Can ----
my %can_loop_range = selectRanges(\%report,'l1','l2','l3','h1','h2');
$profit_code = genProFitCode(\%can_loop_range,\%can_loop_range);
$rms = runProFit($ref_pdb_file,$model_pdb_file,$profit_code);
print "Can Can RMS: $rms\n";

# ---- FR ALL ----
$profit_code = genProFitCode(\%FR_range,\%chain_range);
#print $profit_code;
$rms = runProFit($ref_pdb_file,$model_pdb_file,$profit_code);
print "FR ALL RMS: $rms\n";



# ---- Subs ----

sub runProFit {

	my $ref_pdb = shift;
	my $mobile_pdb = shift;
	my $profit_cmd = shift;

	#print $profit_cmd;

	my $profit_exe = shift || '/home/dla/bin/profit';

	my $random_int = int(rand 10000000000);
	my $tmp_file = ".tmp.profit.$$.$random_int.txt";

	open(FILE,">$tmp_file") or die "Cannot open file $tmp_file\n";
	print FILE "$profit_cmd\n";
	close(FILE);

	my $out = `$profit_exe -h -f $tmp_file $ref_pdb $mobile_pdb`;
	#print $out;

	my $rms;
	($rms) = $out =~ /RMS: (\S+)\s+Finished/;

	# Clean up
	`rm $tmp_file`;

	return $rms;

}

sub selectRanges {
	my $report_ref = shift;
	my @list_range = @_;

	my %report = %{$report_ref};
	my %report_select;

	for my $type (@list_range) {
		#print "COMBINE: $type : $report{$type}{'begin'} : $report{$type}{'end'}\n";
		$report_select{$type}{'begin'} = $report{$type}{'begin'};
		$report_select{$type}{'end'} = $report{$type}{'end'};

	}

	return %report_select;
}

sub inverseRanges {
	my $full_range_ref = shift;
	my $range_ref = shift;

	my %full_range = %{$full_range_ref};
	my %range = %{$range_ref};

	my %inverse_range;
	my $count;
	my $type;
	my $type_rem;
	my $last_pos_rem;

	for my $full_type (sort {$a <=> $b} keys %full_range) {

		my $start = 1;

		#print "CHAIN: $full_type $full_range{$full_type}{'begin'} $full_range{$full_type}{'end'}\n";
		
		for $type (sort {$a cmp $b} keys %range) {
			next unless $type =~ /^$full_type/i;
			#print "LOOP: $type $range{$type}{'begin'} $range{$type}{'end'}\n";

			# Generate Inverse Selection
			if ($start) {
				$count++;
				#print "INVERSE: $full_range{$full_type}{'begin'} ", $range{$type}{'begin'} - 1;
				#print "\n";
				#print "INVERSE: ", $range{$type}{'end'} + 1, " ";

				$inverse_range{$full_type . "_inverse_" . $count}{'begin'} = $full_range{$full_type}{'begin'};
				$inverse_range{$full_type . "_inverse_" . $count}{'end'} = $range{$type}{'begin'} - 1;

				$count++;
				$inverse_range{$full_type . "_inverse_" . $count}{'begin'} = $range{$type}{'end'} + 1;

				$start = 0;
			}
			else {
				#print " ", $range{$type}{'begin'} - 1;
				#print "\n";
				#print "INVERSE: ", $range{$type}{'end'} + 1;

				$inverse_range{$full_type . "_inverse_" . $count}{'end'} = $range{$type}{'begin'} - 1;
				$count++;

				$last_pos_rem = $range{$type}{'end'} + 1;
				$inverse_range{$full_type . "_inverse_" . $count}{'begin'} = $last_pos_rem if ($full_range{$full_type}{'end'} > $last_pos_rem);
			}

			$type_rem = $type;

		}

		if ($type_rem =~ /^$full_type/i) {
			#print " ", $full_range{$full_type}{'end'};
			#print "\n\n";

			# TEST
			#$inverse_range{$full_type . "_inverse_" . $count}{'end'} = $full_range{$full_type}{'end'} if $full_range{$full_type}{'end'} < $range{$type}{'end'} + 1;
			$inverse_range{$full_type . "_inverse_" . $count}{'end'} = $full_range{$full_type}{'end'} if ($full_range{$full_type}{'end'} > $last_pos_rem);
		}

		$count = 0;
	}

	return %inverse_range;
}

sub genProFitCode {

	my $superimpose_range_ref = shift;
	my $rmsd_range_ref = shift;

	my $profit_delzone_all = "DELZONE ALL\n";
	my $profit_atoms = "ATOMS N,CA,C,O\n";
	my $profit_fit = "FIT\n";

	my $profit_code;

	$profit_code .= $profit_delzone_all;

	my %superimpose_range = %{$superimpose_range_ref};
	my %rmsd_range = %{$rmsd_range_ref};

	my $chain;

	# Set ZONES for superimposition
	for my $type (sort {$a cmp $b} keys %superimpose_range) {
		($chain) = $type =~ /^(\w)/;
		$chain =~ tr/a-z/A-Z/;
		$profit_code .= "ZONE $chain$superimpose_range{$type}{'begin'}-$chain$superimpose_range{$type}{'end'}:$chain$superimpose_range{$type}{'begin'}-$chain$superimpose_range{$type}{'end'}\n";
	}

	# Superimpose
	$profit_code .= $profit_atoms;
	$profit_code .= $profit_fit;

	# Set R ZONES for RMSD calculation
	for my $type (sort {$a cmp $b} keys %rmsd_range) {
		($chain) = $type =~ /^(\w)/;
		$chain =~ tr/a-z/A-Z/;
		$profit_code .= "RZONE $chain$rmsd_range{$type}{'begin'}-$chain$rmsd_range{$type}{'end'}:$chain$rmsd_range{$type}{'begin'}-$chain$rmsd_range{$type}{'end'}\n";
	}

	return $profit_code;

}

sub getChainRange {
	my $pdb_file = shift;

	open(FILE,"$pdb_file") or die "Cannot open file \n";

	my %chain_range;
	my $mer;
	my $res_num;
	my $mer_rem;
	my $res_num_rem;

	while (<FILE>) {

		chomp;

		next unless /^ATOM/;

		$res_num = substr($_,22,5);
		$res_num =~ s/ //g;

		$mer = substr($_,20,2);
		$mer =~ s/ //g;

		if ($mer_rem ne $mer) {
			$chain_range{$mer}{'begin'} = $res_num;
			$chain_range{$mer_rem}{'end'} = $res_num_rem if $res_num_rem;
		}

		$mer_rem = $mer;
		$res_num_rem = $res_num;
	}
	close(FILE);

	# Get the last end residue number for the final chain
	$chain_range{$mer_rem}{'end'} = $res_num_rem;

	# Output all of the change ranges
	#for my $chain (sort {$a cmp $b} keys %chain_range) {
		#print "$chain $chain_range{$chain}{'begin'}-$chain_range{$chain}{'end'}\n"
	#}

	return %chain_range;

}
