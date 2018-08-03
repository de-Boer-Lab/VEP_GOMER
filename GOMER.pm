=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

 Ensembl <http://www.ensembl.org/info/about/contact/index.html>
    
=cut

=head1 NAME

 GOMER

=head1 SYNOPSIS

 mv GeneSplicer.pm ~/.vep/Plugins
 ./vep -i variants.vcf --plugin GeneSplicer,[path_to_genesplicer_bin],[path_to_training_dir],[option1=value],[option2=value]

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 runs GeneSplicer (https://ccb.jhu.edu/software/genesplicer/) to get
 splice site predictions.

 It evaluates a tract of sequence either side of and including the
 variant, both in reference and alternate states. The amount of
 sequence included either side defaults to 100bp, but can be modified
 by passing e.g. "context=50" as a parameter to the plugin.

 Any predicted splicing regions that overlap the variant are reported
 in the output with one of four states: no_change, diff, gain, loss

 There follows a "/"-separated string consisting of the following data:

 1) type (donor, acceptor)
 2) coordinates (start-end)
 3) confidence (Low, Medium, High)
 4) score

 Example: loss/acceptor/727006-727007/High/16.231924

 If multiple sites are predicted, their reports are separated by ",".

 For diff, the confidence and score for both the reference and alternate
 sequences is reported as REF-ALT.

 Example: diff/donor/621915-621914/Medium-Medium/7.020731-6.988368

 Several parameters can be modified by passing them to the plugin string:

 context    : change the amount of sequence added either side of
              the variant (default: 50bp)
 tmpdir     : change the temporary directory used (default: /tmp)
 cache_size : change how many sequences' scores are cached in memory
              (default: 50)

 Example: --plugin GOMER,$GS/bin/linux/genesplicer,$GS/human,context=200,tmpdir=/mytmp

 On some systems the binaries provided will not execute, but can be compiled from source:

   cd $GS/sources
   make
   cd -
   ./vep [options] --plugin GOMER,$GS/sources/genesplicer,$GS/human

 On Mac OSX the make step is known to fail; the genesplicer.cpp file requires modification:

   cd $GS/sources
   perl -pi -e "s/^main  /int main  /" genesplicer.cpp
   make
 

=cut

package GOMER;

use strict;
use warnings;

use Digest::MD5 qw(md5_hex);

use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap);

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

our %DEFAULTS = (
  context => 50,
  pseudocount => 0.001,
  scoreScale => 1,
  scoreCenter => 0,
  scoreCutoff => -999999,
  deltaCutoff => 0,
  tmpdir  => '/tmp',
  cache_size => 50,
);

sub testGomer{
	my $parent =  shift @_;
	my @input = ("M6525_1.02",  "/home/unix/cgdeboer/CIS-BP/1.02/Homo_sapiens/pwms_all_motifs/M6525_1.02.txt", "1", "0", "-999999", "0");
	my $consensus = "CGCCAAGGAGC";
	my $pwm = newPWM(\@input, $parent);
	print(sprintf("min score: %g\n",$pwm->{"conc"}));
	print(sprintf("One binding site: %g\n",scoreSeq($pwm,$consensus)));
	print(sprintf("Two binding sites: %g\n",scoreSeq($pwm,$consensus.$consensus)));
	print(sprintf("Three binding sites: %g\n",scoreSeq($pwm,$consensus.$consensus.$consensus)));
	print(sprintf("Random site: %g\n",scoreSeq($pwm,"TGGTAGCCTTA")));
	print(sprintf("Random site: %g\n",scoreSeq($pwm,"TGTGCAGAGTG")));
	print(sprintf("Random site: %g\n",scoreSeq($pwm,"TGCGCGTAGAT")));
	print(sprintf("Random seq %g\n",scoreSeq($pwm,"TGAGAGACGCGCTGTAGAGTGTGAGCTGAGCGC")));
}


our @bases = ("A","C","G","T");
our $base2i = {A => 0, C => 1, G => 2, T =>3};
sub newPWM{
	my @input = @{shift @_};
	my $parent =  shift @_;
	my $pseudocount = $parent->{'_param_pseudocount'};
	my $ID = shift @input || die "No motif ID in PWM DB file!";
	my $file = shift @input || die "No motif file path in PWM DB file!";
	my $scoreScale;
	my $scoreCenter;
	my $scoreCutoff;
	my $deltaCutoff;
	$scoreScale =  shift @input || $parent->{'_param_scoreScale'};
	$scoreCenter = shift @input || $parent->{'_param_scoreCenter'};
	$scoreCutoff = shift @input || $parent->{'_param_scoreCutoff'};
	$deltaCutoff = shift @input || $parent->{'_param_deltaCutoff'};
	open(PWMFile,$file) || die sprintf("Could not open the PWM for %s",$ID);
	my $pos =0;
	#my (@curFs, $curMin, @curFsRC);
	my @line;
	my $totF=0.0;
	my @PWM = ();
	my @rcPWM = ();
	my $conc = 0.0;
  while(<PWMFile>){
		my (@curFs, $curMin, @curFsRC);
		chomp;
		$totF=0.0;
    @line = split(/\t/,$_);
		#print((join "\t", @line)."\n");
		if ($pos>0){
			@curFs = ($line[1]+$pseudocount, $line[2]+$pseudocount, $line[3]+$pseudocount, $line[4]+$pseudocount);
			$totF += $_ for @curFs;
			@curFs = map {log(0.25/($_/$totF))} @curFs;
			$curMin = $curFs[0];
			for (@curFs) {
    		$curMin = $_ if $_ < $curMin;
			}
			$conc += $curMin;
    	push @PWM, \@curFs;
			@curFsRC = reverse @curFs;
			push @rcPWM, \@curFsRC;
    }elsif($pos==0){
			($bases[0] eq $line[1]) || print sprintf("Bad header field {%s!=%s}\n", $bases[0], $line[1]);
			($bases[1] eq $line[2]) || print sprintf("Bad header field {%s!=%s}\n", $bases[1], $line[2]);
			($bases[2] eq $line[3]) || print sprintf("Bad header field {%s!=%s}\n", $bases[2], $line[3]);
			($bases[3] eq $line[4]) || print sprintf("Bad header field {%s!=%s}\n", $bases[3], $line[4]);
			($bases[0] eq $line[1] && $bases[1] eq $line[2] && $bases[2] eq $line[3] && $bases[3] eq $line[4]) || die sprintf("Bases in %s is not in the correct order ({Pos,%s}): {%s}", $file, (join ",", @bases), (join ",", @line));
		}
		$pos=$pos+1
  }
	@rcPWM = reverse @rcPWM;
	return {conc => $conc, PWM => \@PWM, RCPWM => \@rcPWM, ID=> $ID, scoreScale => $scoreScale, scoreCenter => $scoreCenter, scoreCutoff => $scoreCutoff, deltaCutoff => $deltaCutoff, file => $file};
}

sub scoreSeq{
	my $pwm = shift;
	my $seq = shift;
	my $j;
	my $bindingScore = 0.0;
	my $conc = $pwm->{'conc'};
	my ($curKd, $curKdRC);
	#my $pwmLen = scalar( @{$pwm->{'PWM'}});
	#print ($pwmLen."\n");
	for my $i (0 .. (length($seq) - scalar( @{$pwm->{'PWM'}}))){
		$curKd=0.0;
		$curKdRC=0.0;
		for $j (0 .. scalar( @{$pwm->{'PWM'}} ) - 1){
			#my @temp = @{ $pwm->{'PWM'  } };
			#print( (join ",",@{$temp[$j]})."\n");
			#die("Quit");
			#print(sprintf("%s: %i, %i, %s; score=%g,; rc=%g\n",$seq,$i, $j, (substr $seq, $i+$j, 1), ${${$pwm->{'PWM'  }}[$j]}[$base2i->{substr $seq, $i+$j, 1}], ${${$pwm->{'RCPWM'}}[$j]}[$base2i->{substr $seq, $i+$j, 1}]));
			$curKd   += ${${$pwm->{'PWM'  }}[$j]}[$base2i->{substr $seq, $i+$j, 1}];
			$curKdRC += ${${$pwm->{'RCPWM'}}[$j]}[$base2i->{substr $seq, $i+$j, 1}];
		}
		#$bindingScore += 1.0 - (1.0/(1.0 + exp($conc-$curKd)));
		#$bindingScore += 1.0 - (1.0/(1.0 + exp($conc-$curKdRC)));
		#print(sprintf("conc = %g; kd = %g; kd(rc) = %g\n", $conc, $curKd, $curKdRC));
		#print(sprintf("cs = %g\n", (1.0/(1.0 + exp($curKd-$conc)))));
		#print(sprintf("cs(rc) = %g\n", (1.0/(1.0 + exp($curKdRC-$conc)))));
		$bindingScore += (1.0/(1.0 + exp($curKd-$conc)));
		$bindingScore += (1.0/(1.0 + exp($curKdRC-$conc)));
	}
	#print(sprintf("bs = %g\n", $bindingScore));
	return(log($bindingScore));
}

sub new {
  my $class = shift;
  #my %class = shift;

  my $self = $class->SUPER::new(@_);
  
  # we need sequence, so no offline mode unless we have FASTA
  die("ERROR: cannot function in offline mode without a FASTA file\n") if $self->{config}->{offline} && !$self->{config}->{fasta};

  my $params = $self->params;

  my $db = shift @$params;
  die("ERROR: GOMER PWM database not specified\n") unless $db;
  die("ERROR: GOMER PWM database not found\n") unless -e $db;
  $self->{_dbfile} = $db;
	
  # defaults
  $self->{'_param_'.$_} = $DEFAULTS{$_} for keys %DEFAULTS;

  # REST API passes 1 as first param
  shift @$params if $params->[0] && $params->[0] eq '1';

  # set/override with user params
  foreach my $param(@$params) {
    my ($key, $val) = split('=', $param);
    die("ERROR: Failed to parse parameter $param\n") unless defined($key) && defined($val);

    $self->{'_param_'.$key} = $val;
  }
	#testGomer($self);
	#die("quit");
	#read in DB and all PWMs
	open(DBFile,$db) || die "Could not open the PWM database file";
	my @line;
	my $curPWM;
	my @pwms;
	my @pwmIDs;
	while(<DBFile>){
		chomp;
    @line = split(/\t/,$_);
    push @pwmIDs, $line[0];
		$curPWM = newPWM(\@line, $self);
    push @pwms, $curPWM;
	}
	$self->{_pwms}=\@pwms;
	$self->{_pwmsIDs}=\@pwmIDs;
  print "GOMER found ".($#pwms+1)." PWMs\n";
	$self->{numVar} = 0;
  return $self;
}

sub feature_types {
  return ['Feature', 'Intergenic'];
}

sub get_header_info {
  return {
    GOMER => "GOMER predictions"
  };
}

sub run {
  my ($self, $tva) = @_;
	#print( (ref $tva)."\n");
	use Scalar::Util;
	#print( Scalar::Util::blessed($tva)."\n");
  my $vf = $tva->variation_feature;
  # get up and downstream sequences
  my $up_seq = $vf->{slice}->sub_Slice(
    $vf->start() - $self->{'_param_context'},
    $vf->start() - 1,
    $vf->strand
  )->seq;

  my $down_seq = $vf->{slice}->sub_Slice(
    $vf->{end} + 1,
    $vf->{end} + $self->{'_param_context'},
    $vf->strand
  )->seq;

  # create ref seq by grabbing reference TVA
  my $ref_allele = $tva->variation_feature_overlap->get_reference_VariationFeatureOverlapAllele->variation_feature_seq;
	$ref_allele = ($ref_allele eq "-") ?  "" : $ref_allele;
  my $ref_seq = $up_seq.$ref_allele.$down_seq;
	#print("Current ref seq: ".$ref_seq."\n");
  return {} unless $ref_seq =~ /^[ACGT]+$/;

  # create alt seq
  my $alt_allele = $tva->variation_feature_seq;
  $alt_allele =~ s/\-//g;
  my $alt_seq = $up_seq.$alt_allele.$down_seq;
	
  return {} unless $alt_seq =~ /^[ACGT]+$/;
	
	#scan sequences with all PWMs
	my @pwms = @{$self->{'_pwms'}};
	my @pwmIDs = @{$self->{'_pwmsIDs'}};
	my $curScoreA;
	my $curScoreR;
	my $deltaScore;
	my @results;
	for my $i (0 .. $#pwms){
		# (ID= $ID, scoreScale = $scoreScale, scoreCenter = $scoreCenter, scoreCutoff = $scoreCutoff, deltaCutoff = $deltaCutoff, file=$file)
		#$curScoreA = (`gomerScoreSeq.py -is $alt_seq -i ${pwms[$i]->{'file'}}` - $pwms[$i]->{'scoreCenter'})/$pwms[$i]->{'scoreScale'};
		#$curScoreR = (`gomerScoreSeq.py -is $ref_seq -i ${pwms[$i]->{'file'}}` - $pwms[$i]->{'scoreCenter'})/$pwms[$i]->{'scoreScale'};
		$curScoreA = (scoreSeq(${pwms[$i]},  $alt_seq) - $pwms[$i]->{'scoreCenter'})/$pwms[$i]->{'scoreScale'};
		$curScoreR = (scoreSeq(${pwms[$i]},  $ref_seq) - $pwms[$i]->{'scoreCenter'})/$pwms[$i]->{'scoreScale'};
		$deltaScore = $curScoreA - $curScoreR;
		#print (sprintf("%i, %g, %g, %g\n",$i,$deltaScore,$curScoreA,$curScoreR));
		if (abs($deltaScore) >= $pwms[$i]->{'deltaCutoff'} && ( $curScoreA >= $pwms[$i]->{'scoreCutoff'}  || $curScoreR >= $pwms[$i]->{'scoreCutoff'} )){
			#add the current hit to the results
			push @results, [$pwms[$i]->{'ID'}, $curScoreA, $curScoreR, $deltaScore];
		}
	}
	@results = sort { abs($b->[3]) <=> abs($a->[3]) } @results;
	$self->{numVar} +=1;
	print(sprintf("Tested %i variants\n", $self->{numVar}));
  my $returnVal = join(',', map {join('/', @{$_})} @results );
  return $returnVal ? { GOMER => $returnVal } : {};
}

