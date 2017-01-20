#!/usr/bin/perl

use strict;
use warnings;

use Bio::DB::Taxonomy;
use Bio::SeqIO;
use Getopt::Long;
use Statistics::Descriptive;

my $debug = 0;
my $file = 'name_cache.txt';
my %name2string;
if( -f $file ) {
    open(my $fh => $file) || die $!;
    while(<$fh>) {
	my ($name,$str) = split;
	$name2string{$name} = $str;
    }
}
open(my $cache => ">>$file") || die $!;
#if( keys %name2string ) {
#    for my $name ( keys %name2string ) {
#	print $cache join("\t", $name, $name2string{$name}),"\n";
#    }
#}

my %genus_workarounds = 
    ('Ajellomyces' => 5037,
     'Artolenzites' => 526244,
     'Homoloaphlyctis' => 166479,

    );

my $taxdb = Bio::DB::Taxonomy->new(-source => 'entrez');
my ($genomedir,$gffdir) = qw(DNA GFF);
my $repeatfile = 'repeat_summary.tab';
GetOptions("dna|genome:s" => \$genomedir,
	   "gff:s"      => \$gffdir,
	   'r|repeats:s' => \$repeatfile,
	   'debug!' => \$debug,
    );

my $prefix_file = shift || "prefix.tab";
open(my $fh => $prefix_file) || die "cannot open $prefix_file: $!";

open(my $rpts => $repeatfile) || die "cannot open $repeatfile: $!";
my $rptheader = <$rpts>;
my $i = 0;
my %rptheader = map { $_ => $i++ } split(/\s+/,$rptheader);
my %repeats;
while(<$rpts>) {
    chomp;
    my @row = split(/\t/,$_);
    my ($species,$genomesize,$rptcontent) = @row;
    $species =~ s/\.fasta//;
    $species =~ s/\.\w+\.v\d+(\.\d+)?//;
    $repeats{$species}->{Total} = $rptcontent;
    $repeats{$species}->{DNA} = $row[ $rptheader{'DNA.genomefrac'} ];
    $repeats{$species}->{LTR} = $row[ $rptheader{'LTR.genomefrac'} ];
}

my $header = <$fh>;
my %table;
my @statshdr = qw(mRNA_exons mRNA exon CDS);

print join("\t", qw(PREFIX SPECIES TAX.phyla TAX.subphyla TAX.family GENOME_SIZE GC REPEAT LTRTE DNATE),
	   ( map { my $type = $_;
		   map { sprintf("%s.%s",$type,$_) }
		   qw(mean median SD count) }
	     @statshdr),
	   
    ), "\n";

while(<$fh>) {
    my ($prefix,$species) = split;
    $table{$prefix} = { fullname => $species };
    my $dnafile;
    for my $file ( "$genomedir/$species.fasta.bz2",
		   "$genomedir/$species.fasta" ) {
	if( -f $file ) {
	    $dnafile = $file;
	    last;
	}
    }
    
    my $taxstring = 'UNK;UNK;UNK';
    my ($genus) = split(/_/,$species);
    
    if( $name2string{$genus} ) {
	$taxstring = $name2string{$genus};
    } else {
	my $taxonid;
	if( $genus_workarounds{$genus} ) {
	    $taxonid = $genus_workarounds{$genus};
	} else {
	    $taxonid = $taxdb->get_taxonid($genus);
	}
	if( ! $taxonid ) {
	    warn("cannot find taxonid for $genus ($species)\n");
	} else {
	    my $node = $taxdb->get_Taxonomy_Node(-taxonid => $taxonid);
	    
	    while( my $anc = $node->ancestor ) {
		if( $anc->rank =~ /phylum|family/ ) {
		    $table{$prefix}->{taxonomy}->{$anc->rank} = $anc->scientific_name;
		}
#	print join("\t",$anc->rank, $anc->scientific_name, $anc->id), "\n";
		$node = $anc;
	    }
	    if( ! exists $table{$prefix}->{taxonomy}->{subphylum} ) {
		$table{$prefix}->{taxonomy}->{subphylum} = 
		    $table{$prefix}->{taxonomy}->{phylum}; #copy from level above
	    }
	    $taxstring =  join(";",
			       map { $table{$prefix}->{taxonomy}->{$_} } 
			       qw(phylum subphylum family));
	    print $cache join("\t", $genus, $taxstring), "\n";
	}
    }
        
    if( $dnafile ) {
#	warn "found $dnafile for $prefix\n";
    } else {
	if( $dnafile = `ls $genomedir/$species*.fasta.bz2` ) {
	    chomp($dnafile);	    
#	    warn("found by lookup $dnafile");
	} else { 
	    $dnafile = undef;
	}

	warn "Cannot find DNAfile for $prefix DNA/$species\n",next if ! $dnafile;
    }
    
    my $dnain;
    if( $dnafile =~ /\.bz2$/ ) {
	open($dnain => "bzcat $dnafile |") || die "cannot open bzcat $dnafile: $!";
    } else {
	open($dnain => "< $dnafile") || die "cannot open $dnafile: $!";
    }
    
    my $gc = 0;
    my $len = 0;
    if (! $debug ) 
    {
	my $dnaseq = Bio::SeqIO->new(-format => 'fasta', -fh => $dnain);
	while( my $dna = $dnaseq->next_seq ) {
	    $len += $dna->length;
	    $gc += ( $dna->seq =~ tr/gcGC/gcGC/);
	}
    }

# gene info
    my $genefile = $gffdir . "/".$species.".gff3.gz";
    
    if ( ! -f $genefile ) {
	$genefile = `ls $gffdir/$prefix.*.gff3.gz`;
	chomp($genefile);# remove trailing newline
	if( ! -f $genefile ) {	
	    $genefile = `ls $gffdir/$species.*gff3.gz`;
	    chomp($genefile);
	    if ( ! -f $genefile ) {
		warn("cannot find genefile for $prefix $species\n");
		$genefile = undef;
		next;
	    }
	}
    } 
   
    unless( $genefile ) {
	warn "skipping $species $prefix $species\n";
	next;
    }
    warn("using $genefile for $prefix $species\n");
    my %stats = ( 
	'mRNA_exons'  => Statistics::Descriptive::Full->new,
	'mRNA'        => Statistics::Descriptive::Full->new,
	'exon'        => Statistics::Descriptive::Full->new,
	'CDS'         => Statistics::Descriptive::Full->new );
    my %mrna_exon_counts;
    if ( ! $debug ) {
	open(my $gff => "zcat $genefile |") || die "cannot open $genefile: $!";
    
	while(<$gff>) {
	    next if /^\#/ || /^\s+$/;
	    chomp;
	    my @row = split(/\t/,$_);
	    my $last = pop @row;
	    $last =~ s/;Dbxref=\s*$//;
	    $last =~ s/;([^=]+)=;/;/g;
	    my %grp = map { split(/=/,$_) } split(/;/,$last);
		
	    if( $row[2] =~ /^(exon|CDS|mRNA)$/ ) {
		$stats{$row[2]}->add_data(abs($row[4] - $row[3]));
	    }
	    if( $row[2] eq 'exon' ) {
		if( exists $grp{Parent} ) {
		    $mrna_exon_counts{$grp{Parent}}++;
		}
	    }
	}
	$stats{'mRNA_exons'}->add_data(values %mrna_exon_counts);
	my $taxstring_sep = $taxstring;
	$taxstring_sep =~ s/;/\t/g;
	print join("\t", $prefix,$species, $taxstring,	      
		   sprintf("%.2f",$len/10**6), 
		   sprintf("%.2f",$gc / $len),
		   ( map { $repeats{$species}->{$_} || 0 } qw(Total LTR DNA)),
		   ( map { sprintf("%.3f",$stats{$_}->mean() || 0),
			   $stats{$_}->median() || 0,
			   sprintf("%.3f",$stats{$_}->standard_deviation() || 0),
			   $stats{$_}->count() || 0 } @statshdr), 
	    ),"\n";
    }
#    last if $debug;
}
    

