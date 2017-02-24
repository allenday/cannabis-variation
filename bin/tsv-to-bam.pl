#!/usr/bin/perl
use strict;
use Data::Dumper qw(Dumper);

my %sample = ();

<>; #header
while ( my $line = <> ) {
  chomp $line;

  my @F = split /\t/, $line;
  my $project_nick = $F[0];
  my $project_id = $F[1];
  my $biosample_id = $F[2];
  my $sample_id = $F[3];
  my $run_id = $F[4];
  my $library_type = $F[5];
  my $paired = $F[6] eq 'PAIRED' ? 1 : 0;
  my $sample_name = $F[9];
  my $library_name = $F[10];
  my $fq1 = $F[14];
  my $fq2 = $F[16];
  my $rg = $F[18];
  my $interleaved = $F[19];

  if ( $fq1 =~ m#_2.fastq$# ) {
    $fq2 = $fq1;
    $fq1 = undef;
  }

  my $srg = "$sample_id\__$rg";

  $sample{ $srg }{ paired } = $paired;
  $sample{ $srg }{ interleaved } = $interleaved;
  push @{ $sample{ $srg }{ fq1 } }, $fq1 if $fq1;
  push @{ $sample{ $srg }{ fq2 } }, $fq2 if $fq2;
  $sample{ $srg }{ project_id } = $project_id;
  $sample{ $srg }{ run1_id } = $run_id if $fq1;
  $sample{ $srg }{ run2_id } = $run_id if $fq2;
  $sample{ $srg }{ sample_id } = $sample_id;

#die $srg . "\n" . Dumper($sample{$srg}) if $line =~ m#SRR2584300#;
#SRS1107973__1

}

#SRS1107973__1
#$VAR1 = {
#          'interleaved' => '0',
#          'project_id' => 'SRP064442',
#          'sample_id' => 'SRS1107973',
#          'paired' => 0,
#          'run1_id' => 'SRR2584300',
#          'fq1' => [
#                     'gs://allenday-dev/cannabis/SRP064442/primary/SRR2584300_1.fastq'
#                   ]
#        };

foreach my $srg ( sort keys %sample ) {
  print "### $sample{ $srg }{ project_id }\t$srg\n";

  my $interleaved = $sample{ $srg }{ interleaved };
  
  my $count_fq1 = 0;
  if ( defined( $sample{ $srg }{ fq1 } ) ) {
    $count_fq1 = scalar( @{ $sample{ $srg }{ fq1 } } );
  }
  my $count_fq2 = 0;
  if ( defined( $sample{ $srg }{ fq2 } ) ) {
    $count_fq2 = scalar( @{ $sample{ $srg }{ fq2 } } );
  }

  if ( $count_fq1 > 1 ) {
    die "\tcount(fq1) > 1 $srg";
  }
  if ( $count_fq2 > 1 ) {
    die "\tcount(fq2) > 1 $srg";
  }
  if ( $sample{ $srg }{ paired } == 1 && $count_fq2 == 0 && !$interleaved ) {
    die "\tpaired missing fq2 $srg";
  }
  if ( $sample{ $srg }{ paired } == 1 && $count_fq1  > 1 ) {
    die "\tpaired too many fq1 $srg";
  }
  if ( $sample{ $srg }{ paired } == 1 && $count_fq2  > 1 ) {
    die "\tpaired too many fq2 $srg";
  }
  if ( $sample{ $srg }{ paired } == 0 && $count_fq1  > 1 ) {
    die "\tsingle too many fq1 $srg";
  }
  if ( $sample{ $srg }{ paired } == 0 && $count_fq2  > 0 ) {
    die "\tsingle too many fq2 $srg";
  }

  my ( $fq1 ) = $count_fq1 ? @{ $sample{ $srg }{ fq1 } } : ();
  my ( $fq2 ) = $count_fq2 ? @{ $sample{ $srg }{ fq2 } } : ();

  print "## FQ1: $fq1\n";
  print "## FQ2: $fq2\n";

  my ( $local_fq1 ) = $fq1 =~ m#.+/(.+?)$#; 
  my ( $local_fq2 ) = $fq2 =~ m#.+/(.+?)$#; 

  my $run1_id = $sample{ $srg }{ run1_id };
  my $run2_id = $sample{ $srg }{ run2_id };

  my $sample_id = $sample{ $srg }{ sample_id };
  my ( $read_group ) = $srg =~ m#.+__(\d+)$#;

  print "## RG: $read_group\n";

  if ( $sample{ $srg }{ paired } ) {
    my $interleaved_option = $interleaved ? "-p" : "";
    print qq(gsutil -m cp $fq1 ./$local_fq1;\n);
    print qq(gsutil -m cp $fq2 ./$local_fq2;\n) unless $interleaved;

    print qq(sed -i 's,$run1_id,$sample_id,g' $local_fq1;\n);
    print qq(sed -i 's,$run2_id,$sample_id,g' $local_fq2;\n) if $fq2;

    print qq(bwa mem -R '\@RG\\tID:$srg\\tSM:$sample_id' $interleaved_option -t 30 /data/cannabis/reference-cannatonic/MNPR01.fa $local_fq1 $local_fq2 > $srg.pre.sam;\n);
    print qq(samtools view -Sb -@ 30 $srg.pre.sam > $srg.pre.bam;\n);
    print qq(samtools sort -@ 30 $srg.pre.bam $srg;\n);

    print qq(rm -v $srg.pre.sam;\n);
    print qq(rm -v $srg.pre.bam;\n);
    print qq(rm -v $local_fq1;\n);
    print qq(rm -v $local_fq2;\n) unless $interleaved;
  }
  else {
    print qq(## single\n);
    print qq(gsutil -m cp $fq1 ./$local_fq1;\n);

    print qq(sed -i 's,$run1_id,$sample_id,g' $local_fq1;\n);

    print qq(bwa mem -R '\@RG\\tID:$srg\\tSM:$sample_id' -t 30 /data/cannabis/reference-cannatonic/MNPR01.fa $local_fq1 > $srg.pre.sam;\n);
    print qq(samtools view -Sb -@ 30 $srg.pre.sam > $srg.pre.bam;\n);
    print qq(samtools sort -@ 30 $srg.pre.bam $srg;\n);

    print qq(rm -v $srg.pre.sam;\n);
    print qq(rm -v $srg.pre.bam;\n);
    print qq(rm -v $local_fq1;\n);
  }
}

__DATA__
Project Nick	Project	BioSample_s	SRA_Sample_s	Run_s	LibrarySource_s	LibraryLayout_s	InsertSize_l	MBases_l	Sample_Name_s	Library_Name_s	size_ratio			FQ1	FQ1_size	FQ2	FQ2_size	RG	
