#!/usr/bin/perl
use strict;
use LWP::Simple;

my $study_id = shift @ARGV or die "$0 <study_id e.g. SRP008673>";

my $csv = LWP::Simple::get(qq(http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&rettype=runinfo&db=sra&term=$study_id));
die "no csv" unless $csv;

my @rec = split /\n/, $csv;
my $head = shift @rec;

#col name to number map;
my $i = 0;
my %col = map {$_=>$i++} split /,/, $head;

foreach my $rec ( @rec ) {
  my @F = split /,/, $rec;
  my ( $run, $sample, $library, $layout ) = ( $F[ $col{Run} ], $F[ $col{SampleName} ], $F[ $col{LibraryName} ], $F[ $col{LibraryLayout} ] );

  print "echo $run\n";

  if ( $layout =~ m#pair#ii ) {
    print qq(~/sratoolkit.2.8.1-3-ubuntu64/bin/fastq-dump --split-files -O . $run ;\n);
    print qq(picard-tools FastqToSam F1=${run}_1.fastq F2=${run}_2.fastq O=$run.bam SAMPLE_NAME="$sample" LIBRARY_NAME="$library" READ_GROUP_NAME="${sample}___$library" ;\n);
  }
  elsif ( $layout =~ m#single#i ) {
    print qq(~/sratoolkit.2.8.1-3-ubuntu64/bin/fastq-dump -O . $run ;\n);
    print qq(picard-tools FastqToSam F1=$run.fastq O=$run.bam SAMPLE_NAME="$sample" LIBRARY_NAME="$library" READ_GROUP_NAME="${sample}___$library" ;\n);
  }
  else {
    die "unknown layout: $layout";
  }

  print qq(gsutil -m cp $run* gs://allenday-dev/cannabis/$study_id/primary/ ;\n);
  print qq(rm -v $run* ;\n);
}
