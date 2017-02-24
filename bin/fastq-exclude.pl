#!/usr/bin/perl
use strict;
my $exclude_list = shift @ARGV;
die "no fastq" unless scalar @ARGV;
open E, $exclude_list or die $!;
my %exclude = ();
while ( my $x = <E> ) {
  chomp $x;
  ( $x ) =~ m#^@(\S+)#;
  $x =~ s#/[12]$##;
  $exclude{ $x } = 1;
}
close E;
foreach my $fastq ( @ARGV ) {
  open F, $fastq or die $!;
  open O, ">$fastq.remap.fastq" or die $!;
  while ( 1 ) {
    my $w = <F>; #header
    chomp $w;
    last unless $w;
    my $x = <F>;
    my $y = <F>;
    my $z = <F>;
    my ( $id ) = $w =~ m#^@(\S+)#;
    $id =~ s#/[12]$##;
    next if $exclude{ $id };
    print O "$w\n$x$y$z";
  }
  close O;
  close F;
}
