#! /usr/bin/perl

#USAGE: cat DATA | perl filtra_blast_tab.pl OUTPUT1 OUTPUT2 OUTPUT3
#DATA: tab file 2 columns (hit, title)
#OUTPUT1: file .entry_title2hits
#OUTPUT2: file .entry_hit
#OUTPUT3: file .entry_title

#variabili
my $output_file1 = $ARGV[0];
my $output_file2 = $ARGV[1];
my $output_file3 = $ARGV[2];

#memorizzo il dato in standard input
my (%entry_title2hits, %hit2count, %title2count);
while (my $line = <STDIN>)
{
	$line =~ /^([^\t]+)\t([^\t]+)\n$/;
	my $hit = $1;
	my $title = $2;
	if (exists $hit2count{$hit}) {$hit2count{$hit}++;}
	else {$hit2count{$hit} = 1;}
	if (exists $title2count{$title}) {$title2count{$title}++;}
	else {$title2count{$title} = 1;}
	$entry_title2hits{$title}{$hit} = 1;
}

#scrivo su file i dati raccolti
open (O, ">$output_file1");
foreach my $title (sort(keys(%entry_title2hits)))
{
	my $stampa = $title."\t";
	foreach my $hit (sort(keys(%{$entry_title2hits{$title}}))) {$stampa.= $hit.';';}
	chop $stampa;
	print O $stampa."\n";
}
close O;
open (O, ">$output_file2");
foreach my $hit (sort { $hit2count{$b} <=> $hit2count{$a} } keys %hit2count) {print O $hit2count{$hit}."\t".$hit."\n";}
close O;
open (O, ">$output_file3");
foreach my $title (sort { $title2count{$b} <=> $title2count{$a} } keys %title2count) {print O $title2count{$title}."\t".$title."\n";}
close O;
