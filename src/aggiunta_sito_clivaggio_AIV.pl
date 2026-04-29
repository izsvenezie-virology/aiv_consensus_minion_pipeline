#! /usr/bin/perl
#USO: cat FASTA | perl aggiunta_sito_clivaggio_AIV.pl FILE
#FASTA: file multifasta con la sequenza presente su singola riga
#FILE: file ouptut of controllo_sito_clivaggio_AIV.pl

#variabili
$indel_file = $ARGV[0];

#costanti
my $chr = 'HA';

#memorizzo le modifiche da fare attorno al sito di clivaggio
my (%del, %ins);
open (I, "<$indel_file");
while (my $line = <I>)
{
	chomp $line;
	my @line = split ("\t", $line);
	if ($line[0] eq 'DEL')
	{
		if ( ($line[2] + 1  - $line[1]) % 3 == 0)
		{
			if (exists $del{$line[1]}) {print STDERR "Trovata un'altra delezione consenso a partire dalla posizione ". $line[1].", ignorata\n";}
			else {print STDERR "Inserita delezione tra le posizioni ".$line[1]."-".$line[2]."\n"; $del{$line[1]} = $line[2];}
		}
	}
	elsif ($line[0] eq 'INS')
	{
		if (length($line[2]) % 3 == 0)
		{
			if (exists $ins{$line[1]}) {print STDERR "Trovata un'altra inserzione consenso adiacente alla posizione ". $line[1].", ignorata\n";}
			else {$ins{$line[1]} = $line[2]; print STDERR "Inserita inserzione \"$line[2]\" adiacente alla posizione ".$line[1]."\n";}
		}
	}
}
close I;

#eseguo le modifiche sullo standard input
while (my $line = <STDIN>)
{
	if ($line =~ /_$chr\n$/)
	{
		print STDOUT $line;
		my $seq = <STDIN>;
		chomp $seq;
		my @seq = split (//, $seq);
		my $pos = 1;
		while ($pos <= length ($seq))
		{
			if (exists $del{$pos}) {$pos = $del{$pos};}
			else
			{
				print STDOUT $seq[$pos - 1];
				if (exists $ins{$pos}) {print STDOUT $ins{$pos};}
			}
			$pos++;
		}
		print STDOUT "\n";
	}
	else {print STDOUT $line;}
}
