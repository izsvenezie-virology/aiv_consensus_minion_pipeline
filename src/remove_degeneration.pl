#! /usr/bin/perl

#USAGE: cat FASTA | perl remove_degeneration.pl WINDOW
#FASTA: sequenze in formato fasta dalle quali eliminare le degenerazioni
#WINDOW: ampiezza massima della finestra nella quale contare le basi relative a quella degenerata

#variabili
my $win = $ARGV[0];

#degenerazioni
my %deg;
&degenerazioni (\%deg);

#memorizzo la sequenza in standard input
my (@seq);
my (%fasta);
my $cromo;
while (my $line = <STDIN>)
{
	chomp $line;
	if ($line =~ /^>/)
	{
		$line =~ /^>(.+)$/;
		$cromo = $1;
		push (@seq, $cromo);
		$fasta{$cromo} = "";
	}
	else {$fasta{$cromo}.= &maiuscolo ($line);}
}

#stampo a video la sequenza, eliminando le degenerazioni in modo da massimizzare la diversità entro la finestra
for (my $i = 0; $i <= $#seq; $i++)
{
	print STDOUT ">$seq[$i]\n";
	my @fasta = split (//, $fasta{$seq[$i]});
	for (my $j = 0; $j <= $#fasta; $j++)
	{
		if ( ($fasta[$j] eq 'A') || ($fasta[$j] eq 'C') || ($fasta[$j] eq 'G') || ($fasta[$j] eq 'T') || ($fasta[$j] eq 'N') ) {print STDOUT $fasta[$j];}
		else
		{
			#creo la struttura con i conteggi
			my (%count);
			if (exists $deg{$fasta[$j]})
			{
				for (my $k = 0; $k <=$#{$deg{$fasta[$j]}}; $k++) {$count{${$deg{$fasta[$j]}}[$k]} = 0;}
			}
			else {print STDERR "Non trovo la lista di potenziali basi per il simbolo \"$fasta[$j]\".\n";}

			#partendo dalle basi adiacenti, conto quante volte trovo ogni base potenziale
			my $flag = 0;
			for (my $k = 1; $k <= $win; $k++)
			{
				#base a sinistra
				if ($j - $k >= 0) {if (exists $count{$fasta[$j - $k]}) {$count{$fasta[$j - $k]}++;}}

				#base a destra
				if ($j + $k <= $#fasta) {if (exists $count{$fasta[$j + $k]}) {$count{$fasta[$j + $k]}++;}}

				#se una potenziale base risulta minore delle altre, la scelgo e blocco il ciclo
				my $min_count;
				foreach my $base (sort { $count{$a} <=> $count{$b} } keys %count) {$min_count = $count{$base}; last;}
				my $flag2 = 0;
				my $nuc;
				foreach my $base2 (sort(keys(%count)))
				{
					if ($count{$base2} == $min_count)
					{
						$flag2++;
						$nuc = $base2;
					}
				}
				if ($flag2 == 1)
				{
					print STDOUT $nuc;
					$flag = 1;
					last;
				}
			}
			if ($flag == 0) {print STDOUT ${$deg{$fasta[$j]}}[0];}
		}
	}
	print STDOUT "\n";
}

#accerta che tutti i nucleotidi siano scritti in maiuscolo
sub maiuscolo
{
	#variabili iniziali
	my ($seq);
	($seq) = @_;

	#converte gli eventuali nucleotidi scritti in minuscolo in maiscolo
	my $seq_m = '';
	my @seq = split (//, $seq);
	for (my $i = 0; $i <= $#seq; $i++)
	{
		if ( ($seq[$i] eq 'A') || ($seq[$i] eq 'a') ) {$seq_m.= 'A'; next;}
		if ( ($seq[$i] eq 'C') || ($seq[$i] eq 'c') ) {$seq_m.= 'C'; next;}
		if ( ($seq[$i] eq 'G') || ($seq[$i] eq 'g') ) {$seq_m.= 'G'; next;}
		if ( ($seq[$i] eq 'T') || ($seq[$i] eq 't') ) {$seq_m.= 'T'; next;}
		if ( ($seq[$i] eq 'N') || ($seq[$i] eq 'n') ) {$seq_m.= 'N'; next;}
		if ( ($seq[$i] eq 'R') || ($seq[$i] eq 'r') ) {$seq_m.= 'R'; next;}
		if ( ($seq[$i] eq 'Y') || ($seq[$i] eq 'y') ) {$seq_m.= 'Y'; next;}
		if ( ($seq[$i] eq 'S') || ($seq[$i] eq 's') ) {$seq_m.= 'S'; next;}
		if ( ($seq[$i] eq 'W') || ($seq[$i] eq 'w') ) {$seq_m.= 'W'; next;}
		if ( ($seq[$i] eq 'K') || ($seq[$i] eq 'k') ) {$seq_m.= 'K'; next;}
		if ( ($seq[$i] eq 'M') || ($seq[$i] eq 'm') ) {$seq_m.= 'M'; next;}
		if ( ($seq[$i] eq 'B') || ($seq[$i] eq 'b') ) {$seq_m.= 'B'; next;}
		if ( ($seq[$i] eq 'D') || ($seq[$i] eq 'd') ) {$seq_m.= 'D'; next;}
		if ( ($seq[$i] eq 'H') || ($seq[$i] eq 'h') ) {$seq_m.= 'H'; next;}
		if ( ($seq[$i] eq 'V') || ($seq[$i] eq 'v') ) {$seq_m.= 'V'; next;}
		if ( ($seq[$i] eq '.') || ($seq[$i] eq '-') ) {$seq_m.= $seq[$i]; print STDERR "ATTENZIONE: ho trovato un gap nella reference\n"; next;}
		print STDERR "ATTENZIONE: ho trovato un nucleotide diverso da quelli IUPAC: $seq[$i]\n";
	}
	return $seq_m;
}

#degenerazioni
sub degenerazioni
{
	#variabili iniziali
	my ($hash);
	($hash) = @_;

	$hash->{'R'}->[0] = 'A';
	$hash->{'R'}->[1] = 'G';
	$hash->{'Y'}->[0] = 'C';
	$hash->{'Y'}->[1] = 'T';
	$hash->{'S'}->[0] = 'C';
	$hash->{'S'}->[1] = 'G';
	$hash->{'W'}->[0] = 'A';
	$hash->{'W'}->[1] = 'T';
	$hash->{'K'}->[0] = 'G';
	$hash->{'K'}->[1] = 'T';
	$hash->{'M'}->[0] = 'A';
	$hash->{'M'}->[1] = 'C';
	$hash->{'B'}->[0] = 'C';
	$hash->{'B'}->[1] = 'G';
	$hash->{'B'}->[2] = 'T';
	$hash->{'D'}->[0] = 'A';
	$hash->{'D'}->[1] = 'G';
	$hash->{'D'}->[2] = 'T';
	$hash->{'H'}->[0] = 'A';
	$hash->{'H'}->[1] = 'C';
	$hash->{'H'}->[2] = 'T';
	$hash->{'V'}->[0] = 'A';
	$hash->{'V'}->[1] = 'C';
	$hash->{'V'}->[2] = 'G';
}
