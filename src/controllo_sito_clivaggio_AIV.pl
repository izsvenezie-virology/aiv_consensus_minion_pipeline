#! /usr/bin/perl
#USAGE: cat DATA | perl controllo_sito_clivaggio_AIV.pl REF
#DATA: fields number 3,4,6,10 of a bam alignment
#REF: reference against which the alignment in DATA refers to

#variabili
my $ref_file = $ARGV[0];

#costanti
my $chr = 'HA';
my $start = 800;
my $end = 1200;

#memorizzo la sequenza dell'HA
my $ref = '';
my $flag = 0;
open (I, "<$ref_file");
while (my $line = <I>)
{
	chomp $line;
	if ($flag == 1)
	{
		if ($line =~ /^>/) {$flag = 0;}
		else {$ref.= &maiuscolo ($line);}
	}
	if ($line eq '>'.$chr) {$flag = 1;}
}
close I;

#localizzo con precisione il sito di clivaggio
my (%codice_genetico);
&codice_genetico (\%codice_genetico);
my @ref = split (//, $ref);
my (%frame_length, %count, %pos);
($frame_length{0}, $count{0}, $pos{0}) = &frame_len (0, \@ref, \%codice_genetico);
($frame_length{1}, $count{1}, $pos{1}) = &frame_len (1, \@ref, \%codice_genetico);
($frame_length{2}, $count{2}, $pos{2}) = &frame_len (2, \@ref, \%codice_genetico);
my $chosen_frame = -1;
foreach my $key (sort {$frame_length{$b}<=>$frame_length{$a}} keys %frame_length)
{
	if ($count{$key} == 1)
	{
		$chosen_frame = $key;
		last;
	}
}

#inizializzo e riempio le strutture per memorizzare le indel nella regione del sito di clivaggio
my (@match, @del, %ins);
for (my $i = 0; $i < $end - $start; $i++)
{
	$ins{$i}{"sum"} = 0;
	for (my $j = $i + 1; $j <= $end - $start; $j++)
	{
		$match[$i][$j] = 0;
		$del[$i][$j] = 0;
	}
}
while (my $line = <STDIN>)
{
	chomp $line;
	if ($line =~ /^$chr/)
	{
		#memorizzo singolarmente i valori nei campi
		my @line = split ("\t", $line);
		my $pos = $line[1];
		my $cigar = $line[2];
		my $seq = $line[3];
		if ($pos >= $end) {next;}

		#il cigar inizia o finisce con una H: non dovrebbe succedere, in quanto l'hard clipping viene usato negli allineamenti secondari al posto del soft clipping per evitare di rappresentare più volte la stessa sequenza
		#il cigar contiene N: non dovrebbe succedere, in quanto le skipped region appaiano nell'allineamento del cDNA contro il genoma e rappresentano gli introni
		#il cigar contiene P: non dovrebbe succedere, in quanto per avere padding è necessario avere una reference padded che, da quanto capisco, è una reference con esplicitata la posizione di eventuali inserzioni
		if ( ($cigar =~ /^\d+H/) || ($cigar =~ /\d+H$/) || ($cigar =~ /\d+N/) || ($cigar =~ /\d+P/) )
		{
			print STDERR "Found forbidden operator at \"$line\"\n";
			next;
		}

		#il cigar contiene S all'inizio
		if ($cigar =~ /^\d+S/)
		{
			$cigar =~ /^(\d+)S/;
			my $len = $1;
			$cigar =~ s/^\d+S//;
			$seq = substr ($seq, $len);
		}

		#itero dinamicamente il resto del cigar fino a identificare tutti gli operatori
		while (1)
		{
			#controllo se bloccare il ciclo
			if ( (length ($cigar) == 0) || ($cigar =~ /^\d+S$/) ) {last;}

			#estraggo le informazioni relative al primo operatore
			$cigar =~ /^(\d+)([MID=X])/;
			my $len = $1;
			my $op = $2;
			$cigar =~ s/^$len$op//;

			#trovo M o = o X
			if ( ($op eq 'M') || ($op eq '=') || ($op eq 'X') )
			{
				for (my $i = $pos - $start; $i < $pos + $len - 1 - $start; $i++)
				{
					if ( (0 <= $i) && ($i <= $end - $start) )
					{
						for (my $j = $i + 1; $j <= $pos + $len - 1 - $start; $j++)
						{
							if ( (0 <= $j) && ($j <= $end - $start) ) {$match[$i][$j]++;}
						}
					}
				}
				$pos = $pos + $len;
				$seq = substr ($seq, $len);
			}

			#trovo D
			if ($op eq 'D')
			{
				if ( ($start <= $pos - 1) && ($pos + $len <= $end) ) {$del[$pos - 1 - $start][$pos + $len - $start]++}
				$pos = $pos + $len;
			}

			#trovo I
			if ($op eq 'I')
			{
				if ( ($start <= $pos - 1) && ($pos <= $end) )
				{
					my $stringa = substr ($seq, 0, $len);
					$ins{$pos - 1 - $start}{"sum"}++;
					if (exists $ins{$pos - 1 - $start}{$stringa}) {$ins{$pos - 1 - $start}{$stringa}++;}
					else {$ins{$pos - 1 - $start}{$stringa} = 1;}
				}
				$seq = substr ($seq, $len);
			}
		}
	}
}

#cerco delezioni consenso
for (my $i = 0; $i < $end - $start; $i++)
{
	for (my $j = $i + 1; $j <= $end - $start; $j++)
	{
		if ($del[$i][$j] > $match[$i][$j])
		{
			my $sinistra = $start + $i;
			my $destra = $start + $j;
			print STDOUT "DEL\t".$sinistra."\t".$destra."\n";
		}
	}
}

#cerco inserzioni consenso
for (my $i = 0; $i < $end - $start; $i++)
{
	foreach $seq (sort {$ins{$i}{$b}<=>$ins{$i}{$a}} keys %{$ins{$i}})
	{
		if ( ($seq ne 'sum') && ($ins{$i}{$seq} > $match[$i][$i+1]) )
		{
			my $sinistra = $start + $i;
			print STDOUT "INS\t".$sinistra."\t".$seq."\n";
		}
	}
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

#codice genetico
sub codice_genetico
{
	#variabili iniziali
	my ($hash);
	($hash) = @_;

	$hash->{'TTT'} = 'F';
	$hash->{'TTC'} = 'F';
	$hash->{'TTA'} = 'L';
	$hash->{'TTG'} = 'L';
	$hash->{'TCT'} = 'S';
	$hash->{'TCC'} = 'S';
	$hash->{'TCA'} = 'S';
	$hash->{'TCG'} = 'S';
	$hash->{'TAT'} = 'Y';
	$hash->{'TAC'} = 'Y';
	$hash->{'TAA'} = '*';
	$hash->{'TAG'} = '*';
	$hash->{'TGT'} = 'C';
	$hash->{'TGC'} = 'C';
	$hash->{'TGA'} = '*';
	$hash->{'TGG'} = 'W';
	$hash->{'CTT'} = 'L';
	$hash->{'CTC'} = 'L';
	$hash->{'CTA'} = 'L';
	$hash->{'CTG'} = 'L';
	$hash->{'CCT'} = 'P';
	$hash->{'CCC'} = 'P';
	$hash->{'CCA'} = 'P';
	$hash->{'CCG'} = 'P';
	$hash->{'CAT'} = 'H';
	$hash->{'CAC'} = 'H';
	$hash->{'CAA'} = 'Q';
	$hash->{'CAG'} = 'Q';
	$hash->{'CGT'} = 'R';
	$hash->{'CGC'} = 'R';
	$hash->{'CGA'} = 'R';
	$hash->{'CGG'} = 'R';
	$hash->{'ATT'} = 'I';
	$hash->{'ATC'} = 'I';
	$hash->{'ATA'} = 'I';
	$hash->{'ATG'} = 'M';
	$hash->{'ACT'} = 'T';
	$hash->{'ACC'} = 'T';
	$hash->{'ACA'} = 'T';
	$hash->{'ACG'} = 'T';
	$hash->{'AAT'} = 'N';
	$hash->{'AAC'} = 'N';
	$hash->{'AAA'} = 'K';
	$hash->{'AAG'} = 'K';
	$hash->{'AGT'} = 'S';
	$hash->{'AGC'} = 'S';
	$hash->{'AGA'} = 'R';
	$hash->{'AGG'} = 'R';
	$hash->{'GTT'} = 'V';
	$hash->{'GTC'} = 'V';
	$hash->{'GTA'} = 'V';
	$hash->{'GTG'} = 'V';
	$hash->{'GCT'} = 'A';
	$hash->{'GCC'} = 'A';
	$hash->{'GCA'} = 'A';
	$hash->{'GCG'} = 'A';
	$hash->{'GAT'} = 'D';
	$hash->{'GAC'} = 'D';
	$hash->{'GAA'} = 'E';
	$hash->{'GAG'} = 'E';
	$hash->{'GGT'} = 'G';
	$hash->{'GGC'} = 'G';
	$hash->{'GGA'} = 'G';
	$hash->{'GGG'} = 'G';
	$hash->{'AAN'} = '-';
	$hash->{'ACN'} = '-';
	$hash->{'AGN'} = '-';
	$hash->{'ATN'} = '-';
	$hash->{'ANA'} = '-';
	$hash->{'ANC'} = '-';
	$hash->{'ANG'} = '-';
	$hash->{'ANT'} = '-';
	$hash->{'ANN'} = '-';
	$hash->{'CAN'} = '-';
	$hash->{'CCN'} = '-';
	$hash->{'CGN'} = '-';
	$hash->{'CTN'} = '-';
	$hash->{'CNA'} = '-';
	$hash->{'CNC'} = '-';
	$hash->{'CNG'} = '-';
	$hash->{'CNT'} = '-';
	$hash->{'CNN'} = '-';
	$hash->{'GAN'} = '-';
	$hash->{'GCN'} = '-';
	$hash->{'GGN'} = '-';
	$hash->{'GTN'} = '-';
	$hash->{'GNA'} = '-';
	$hash->{'GNC'} = '-';
	$hash->{'GNG'} = '-';
	$hash->{'GNT'} = '-';
	$hash->{'GNN'} = '-';
	$hash->{'TAN'} = '-';
	$hash->{'TCN'} = '-';
	$hash->{'TGN'} = '-';
	$hash->{'TTN'} = '-';
	$hash->{'TNA'} = '-';
	$hash->{'TNC'} = '-';
	$hash->{'TNG'} = '-';
	$hash->{'TNT'} = '-';
	$hash->{'TNN'} = '-';
	$hash->{'NAA'} = '-';
	$hash->{'NAC'} = '-';
	$hash->{'NAG'} = '-';
	$hash->{'NAT'} = '-';
	$hash->{'NAN'} = '-';
	$hash->{'NCA'} = '-';
	$hash->{'NCC'} = '-';
	$hash->{'NCG'} = '-';
	$hash->{'NCT'} = '-';
	$hash->{'NCN'} = '-';
	$hash->{'NGA'} = '-';
	$hash->{'NGC'} = '_';
	$hash->{'NGG'} = '_';
	$hash->{'NGT'} = '_';
	$hash->{'NGN'} = '_';
	$hash->{'NTA'} = '_';
	$hash->{'NTC'} = '_';
	$hash->{'NTG'} = '_';
	$hash->{'NTT'} = '_';
	$hash->{'NTN'} = '_';
	$hash->{'NNA'} = '_';
	$hash->{'NNC'} = '_';
	$hash->{'NNG'} = '_';
	$hash->{'NNT'} = '_';
	$hash->{'NNN'} = '_';
}

#trova la lunghezza del frame e la presenza di esattamente un GLFGAIA
sub frame_len
{
	#variabili iniziali
	my ($start, $ref, $code);
	($start, $ref, $code) = @_;

	#cerco il frame più lungo
	my $frame_length = 0;
	my $flag_start = -1;
	my $flag_end = 0;
	my @prot;
	for (my $i = $start; $i <= $#{$ref}; $i = $i + 3)
	{
		my $codone = $ref->[$i].$ref->[$i + 1].$ref->[$i + 2];
		if ( ($flag_start < 0) && ($codone eq 'ATG') ) {$flag_start = $i;}
		if ( ($flag_start >= 0) && ($flag_end == 0) )
		{
			$frame_length+= 3;
			if (exists $code->{$codone})
			{
				push (@prot, $code->{$codone});
				if ($code->{$codone} eq '*') {$flag_end = 1;}
			}
			else {push (@prot, 'X');}
		}
	}

	#conto quante volte è presente il pattern "GLFGAIA" e in che posizione
	my $count = 0;
	my $pos = 0;
	if ( ($flag_start >= 0) && ($flag_end == 1) )
	{
		for (my $i = 0; $i <= $#prot - 6; $i++)
		{
			my $pattern = $prot[$i].$prot[$i + 1].$prot[$i + 2].$prot[$i + 3].$prot[$i + 4].$prot[$i + 5].$prot[$i + 6];
			if ($pattern eq 'GLFGAIA')
			{
				$count++;
				$pos = $flag_start + 3 * $i + 1;
			}
		}
	}
	else {$frame_length = 0;}
	return ($frame_length, $count, $pos);
}
