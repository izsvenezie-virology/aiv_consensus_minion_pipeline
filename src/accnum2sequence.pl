#variabili
my $root = $ARGV[0];
my $fasta = $ARGV[1];

#leggo lo standard input: collego il segmento di influenza al best title
my $special_word = 'NON-SENSE_WORD_LIKE_ASH_NAZG_DURBATULUK,_ASH_NAZG_GIMBATUL,_ASH_NAZG_THRAKATULUK,_AGH_BURZUM-ISHI_KRIMPATUL';
my %seg2title;
while (my $line = <STDIN>)
{
	$line =~ /^(\S+)\t(\S+)\n$/;
	my $seg = $1;
	my $title = $2;
	if ($title ne $special_word) {$seg2title{$seg} = $title;}
}

#collego il title alle hits
my %title2hits;
&title2hits($root, \%title2hits);

#collego la hit alla sua conta
my %hit2count;
&hit2count($root, \%hit2count);

#collego il segmento al best hit
my (%seg2hit, %hit2seg);
foreach my $seg (sort(keys(%seg2title)))
{
	my $title = $seg2title{$seg};
	if (exists $title2hits{$title})
	{
		my @hits = split (/;/, $title2hits{$title});
		my %hits;
		for (my $i = 0; $i <= $#hits; $i++)
		{
			if (exists $hit2count{$hits[$i]}) {$hits{$hits[$i]} = $hit2count{$hits[$i]};}
			else {print STDERR "Hit \"$hits[$i]\" is not present in file .hit2count\n";}
		}
		foreach my $best_hit (sort { $hits{$b} <=> $hits{$a} } keys %hits)
		{
			$seg2hit{$seg} = $best_hit;
			if (exists $hit2seg{$best_hit}) {print STDERR "ATTENZIONE: l'accession number \"$best_hit\" è associato a più di un segmento di influenza!\n";}
			else {$hit2seg{$best_hit} = $seg;}
			last;
		}
	}
	else {print STDERR "Title \"$title\" is not present in file .entry_title2hits\n";}
}

#scarico le sequenze dal fasta del database e creo la reference
my $acc_num;
my %stampa;
my $flag;
open (I, "<$fasta");
while (my $line = <I>)
{
	if ($line =~ /^>/)
	{
		$line =~ /^>([^ ]+) /;
		$acc_num = $1;
		if (exists $hit2seg{$acc_num})
		{
			$stampa{$hit2seg{$acc_num}} = '';
			$flag = 1;			
		}
		else {$flag = 0;}
	}
	else
	{
		if ($flag == 1)
		{
			chomp $line;
			$stampa{$hit2seg{$acc_num}}.= $line;
		}
	}
}
close I;
foreach my $seg (sort(keys(%stampa))) {print STDOUT '>'.$seg."\n".$stampa{$seg}."\n";}

#collega il title al hit
sub title2hits
{
	#variabili iniziali
	my ($root, $title2hits);
	($root, $title2hits) = @_;

	my $file = $root.'.entry_title2hits';
	open (I, "<$file");
	while (my $line = <I>)
	{
		$line =~ /^(\S+)\t(\S+)+\n$/;
		my $title = $1;
		my $hits = $2;
		$title2hits->{$title} = $hits;
	}
	close I;
}

#collega la hit alla sua conta
sub hit2count
{
	#variabili iniziali
	my ($root, $hit2count);
	($root, $hit2count) = @_;

	my $file = $root.'.entry_hit';
	open (I, "<$file");
	while (my $line = <I>)
	{
		$line =~ /^(\d+)\t(\S+)\n$/;
		my $count = $1;
		my $hit = $2;
		$hit2count->{$hit} = $count;
	}
	close I;
}
