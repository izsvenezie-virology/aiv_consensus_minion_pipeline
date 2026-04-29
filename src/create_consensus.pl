#USAGE: cat DATA | perl create_consensus.pl REF QUAL COV CALIB DEG NAME
#DATA: fields number 3,4,6,10,11 of a bam alignment
#REF: reference against which the alignment in DATA refers to
#QUAL: minimum quality to consider a nucleotide
#COV: minimum coverage to set a defined nucleotide
#CALIB: manually replace in CIGARs a 1-bp insertion followed, after maximum CALIB match, by a 1-bp deletion (or viceversa) with a match equal to the middle one plus 1
#DEG: file in which write degeneration details, if necessary
#NAME: add this option if you want to add a common string to each segment name beside the default one

#variables
my $ref_file = $ARGV[0];
my $min_qual = $ARGV[1];
my $min_cov = $ARGV[2];
my $calib = $ARGV[3];
my $deg_file = $ARGV[4];
my $name = $ARGV[5];

#quality values
my %qual_values;
&quality_values (\%qual_values);

#reference length
my (@ref_order, %ref_length);
my $seg;
open (I, "<$ref_file");
while (my $line = <I>)
{
	chomp $line;
	if ($line =~ /^>/)
	{
		$line =~ /^>(.+)$/;
		$seg = $1;
		push (@ref_order, $seg);
		$ref_length{$seg} = 0;
	}
	else {$ref_length{$seg}+= length ($line);}
}
close I;

#initial frame
my %data;
foreach my $seg (sort(keys(%ref_length)))
{
	for (my $i = 0; $i < $ref_length{$seg}; $i++)
	{
		$data{$seg}[$i]{"A"} = 0;
		$data{$seg}[$i]{"C"} = 0;
		$data{$seg}[$i]{"G"} = 0;
		$data{$seg}[$i]{"T"} = 0;
	}
}

#putting alignment info into the frame
while (my $line = <STDIN>)
{
	chomp $line;
	my @line = split ("\t", $line);
	my $seg = $line[0];
	my $pos	= $line[1];
	my $cigar = $line[2];
	my @seq	= split (//, $line[3]);
	my @qual = split (//, $line[4]);
	my $index = 0;
	if (exists $data{$seg})
	{
		if ( ($cigar !~ /^\d+H/) && ($cigar !~ /\d+H$/) && ($cigar !~ /\d+N/) && ($cigar !~ /\d+P/) )
		{
			$cigar = &correct_cigar ($cigar, $calib);
			while (length ($cigar) > 0)
			{
				if ($cigar =~ /^\d+[SIDM=X]/)
				{
					$cigar =~ /^(\d+)([SIDM=X])/;
					my $len = $1;
					my $op = $2;
					$cigar =~ s/^$len$op//;
					if ( ($op eq 'S') || ($op eq 'I') ) {$index+= $len;}
					elsif ($op eq 'D') {$pos+= $len;}
					else
					{
						for (my $i = 0; $i < $len; $i++)
						{
							if ($pos + $i > $ref_length{$seg})
							{
								print STDERR "Operator \"$op\" goes outside reference \"$seg\", this should not happen.\n";
								last;
							}
							if (exists $qual_values{$qual[$index + $i]})
							{
								if ($qual_values{$qual[$index + $i]} >= $min_qual)
								{
									if ( ($seq[$index + $i] eq 'A') || ($seq[$index + $i] eq 'C') || ($seq[$index + $i] eq 'G') || ($seq[$index + $i] eq 'T') ) {$data{$seg}[$pos - 1 + $i]{$seq[$index + $i]}++;}
									else {print STDERR "Undefined nucleotide in the alignment \"$seq[$index + $i]\"\n";}
								}
							}
							else {print STDERR "Unknown quality value \"$qual[$index + $i]\"\n";}

						}
						$index+= $len;
						$pos+= $len;
					}
				}
			}
		}
	}
	else {print STDERR "I cannot find segment \"$seg\" in the internal map, this should not happen.\n"}
}

#consensus creation
my $stampa = '';
my %nucs;
&nucleotide (\%nucs);
for (my $i = 0; $i <= $#ref_order; $i++)
{
	print STDOUT '>'.$name.$ref_order[$i]."\n";
	for ($j = 0; $j < $ref_length{$ref_order[$i]}; $j++)
	{
		my $cov = $data{$ref_order[$i]}[$j]{'A'} + $data{$ref_order[$i]}[$j]{'C'} + $data{$ref_order[$i]}[$j]{'G'} + $data{$ref_order[$i]}[$j]{'T'};
		if ($cov >= $min_cov)
		{
			my $bases = '';
			my $count = 0;
			if (4 * $data{$ref_order[$i]}[$j]{'A'} >= $cov) {$bases.= 'A'; $count++;}
			if (4 * $data{$ref_order[$i]}[$j]{'C'} >= $cov) {$bases.= 'C'; $count++;}
			if (4 * $data{$ref_order[$i]}[$j]{'G'} >= $cov) {$bases.= 'G'; $count++;}
			if (4 * $data{$ref_order[$i]}[$j]{'T'} >= $cov) {$bases.= 'T'; $count++;}
			if (exists $nucs{$bases}) {print STDOUT $nucs{$bases};}
			else
			{
				print STDOUT 'N';
				print STDERR "I cannot find bases composition \"$bases\" in degeneration at 0-index position $j of segment \"$name.$ref_order[$i]\", this should not happen.\n";
			}
			if ($count > 1)
			{
				my $pos = $j + 1;
				$stampa.= $name.$ref_order[$i]."\t".$pos."\t".$nucs{$bases};
				$stampa.= "\t".sprintf("%.1f", $data{$ref_order[$i]}[$j]{'A'} / $cov);
				$stampa.= "\t".sprintf("%.1f", $data{$ref_order[$i]}[$j]{'C'} / $cov);
				$stampa.= "\t".sprintf("%.1f", $data{$ref_order[$i]}[$j]{'G'} / $cov);
				$stampa.= "\t".sprintf("%.1f", $data{$ref_order[$i]}[$j]{'T'} / $cov)."\n";
			}
		}
		else {print STDOUT 'N';}
	}
	print STDOUT "\n";
}
if ($stampa)
{
	open (O, ">$deg_file");
	print O "Segment\tPosition\tDegeneration\t%A\t%C\t%G\t%T\n".$stampa;
	close O;
}

#quality values
sub quality_values
{
	my ($hash);
	($hash) = @_;
	$hash->{'!'} = '0';
	$hash->{'"'} = '1';
	$hash->{'#'} = '2';
	$hash->{'$'} = '3';
	$hash->{'%'} = '4';
	$hash->{'&'} = '5';
	$hash->{'\''} = '6';
	$hash->{'('} = '7';
	$hash->{')'} = '8';
	$hash->{'*'} = '9';
	$hash->{'+'} = '10';
	$hash->{','} = '11';
	$hash->{'-'} = '12';
	$hash->{'.'} = '13';
	$hash->{'/'} = '14';
	$hash->{'0'} = '15';
	$hash->{'1'} = '16';
	$hash->{'2'} = '17';
	$hash->{'3'} = '18';
	$hash->{'4'} = '19';
	$hash->{'5'} = '20';
	$hash->{'6'} = '21';
	$hash->{'7'} = '22';
	$hash->{'8'} = '23';
	$hash->{'9'} = '24';
	$hash->{':'} = '25';
	$hash->{';'} = '26';
	$hash->{'<'} = '27';
	$hash->{'='} = '28';
	$hash->{'>'} = '29';
	$hash->{'?'} = '30';
	$hash->{'@'} = '31';
	$hash->{'A'} = '32';
	$hash->{'B'} = '33';
	$hash->{'C'} = '34';
	$hash->{'D'} = '35';
	$hash->{'E'} = '36';
	$hash->{'F'} = '37';
	$hash->{'G'} = '38';
	$hash->{'H'} = '39';
	$hash->{'I'} = '40';
	$hash->{'J'} = '41';
	$hash->{'K'} = '42';
	$hash->{'L'} = '43';
	$hash->{'M'} = '44';
	$hash->{'N'} = '45';
	$hash->{'O'} = '46';
	$hash->{'P'} = '47';
	$hash->{'Q'} = '48';
	$hash->{'R'} = '49';
	$hash->{'S'} = '50';
	$hash->{'T'} = '51';
	$hash->{'U'} = '52';
	$hash->{'V'} = '53';
	$hash->{'W'} = '54';
	$hash->{'X'} = '55';
	$hash->{'Y'} = '56';
	$hash->{'Z'} = '57';
	$hash->{'['} = '58';
	$hash->{'\\'} = '59';
	$hash->{']'} = '60';
	$hash->{'^'} = '61';
	$hash->{'_'} = '62';
	$hash->{'`'} = '63';
	$hash->{'a'} = '64';
	$hash->{'b'} = '65';
	$hash->{'c'} = '66';
	$hash->{'d'} = '67';
	$hash->{'e'} = '68';
	$hash->{'f'} = '69';
	$hash->{'g'} = '70';
	$hash->{'h'} = '71';
	$hash->{'i'} = '72';
	$hash->{'j'} = '73';
	$hash->{'k'} = '74';
	$hash->{'l'} = '75';
	$hash->{'m'} = '76';
	$hash->{'n'} = '77';
	$hash->{'o'} = '78';
	$hash->{'p'} = '79';
	$hash->{'q'} = '80';
	$hash->{'r'} = '81';
	$hash->{'s'} = '82';
	$hash->{'t'} = '83';
	$hash->{'u'} = '84';
	$hash->{'v'} = '85';
	$hash->{'w'} = '86';
	$hash->{'x'} = '87';
	$hash->{'y'} = '88';
	$hash->{'z'} = '89';
	$hash->{'{'} = '90';
	$hash->{'|'} = '91';
	$hash->{'}'} = '92';
	$hash->{'~'} = '93';
}

#correct cigar
sub correct_cigar
{
	my ($cigar, $calib);
	($cigar, $calib) = @_;
	my (@len, @op);
	while ($cigar)
	{
		$cigar =~ /^(\d+)([SIDM=X])/;
		my $len = $1;
		my $op = $2;
		push (@len, $len);
		push (@op, $op);
		$cigar =~ s/^$len$op//;
	}
	my $index = 0;
	while ($index <= $#len)
	{
		my $flag = 0;
		if ( ( ($op[$index] eq 'I') || ($op[$index] eq 'D') ) && ($len[$index] == 1) && ($index <= $#len - 2) )
		{
			if ( ( ($op[$index + 1] eq 'M') || ($op[$index + 1] eq '=') || ($op[$index + 1] eq 'X') ) && ($len[$index + 1] <= $calib) )
			{
				if ( ( ( ($op[$index] eq 'I') && ($op[$index + 2] eq 'D') ) || ( ($op[$index] eq 'D') && ($op[$index + 2] eq 'I') ) ) && ($len[$index + 2] == 1) )
				{
					$flag = 1;
					my $new_len = $len[$index + 1] + 1;
					$cigar.= $new_len.$op[$index + 1];
					$index+= 2;
				}
			}
		}
		if ($flag == 0) {$cigar.= $len[$index].$op[$index];}
		$index++;
	}
	return ($cigar);
}

#base composition
sub nucleotide
{
	my ($hash);
	($hash) = @_;
	$hash->{'A'} = 'A';
	$hash->{'C'} = 'C';
	$hash->{'G'} = 'G';
	$hash->{'T'} = 'T';
	$hash->{'AC'} = 'M';
	$hash->{'AG'} = 'R';
	$hash->{'AT'} = 'W';
	$hash->{'CG'} = 'S';
	$hash->{'CT'} = 'Y';
	$hash->{'GT'} = 'K';
	$hash->{'ACG'} = 'V';
	$hash->{'ACT'} = 'H';
	$hash->{'AGT'} = 'D';
	$hash->{'CGT'} = 'B';
	$hash->{'ACGT'} = 'N';
}
