use strict;

my $root_dir = "/Users/ashishkabra/Dropbox/Ashish_Karan/Sustainable\ Transportation\ System/Data/Code/R/Model_v22";


my @in_file_list = ("r70.5_noRD_bytw_metrodummies_logmetrocorrected_moments_1:10_full_nldis_pr3.txt",
	"r70.5_noRD_bytw_metrodummies_logmetrocorrected_moments_1:10_full_nldis_pr6.txt",
	"r70.5_noRD_bytw_metrodummies_logmetrocorrected_moments_1:10_full_stepf_2_pr3.txt",
	"r70.5_noRD_bytw_metrodummies_logmetrocorrected_moments_1:10_full_stepf_2_pr6.txt",
	"r70.5_noRD_bytw_metrodummies_logmetrocorrected_moments_1:10_full_stepf_pr3.txt",
	"r70.5_noRD_bytw_metrodummies_logmetrocorrected_moments_1:10_full_stepf_pr6.txt",
	"r70.5_noRD_bytw_metrodummies_logmetrocorrected_moments_1:20_full_nldis_pr3.txt",
	"r70.5_noRD_bytw_metrodummies_logmetrocorrected_moments_1:20_full_nldis_pr6.txt",
	"r70.5_noRD_bytw_metrodummies_logmetrocorrected_moments_1:20_full_stepf_2_pr3.txt",
	"r70.5_noRD_bytw_metrodummies_logmetrocorrected_moments_1:20_full_stepf_2_pr6.txt",
	"r70.5_noRD_bytw_metrodummies_logmetrocorrected_moments_1:20_full_stepf_pr3.txt",
	"r70.5_noRD_bytw_metrodummies_logmetrocorrected_moments_1:20_full_stepf_pr6.txt");

open(OUTFILE, '>temp.csv') or die "no file";

for my $in_file (@in_file_list) {
	#my $in_file = "r78.1_noRD_bytw_metrodummies_logmetrocorrected_moments_1:20_full_stepf_pr3.txt";
	print OUTFILE $in_file,"\n";
	my $file_name = "$root_dir/"."$in_file";
	open(INFILE, "<$file_name") or die "no file";
	my @lines = <INFILE>;

	my $first_coef = undef;

	for(my $line_no=0; $line_no<=$#lines; $line_no++) {
		if($lines[$line_no] =~ /Step 1 output:/) {
			#print $line_no;
			#get theta_nonden string from $line_no+2		
			$line_no = $line_no + 2;
			$first_coef = get_first_coef($lines[$line_no]);
			#print $first_coef;	
			$line_no++; #goto next line
		}
		#find line for $first_coef
		if(defined($first_coef)) {
			if($lines[$line_no] =~ /$first_coef/) {
				print $lines[$line_no],"\n";				
				print OUTFILE "Step 1,t-stat,Step 2,t-stat","\n";
				while($line_no <= $#lines) {
					print OUTFILE $lines[$line_no];
					$line_no++;
				}
				
			}
		}
	}
	close(INFILE);
	print OUTFILE "\n\n";
}

sub get_first_coef {
	my $text = shift;
	my (undef, $theta_nonden) = split('\s', $text, 2);
	#remove quotes
	$theta_nonden =~ s/"//;
	#print join(':', split(/, /, $theta_nonden,2)), "\n";
	my ($theta_nonden_1, undef) = split(', ', $theta_nonden, 2);
	#print $theta_nonden_1,"\n";
	return($theta_nonden_1);
}

1;


