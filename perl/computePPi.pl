#!/usr/bin/perl

## Macd location (get from my github)
use lib '/home/oliver/GIT_REPOS/macd/lib';
use Macd;
use FindBin qw($Bin);
use File::Path qw(make_path);
use File::Find;
use strict;
use Data::Dumper;

## get these from environment variables

my $HOME = $ENV{HOME} || die "NO \$HOME variable set\n";
my $GRPATH = $ENV{GRPATH} || die "NO \$GRPATH variable set\n";
my $R_LIBS = $ENV{R_LIBS} || die "NO \$R_LIBS variable\n";
my $CHIGP_PATH = "$ENV{GRPATH}/CHIGP/";

## site specific conf file for macd - see http://github.com/ollyburren/macd
## for more details
my $grid_cnf = "$GRPATH/macd/example/ini/example.cnf";

my $r_lib_dir = $R_LIBS;
	
my $DRIVER = Macd::GRIDDriverFactory->instantiate('SGE',inifile => $grid_cnf);


###############
#CONFIGURATION#
###############
## where to mail to once finished.
my $BASE_DIR = "$GRPATH/cogs_paper/DATA/";
my $DOTEST = 0; #IF set to true allows us to test the script by running only a few jobs
my $LD_REGION_DIR = "$BASE_DIR/support/0.1cM_shuffled_regions/";
my $MANIFEST_FILE = "$BASE_DIR/support/gwas_manifest.csv";
my $GWAS_DIR = "$BASE_DIR/gwas/";
my $OUTDIR = "$BASE_DIR/out/ppi/";
my $FILEPATTERN='^0\.1';
my $RSCRIPT="/home/oliver/bin/Rscript --vanilla $CHIGP_PATH/R/computePPi.R";
my $THOU_GENOMES_PREFIX="$BASE_DIR/1kgenome/VCF/EUR/by.chr.phase3/ALL.";
my $THOU_GENOMES_SUFFIX=".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR.vcf.gz";
my $TABIX_BIN="/home/oliver/bin/tabix_0.2.5";
## prior
my $PI_I=1e-4;
my $LOG_DIR = "$GRPATH/cogs_paper/log/cogs_ppi/";

my $step = Macd::Step->new(
	logdir=> $LOG_DIR,
	driver=>$DRIVER,
	env_settings=>"R_LIBS=$r_lib_dir,GRPATH=$GRPATH"
);







my @GWAS;
## open and parse the manifest file
open(GWAS,$MANIFEST_FILE) || die "Cannot open $MANIFEST_FILE\n";
my @header=();
my @TRAIT;
while(<GWAS>){
	chomp;
	if(!@header){
		@header=split(",",$_);
		next;
	}else{
		my %tmp;
		my @vals=split(",",$_);
		for(my $x=0;$x<@header;$x++){
			$tmp{$header[$x]}=$vals[$x];
		}
		$tmp{n_samples}=$tmp{cases} + $tmp{controls};
		if($tmp{type} eq 'CC'){
			$tmp{prop_cases}=sprintf("%.3f",($tmp{cases}/$tmp{n_samples}));
		}else{
			$tmp{prop_cases}=1;
		}
		push @TRAIT,\%tmp;
	}
}

#print Dumper(\@TRAIT);


foreach my $t(@TRAIT){
	my %th = %{$t};
	my $TEST=$DOTEST;
	my $addcount;
	my $TABIX_GWAS_FILE="${GWAS_DIR}$th{filename}";
	$TABIX_GWAS_FILE.='.gz' if $TABIX_GWAS_FILE!~/gz$/;
	if(! -e $TABIX_GWAS_FILE){
		print "ERROR: Cannot find $TABIX_GWAS_FILE skipping\n";
		next;
	}
	my $ODIR = "${OUTDIR}$th{label}";;
	if(! -e $ODIR){
		#die "Making $ODIR\n";
		make_path($ODIR) || die "Cannot make dir $ODIR\n";
	}
	find(sub {
		if(/$FILEPATTERN/){
			return if -e "$ODIR/$_.ppi";
			return if $TEST > 1;
			my @param = "$RSCRIPT";
      push @param, "--region_file=$File::Find::name";
      push @param, "--out_dir=$ODIR/";
      push @param, "--gwas_tbx=$TABIX_GWAS_FILE";
      push @param, "--gwas_type=$th{type}";
      push @param, "--n_samples=$th{n_samples}";
      push @param, "--prop_cases=$th{prop_cases}";
      push @param, "--kg_compress_dir=$THOU_GENOMES_PREFIX";
      push @param, "--kg_compress_suffix=$THOU_GENOMES_SUFFIX";
      push @param, "--tabix_bin=$TABIX_BIN";
      push @param, "--pi_i=$PI_I";
      my $cmd = join(" ",@param);
      #die "$cmd\n";
			my $job =  Macd::Step::Job->new(command=>$cmd);
			$step->add_job($job);
			$addcount++;
			$TEST++ if $TEST;
		}
	}, ($LD_REGION_DIR));
	##last;
	print "Added $addcount jobs for $th{label}\n";
}


if($step->execute()){
	print "Step submitted successfully\n";
	## we can hold up prog execution as follows
	## in case we have a step below that requires the output
	$step->wait_on_complete();
	print "Step completed successfully\n";
}else{
	print "Step did not complete successfully\n";
}

