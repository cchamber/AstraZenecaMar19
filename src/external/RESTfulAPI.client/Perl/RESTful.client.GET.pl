# http://search.cpan.org/~mcrawfor/REST-Client/lib/REST/Client.pm
# Example install using cpanm:
#   sudo cpanm -i REST::Client
require 'lib/REST/Client.pm';
use LWP;
my $browser = LWP::UserAgent->new;

sub REST_Get
{
	my $Bioconcept = $_[0];
	my $InputFile = $_[1];
	my $Format = $_[2];	
	open input,"<$InputFile";
	while(<input>)
	{
		my $pmid=$_;
		$pmid=~s/[\n\r]//g;
		my $host = "https://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/RESTful/tmTool.cgi";
		my $client = REST::Client->new(host => $host);
		my $Query = "/".$Bioconcept."/".$pmid."/".$Format."/";
		print $client->GET($Query)->responseContent()."\n";
	}
	close input;
}

sub main
{
	my $BioConcept=$ARGV[0];
	my $InputFile=$ARGV[1];
	my $Format=$ARGV[2];
	if($BioConcept eq "" || $InputFile eq "" || $Format eq "")
	{
		print "
perl RESTful.client.GET.pl [BioConcept] [InputFile] [Format]

	Bioconcept: We support five kinds of bioconcepts, i.e., Gene, Disease, Chemical, Species, Mutation. When 'BioConcept' is used, all five are included.
	Inputfile: a file with a pmid list
	Format: PubTator (tab-delimited text file), BioC (xml), and JSON

Eg., perl RESTful.client.GET.pl BioConcept ../examples/ex.pmid PubTator
";	
	}
	else
	{
		&REST_Get($BioConcept,$InputFile,$Format);
	}
}
&main();