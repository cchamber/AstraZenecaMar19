# http://search.cpan.org/~mcrawfor/REST-Client/lib/REST/Client.pm
# Example install using cpanm:
#   sudo cpanm -i REST::Client
require 'lib/REST/Client.pm';
use LWP;
my $browser = LWP::UserAgent->new;

sub REST_Post
{
	my $Inputfile =	@_[0];
	my $trigger = @_[1];
	my $tax_id = @_[2];
	
	# Set the request parameters
	my $host = 'https://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/RESTful/tmTool.cgi';
	my $client = REST::Client->new(host => $host);

	my $STR="";
	open input,"<$Inputfile";
	while(<input>)
	{
		$STR=$STR.$_;
	}
	close input;
	
	my $SessionNumber = "";
	if($tax_id ne "")
	{
		$SessionNumber = $client->POST("/".$trigger."/".$tax_id."/",$STR)->responseContent;
	}
	else
	{
		$SessionNumber = $client->POST("/".$trigger."/Submit/",$STR)->responseContent;
	}
	
	if($tax_id=~/^Submit:(.+)$/)
	{
		print "Thanks for your submission (Session number: ".$SessionNumber.").\nThe result will be sent to your E-mail: ".$1.".\n";
		exit;
	}
	else
	{
		print "Thanks for your submission. The session number is : ".$SessionNumber."\n";
		print "\nThe request is received and processing....\n\n";
		my $receive_result = $client->GET("/".$SessionNumber."/Receive/");
		#print $receive_result->responseCode()."\n";
		while($receive_result->responseCode()=~/^(501|404|500)$/) # Result is not ready
		{
			sleep(5); # Don't decrease the sleep time. Too frequent accessing in NCBI may by blocked.
			$receive_result = $client->GET("/".$SessionNumber."/Receive/");
			#print $receive_result->responseCode()."\n";
		}
		my $result = $receive_result->responseContent();
		open output,">$Inputfile.out";
		print output $result;
		close output;
		print "\nThe results are saved in $Inputfile.out\n\n";
	}

}

sub main
{
	my $InputfileQuery=$ARGV[0];
	my $trigger=$ARGV[1];
	my $tax_id=$ARGV[2];
	if($InputfileQuery ne "" && $trigger ne "")
	{
		my ($sec1,$min1,$hour1,$day1,$mon,$year)=localtime(time);
		
		&REST_Post($InputfileQuery,$trigger,$tax_id);
		
		my ($sec2,$min2,$hour2,$day2,$mon,$year)=localtime(time);
		$hour1=$hour1+$day1*24;
		$hour2=$hour2+$day2*24;
		my $timecost=((($hour2-$hour1)*60)+($min2-$min1)*60)+($sec2-$sec1);
		print "Running RESTful in $timecost seconds. \n";
	}
	else
	{
		print "
perl RESTful.client.POST.pl [Input file] [Trigger:tmVar|tmChem|DNorm|GNormPlus]] Submit:[PubTator username](optional)
	Eg., perl RESTful.client.pl input.txt tmChem
	Eg., perl RESTful.client.pl input.txt tmChem Submit:[E-mail]
	
perl RESTful.client.POST.pl [Input file] [GNormPlus] [Taxonomy ID]
	Eg., perl RESTful.client.pl input.txt GNormPlus 9606
";	
	}
}
&main();