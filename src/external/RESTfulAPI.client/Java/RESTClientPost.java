import java.io.*;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URL;

public class RESTClientPost {
    public static void main(String[] args)
	{
		if(args.length<2)
		{
			System.out.println("\n$ java RESTClientPost [Inputfile] [Trigger] Submit:[E-mail](optional)\n$ java RESTClientPost [Inputfile] GNormPlus [Taxonomy ID]\n		e.g., java RESTClientPost input.PubTator tmChem Submit:[PubTator username](optional)\n		e.g., java RESTClientPost input.PubTator GNormPlus 10090\n\nParameters:\n\n	[Inputfile]:The file you would like to process.\n	[Trigger]:tmChem|DNorm|tmVar|GNormPlus\n	[Taxonomy ID]: NCBI Taxonomy identifier (e.g., 10090 for mouse). The species you would like to focus on. Only avaliable for GNormPlus.\n\n");
		}
		else
		{
			String Inputfile=args[0];
			String Trigger=args[1];
			String Taxonomy="";
			if(args.length > 2)
			{
				Taxonomy=args[2];
			}
			
			try {
				//Submit
				URL url_Submit;
				if(Taxonomy != "")
				{
					url_Submit = new URL("https://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/RESTful/tmTool.cgi/" + Trigger + "/" + Taxonomy + "/");
				}
				else
				{
					url_Submit = new URL("https://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/RESTful/tmTool.cgi/" + Trigger + "/Submit/");
				}
				HttpURLConnection conn_Submit = (HttpURLConnection) url_Submit.openConnection();
				conn_Submit.setDoOutput(true);
				conn_Submit.setRequestMethod("POST");
				
				BufferedReader fr= new BufferedReader(new FileReader(Inputfile));
				String input="";
				String str = null;
				while((str = fr.readLine()) != null)
				{
					input=input+str+"\n";
				}
				OutputStream os = conn_Submit.getOutputStream();
				os.write(input.getBytes());
				os.flush();
				BufferedReader br_Sumbit = new BufferedReader(new InputStreamReader(conn_Submit.getInputStream()));
				String SessionNumber="";
				String line="";
				while((line = br_Sumbit.readLine()) != null)
				{
					SessionNumber = SessionNumber + line;
				}
				conn_Submit.disconnect();

				String sub="";
				String email="";
				if(args.length > 2)
				{
					if(Taxonomy.length()>7)
					{
						sub=Taxonomy.substring(0,7);
						email=Taxonomy.substring(8);
					}
				}
				if(sub.equals("Submit:"))
				{
					System.out.println("Thanks for your submission (Session number: " + SessionNumber + ").\nThe result will be sent to your E-mail: " + email + ".\n");
				}
				else
				{
					try {
						System.out.println("Thanks for your submission. The session number is : "+SessionNumber+"\n");
						
						System.out.println("The request is received and processing....\n");
						
						//Receive
						URL url_Receive = new URL("https://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/RESTful/tmTool.cgi/" + SessionNumber + "/Receive/");
						HttpURLConnection conn_Receive = (HttpURLConnection) url_Receive.openConnection();
						BufferedReader br_Receive;
						String outputSTR="";
						line="";
						int code=404;
						while(code == 404 || code == 501)
						{
							try {Thread.sleep(5000);} catch(InterruptedException e) {}
							conn_Receive = (HttpURLConnection) url_Receive.openConnection();
							conn_Receive.setDoOutput(true);
							conn_Receive.setRequestMethod("GET");
							code = conn_Receive.getResponseCode();
						} 
						if(code == 200)
						{
							br_Receive = new BufferedReader(new InputStreamReader(conn_Receive.getInputStream()));
							while((line = br_Receive.readLine()) != null)
							{
								System.out.println(line);
							}
						}
						conn_Receive.disconnect();
					}
					catch(NullPointerException e){}
					
					
				}
			}
			catch (MalformedURLException e) 
			{
				e.printStackTrace();
			}
			catch (IOException e) 
			{
				e.printStackTrace();
			}
		}
    }
}