package apply;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.rmi.RemoteException;

import javax.xml.rpc.ServiceException;

import ochem.eadmet.wsapi.ModelResponse;
import ochem.eadmet.wsapi.ModelServiceLocator;
import ochem.eadmet.wsapi.ModelServicePortType;
import ochem.eadmet.wsapi.Prediction;
import ochem.eadmet.wsapi.PropertyPrediction;

import org.apache.axis.AxisFault;

public class TestOchemWS 
{
	public static void main(String[] args) throws ServiceException{
		long modelID = 0;
		long taskID = 0;
		String sdf[] = null;
		
		if (args.length < 2) {
			printError("No sufficient input parameters");
			System.exit(0);
		}
		
		if (args[0].equals("fetch")) 
		{
			try {
				taskID = Long.parseLong(args[1]);
			}
			catch (Exception e) {
				printError("Task-ID is not valid (only numbers)");
				System.exit(0);
			}
			try {
				ModelServicePortType stub  = getService();
				ModelResponse modelResponse = stub.fetchModel(null,taskID);
				printResult(modelResponse, true);
			} catch (AxisFault e) {
				System.out.println(e);
			} catch (RemoteException e) {
				System.out.println(e);
			}
		}
		else if (args[0].equals("post") || args[0].equals("postm")) 
		{
			try 
			{
				modelID = Long.parseLong(args[1]);
			}
			catch (Exception e) {
				printError("Model-ID is not valid (only numbers)");
				System.exit(0);
			}
			if (args.length < 3) {
				printError("No sufficient input parameters");
				System.exit(0);
			}
			sdf = new String[args.length -2];
			for (int i = 2; i < args.length; i++) {
				String checked = checkSDF(args[i]);
				if (null != checked) {
					printError(checked);
					System.exit(0);
				}
				try {
					sdf[i -2] = readSDF(args[i]);
				} catch (Exception e) {
					printError("Provided file can not be found or read: " + args[2]);
					System.exit(0);
				}
			}
			try 
			{
				ModelServicePortType stub  = getService();
				ModelResponse modelResponse = null;
				if (args[0].equals("post")) {
					modelResponse = stub.postModelWithSession(null, modelID, sdf);
				}
				else if (args[0].equals("postm")) {
					modelResponse = stub.applyModelSingleSDF(null,modelID, sdf[0]);
				}
				printResult(modelResponse, false);
			} catch (AxisFault e) {
				System.out.println(e);
			} catch (RemoteException e) {
				System.out.println(e);
			}
		}
		else {
			printError("No action specified");
			System.exit(0);
		}
	}
	
	public static String readSDF(String fileLocation) throws Exception {
		String content = "";
		BufferedReader br;
		br = new BufferedReader(new FileReader(fileLocation));
		String line = null; 
		while ((line = br.readLine()) != null) {
			content = content + line + "\n";
		}
		br.close();
		return content;
	}
	
	public static void printError(String error) {
		System.out.println("ERROR: " + error);
		System.out.println("Usage:");
		System.out.println("Posting a request:                   java -jar Thesaurus.jar post [public_model_id] [sdf_file1] [sdf_file2] ...");
		System.out.println("Posting a request (multi sdf-file):  java -jar Thesaurus.jar postm [public_model_id] [sdf_file1]");
		System.out.println("Fetching a task:                     java -jar Thesaurus.jar fetch [public_model_id]");
	}
	
	private static ModelServicePortType getService() throws ServiceException
	{
		ModelServiceLocator locator = new ModelServiceLocator();
		ModelServicePortType service = locator.getModelServiceHttpSoap11Endpoint();
		return service;
	}
	
	public static String checkSDF(String filename) {
		File file = new File(filename);
		String error = null;
		if (!file.exists() || !file.isFile()) {
			error = "Provided file can not be found or read: " + filename;
		}
		else if (file.length() > 10000000) {
			error = "Input file is too big: " + filename + "\n       (input files habe to smaller than 10MB)";
		}
		return error;
	}
	
	public static void printResult(ModelResponse modelResponse, boolean fetch) 
	{
		if (fetch)
			System.out.println("Fetching task ...");
		else
			System.out.println("Posting request ...");
		System.out.println();
		System.out.println("STATUS:                " + modelResponse.getStatus());
		
		if (fetch) 
		{
			if (modelResponse.getStatus().equals("success")) 
			{
				System.out.println("MODEL-URL:             " + modelResponse.getModelDescriptionUrl());
				if (null != modelResponse.getPredictions()) 
				{
					System.out.println();
					System.out.println("MODEL RESPONSE:");
					Prediction preds[] = modelResponse.getPredictions();
					System.out.println("| Property                     | Predicted value              | Unit                         | Accuracy                     |");
					for (int i = 0; i < preds.length; i++) 
					{
						PropertyPrediction proppred[] = preds[i].getPredictions();
						if (null == preds[i].getError() && null != proppred) 
						{
							for (int j = 0; j < proppred.length; j++) 
							{
								String property = addWhiteSpace(proppred[j].getProperty());
								String value = addWhiteSpace(proppred[j].getValue() + "");
								String unit = addWhiteSpace(proppred[j].getUnit());
								String accuracy = addWhiteSpace(proppred[j].getAccuracy() + "");
								System.out.println("|" + property + "|" + value + "|" + unit + "|" + accuracy + "|");
							}
						}
						else if (null != preds[i].getError()) 
						{
							System.out.println("\t" + preds[i].getError());
						}
						else 
						{
							System.out.println("*** UNEXPECTED ERROR ***");
						}
					}
				}
			}
		}
		else
			System.out.println("TASK-ID:               " + modelResponse.getTaskId());
	}
	
	public static String addWhiteSpace(String str) {
		str = str + " ";
		while (str.length() < 30)
			str = " " + str;
		return str;
	}
}
