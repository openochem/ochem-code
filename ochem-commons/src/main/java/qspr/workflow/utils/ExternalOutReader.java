/* Copyright (C) 2022 BIGCHEM GmbH <info@bigchem.de>
 *
 * Contact: info@bigchem.de
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License (AGPL)
 * as published by the Free Software Foundation; either version 3.0
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the Affero GNU General Public License for more details.
 *
 * You should have received a copy of the Affero GNU Lesser General Public License
 * along with this program; If not, see <https://www.gnu.org/licenses/>. 
 */

package qspr.workflow.utils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.InputStream;
import java.io.PrintWriter;

import qspr.metaserver.CalculationServer;

import com.eadmet.utils.exe.OutReader;

/**
 * Extends OutReader to report status and keep track of errors
 * */

public class ExternalOutReader extends OutReader {
	public static final String messageTemplate = "MESSAGE:";
	public String errorMessage = null; 
	private static final String errorTemplate = "ERROR:";
	private static final String completedTemplate = "COMPLETED:";
	protected CalculationServer parentServer = null;
	PrintWriter pwOutputFile = null;
	private String lastMessages[] = new String[10];

	InterruptTimerTask task;

	public ExternalOutReader(InputStream is, Thread parentThread, CalculationServer parentServer) 
	{
		super(is, parentThread);
		this.parentServer = parentServer; 
	}

	/**
	 * Controlling whether task is stopped and should be restarted
	 * @param task
	 */

	public void setExternalMonitor(InterruptTimerTask task){
		this.task = task;
	}

	public String getLastMessage() {
		String s = "";
		for(int i=lastMessages.length-1;i>=0;i--) {
			s +=lastMessages[i] == null?"":lastMessages[i]+"\n";
			lastMessages[i] = null;
		}
		return s;
	}

	public void onLineRead(String line) throws Exception
	{
		super.onLineRead(line);

		for(int i=lastMessages.length-1;i>0;i--)
			lastMessages[i]=lastMessages[i-1];
		lastMessages[0]=line;


		if (line.toUpperCase().indexOf(messageTemplate)!=-1){

			if(task != null)
				task.resetInterruption();

			int n=line.toUpperCase().indexOf(completedTemplate);
			if(n!=-1){
				double val=Double.valueOf(line.substring(n+completedTemplate.length()+1));
				if(parentServer != null) {
					parentServer.setPercentageCompleted(val);
					line=line.substring(0,n)+" to finish in "+parentServer.getTimeToComplete(); 
				}
			}
			if(parentServer != null)
				parentServer.setStatus(line.substring(line.toUpperCase().indexOf(messageTemplate)+messageTemplate.length()));
		}

		if(parentServer != null)
			parentServer.scanStatus(line);

		if (pwOutputFile != null)
			pwOutputFile.println(line);

		if (isErrorMessage(line))
			errorMessage =  errorMessage == null? line : errorMessage +" " + line;
	}

	public ExternalOutReader saveOutTo(String fileName) throws FileNotFoundException
	{
		pwOutputFile = new PrintWriter(new File(fileName));
		return this;
	}

	public void flush() {
		if(pwOutputFile != null)
			pwOutputFile.flush();
	}

	public void close()
	{
		super.close();
		if (pwOutputFile != null)
			pwOutputFile.close();
	}

	static public boolean isErrorMessage(String line){
		line = line.toUpperCase();
		return line.startsWith(errorTemplate)
				|| line.contains("CUBLASError".toUpperCase())
				|| line.contains("CUDARuntimeError".toUpperCase())
				|| line.contains("CURANDError".toUpperCase())
				|| line.contains("ValueError:".toUpperCase())
				;
	}
}
