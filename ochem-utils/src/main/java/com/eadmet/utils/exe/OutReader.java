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

package com.eadmet.utils.exe;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;

// Midnighter

// Since it is heavily used inside of servers,
// it is extended to report messages of CalculationServer
// in order to do it, the statusServer  should be set by process which starts OutReader

public class OutReader extends Thread
{
	InputStream inputStream = null;
	BufferedReader reader;
	Thread parentThread;
	String preffix = null;
	PrintWriter out = new PrintWriter(System.out, true);

	public void run()
	{
		try
		{
			String line;
			while ((line = reader.readLine()) != null && parentThread.isAlive())
			{
				try
				{
					onLineRead(line);
				} catch (Exception e)
				{
					e.printStackTrace();
				}
			}
		} catch (IOException e)
		{
			// Stream closed.. its ok.
		} catch (Exception e)
		{
			e.printStackTrace();
		} finally
		{
			close(inputStream);
		}
	}

	public void onLineRead(String line) throws Exception
	{
		if (out != null)
		{
			if (preffix != null)
				out.println(preffix + line);
			else
				out.println(line);
		}
	}

	public OutReader setPreffix(String preffix)
	{
		this.preffix = preffix;
		return this;
	}

	public OutReader(InputStream is, Thread parentThread)
	{
		inputStream = is;
		this.reader = new BufferedReader(new InputStreamReader(is));
		this.parentThread = parentThread;
	}

	public OutReader setOutputWriter(PrintWriter out)
	{
		this.out = out;
		return this;
	}

	private void close(Closeable c)
	{
		if (c != null)
			try
		{
				c.close();
		} catch (IOException e)
		{
			// ignored
		}
	}

	public void close() {
		try {
			if(reader!=null)reader.close();
			if(inputStream !=null)inputStream.close();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}		
	}
}
