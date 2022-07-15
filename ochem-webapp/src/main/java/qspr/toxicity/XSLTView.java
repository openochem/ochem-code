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

package qspr.toxicity;

import java.io.File;
import java.io.StringReader;
import java.io.StringWriter;
import java.util.Calendar;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Pattern;

import javax.servlet.ServletContext;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import javax.xml.bind.Marshaller;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.stream.StreamResult;
import javax.xml.transform.stream.StreamSource;

import net.sf.saxon.trans.XPathException;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.json.JSONArray;
import org.json.JSONObject;
import org.json.XML;
import org.springframework.web.servlet.view.AbstractView;

import qspr.Globals;
import qspr.ThreadScope;
import qspr.VirtualHostConfiguration;
import qspr.frontend.WebModel;

import com.eadmet.exceptions.UserFriendlyException;

public class XSLTView extends AbstractView
{
	private static transient final Logger logger = LogManager.getLogger(XSLTView.class);

	@SuppressWarnings("rawtypes")
	@Override
	public void renderMergedOutputModel(Map _model, HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		String out = (req.getParameter("out") != null) ? req.getParameter("out") : "html";
		if ("json".equals(_model.get("view")))
			out = "json";
		transformAndRender(req, res, _model.get("object"), out);
	}


	public void transformAndRender(HttpServletRequest request, HttpServletResponse response, Object object, String out) throws Exception
	{
		ThreadScope.setStatus("Marshalling XML data...");


		StringWriter outWriter;
		String content;
		long time;
		//StringWriter xmlWriter2 = new StringWriter();
		//StringWriter xmlWriter3 = new StringWriter();

		try
		{
			System.setProperty("jaxb.debug", "true");
			Marshaller marshaller = Globals.jaxbContext.createMarshaller();

			if (Globals.debugJaxb)  // Uncomment to profile marshalling 
				marshaller.setListener(new Marshaller.Listener()
				{
					// Marshaller profiling by Midnighter
					public Map<Object, Long> time = new HashMap<Object, Long>();
					public String prefix = "";
					public void beforeMarshal(Object source)
					{
						time.put(source, Calendar.getInstance().getTimeInMillis());
						logger.info(prefix + source.getClass());
						prefix += "   ";
					}

					public void afterMarshal(Object source)
					{
						prefix = prefix.substring(3);
						logger.info(prefix + "/" + source.getClass() + " in " + (Calendar.getInstance().getTimeInMillis() - time.get(source)) + "ms.");
					}
				});

			//CharacterEscapeHandler escapeHandler = MidnightersSpecialEscapeHandler.theInstance;
			//marshaller.setProperty("com.sun.xml.bind.characterEscapeHandler", escapeHandler);
			marshaller.setProperty(Marshaller.JAXB_ENCODING, "utf-8");

			time = Calendar.getInstance().getTimeInMillis();
			outWriter = new StringWriter();

			marshaller.marshal(object, outWriter);
			content = outWriter.toString();
			time = Calendar.getInstance().getTimeInMillis() - time;
			if (time > 50)
				logger.info("Marshalled XML of " + content.length() + " bytes in " + time + "ms.");		

			if (out.equals("markup") || out.equals("html"))
			{
				WebModel webModel = (WebModel) object;
				if (request.getParameter("render-mode") != null)
					webModel.renderMode = request.getParameter("render-mode");
				time = Calendar.getInstance().getTimeInMillis();
//				TransformerFactory transFact = TransformerFactory.newInstance();
				TransformerFactory transFact = TransformerFactory.newInstance("net.sf.saxon.TransformerFactoryImpl",null);
//				System.out.println(transFact.getClass());
				ServletContext context = request.getSession().getServletContext();
				String realContextPath = context.getRealPath("/");

				Transformer contentTransformer = transFact.newTransformer(new StreamSource(new File(realContextPath+"/xslt/templates/" + webModel.templateName + ".xslt")));
				if (contentTransformer == null)
					throw new UserFriendlyException("Cannot find or compile template xslt/templates/" + webModel.templateName + ".xslt");
				contentTransformer.setOutputProperty("omit-xml-declaration", "no");
				contentTransformer.setOutputProperty("encoding", "utf-8");

				outWriter = new StringWriter();
				try
				{
					contentTransformer.transform(new StreamSource(new StringReader(content)), new StreamResult(outWriter));
				} catch (XPathException e) {
					if (!e.getMessage().contains("Unicode"))
						throw e;
					e.printStackTrace();
					// Replace weird unicode characters and try again
					Pattern p = Pattern.compile("[^\\u0009\\u000A\\u000D\u0020-\\uD7FF\\uE000-\\uFFFD\\u10000-\\u10FFF]+");
					//					Matcher m = p.matcher(content);
					//					while (m.find())
					//					{
					//						System.out.println(m.start()+"-"+m.end()+":"+m.group());
					//					}
					content = p.matcher(content).replaceAll("");
					outWriter = new StringWriter();
					contentTransformer.transform(new StreamSource(new StringReader(content)), new StreamResult(outWriter)); 
				}
				content = outWriter.toString();
				
				logger.info("SIze of content is:"+content.length());

				if (out.equals("html"))
				{
					if (webModel.outerTemplate == null)
					{
						String overridenTemplate = VirtualHostConfiguration.getOuterTemplate();
						if (overridenTemplate != null)
							webModel.outerTemplate = overridenTemplate;
						else
							webModel.outerTemplate = "global"; // The default theme
					}

					Transformer globalTransformer = transFact.newTransformer(new StreamSource(new File(realContextPath+"/xslt/" + webModel.outerTemplate + ".xslt")));
					Transformer markupTransformer = transFact.newTransformer(new StreamSource(new File(realContextPath+"/xslt/markup.xslt")));
					globalTransformer.setOutputProperty("encoding", "utf-8");
					globalTransformer.setOutputProperty("omit-xml-declaration", "no");
					markupTransformer.setOutputProperty("encoding", "utf-8");	
					markupTransformer.setOutputProperty("omit-xml-declaration", "yes");

					outWriter = new StringWriter();
					globalTransformer.transform(new StreamSource(new StringReader(content)), new StreamResult(outWriter));
					content = outWriter.toString();

					outWriter = new StringWriter();
					markupTransformer.transform(new StreamSource(new StringReader(content)), new StreamResult(outWriter));
					content = outWriter.toString();

					outWriter.flush();

					response.setHeader("Content-type", "text/html");
					setContentType("text/html");

				} else
				{
					response.setHeader("Content-type", "text/xml");
					setContentType("text/xml");	
				}
				time = Calendar.getInstance().getTimeInMillis() - time;
				if (time > 100)
					logger.info("XSLT time: "+time+"ms.");		
			} else
				if (out.equals("json"))
				{
					time = Calendar.getInstance().getTimeInMillis();
					JSONObject o = XML.toJSONObject(content);
					if (request.getParameter("list-only") != null)
					{
						try
						{
							JSONArray array = o.getJSONObject("model").getJSONObject("others").getJSONArray("labeledValue");
							content = array.toString();
						}
						catch (org.json.JSONException e)
						{
							try
							{
								content = "[" + o.getJSONObject("model").getJSONObject("others").getJSONObject("labeledValue").toString() + "]";
							}
							catch (org.json.JSONException e2)
							{
								content = "[]";
							}
						}
					}
					else
					{
						o = o.getJSONObject("model");
						content = o.toString(1);
					}
					time = Calendar.getInstance().getTimeInMillis() - time;
					if (time > 50)
						logger.info("XML2JSON time: "+time+"ms.");		

					response.setHeader("Content-type", "text/html");
					setContentType("text/html");

				} else
				{ //out = xml
					response.setHeader("Content-type", "text/xml");
					setContentType("text/xml");
				}
			response.setCharacterEncoding("utf-8");
			response.getWriter().write(content);		
		}
		catch (Exception e)
		{
			e.printStackTrace();
			throw e;
		}
	}
}

