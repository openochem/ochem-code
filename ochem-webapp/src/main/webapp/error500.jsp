<%@page import="java.io.*,java.net.URLDecoder,qspr.*,org.hibernate.*,qspr.entities.*,qspr.util.*"%>
<%
	response.setContentType("text/html");
	
	String code = null, message = null, type = null;
	Object codeObj, messageObj, typeObj;
	
	// Retrieve the three possible error attributes, some may be null
	Exception exception = (Exception)request.getAttribute("javax.servlet.error.exception");
	String trace;
	if (exception != null)
	{
		StringWriter sw = new StringWriter();
	    PrintWriter pw = new PrintWriter(sw);
	    exception.printStackTrace(pw);
	    pw.flush();
	    trace = sw.toString();
	} else
		trace = "Tace unknown";
	
	String uri = (String)request.getAttribute("javax.servlet.error.request_uri");
	if (uri == null)
		uri = request.getRequestURI();
%>
<html lang="en-US"><head>
    <meta http-equiv="Content-Type" content="text/html; charset=UTF-8">
    <title>OCHEM - Online Chemical Database with Modeling Environment</title>
    <BASE href="<%=OCHEMConfiguration.rootHost+OCHEMConfiguration.rootDir %>">
    <link rel="stylesheet" type="text/css" href="css/reset-fonts-grids.css">
    <link rel="stylesheet" type="text/css" href="css/template.css">
    <link rel="stylesheet" type="text/css" href="css/content.css">
    <link rel="stylesheet" type="text/css" href="css/details.css">
</head><body>
<div id="wrapper">
<div id="doc">
    <!-- start #header -->
    <div id="main">
    <div id="header">
    </div>
    <!-- end #header -->
    <div id="main-feature">
        <img src="img/molecule.png"/>
        <h2>Error 500</h2>
        <p>Error while accessing page <b><%= uri %></b></p>
    </div>
      
    <div id="main-content">
        <div id="sub-features">
            <div class="sub-feature-wide">
            	<div style="font-size: 130%">
            	Dear user,<br/><br/>
            	we are experiencing a problem in serving your request. <br/><br/>We would greatly appreciate it, if you could provide us with a brief description of the problem 
            	(a screenshot or copy of the failure, what you have tried to do, which files or sets, you have used, etc - any hint is helpful!).<br/>
            	<br/>
            	You can <a href="http://ochem.eu/support/open.php" target="_blank">submit your problem description at the support center.</a>
            	<br/><br/>
            	We appreciate your feedback!
            	</div>
                <h3>Stacktrace of the problem:</h3>
                <p>
                <code>
                <%= trace %>
                </code>
                </p>
            </div>
            <div class="clear"></div>   
        </div>
    </div>
    <div id="placeholder"></div>
    </div>
    <div id="footer"></div>
</div><!-- end #doc -->
</div><!-- end #wrapper -->
</body></html>