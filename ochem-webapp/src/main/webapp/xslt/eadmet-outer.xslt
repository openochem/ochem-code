<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"  version="2.0">
	<xsl:output method="xml" omit-xml-declaration="no" indent="yes"/>
	
	<xsl:template match="model">
		<html>
			<head>
				<meta name="Online Chemical Modeling Environment by QSPR Team" http-equiv="content-type" content="text/html; charset=UTF-8" />
				<style type="text/css">
					H1 {color: #093A7B;
						    font-family: "Century Gothic","Trebuchet MS","Lucida Grande","Lucida Sans Unicode",sans-serif;
						    font-size: 36px;
						    font-weight: normal;
						    line-height: 36px;
						    margin: 0;
						   }
						   
					H2 {
						color: #456fa7;
					    font-family: "Century Gothic","Trebuchet MS","Lucida Grande","Lucida Sans Unicode",sans-serif;
					    font-size: 22px;
					    font-weight: normal;
					    line-height: 22px;
					    margin-top: 30px;
					}
					
					P {font-family: Arial;}
					
					HTML, BODY {
						width: 1024px;
						margin: 20px auto;
					}
				</style>
			</head>
			<body>
				<xsl:apply-templates select="content" />
			</body>
		</html>
		
	</xsl:template>
	
	<xsl:template match="@*|node()" priority="-1">
			<xsl:copy>
				<xsl:apply-templates select="@*|node()" />
			</xsl:copy>
		</xsl:template>
</xsl:stylesheet>