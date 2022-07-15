<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:xs="http://www.w3.org/2001/XMLSchema" version="2.0">
	<xsl:include href="../helper.xslt" />
	<xsl:include href="inc/data-report.xslt" />
	<xsl:template name="content">
		<title>Editor's corner</title>
		
		<table height="100%" width="100%">
			<tr>
				<td class="itunes-up" colspan="2">
					<img src="img/icons/editor.png"/>
					<h1>Editor's corner</h1>
				</td>
			</tr>
			<tr>
				<td class="itunes-right">
					<h1>OCHEM data: analytical report</h1>
					<xsl:call-template name="data-report"/>
				</td>
			</tr>
		</table>
	</xsl:template>
	
	
</xsl:stylesheet>