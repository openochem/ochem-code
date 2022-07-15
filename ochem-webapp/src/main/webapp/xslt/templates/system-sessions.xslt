<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />
	<xsl:template name="content">
		<style type="text/css">
			tr.header td {font-size: 80%; text-align: center; padding: 5px 10px 5px 10px;}
			.session
			{
				border-collapse: separate;
				border-spacing: 1px;
			}
			
			.session td
			{
				padding: 5px 5px 5px 5px;			
				background-color: #CCCCFF;		
			}
			
			.session th
			{
				padding: 2px 2px 2px 2px;			
				background-color: #AAAAEE;	
				vertical-align: top;	
			}
		</style>
		<table width="100%">
	  		<tr>
	  			<td class="itunes-up">
	  				<h1>System sessions</h1>
	  				Sessions with activity during last 24 hours
	  			</td>
	  		</tr>
			<tr>
			<td style="padding: 10px 10px 10px 10px; vertical-align:top;">
				<table class="session">
				<tr><th><b>Login</b></th><th><b>IP-address</b></th><th><b>Activity time</b></th></tr>
				<xsl:for-each select="//model/others/session">
					<tr>
					<td><xsl:value-of select="user/@login"/></td>
					<td><xsl:value-of select="ip-address"/></td>
					<td><xsl:value-of select="activity-time"/></td>
					</tr>
				</xsl:for-each>
				</table>
			</td>
			</tr>
		</table>
	</xsl:template>
	
</xsl:stylesheet>
