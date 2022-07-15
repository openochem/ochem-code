<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />

	<xsl:template name="content">
		<style type="text/css">
			TABLE.torefactor TD, TABLE.torefactor TH {background-color: #FFD; padding: 10px; text-align: left; border-bottom: 1px solid white;}
			TABLE.torefactor TH {font-style: bold; border-bottom: 1px solid black;}
			TR.ready TD {background-color: #c0e793;}
		</style>
		<table width="100%">
			<tr><td class="itunes-up"><h1>Running operation</h1></td></tr>
			<tr><td class="itunes-right">
				<table class="torefactor">
					<tr>
						<th>Time started</th>
						<th>ID</th>
						<th>User</th>
						<th>Status</th>
					</tr>
					<xsl:for-each select="//operation">
						<tr>
							<td><xsl:value-of select="time-started"/></td>
							<td><xsl:value-of select="operationId"/></td>
							<td><xsl:value-of select="userLogin"/></td>
							<td><xsl:value-of select="status"/></td>
						</tr>
					</xsl:for-each>
				</table>			
			</td></tr>
		</table>
		
	</xsl:template>
</xsl:stylesheet>