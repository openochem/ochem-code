<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="configurators-common.xslt" />
	
	<xsl:template name="configuration-content">
		<title>Model builder - MLRA</title>
		<h1>Configure MLRA method</h1>
		<table class="configuration">
			<tr>
				<td>ALPHA</td>
				<td><input type="text" class="small" name="alpha" value="0.05"/> </td>
				<td>Maximal variables in equation</td>
				<td><input type="text" class="small" name="nvariables" value="100"/> </td>
			</tr>
		</table>
		<br/>
		<input type="checkbox" name="leverage"/> Calculate leverage for AD assessment<br></br>
		<input type="checkbox" name="limitrange"/> Limit predicted values to the training set range
			
		<script language="javascript">
			<xsl:if test="//method = 'MLRA'">
				setValue("alpha", '<xsl:value-of select="//attachment/configuration/modelConfiguration/alpha"/>');
				setValue("nvariables", '<xsl:value-of select="//attachment/configuration/modelConfiguration/nvariables"/>');
			</xsl:if>
		</script>
	</xsl:template>
</xsl:stylesheet>