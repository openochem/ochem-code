<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="configurators-common.xslt" />
	
	<xsl:template name="configuration-content">
		<title>MACAU parameters</title>
		<h1>Configure MACAU method</h1>
		<table class="configuration">
			<tr>
				<td>Test ratio (fraction of samples for internal validation)</td>
				<td><input type="text" class="small" name="test_percent" value="0.2"/> </td>
			</tr>
			<tr>
				<td>Number of latent variables (int)</td>
				<td><input type="text" class="small" name="num_latent" value="32"/> </td>
			</tr>
			<tr>
				<td>Burnin</td>
				<td><input type="text" class="small" name="burnin" value="300"/> </td>
			</tr>
			<tr>
				<td>Number of generated models (large models will not be saved)</td>
				<td><input type="text" class="small" name="samples" value="512"/> </td>
			</tr>			
			<tr>
				<td>Experimental error of data (variance)</td>
				<td><input type="text" class="small" name="accuracy" value="0.5"/> </td>
			</tr>
		</table>
		<table class="configuration">
			<tr>
				OR  <input type="checkbox" name="adaptive"/><label> Determine data precision automatically</label>
			</tr>
		</table>
		<br/>
			
		<script language="javascript">
			<xsl:if test="//method = 'MACAU'">
				setCheckbox("adaptive", '<xsl:value-of select="//attachment/configuration/descriptors/types/configuration/adaptive"/>');
				setValue("test_percent", '<xsl:value-of select="//attachment/configuration/modelConfiguration/test_percent"/>');
				setValue("num_latent", '<xsl:value-of select="//attachment/configuration/modelConfiguration/num_latent"/>');
				setValue("accuracy", '<xsl:value-of select="//attachment/configuration/modelConfiguration/accuracy"/>');
				setValue("burnin", '<xsl:value-of select="//attachment/configuration/modelConfiguration/burnin"/>');
				setValue("samples", '<xsl:value-of select="//attachment/configuration/modelConfiguration/samples"/>');
			</xsl:if>
		</script>
	</xsl:template>
</xsl:stylesheet>
