<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="configurators-common.xslt" />
	
	<xsl:template name="configuration-content">
		<title>Model builder - ChemProp</title>
		<h1>ChemProp</h1>
		
		<table class="configuration">
			<tr>
			<td>Epochs (max: 1000): <input type="text" class="small" name="nepochs" value="200"/> </td>
			</tr><tr>
			<td>Batch size [8,128]: <input type="text" class="small" name="batch" value="32"/> </td>
			</tr><tr>
			<td>Depth size: <input type="text" class="small" name="depth" value="5"/> - number of message passing steps</td>
			</tr><tr>
			<td>Hidden neurones size: <input type="text" class="small" name="hidden" value="200"/></td>
			</tr>
		</table>
		<br/>

		<script language="javascript">
			<xsl:if test="//method = 'ChemProp'">
				setValue("nepochs", '<xsl:value-of select="//attachment/configuration/modelConfiguration/nepochs"/>');
				setValue("batch", '<xsl:value-of select="//attachment/configuration/modelConfiguration/batch"/>');
				setValue("depth", '<xsl:value-of select="//attachment/configuration/modelConfiguration/depth"/>');
				setValue("hidden", '<xsl:value-of select="//attachment/configuration/modelConfiguration/hidden"/>');
			</xsl:if>
		</script>
	</xsl:template>
</xsl:stylesheet>
