<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="configurators-common.xslt" />
	
	<xsl:template name="configuration-content">
		<title>Model builder - Transformer</title>
		<h1>Transformer Neural Networks</h1>
		
		<table class="configuration">
			<tr>
			<td>Epochs (max: 1000): <input type="text" class="small" name="nepochs" value="100"/> </td>
			</tr><tr>
			<td>Batch size [8,128]: <input type="text" class="small" name="batch" value="64"/> </td>
			</tr><tr>
			<td>Depth size: <input type="text" class="small" name="depth" value="3"/> - number of message passing steps</td>
			</tr><tr>
			<td>Hidden neurones size: <input type="text" class="small" name="hidden" value="300"/></td>
			</tr>
		</table>
		<br/>

		<script language="javascript">
			<xsl:if test="//method = 'ChemProp'">
				setValue("nepochs", '<xsl:value-of select="//attachment/configuration/modelConfiguration/nepochs"/>');
				setValue("batch", '<xsl:value-of select="//attachment/configuration/modelConfiguration/batch"/>');
				setValue("depth", '<xsl:value-of select="//attachment/configuration/modelConfiguration/depth"/>');
				setValue("hidden", '<xsl:value-of select="//attachment/configuration/modelConfiguration/hidden"/>');
				setCheckbox("edges", '<xsl:value-of select="//attachment/configuration/modelConfiguration/edges"/>');
			</xsl:if>
		</script>
	</xsl:template>
</xsl:stylesheet>
