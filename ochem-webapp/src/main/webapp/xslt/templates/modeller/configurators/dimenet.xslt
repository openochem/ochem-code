<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="configurators-common.xslt" />
	<xsl:include href="../../model/inc/structure-optimisation-commons.xslt" />
	
	<xsl:template name="configuration-content">
		<title>Model builder - DIMENET</title>
		<h1>Configure Directional Message Passing Neural Network</h1>
		
		<table class="configuration">
			<tr>
				<td>Epochs<br/></td>
				<td><input type="text" class="small" name="nbepochs" value="1000"/> </td>
			</tr>
			<tr>
				<td>Batch size<br/></td>
				<td><input type="text" class="small" name="batch" value="1024"/> </td>
			</tr>
 			<tr>
				<td>Early stopping<br/></td>
				<td><input type="checkbox" name="early" checked="checked"/></td>
			</tr>			
		</table>
		<br/>
		<table class="configuration">
 			<h1>By default RDKIT is used to generate 3D (No optimisation). You can also select another method</h1>
 			<tr>
				<xsl:call-template name="structure-optimisation"/>
			</tr>			
		</table>	
		<br/>
		<br/>
		
		<script language="javascript">
			<xsl:if test="//method = 'DIMENET'">
				setValue("batch", '<xsl:value-of select="//attachment/configuration/modelConfiguration/batch"/>');
				setValue("nbepochs", '<xsl:value-of select="//attachment/configuration/modelConfiguration/nbepochs"/>');
				setValue("early", '<xsl:value-of select="//attachment/configuration/modelConfiguration/early"/>');
			</xsl:if>
		</script>
	</xsl:template>
</xsl:stylesheet>

	