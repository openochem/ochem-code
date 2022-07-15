<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="configurators-common.xslt" />
	
	<xsl:template name="configuration-content">
		<title>Model builder - Random Forest</title>
		<h1>Configure Random Forest method</h1>
		<table class="configuration">
			<tr>
				<td>Number of trees:</td>
				<td><input type="text" class="small" name="trees" value="128"/></td>
			</tr>
			<tr>
				<td>Number of features to consider:<br/><small>(0 = log(number of inputs) + 1)</small></td>
				<td><input type="text" class="small" name="features" value="0"/></td>
			</tr>
			<tr>
				<td>Maximum tree depth:<br/><small>(0 = infinite)</small></td>
				<td><input type="text" class="small" name="depth" value="0"/></td>
			</tr>			
		</table>
		
		<script language="javascript">
			<xsl:if test="//attachment/configuration/modelConfiguration/numTrees">
				setValue("trees", '<xsl:value-of select="//attachment/configuration/modelConfiguration/numTrees"/>');
				setValue("features", '<xsl:value-of select="//attachment/configuration/modelConfiguration/numFeatures"/>');
				setValue("depth", '<xsl:value-of select="//attachment/configuration/modelConfiguration/maxDepth"/>');
			</xsl:if>
		</script>
			
	</xsl:template>
</xsl:stylesheet>