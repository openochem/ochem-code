<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="configurators-common.xslt" />
	
	<xsl:template name="configuration-content">
		<title>Model builder - KNN</title>
		<h1>Configure KNN method</h1>
		<table class="configuration">
			<tr>
				<td>Distance:</td>
				<td>
				<select name="distance">
					<option value="0">Euclidean distance</option>
					<option value="1">Pearson's correlation</option>
				</select></td>
			</tr>
			<tr>
				<td>KNN neighbors<br/><small>(leave 0 for automatic determination)</small></td>
				<td><input type="text" class="small" name="knn" value="0"/> </td>
			</tr>
			<tr>
				<td>Max KNN neighbors<br/><small>(in case of automatic determination)</small></td>
				<td><input type="text" class="small" name="max-knn" value="100"/> </td>
			</tr>
		</table>	
		
		<script language="javascript">
			<xsl:if test="//method = 'KNN'">
				setValue("max-knn", '<xsl:value-of select="//attachment/configuration/modelConfiguration/maxKnn"/>');
				setValue("knn", '<xsl:value-of select="//attachment/configuration/modelConfiguration/knn"/>');
				setValue("distance", '<xsl:value-of select="//attachment/configuration/modelConfiguration/distance"/>');
			</xsl:if>
		</script>
	</xsl:template>
</xsl:stylesheet>