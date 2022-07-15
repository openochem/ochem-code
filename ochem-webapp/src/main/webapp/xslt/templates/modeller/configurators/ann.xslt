<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="configurators-common.xslt" />
	
	<xsl:template name="configuration-content">
		<title>Model builder - ANN</title>
		<h1>Configure ANN method<a class="infolink" href="https://docs.ochem.eu/pages/viewpage.action?pageId=7012845" target="_blank"></a></h1>
		<table class="configuration">
			<tr>
				<td>Training method:</td>
				<td>
				<select name="training-method">
					<option value="0">Momentum</option>
					<option value="1" selected="1">SuperSAB</option>
					<option value="2">RPROP</option>
					<option value="3">QuickProp</option>
					<option value="4">Differential equations</option>
					<option value="5">QuickProp II</option>
					<option value="6">Levenberg-Marquardt</option>
				</select></td>
			</tr>
			<tr>
				<td>Number of neurons<br/><small>in hidden layer</small><br/></td>
				<td><input type="text" class="small" name="neurons" value="3"/> </td>
			</tr>
			<tr>
				<td>Learning iterations<br/><small>(learning iterations)</small></td>
				<td><input type="text" class="small" name="iterations" value="1000"/> </td>
			</tr>
			<tr>
				<td>Ensemble</td>
				<td><input type="text" class="small" name="ensemble" value="64"/></td>
			</tr>
			<tr>
				<td>Enable ASNN<br/></td>
				<td><input type="checkbox" checked="checked" name="asnn"/></td>
			</tr>			
			<tr>
				<td>Additional Parameters<br/><small>(separated by comma)</small></td>
				<td><input type="text" name="additionalParam" value=""/></td>
			</tr>
			<tr>
				<td>Experimental Parameters<br/><small>(conditions merging)</small></td>
				<td><input type="text" name="experimentalParam" value=""/></td>
			</tr>
		</table>	
		
		<script language="javascript">
			<xsl:if test="//method = 'ASNN'">
				setValue("training-method", '<xsl:value-of select="//attachment/configuration/modelConfiguration/training"/>');
				setValue("neurons", '<xsl:value-of select="//attachment/configuration/modelConfiguration/neurons"/>');
				setValue("iterations", '<xsl:value-of select="//attachment/configuration/modelConfiguration/iterations"/>');
				setValue("ensemble", '<xsl:value-of select="//attachment/configuration/modelConfiguration/ensemble"/>');
				setCheckbox("asnn", '<xsl:value-of select="//attachment/configuration/modelConfiguration/asnn"/>');
				setValue("additionalParam", '<xsl:value-of select="//attachment/configuration/modelConfiguration/additionalParam"/>');
				setValue("experimentalParam", '<xsl:value-of select="//attachment/configuration/modelConfiguration/experimentalParam"/>');
			</xsl:if>
		</script>
	</xsl:template>
</xsl:stylesheet>