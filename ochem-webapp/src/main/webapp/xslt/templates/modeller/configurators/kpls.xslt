<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="configurators-common.xslt" />
	
	<xsl:template name="configuration-content">
		<title>Model builder - KPLS</title>
		<h1>Configure KPLS (Kernel Partial Least Squares) method</h1>
		<table class="configuration">
			<tr>			
			<td>Number of latent variables:<br/></td>
				<td>
				<select name="latentVariables">
					<option value="0">select automatically (longer calculations)</option>
					<option value="5" selected="1">5</option>
					<option value="12" >12</option>
				</select>
				</td>
			</tr>
			<tr>
				<td>Use alternative Lambda Tune<br/></td>
				<td><input type="checkbox" checked="checked" name="useLambda"/></td>
			</tr>
		</table>	
	<script language="javascript">
			<xsl:if test="//method = 'KPLS'">
				setValue("latentVariables", '<xsl:value-of select="//attachment/configuration/modelConfiguration/latentVariables"/>');
				setCheckbox("useLambda", '<xsl:value-of select="//attachment/configuration/modelConfiguration/useLambda"/>');
			</xsl:if>
			$(document).ready(function(){
				updateVisibility();	
			});
		</script>
	</xsl:template>
</xsl:stylesheet>