<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="configurators-common.xslt" />
	
	<xsl:template name="configuration-content">
		<title>Least Squares Support Vector Machine builder - LSSVMG (GPU)</title>
		<h1>Configure LSSVMG method</h1>
		<table class="configuration">
			<tr>
				<td>Kernel</td>
				<td>
				<select name="kernel-lssvmg">
					<option value="rbf" selected="1">RBF kernel</option>
					<option value="linear">Linear kernel</option>
				</select>
                </td>
			</tr>
			<tr>
				<td>Internal CV folds</td>
				<td><input type="text" class="small" name="cv-lssvmg" value="5"/> </td>
			</tr>
			<tr>
				<td>  Use global optimisation of parameters</td>
				<td><input type="checkbox" name="use-global"/></td>
			</tr>
			<tr>
				<td>Additional Parameters<br/><small>(separated by comma)</small></td>
				<td><input type="text" name="additionalParam" value=""/></td>
			</tr>
		</table>
		<br/>
			
		<script language="javascript">
			<xsl:if test="//method = 'LSSVMG'">
				setValue("cv-lssvm", '<xsl:value-of select="//attachment/configuration/modelConfiguration/cv"/>');
				setValue("kernel-lssvm", '<xsl:value-of select="//attachment/configuration/modelConfiguration/kernel"/>');
				setCheckbox("use-global", '<xsl:value-of select="//attachment/configuration/modelConfiguration/useGlobal"/>');
				setValue("additionalParam", '<xsl:value-of select="//attachment/configuration/modelConfiguration/additionalParam"/>');
			</xsl:if>
		</script>
	</xsl:template>
</xsl:stylesheet>