<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="configurators-common.xslt" />
	
	<xsl:template name="configuration-content">
		<title>Model builder - XGBOOST</title>
		<h1>Configure XGBOOST method</h1>
		<table class="configuration">
			<tr>
				<td>Max - Maximum depth of a tree, increase this value will make the model more complex<br/></td>
				<td><input type="text" class="small" name="depth" value="8"/> </td>
			</tr>
			<tr>
				<td>ETA - Step size shrinkage<br/></td>
				<td><input type="text" class="small" name="eta" value="0.2"/> </td>
			</tr>
			<tr>
				<td>Lambda - L2 regularization term on weights, increase this value will make model more conservative. <br/></td>
				<td><input type="text" class="small" name="lambda" value="1"/> </td>
			</tr>			
			<tr>
				<td>Rounds - The number of rounds for boosting<br/></td>
				<td><input type="text" class="small" name="rounds" value="1024"/> </td>
			</tr>
			<tr>
				<td>Objective:</td>
				<td>
				<select name="objective">
					<option value="reg:linear">Linear regression</option>
					<option value="reg:logistic">Logistic regression</option>
					<option value="binary:logistic">logistic regression for binary classification</option>
				</select></td>
			</tr>
		</table>	
		
		<script language="javascript">
			<xsl:if test="//method = 'XGBOOST'">
				setValue("depth", '<xsl:value-of select="//attachment/configuration/modelConfiguration/lambda"/>');
				setValue("depth", '<xsl:value-of select="//attachment/configuration/modelConfiguration/depth"/>');
				setValue("eta", '<xsl:value-of select="//attachment/configuration/modelConfiguration/eta"/>');
				setValue("rounds", '<xsl:value-of select="//attachment/configuration/modelConfiguration/rounds"/>');
				setValue("objective", '<xsl:value-of select="//attachment/configuration/modelConfiguration/objective"/>');
			</xsl:if>
		</script>
	</xsl:template>
</xsl:stylesheet>