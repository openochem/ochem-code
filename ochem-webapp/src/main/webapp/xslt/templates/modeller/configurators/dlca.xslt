<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="configurators-common.xslt" />
	
	<xsl:template name="configuration-content">
		<title>Model builder - Deep Learning Consensus Architecture</title>
		<h2>DLCA: https://github.com/ncats/ld50-multitask</h2>

		<h1>Method specific hyperparameters</h1>
		<div class="configuration">
		<table>
		<tr>
			<td>Epochs<br/></td>
				<td><input type="text" class="small" name="epochs" value="20"/> </td>
		</tr>
		<tr>
			<td>Batch size<br/></td>
				<td><input type="text" class="small" name="batch" value="32"/> </td>
		</tr>
		</table>
		</div>
		
		<script language="javascript">
			<xsl:if test="//method = 'DLCA'">
				setValue("epochs", '<xsl:value-of select="//attachment/configuration/modelConfiguration/epochs"/>');
				setValue("batch", '<xsl:value-of select="//attachment/configuration/modelConfiguration/batch"/>');
			</xsl:if>
		</script>
	</xsl:template>
</xsl:stylesheet>

	