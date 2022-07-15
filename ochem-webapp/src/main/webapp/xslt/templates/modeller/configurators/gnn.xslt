<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="configurators-common.xslt" />
	
	<xsl:template name="configuration-content">
		<title>Model builder - GNN</title>
		<h1>Graph Neural Networks</h1>
		
		<table class="configuration">
			<tr>
			<td>Model:</td>
			<td>
				<select name="gnn">
					<option value="GIN" selected="1">GIN</option>
					<option value="GAIN">GAIN</option>
					<option value="GGRNet">GGRNet</option>
				</select>
			</td>
			</tr>
		</table>
		<br/>
		<table class="configuration">
			<tr>
			<td>Epochs (max: 1000):<input type="text" class="small" name="nepochs" value="200"/> </td>
			</tr><tr>
			<td>Batch size:<input type="text" class="small" name="batch" value="32"/> </td>
			</tr><tr>
			<td>Patience:<input type="text" class="small" name="patience" value="10"/> </td>
			</tr><tr>
			<td>Patience early:<input type="text" class="small" name="patienceearly" value="40"/> </td>
			</tr><tr>
			<td>LR decay:<input type="text" class="small" name="lr_decay" value="0.5"/> </td>
			</tr><tr>
			<td>Dimension:<input type="text" class="small" name="dim" value="95"/> </td>
			</tr>
		</table>
		<br/>

		<script language="javascript">
			<xsl:if test="//method = 'GNN'">
				setValue("gnn", '<xsl:value-of select="//attachment/configuration/modelConfiguration/gnn"/>');
				setValue("batch", '<xsl:value-of select="//attachment/configuration/modelConfiguration/batch"/>');
				setValue("patienceearly", '<xsl:value-of select="//attachment/configuration/modelConfiguration/patienceearly"/>');
				setValue("patience", '<xsl:value-of select="//attachment/configuration/modelConfiguration/patience"/>');
				setValue("lr_decay", '<xsl:value-of select="//attachment/configuration/modelConfiguration/lr_decay"/>');
				setValue("nepochs", '<xsl:value-of select="//attachment/configuration/modelConfiguration/nepochs"/>');
				setValue("dim", '<xsl:value-of select="//attachment/configuration/modelConfiguration/dim"/>');
			</xsl:if>
		</script>
	</xsl:template>
</xsl:stylesheet>
