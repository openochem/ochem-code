<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="configurators-common.xslt" />
	
	<xsl:template name="configuration-content">
		<title>Model builder - KGCNN configurations</title>
		<h1>Graph Neural Networks: KGCNN configurations </h1>
		
		<br/><table class="configuration">
			<tr>
				<td>Model:</td>
				<td>
				<select name="method">
					<option value="AttFP" selected="1">AttFP</option>
					<option value="ChemProp">ChemProp</option>
					<option value="GIN">GIN</option>
					<option value="GCN">GCN</option>
					<option value="GINE">GINE</option>
					<option value="PAiNN">PAiNN</option>
					<option value="Schnet">Schnet</option>
					<option value="GATv2">GATv2</option>
					<option value="GAT">GAT</option>
					<option value="GraphSAGE">GraphSAGE</option>
					<option value="GCN">GCN</option>
					<option value="DimeNetPP">DimeNetPP</option>
					<option value="HamNet">HamNet</option>
				</select></td>
			</tr>
		</table>

		<br/>
		<table class="configuration">
			<tr>
			<td>Epochs (max: 1000):<input type="text" class="small" name="nepochs" value="200"/> </td>
			</tr><tr>
			<td>Batch size:<input type="text" class="small" name="batch" value="32"/> </td>
			</tr>
			<tr>
				<td>Sanitize data<br/></td>
				<td><input type="checkbox" name="sanitize"/></td>
			</tr>	
		</table>
		<br/>


		<script language="javascript">
			<xsl:if test="//method = 'KERAS'">
				setCheckbox("sanitize", '<xsl:value-of select="//attachment/configuration/modelConfiguration/sanitize"/>');
				setValue("method", '<xsl:value-of select="//attachment/configuration/modelConfiguration/method"/>');
				setValue("x", '<xsl:value-of select="//attachment/configuration/modelConfiguration/x"/>');
				setValue("batch", '<xsl:value-of select="//attachment/configuration/modelConfiguration/batch"/>');
				setValue("nepochs", '<xsl:value-of select="//attachment/configuration/modelConfiguration/nepochs"/>');
			</xsl:if>
		</script>
	</xsl:template>
</xsl:stylesheet>
