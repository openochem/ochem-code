<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="configurators-common.xslt" />
	<xsl:include href="../../model/inc/structure-optimisation-commons.xslt" />
	
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
					<option value="DimeNetPP">DimeNetPP</option>
					<option value="HamNet">HamNet</option>

					<option value="NMPN">NMPN</option>
					<option value="Megnet">Megnet</option>
					<option value="MoGAT">MoGAT</option>
					<option value="MAT">MAT</option>
<!-- 				<option value="MEGAN">MEGAN</option>  -->	
					<option value="RGCN">RGCN</option>
					<option value="GNNFilm">GNNFilm</option>
					<option value="HDNNP2nd">HDNNP2nd</option>
					<option value="rGIN">rGIN</option>
					<option value="CMPNN">CMPNN</option>
					
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
		</table>
		<br/>

		<br/>
		<table class="configuration">
 			<h1>By default RDKIT is used to generate 3D for PAiNN and DimeNetPP. You can also select another method</h1>
 			<tr>
				<xsl:call-template name="structure-optimisation"/>
			</tr>			
		</table>	
		<br/>

		<script language="javascript">
				setValue("method", '<xsl:value-of select="//attachment/configuration/modelConfiguration/method"/>');
				setValue("batch", '<xsl:value-of select="//attachment/configuration/modelConfiguration/batch"/>');
				setValue("nepochs", '<xsl:value-of select="//attachment/configuration/modelConfiguration/nepochs"/>');
		</script>
	</xsl:template>
</xsl:stylesheet>
