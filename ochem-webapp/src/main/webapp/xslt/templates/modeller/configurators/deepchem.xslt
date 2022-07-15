<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="configurators-common.xslt" />
	
	<xsl:template name="configuration-content">
		<title>Model builder - DEEPCHEM</title>
		<h1>Select DEEPCHEM method</h1>
		<h2>DEEPCHEM: https://github.com/deepchem/deepchem</h2>
		
		<br/><table class="configuration">
			<tr>
				<td>Epochs (max: 100):<input type="text" class="small" name="nepochs" value="25"/> </td>
			</tr>
			<tr>
				<td>Early stopping fraction [0.05,1]:<input type="text" class="small" name="early" value="0.2"/> </td>
			</tr>
			<tr>
				<td>Use stereochemistry<br/></td>
				<td><input type="checkbox" checked="checked" name="chirality"/></td>
			</tr>
			<tr>
				<td>Train with all data (can overfit)<br/></td>
				<td><input type="checkbox" name="shuffle"/></td>
			</tr>
			<tr>			
			<td>SMILES augmentation for training (only TEXTCNN): 1 - no augmentation:<br/></td>
				<td>
				<select name="augmentation">
					<option value="1">no augmentation</option>
					<option value="5">5</option>
					<option value="10" selected="1">10</option>
					<option value="25">25</option>
				</select>
				</td>
			</tr>
			<tr>
			<td>SMILES augmentation for application (only TEXTCNN): 1 - no augmentation:<br/></td>
				<td>
				<select name="augmentApplySet">
					<option value="1">no augmentation</option>
					<option value="5">5</option>
					<option value="10" selected="1">10</option>
					<option value="25">25</option>
					<option value="50">50</option>
				</select>
				</td>
			</tr>
		</table>
	
		
		<br/><table class="configuration">
			<tr>
				<td>Model:</td>
				<td>
				<select name="method">
					<option value="DAG">Directed acyclic graph</option>
					<option value="GRAPH_CONV">Graph Convolution</option>
					<option value="MPNN">Message Passing Neural Networks (MPNN)</option>
					<option value="TEXTCNN"  selected="1">TEXTCNN</option>
				<!-- 	<option value="WEAVE">WEAVE</option>  -->
				</select></td>
			</tr>
		</table>	
<!-- 
		<h1>Method specific hyperparameters</h1>
		<div class="configuration">
		<table>
			<tr>
				<td>Learning Rate<br/></td>
				<td><input type="text" class="small" name="learning_rate" value="0.001"/> </td>
			</tr>
			<tr>
				<td>Dropout <br/></td>
				<td><input type="text" class="small" name="dropout" value="0.25"/> </td>
			</tr>
			<tr>
				<td>Graph Convolution: Dense Layer Size<br/></td>
				<td><input type="text" class="small" name="dense_layer_size" value="128"/> </td>
			</tr>
			<tr>
				<td>Graph Convolution: Convolution layers <br/></td>
				<td><input type="text" class="small" name="graph_conv_layers" value="64,64"/> </td>
			</tr>
			<tr>
				<td>MPNN and WEAVE: Number of Hidden neurones<br/></td>
				<td><input type="text" class="small" name="n_hidden" value="100"/> </td>
			</tr>
			<tr>
				<td>MPNN: M<br/></td>
				<td><input type="text" class="small" name="M" value="5"/> </td>
			</tr>
			<tr>
				<td>MPNN: T<br/></td>
				<td><input type="text" class="small" name="T" value="3"/> </td>
			</tr>
			<tr>
				<td>TEXTCNN: Embedding<br/></td>
				<td><input type="text" class="small" name="n_embedding" value="75"/> </td>
			</tr>
		</table>
		</div>
 -->	
		<div class="other parameters">
		<table>
			<tr>
				<td>Balance imbalanced data with augmentation<br/></td>
				<td><input type="checkbox" name="balance"/></td>
		</tr><tr>
				<td>Use CRS for reactions<br/></td>
				<td><input type="checkbox" name="crs"/></td>
			</tr>	

		</table>
		</div>
		
		
		<script language="javascript">
			<xsl:if test="//method = 'DEEPCHEM'">
				setCheckbox("crs", '<xsl:value-of select="//attachment/configuration/modelConfiguration/crs"/>');
				setCheckbox("balance", '<xsl:value-of select="//attachment/configuration/modelConfiguration/balance"/>');
				setCheckbox("sanitize", '<xsl:value-of select="//attachment/configuration/modelConfiguration/sanitize"/>');
				setValue("nepochs", '<xsl:value-of select="//attachment/configuration/modelConfiguration/nepochs"/>');
				setValue("early", '<xsl:value-of select="//attachment/configuration/modelConfiguration/early"/>');
				setCheckbox("chirality", '<xsl:value-of select="//attachment/configuration/modelConfiguration/chirality"/>');
				setCheckbox("shuffle", '<xsl:value-of select="//attachment/configuration/modelConfiguration/shuffle"/>');
				setValue("augmentation", '<xsl:value-of select="//attachment/configuration/modelConfiguration/augmentation"/>');
				setValue("augmentApplySet", '<xsl:value-of select="//attachment/configuration/modelConfiguration/augmentApplySet"/>');
				
				setValue("method", '<xsl:value-of select="//attachment/configuration/modelConfiguration/method"/>');
				setValue("learning_rate", '<xsl:value-of select="//attachment/configuration/modelConfiguration/learning_rate"/>');
				setValue("dropout", '<xsl:value-of select="//attachment/configuration/modelConfiguration/dropout"/>');
				setValue("n_hidden", '<xsl:value-of select="//attachment/configuration/modelConfiguration/n_hidden"/>');
				setValue("dense_layer_size", '<xsl:value-of select="//attachment/configuration/modelConfiguration/dense_layer_size"/>');
				setValue("graph_conv_layers", '<xsl:value-of select="//attachment/configuration/modelConfiguration/graph_conv_layers"/>');
				setValue("M", '<xsl:value-of select="//attachment/configuration/modelConfiguration/M"/>');
				setValue("T", '<xsl:value-of select="//attachment/configuration/modelConfiguration/T"/>');
				setValue("n_embedding", '<xsl:value-of select="//attachment/configuration/modelConfiguration/n_embedding"/>');
			</xsl:if>
		</script>
	</xsl:template>
</xsl:stylesheet>

	