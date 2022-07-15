<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="configurators-common.xslt" />
	
	<xsl:template name="configuration-content">
		<title>Model builder - CNF</title>
		<h1>Configure Convolutional Neural Network Fingerprint method</h1>
		
		<br/><table class="configuration">
			<tr>
				<td>Epochs (max: 10000):<input type="text" class="small" name="nepochs" value="1000"/> </td>
			</tr><tr>
				<td>Use stereochemistry<br/></td>
				<td><input type="checkbox" checked="checked" name="chirality"/></td>
			</tr><tr>
				<td>Use sanitization of molecules<br/></td>
				<td><input type="checkbox" name="sanitize"/></td>
			</tr><tr>
				<td>Train with all data (can overfit)<br/></td>
				<td><input type="checkbox" name="shuffle"/></td>
			</tr><tr>			
				<td>Early stopping fraction [0.05,1]:<input type="text" class="small" name="early" value="0.1"/> </td>
			</tr><tr>			
			<td>SMILES augmentation for training:<br/></td>
				<td>
				<select name="augmentation">
					<option value="-1">online augmentation</option>
					<option value="1">no augmentation</option>
					<option value="5">5</option>
					<option value="10" selected="1">10</option>
					<option value="25">25</option>
				</select>
				</td>
			</tr><tr>
			<td>SMILES augmentation for application:<br/></td>
				<td>
				<select name="augmentApplySet">
					<option value="1">no augmentation</option>
					<option value="5">5</option>
					<option value="10" selected="1">10</option>
					<option value="25">25</option>
				</select>
				</td>
			</tr>
		</table>
		
		
		<br/><table class="configuration">
			<tr>
				<td>FPDim<br/></td>
				<td><input type="text" class="small" name="FPDim" value="32"/> </td>
			</tr>
			<tr>
				<td>nLayers<br/></td>
				<td><input type="text" class="small" name="nLayers" value="8"/> </td>
			</tr>
			<tr>
				<td>FFNET_dim<br/></td>
				<td><input type="text" class="small" name="FFNET_dim" value="128,64"/> </td> 
			</tr>
			<tr>
				<td>Batch size<br/></td>
				<td><input type="text" class="small" name="batch" value="32"/> </td>
			</tr>
			<tr>
				<td>Learning rate<br/></td>
				<td><input type="text" class="small" name="rate" value="0.001"/> </td>
			</tr>
<!-- 			<tr>
				<td>Kernel width (only for v0, v1.1 and v2.1)</td>
				<td><input type="text" class="small" name="filter_size" value="2"/> </td>
			</tr>
			<tr>
				<td>Alpha: Kernel averaging (only for v2, alpha is in [0,1] range)</td>
				<td><input type="text" class="small" name="nfp_alpha" value="0.5"/> </td>
			</tr>
 -->		<tr>
			<td>Normalise:</td>
			<td>
				<select name="normalisation">
					<option value="NONE">No normalisation</option>
					<option value="LAYER" selected="1">Layer normalization</option>
					<option value="BATCH">Batch  normalization</option>
				</select>
			</td>
			</tr><tr>
<!--			<td>Loss function:</td>
			<td>
				<select name="error_function">
					<option value="RMSE" selected="1">RMSE - regression</option>
					<option value="ENTROPY">ENTROPY - classification</option>
					<option value="TAXONOMY">TAXONOMY  error function</option> 
					<option value="MULTILABELNN">BP-MLL error function</option>
				</select>
			</td>
			</tr><tr>
-->			<td>Activation function:</td>
			<td>
				<select name="activation_function">
					<option value="ELU" selected="1">ELU</option>
					<option value="RELU">RELU</option>
					<option value="CRELU">CRELU</option>
					<option value="LRELU">Leaky Relu</option>
					<option value="SWISH">SWISH</option>
				</select>
			</td>
			</tr><tr>
			<td>Convolutional Neural Network Fingerprint (CNF) function:</td>
			<td>
				<select name="fingerprint_function">
					<option value="V0">cascade processing with fixed Kernel width (v0)</option>
					<option value="V11">kernel concatenate with fixed Kernel width (v1.1)</option>
					<option value="V12">kernel concatenate with variable Kernel width: 1,...,nLayers (v1.2)</option>
					<option value="V21">combined cascade/concatenate kernel (by alpha) with fixed Kernel width (v2.1)</option>
					<option value="V22" selected="1">combined cascade/concatenate kernel (by alpha) with variable Kernel width: 1,...,nLayers (v2.2)</option>
				</select>
			</td>
			</tr><tr>
				<td>Select tokenizer<br/></td>
				<td>
				<select name="tokenizer">
					<option value="ATOM">atom</option>
					<option value="CHAR" selected="1">character</option>
					<option value="OLD">fixed set of characters</option>
				</select>
				</td>
			</tr><tr>
				<td>Dropout (-1, [0-0.5])<br/></td>
				<td><input type="text" class="small" name="dropout" value="0.2"/> </td>
			</tr><tr>
				<td>Use CRS for reactions<br/></td>
				<td><input type="checkbox" name="crs"/></td>
			</tr>	
		</table>	
		<br/>
		<br/>
<!--		
		<table class="configuration">
			Provide taxonomy information (parent1,child1;parent2,child2): <input type="text" name="taxonomy" value="" style="width: 300px;"/><br/>
		</table>
-->		
	
		<script language="javascript">
			<xsl:if test="//method = 'CNF'">
				setValue("tokenizer", '<xsl:value-of select="//attachment/configuration/modelConfiguration/tokenizer"/>');
				setValue("fingerprint_function", '<xsl:value-of select="//attachment/configuration/modelConfiguration/type"/>');
				setValue("activation_function", '<xsl:value-of select="//attachment/configuration/modelConfiguration/activationFunction"/>');
				setValue("dropout", '<xsl:value-of select="//attachment/configuration/modelConfiguration/dropout"/>');
				setValue("batch", '<xsl:value-of select="//attachment/configuration/modelConfiguration/batch"/>');
				setValue("rate", '<xsl:value-of select="//attachment/configuration/modelConfiguration/rate"/>');
				setValue("FPDim", '<xsl:value-of select="//attachment/configuration/modelConfiguration/FPDim"/>');
				setValue("nLayers", '<xsl:value-of select="//attachment/configuration/modelConfiguration/nLayers"/>');
				setValue("FFNET_dim", '<xsl:value-of select="//attachment/configuration/modelConfiguration/FFNET_dim"/>');
				setCheckbox("crs", '<xsl:value-of select="//attachment/configuration/modelConfiguration/crs"/>');

				setValue("early", '<xsl:value-of select="//attachment/configuration/modelConfiguration/early"/>');
				setValue("nepochs", '<xsl:value-of select="//attachment/configuration/modelConfiguration/nepochs"/>');
				setCheckbox("chirality", '<xsl:value-of select="//attachment/configuration/modelConfiguration/chirality"/>');
				setCheckbox("shuffle", '<xsl:value-of select="//attachment/configuration/modelConfiguration/shuffle"/>');
				setValue("augmentation", '<xsl:value-of select="//attachment/configuration/modelConfiguration/augmentation"/>');
				setValue("augmentApplySet", '<xsl:value-of select="//attachment/configuration/modelConfiguration/augmentApplySet"/>');
				setCheckbox("sanitize", '<xsl:value-of select="//attachment/configuration/modelConfiguration/sanitize"/>');

				setValue("error_function", '<xsl:value-of select="//attachment/configuration/modelConfiguration/errorFunction"/>');
				setValue("nfp_alpha", '<xsl:value-of select="//attachment/configuration/modelConfiguration/nfp_alpha"/>');
				setValue("highway", '<xsl:value-of select="//attachment/configuration/modelConfiguration/highway"/>');
				setValue("filter_size", '<xsl:value-of select="//attachment/configuration/modelConfiguration/filterSize"/>');
				setValue("normalisation", '<xsl:value-of select="//attachment/configuration/modelConfiguration/normalisation"/>');
				setValue("taxonomy", '<xsl:value-of select="//attachment/configuration/modelConfiguration/originalTaxonomy"/>');

			</xsl:if>
		</script>
	</xsl:template>
</xsl:stylesheet>
