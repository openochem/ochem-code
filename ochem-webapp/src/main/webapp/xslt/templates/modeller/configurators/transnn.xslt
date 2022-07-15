<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="configurators-common.xslt" />
	
	<xsl:template name="configuration-content">
		<title>Model builder - Transformer</title>
		<h1>Transformer Neural Networks (TRANS-CNN)</h1>
		
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
			<td>SMILES augmentation for training:<br/></td>
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
			<td>SMILES augmentation for application:<br/></td>
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
				<td>Batch size [8,128]:<input type="text" class="small" name="batch" value="64"/> </td>
			</tr><tr>
			<td>Fixed learning rate<br/></td>
				<td><input type="checkbox" checked="checked" name="fixedrate"/></td>
			</tr><tr>
				<td>Balance imbalanced data with augmentation<br/></td>
				<td><input type="checkbox" name="balance"/></td>
			</tr>
		</table>
		<br/>


		<script language="javascript">
			<xsl:if test="//method = 'TRANSNN'">
				setCheckbox("fixedrate", '<xsl:value-of select="//attachment/configuration/modelConfiguration/fixedrate"/>');
				setValue("batch", '<xsl:value-of select="//attachment/configuration/modelConfiguration/batch"/>');
				setCheckbox("balance", '<xsl:value-of select="//attachment/configuration/modelConfiguration/balance"/>');
				
				setValue("nepochs", '<xsl:value-of select="//attachment/configuration/modelConfiguration/nepochs"/>');
				setValue("early", '<xsl:value-of select="//attachment/configuration/modelConfiguration/early"/>');
				setCheckbox("chirality", '<xsl:value-of select="//attachment/configuration/modelConfiguration/chirality"/>');
				setCheckbox("shuffle", '<xsl:value-of select="//attachment/configuration/modelConfiguration/shuffle"/>');
				setValue("augmentation", '<xsl:value-of select="//attachment/configuration/modelConfiguration/augmentation"/>');
				setValue("augmentApplySet", '<xsl:value-of select="//attachment/configuration/modelConfiguration/augmentApplySet"/>');				
			</xsl:if>
		</script>
	</xsl:template>
</xsl:stylesheet>
