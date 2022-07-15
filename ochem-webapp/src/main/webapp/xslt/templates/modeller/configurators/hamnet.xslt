<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="configurators-common.xslt" />
	
	<xsl:template name="configuration-content">
		<title>Model builder - HamNet</title>
		<h1>Conformation-Guided Molecular Representation with Hamiltonian Neural Networks (HamNet)</h1>
		
		
		<br/><table class="configuration">
			<tr>
				<td>Epochs (max: 100):<input type="text" class="small" name="nepochs" value="25"/> </td>
			</tr><tr>
				<td>Early stopping fraction [0.05,1]:<input type="text" class="small" name="early" value="0.2"/> </td>
			</tr><tr>
				<td>Use stereochemistry<br/></td>
				<td><input type="checkbox" checked="checked" name="chirality"/></td>
			</tr><tr>
				<td>Train with all data (can overfit)<br/></td>
				<td><input type="checkbox" name="shuffle"/></td>
			</tr>
		</table>
	
		
		<br/><table class="configuration">
			<tr>
				<td>Batch size:<input type="text" class="small" name="batch_size" value="64"/> </td>
			</tr><tr>
				<td>Dropout:<input type="text" class="small" name="dropout" value="0.3"/> </td>
			</tr><tr>
				<td>Learning rate <br/></td>
				<td><input type="text" class="small" name="learningRate" value="0.001"/> </td>
			</tr><tr>
			</tr><tr>
				<td>Use REFINE mode<br/></td>
				<td><input type="checkbox" name="refine"/></td>
			</tr>	
		</table>

	
		<script language="javascript">
			<xsl:if test="//method = 'TRANSNNI'">
				setValue("learningRate", '<xsl:value-of select="//attachment/configuration/modelConfiguration/learningRate"/>');
				setValue("dropout", '<xsl:value-of select="//attachment/configuration/modelConfiguration/dropout"/>');
				setValue("batch_size", '<xsl:value-of select="//attachment/configuration/modelConfiguration/batch_size"/>');
				setValue("optimizer", '<xsl:value-of select="//attachment/configuration/modelConfiguration/optimizer"/>');

				setCheckbox("refine", '<xsl:value-of select="//attachment/configuration/modelConfiguration/refine"/>');
				setCheckbox("sanitize", '<xsl:value-of select="//attachment/configuration/modelConfiguration/sanitize"/>');
				setValue("nepochs", '<xsl:value-of select="//attachment/configuration/modelConfiguration/nepochs"/>');
				setValue("early", '<xsl:value-of select="//attachment/configuration/modelConfiguration/early"/>');
				setCheckbox("chirality", '<xsl:value-of select="//attachment/configuration/modelConfiguration/chirality"/>');
				setCheckbox("shuffle", '<xsl:value-of select="//attachment/configuration/modelConfiguration/shuffle"/>');
			</xsl:if>
		</script>
	</xsl:template>
</xsl:stylesheet>
