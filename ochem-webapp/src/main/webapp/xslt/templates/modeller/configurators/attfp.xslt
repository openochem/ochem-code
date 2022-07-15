<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="configurators-common.xslt" />
	
	<xsl:template name="configuration-content">
		<title>Model builder - ATTFP</title>
		<h1>AttentiveFP</h1>
		
		<br/>
		<table class="configuration">
			<tr>
			<td>Epochs (max: 1000):<input type="text" class="small" name="nepochs" value="200"/> </td>
			</tr><tr>
			<td>Batch size:<input type="text" class="small" name="batch" value="64"/> </td>
			</tr><tr>
			<td>patience_reduce:<input type="text" class="small" name="patience_reduce" value="10"/> </td>
			</tr><tr>
			<td>patience_early:<input type="text" class="small" name="patience_early" value="40"/> </td>
			</tr><tr>
			<td>radius:<input type="text" class="small" name="radius" value="2"/> </td>
			</tr><tr>
			<td>T:<input type="text" class="small" name="T" value="2"/> </td>
			</tr><tr>
			<td>fp_dim:<input type="text" class="small" name="fp_dim" value="200"/> </td>
			</tr><tr>
			<td>cosineT:<input type="text" class="small" name="cosineT" value="14"/> </td>
			</tr><tr>
			<td>Dropout:<input type="text" class="small" name="dropout" value="0.2"/> </td>
			</tr><tr>
			<td>lr:<input type="text" class="small" name="lr" value=" 0.0032"/> </td>
			</tr><tr>
			<td>weight_decay:<input type="text" class="small" name="weight_decay" value="0.00001"/> </td>
			</tr>
			<tr>
				<td>simpleO<br/></td>
				<td><input type="checkbox" name="simpleO"/></td>
			</tr>			
			<tr>
				<td>lngru<br/></td>
				<td><input type="checkbox" name="lngru"/></td>
			</tr>			
			<tr>
				<td>singleT<br/></td>
				<td><input type="checkbox" name="singleT"/></td>
			</tr>			
			<tr>
				<td>Cosine<br/></td>
				<td><input type="checkbox" name="cosine" checked="checked"/></td>
			</tr>			
 			<tr>
				<td>Early stopping<br/></td>
				<td><input type="checkbox" name="early" checked="checked"/></td>
			</tr>			
			<tr>
				<td>Maximal molecule size: <br/></td>
				<td><input type="text" class="small" name="molsize" value="200"/> (GPU memory increases with N and can fail due to memory size) </td>
			</tr>		
			</table>
		<br/>

		<script language="javascript">
			<xsl:if test="//method = 'ATTFP'">
				setCheckbox("simpleO", '<xsl:value-of select="//attachment/configuration/modelConfiguration/simpleO"/>');
				setCheckbox("singleT", '<xsl:value-of select="//attachment/configuration/modelConfiguration/singleT"/>');
				setValue("batch", '<xsl:value-of select="//attachment/configuration/modelConfiguration/batch"/>');
				setValue("nepochs", '<xsl:value-of select="//attachment/configuration/modelConfiguration/nepochs"/>');
				setValue("patience_reduce", '<xsl:value-of select="//attachment/configuration/modelConfiguration/patience_reduce"/>');
				setValue("patience_early", '<xsl:value-of select="//attachment/configuration/modelConfiguration/patience_early"/>');
				setValue("radius", '<xsl:value-of select="//attachment/configuration/modelConfiguration/radius"/>');
				setValue("T", '<xsl:value-of select="//attachment/configuration/modelConfiguration/T"/>');
				setValue("fp_dim", '<xsl:value-of select="//attachment/configuration/modelConfiguration/fp_dim"/>');
				setValue("cosineT", '<xsl:value-of select="//attachment/configuration/modelConfiguration/cosineT"/>');
				setValue("dropout", '<xsl:value-of select="//attachment/configuration/modelConfiguration/dropout"/>');
				setValue("lr", '<xsl:value-of select="//attachment/configuration/modelConfiguration/lr"/>');
				setValue("weight_decay", '<xsl:value-of select="//attachment/configuration/modelConfiguration/weight_decay"/>');
				setCheckbox("cosine", '<xsl:value-of select="//attachment/configuration/modelConfiguration/cosine"/>');
				setCheckbox("early", '<xsl:value-of select="//attachment/configuration/modelConfiguration/early"/>');
				setCheckbox("lngru", '<xsl:value-of select="//attachment/configuration/modelConfiguration/lngru"/>');
				setValue("molsize", '<xsl:value-of select="//attachment/configuration/modelConfiguration/molsize"/>');
			</xsl:if>
		</script>
	</xsl:template>
</xsl:stylesheet>


