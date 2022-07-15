<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="configurators-common.xslt" />
	
	<xsl:template name="configuration-content">
		<title>Model builder - Deep Neural Networks (under development)</title>
		<h1>Configure DNN method</h1>
		<table class="configuration">
		<tr>
				<td>Batch size <br/></td>
				<td><input type="text" class="small" name="batchSize" value="1024"/> </td>
			</tr>
			<tr>
				<td>Epochs <br/></td>
				<td><input type="text" class="small" name="epochs" value="10000"/> </td>
			</tr>
			<tr>
				<td>Training/internal validation set ratio <br/></td>
				<td><input type="text" class="small" name="ratio" value="0.8"/> </td>
			</tr>
			<tr>
				<td>Model:</td>
				<td>
				<select name="modeltype">
					<option value="dense7new" selected="1">Dense7 (seven layers) Net</option>
					<option value="dense_log2">Number of neurons decreases by half at each layer</option>
					<option value="dense_exp">Number of neurons exponentially decreases at each layer</option>
<!-- 				<option value="merck">Merk's model from https://github.com/Merck</option>  -->
				</select></td>
			</tr>
			<tr>
				<td>Optimisation method:</td>
				<td>
				<select name="optimizertype">
					<option value="adam" selected="1">Adam</option>
					<option value="sgd">MomentumSGD</option>
					<option value="rmsprop">RMSprop</option>
					<option value="smorms3">SMORMS3</option>
				</select></td>
			</tr>
			<tr>
				<td>Activation function:</td>
				<td>
				<select name="activation">
					<option value="relu" selected="1">relu</option>
					<option value="elu">elu</option>
					<option value="prelu">prelu</option>
					<option value="sigmoid">sigmoid</option>
				</select></td>
			</tr>
			<xsl:if test="//user/superuser = 'true'"> 
				<tr>
					<td>Additional Parameters<br/><small>(separated by spaces)</small></td>
					<td><input type="text" name="additionalParam" value=""/></td>
				</tr>
				<tr>
					<td>Model Specific Parameters<br/><small>(separated by spaces)</small></td>
					<td><input type="text" name="additionalParamModel" value=""/></td>
				</tr>
				<tr>
					<td>Experimental Parameters<br/><small>(conditions merging)</small></td>
					<td><input type="text" name="experimentalParam" value=""/></td>
				</tr>
			</xsl:if>
		</table>	
		<script language="javascript">
			<xsl:if test="//method = 'DNN'">
				setValue("activation", '<xsl:value-of select="//attachment/configuration/modelConfiguration/activation"/>');
				setValue("optimizertype", '<xsl:value-of select="//attachment/configuration/modelConfiguration/optimizertype"/>');
				setValue("modeltype", '<xsl:value-of select="//attachment/configuration/modelConfiguration/modeltype"/>');
				setValue("batchSize", '<xsl:value-of select="//attachment/configuration/modelConfiguration/batchSize"/>');
				setValue("epochs", '<xsl:value-of select="//attachment/configuration/modelConfiguration/epochs"/>');
				setValue("ratio", '<xsl:value-of select="//attachment/configuration/modelConfiguration/ratio"/>');
				setValue("additionalParam", '<xsl:value-of select="//attachment/configuration/modelConfiguration/additionalParam"/>');
				setValue("additionalParamModel", '<xsl:value-of select="//attachment/configuration/modelConfiguration/additionalParam"/>');
				setValue("experimentalParam", '<xsl:value-of select="//attachment/configuration/modelConfiguration/experimentalParam"/>');
			</xsl:if>
		</script>
	</xsl:template>
</xsl:stylesheet>