<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="configurators-common.xslt" />
	
	<xsl:template name="configuration-content">
		<title>Edge Attention based Multi-relational Graph Convolutional Networks (EAGCNG)</title>
		<h1>Configure EAGCNG method</h1>
		
		<table class="configuration">
			<tr>
				<td>Epochs <br/></td>
				<td><input type="text" class="small" name="nepochs" value="1000"/> </td>
			</tr>
			<tr>
				<td>Batch Size <br/></td>
				<td><input type="text" class="small" name="batchsize" value="512"/> </td>
			</tr>
			<tr>
				<td>Dropout <br/></td>
				<td><input type="text" class="small" name="dropout" value="0.2"/> </td>
			</tr>
			<tr>
				<td>Learning Rate <br/></td>
				<td><input type="text" class="small" name="learningRate" value="0.01"/> </td>
			</tr>
			<tr>
				<td>Weight Decay <br/></td>
				<td><input type="text" class="small" name="weightDecay" value="0.0001"/> </td>
			</tr>
			<tr><td>Method</td><td>
				<select name="method">
					<option value="concate" checked="true">Concatenate</option>
					<option value="weighted">Weighted</option>
				</select>
				</td>
			</tr>
			<tr><td>Normalization</td><td>
				<select name="normalisation">
					<option value="BATCH">batch</option>
					<option value="LAYER" checked="true">layer</option>
				</select>
				</td>
			</tr>
			<tr><td>Activation function</td><td>
				<select name="activation_function">
					<option value="RELU" checked="true">relu</option>
					<option value="LRELU">leaky relu</option>
					<option value="ELU">elu</option>
				</select>
				</td>
			</tr>
			<tr>
				<td>Use GATE mode</td>
				<td><input type="checkbox" name="gate"/></td>
			</tr>
		</table>
		<br/>
		<table>
		<tr>
				<td>n_sgc1<br/></td>
				<td><input type="text" class="medium" name="n_sgc1" value="30,15,15,15,15"/> </td>
			</tr>

			<tr>
				<td>n_sgc2<br/></td>
				<td><input type="text" class="medium" name="n_sgc2" value="60,30,30,30,30"/> </td>
			</tr>

			<tr>
				<td>n_sgc3<br/></td>
				<td><input type="text" class="medium" name="n_sgc3" value=""/> </td>
			</tr>

			<tr>
				<td>n_den<br/></td>
				<td><input type="text" class="medium" name="n_den" value="64,32,12"/> </td>
			</tr>		
			<tr/><br/>	
			<tr>
				<td>Maximal molecule size: <br/></td>
				<td><input type="text" class="small" name="molsize" value="100"/> (GPU memory increases as N*N; N ~ 100 requires ca 4GB) </td>
			</tr>
			</table>
		<br/>
			
		<script language="javascript">
			<xsl:if test="//method = 'EAGCNG'">
				setCheckbox("sanitize", '<xsl:value-of select="//attachment/configuration/modelConfiguration/sanitize"/>');
				setValue("activation_function", '<xsl:value-of select="//attachment/configuration/modelConfiguration/activationFunction"/>');
				setValue("normalisation", '<xsl:value-of select="//attachment/configuration/modelConfiguration/normalisation"/>');
				setCheckbox("gate", '<xsl:value-of select="//attachment/configuration/modelConfiguration/gate"/>');
				setValue("nepochs", '<xsl:value-of select="//attachment/configuration/modelConfiguration/nepochs"/>');
				setValue("batchsize", '<xsl:value-of select="//attachment/configuration/modelConfiguration/batchsize"/>');
				setValue("dropout", '<xsl:value-of select="//attachment/configuration/modelConfiguration/dropout"/>');
				setValue("learningRate", '<xsl:value-of select="//attachment/configuration/modelConfiguration/learningRate"/>');
				setValue("weightDecay", '<xsl:value-of select="//attachment/configuration/modelConfiguration/weightDecay"/>');
				setValue("n_sgc1", '<xsl:value-of select="//attachment/configuration/modelConfiguration/n_sgc1"/>');
				setValue("n_sgc2", '<xsl:value-of select="//attachment/configuration/modelConfiguration/n_sgc2"/>');
				setValue("n_sgc3", '<xsl:value-of select="//attachment/configuration/modelConfiguration/n_sgc3"/>');
				setValue("n_den", '<xsl:value-of select="//attachment/configuration/modelConfiguration/n_den"/>');
				setValue("method", '<xsl:value-of select="//attachment/configuration/modelConfiguration/method"/>');
				setValue("seed", '<xsl:value-of select="//attachment/configuration/modelConfiguration/seed"/>');
				setValue("molsize", '<xsl:value-of select="//attachment/configuration/modelConfiguration/molsize"/>');
			</xsl:if>
		</script>
	</xsl:template>
</xsl:stylesheet>