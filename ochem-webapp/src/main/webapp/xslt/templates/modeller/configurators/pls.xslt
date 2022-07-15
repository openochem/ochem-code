<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="configurators-common.xslt" />
	
	<xsl:template name="configuration-content">
		<title>Model builder - PLS</title>
		<h1>Configure PLS method</h1>
		<table class="configuration">
			<tr>
				<td>Number of latent variables:</td><td><input type="text" name="num-latent-vars" value="0" id="num"/></td>
			</tr>
			<tr>
				<td><input type="checkbox" name="optimize" id="optimize" checked="checked"/>Optimize the number of latent variables automatically</td>
			</tr>			
			<tr>
				<td><br></br><input type="checkbox" name="limitrange"/> Limit predicted values to the training set range</td>
			</tr>
		</table>	
		<script language="javascript">
			$("#optimize").change(function(){
				if ($(this).is(":checked"))
					$("#num").attr("disabled", "disabled").val("0");
				else
					$("#num").removeAttr("disabled").val("5");
			});
			$("#optimize").change();
			
			<xsl:if test="//method = 'PLS'">
				var numVars = '<xsl:value-of select="//attachment/configuration/modelConfiguration/numLatentVariables"/>';
				if (numVars != '')
				{
					$("#optimize").setChecked(1*numVars == 0).change();
					if (1*numVars > 0)
						$("#num").val(numVars);
				}
			</xsl:if>
		</script>
	</xsl:template>
</xsl:stylesheet>