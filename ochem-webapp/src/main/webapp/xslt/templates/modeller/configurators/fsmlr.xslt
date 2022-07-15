<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="configurators-common.xslt" />
	
	<xsl:template name="configuration-content">
		<style type="text/css">
			.left{padding-left:15px;font-size: 10.5pt; padding-top:5px;}
			.left input{padding-right:5px;}
			.left label{text-indent: 10px; margin-left: 1.5em}
		</style>
		<title>Model builder - FSMLR</title>
		<h1><doc term="FSMLR">Configure Fast Stepwise Miltiple Linear Regression (FSMLR)  method</doc></h1>
		<table class="configuration">
			<tr>
				<td class="left">
					<label>Shrinkage<small>(Lies in interval [0,1])</small></label><input type="text" class="small" name="shrinkage" value="1"/>
					<label>Max Delta Descriptors</label><input type="text" class="small" name="delta" value="20"/>
				</td>
			</tr>
			<tr>
				<td class="left">
					<label>nfolds<small>(number of folds)</small></label><input type="text" class="small" name="nfolds" value="10"/>
					<label>Descriptor Num Factor</label><input type="text" class="small" name="num-factor" value="1.0"/>			
				</td>		
			</tr>
			<tr>
				<td class="left">
					<label>Descriptor Extra Factor</label><input type="text" class="small" name="ex-factor" value="5.0"/>
					<label>Start</label><input type="text" class="small" name="start" value="0"/>
				</td>
			</tr>
		</table>	
		
		<script language="javascript">
			<xsl:if test="//method = 'FSMLR'">
				setValue("shrinkage", '<xsl:value-of select="//attachment/configuration/modelConfiguration/shrinkage"/>');
				setValue("delta", '<xsl:value-of select="//attachment/configuration/modelConfiguration/delta"/>');
				setValue("nfolds", '<xsl:value-of select="//attachment/configuration/modelConfiguration/nfolds"/>');
				setValue("num-factor", '<xsl:value-of select="//attachment/configuration/modelConfiguration/numFactor"/>');
				setValue("ex-factor", '<xsl:value-of select="//attachment/configuration/modelConfiguration/extraFactor"/>');
				setValue("start", '<xsl:value-of select="//attachment/configuration/modelConfiguration/start"/>');
			</xsl:if>
		</script>
		
	</xsl:template>
</xsl:stylesheet>