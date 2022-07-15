<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />
	<xsl:include href="../model/inc/select-compounds.xslt" />
	
	<xsl:template name="content">
		<style type="text/css">
			.step {borsder: 1px solid #CCC; width: 650px; padding: 10px; margin-bottom: 5px;}
			.step TABLE TD {padding-bottom: 15px;}
			.step INPUT[type=checkbox] {margin-right: 5px;}
			.step {background-color: #FFF;}
			.step H1 {border-bottom: 1px solid black;}
			.submitform {background-color: #EEE; padding: 5px; cursor: hand; cursor: pointer;}
			.submitform:hover {background-color: #ADA;}
		</style>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script language="javascript" src="js/commons/dynamic-select.js"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<title>OCHEM predictor</title>
		<table width="100%">
			<tr><td class="itunes-up silver">
				<h1>OCHEM Predictor</h1>
				Predict a number of supported properties/activities for your chemical compounds
			</td></tr>
			<tr><td class="itunes-right">
			
			<form method="post" enctype="multipart/form-data" action="predictor/start.do">
				<div class="step">
					<h1>What would you like to predict?</h1>
					<div style="margin-left: 20px;">
						<xsl:for-each select="//others/model"> 
							<input type="checkbox" name="model" value="{publicId}"/><xsl:value-of select="@featuredName"/><br/>
						</xsl:for-each>
						<br/>
						<small><a href="model/select.do?render-mode=popup">Browse the full list of public models</a></small>
					</div>
				</div>
					
				<div class="step">
					<h1>Select the compounds to be predicted</h1>
					<div style="margin-left: 20px;">
					<xsl:call-template name="select-compounds"/>
					</div>
				</div>
				<input type="checkbox" name="disable-cache" /> Disable prediction cache (predictions but not descriptors will be recalculated)<br/>
				<xsl:if test="//session/user/superuser='true'">
					<input type="checkbox" id="force_cache" name="force_cache"/><label for="force_cache">Force recalculation of descriptors in cache </label><br/>
				</xsl:if>
				<input type="submit" name="submit" class="fancy-button" value="Run predictions!" onclick="return onSubmit();"/>
			</form>
			</td></tr>
		</table>
		<script language="javascript">
			function onSubmit()
			{
				if ($("input[name=model]:checked").length == 0)
				{
					window.alert("Please, select at least one property to predict.")
					return false;
				}
				
				return true;
			}
		</script>
	</xsl:template>
</xsl:stylesheet>