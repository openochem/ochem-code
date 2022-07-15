<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="configurators-common.xslt" />
	
	<xsl:template name="configuration-content">
		<title>Model builder - Deep Learning neural network</title>
		<h1>Configure DEEP (eep Learning neural network) method</h1>
	<script language="javascript">
			<xsl:if test="//method = 'DEEP'">
				setValue("epochs", '<xsl:value-of select="//attachment/configuration/modelConfiguration/epochs"/>');
				setValue("ratio", '<xsl:value-of select="//attachment/configuration/modelConfiguration/ratio"/>');
			</xsl:if>
			$(document).ready(function(){
				updateVisibility();	
			});
		</script>
	</xsl:template>
</xsl:stylesheet>