<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="configurators-common.xslt" />
	<xsl:include href="../../model/inc/data-preprocessing-commons.xslt" />
	
	<xsl:template name="configuration-content">
		<title>Model builder - data preprocessing options</title>
		<h1>Select the preferred data preprocessing options</h1>
		<xsl:call-template name="data-preprocessing"/>
	</xsl:template>
	
</xsl:stylesheet>