<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="configurators-common.xslt" />
	
	<xsl:template name="configuration-content">
		<title>Upload model file</title>

		<h1>Upload model file</h1>
		<xsl:if test="//message/message">
		<div class="warning">
			<xsl:value-of select="//message/message"/>
		</div>
		</xsl:if>
<!--    <input type="checkbox" name="xls"/> -->
		<br/>
		<label> Excel sheet with model parameters : </label>
		<input type="file" name="xls-file"/><br/><br/>		
	</xsl:template>
</xsl:stylesheet>

<!-- Toxicity/xslt/templates/modeller/configurators/upload-file.xslt tee 20100420 -->
