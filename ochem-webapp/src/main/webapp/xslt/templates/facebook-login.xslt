<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />
	
	<xsl:template name="content">
		<div id="area">
			<img style="float: left; margin-right: 30px;" src="img/facebook-logo.png"/>
			<h1>Successful login!</h1>
			You have been successfully logged in with your facebook account.<br/><br/>
			Your OCHEM user ID <b><xsl:value-of select="//user/@login"/></b> was identified automatically using your email <b><xsl:value-of select="//user/E-mail"/></b>
			<br/>
			<a href="/home/show.do?render-mode=full" target="_top">Proceed to the home page</a>
		</div>
	
	</xsl:template>
</xsl:stylesheet>