<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />
	
	<xsl:template name="content">
		<div id="area">
			<img style="float: left; margin-right: 30px;" src="img/facebook-logo.png"/>
			<h1>Your OCHEM account created!</h1>
			We have created an OCHEM account based on your facebook account and sent you a confirmation email.
			Your new user ID on OCHEM is <b><xsl:value-of select="//user/@login"/></b><br/>
			You can always login into OCHEM either your facebook account or your new user ID.
			
			<br/><br/>
			You can edit your personal data in <a href="user/show.do">your profile</a>.
			
			<br/><br/>
			<a href="/home/show.do?render-mode=full" target="_top">Enjoy and proceed to the home page!</a>
		</div>
	
	</xsl:template>
</xsl:stylesheet>