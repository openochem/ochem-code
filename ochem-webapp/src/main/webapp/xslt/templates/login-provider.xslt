<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />
	<xsl:template name="content">
	
	<script type="text/javascript">

	var current_url = new URL(window.location.href);
	var login = current_url.searchParams.get("login");
	if (login) {
		window.opener.location.replace("user/profile.do?login=" + login + "&amp;first-login=true");
		window.close();
	}
	
	</script>
	
	<style type="text/css">
		DIV.comp-block-login
		{
			padding: 30px;
			width: 770px;
		}
		DIV.comp-block-login DIV {padding: 20px; border-bottom: 1px solid #999; margin: 5px;}
		DIV.comp-block-login TABLE TD {padding: 5px 15px 5px 0px;}
		TABLE.comp-block-login
		{
			border: outset 1px #AAA;
		  	border-spacing: 30pt 10pt;
			background-color:#FFF;
			margin: 10px 10px 10px 10px;
			position: relative;
			
		}
		.button-link {display: block; margin-top: 15px;}
		.comp-block-login INPUT {border: 1px solid #555; width: 200px;}
		TABLE.comp-block-login TD {padding: 15px 10px 15px 10px;}
		.comp-block-login H1 {font-size: 200%; margin-bottom: 13px; border-bottom: 1px solid black;}
		DIV.comp-block-login H2 {font-weight: bold; color: #000; margin-bottom: 6px; font-size: 130%;}
	</style>
	
	<title>Login Via Provider</title>
	<div class="comp-block-login">
	<xsl:if test="not(/model/@nomoreusers = 'true')">
	   <div>
				<h2>Please, choose the identity provider for your login.</h2>
				
	      <xsl:for-each select="/model/others/handlerInfo">
	      	<br class="cb"/>
				<a class="button-link">
					<xsl:attribute name="href">login/withProvider.do?provider=<xsl:value-of select="@name" /></xsl:attribute>
					<xsl:value-of select="@name" />
				</a>
	      	<br class="cb"/>
	      </xsl:for-each>
	      
	   </div>
	</xsl:if>
	</div>
	</xsl:template>
</xsl:stylesheet>