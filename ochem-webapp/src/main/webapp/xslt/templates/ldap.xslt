<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />
	<xsl:template name="content">
	
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
	
	<title>Login with LDAP</title>
	<div class="comp-block-login">
	<xsl:if test="not(/model/@nomoreusers = 'true')">
	   <div>
				
	      <form action="login/ldap.do" method="post" name="mainform" target="_parent" accept-charset="UTF-8,ISO-8859-1">
		
			<h2>Authenticate with a remote LDAP server.</h2>
			Select your server from the list below and fill in your username and password.
		<table>
			<tr>
			<td>
				LDAP Server
			</td>
			<td>
				<select name="ldapServer">
					<xsl:for-each select="/model/others/ldapServer">
						<option value="{name}"><xsl:value-of select="name"/></option>
					</xsl:for-each>
				</select>
			</td>
			</tr>
			<tr>
				<td>Username</td>
				<td>
					<input class="text2" type="text" name="username" />
				</td>
				<xsl:if test="starts-with(@web-root, 'https')">
				<td rowspan="3">
				<img src="img/PossitiveSSL_tl_trans.gif"/>
				</td>
				</xsl:if>
			</tr>
			<tr>
				<td>Password</td>
				<td>
					<input class="text2" type="password" name="password" />
				</td>
			</tr>
			<tr>
				<td colspan="2" height="30"	valign="bottom">
					<a onclick="document.forms[0].submit(); return false;" href="#" class="button-link">login</a>
					<input type="image" value="submitname" src="submit-button.gif" width="1" height="1" border="0" alt="SUBMIT!" name="image" class="invisible"/>
				</td>
			</tr>
		</table>
	</form>
	   </div>
	</xsl:if>
	</div>
	</xsl:template>
</xsl:stylesheet>