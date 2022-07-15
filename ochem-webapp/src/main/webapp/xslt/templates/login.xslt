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
	
	<title>Login</title>
	<div class="comp-block-login">
	<h1>Please, login</h1>
	<xsl:if test="not(/model/@noanonymous = 'true') and /model/loginSettings/guestLogin = 'true'">
	  <div>
		<h2>Instant login</h2>
		In order to access OCHEM, you must login. 
		If you do not wish to register now, you can login as a guest. Guest users have access to less features than registered users.<br/> 
		<a href="login/login.do?anonymous=1&amp;render-mode=popup" class="button-link">Login as a guest</a>
		<br class="cb"/>
	  </div>
	</xsl:if>
	
	<xsl:if test="not(/model/@noanonymous = 'true') and /model/loginSettings/providerLogin = 'true'">
	  <div>
		<h2>Instant login with an Identity Provider</h2>
		You can also login with an external identity provider such as a social network or an LDAP server. 
		<a class="button-link" target="popup" href="login/chooseProvider.do?render-mode=popup">Choose Identity Provider</a>
		<!--  onclick="window.open('login/chooseProvider.do?render-mode=popup','name','width=600,height=400')" -->
		<br class="cb"/>
	  </div>
	</xsl:if>
	
	<xsl:if test="not(/model/@noanonymous = 'true') and /model/loginSettings/registerLogin = 'true'">
	
		<div>
		
		<form action="login/login.do" method="post" name="mainform" target="_parent" accept-charset="UTF-8,ISO-8859-1">
			
			<input type="hidden" name="render-mode" value="redirect"/>
			<xsl:if test="/model/loginSettings/oauth">
				<input type="hidden" name="oauth" value="true"/>
				<input type="hidden">
					<xsl:attribute name="name">clientID</xsl:attribute> 
					<xsl:attribute name="value">
					  <xsl:value-of select="/model/loginSettings/oauth/clientID"/>
					</xsl:attribute> 
				</input>
				<input type="hidden">
					<xsl:attribute name="name">redirectURI</xsl:attribute> 
					<xsl:attribute name="value">
					  <xsl:value-of select="/model/loginSettings/oauth/redirectURI"/>
					</xsl:attribute> 
				</input>
			</xsl:if>
			
			
				<h2>Already have an account?</h2>
					If you already have an account, please enter you login and password below:
			<table>
				<tr>
					<td>Login ID</td>
					<td>
						<input class="text2" type="text" name="login" />
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
						<input class="text2" type="password" name="pwd" />
					</td>
				</tr>
				<tr>
					<td colspan="2" height="30"	valign="bottom">
						<a onclick="document.forms[0].submit(); return false;" href="#" class="button-link">login</a>
						<a href="static/pwd-reminder.do" class="button-link">password reminder</a>
						<input type="image" value="submitname" src="submit-button.gif" width="1" height="1" border="0" alt="SUBMIT!" name="image" class="invisible"/>
					</td>
				</tr>
			</table>
		</form>
		</div>
		<xsl:if test="not(/model/@nomoreusers = 'true')">
		   <div>
					<h2>Join OCHEM - register a new user!</h2>
					Create a free account to upload data, create and apply QSAR models, screen chemical libraries and many more.
					Registered users can correct data uploaded by other registered users publish models. As a registered user, you can
		      configure flexible access policies for your data and models.
		      <br class="cb"/>
					<a class="button-link" href="user/newuser.do">Register a New User</a>
		      <br class="cb"/>
		   </div>
		</xsl:if>
	
	</xsl:if>
	
	</div>
	
	</xsl:template>
</xsl:stylesheet>