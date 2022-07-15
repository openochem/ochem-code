<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="2.0">
	<xsl:include href="../helper.xslt" />
	<xsl:include href="inc/license.xslt" />

	<xsl:template name="content">
		<style type="text/css">
		.disable{
pointer-events:none;
background:lightgrey;
}
			.style {padding: 10px; width: 900px; margin: 10px 10px 10px 50px;}
			.trStyle {border: 1px none #999;}
			.trStyle > TD {padding-top: 10px; padding-bottom: 20px;}
			.trStyle.header {font-weight: bold; border-bottom: 1px solid black; text-align: left !important; }
			.trStyle.header TD {padding-top: 15px; padding-bottom: 0px;}
			TD {padding: 3px;}
			LABEL {width: 200px}
			INPUT, TEXTAREA {border: 1px solid black; padding: 2px 2px 2px 2px; margin: 1px 1px 1px 1px;}
			.warning {color: red !important;}
			.warn {color: red !important; font-size: 11px; font-style: italic;}
			.redirect {background-color: #FFDDDD; font-size: 90%; padding: 10px; color: red;}
			.message {font-size: 11px; font-style: italic;}
			<xsl:if test="not(@id)">
				.required {background-color: #FFA;}
			</xsl:if>
					
		</style>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script language="javascript" src="js/blocks/user.js"></script>
		<xsl:apply-templates select="/model/authError"/>
		<xsl:apply-templates select="/model/loginSettings/oauth"/>
	</xsl:template>
	<xsl:template match="/model/loginSettings/oauth">
		<title>Grant Permissions to External Application</title><table height="100%" width="100%">
			<tr>
				<td class="itunes-up" colspan="2">
				<img src="img/icons/article.png"/>
					<h1>External application is requesting permissions to access your OCHEM account data.</h1>
				</td>
			</tr>
			
			<tr>
				<td class="itunes-up" colspan="2">
					<form action="oauth/authorize.do" method="post" name="mainform" target="_parent" accept-charset="UTF-8,ISO-8859-1">
					
					<input type="hidden" name="render-mode" value="redirect"/>
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
					
					
						<h2>The following data will be provided:</h2>
						<table>
							<tr>
								<td>Login ID</td>
								<td>
									<input class="text2 disable" type="text" name="login">
										<xsl:attribute name="value">
											<xsl:value-of select="/model/session/user/@login"/>
										</xsl:attribute>
									</input>
								</td>
							</tr>
							<tr>
								<td>First Name</td>
								<td>
									<input class="text2 disable" type="text" name="firstName">
										<xsl:attribute name="value">
											<xsl:value-of select="/model/session/user/FirstName"/>
										</xsl:attribute>
									</input>
								</td>
							</tr>
							<tr>
								<td>Last Name</td>
								<td>
									<input class="text2 disable" type="text" name="lastName">
										<xsl:attribute name="value">
											<xsl:value-of select="/model/session/user/LastName"/>
										</xsl:attribute>
									</input>
								</td>
							</tr>
							<tr>
								<td colspan="2" height="30"	valign="bottom">
									<a onclick="document.forms[0].submit(); return false;" href="#" class="button-link">Grant Access</a>
									<input type="image" value="submitname" src="submit-button.gif" width="1" height="1" border="0" alt="SUBMIT!" name="image" class="invisible"/>
								</td>
							</tr>
						</table>
					</form>			
				</td>
			</tr>

		</table>
		
	</xsl:template>	
	<xsl:template match="/model/authError">
		<title>Unauthorized</title><table height="100%" width="100%">
			<tr>
				<td class="itunes-up" colspan="2">
				<img src="img/icons/article.png"/>
					<h1>You cannot access this page if not logged in.</h1>
					
					<xsl:value-of select="message"/>
				</td>
			</tr>

		</table>
		
	</xsl:template>	
</xsl:stylesheet>