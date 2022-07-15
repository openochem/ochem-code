<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />
	
	
	<xsl:template name="content">
	<div class="single-page">
	<h1>Forgot your password? We can help</h1>
	Please enter your username or email and we will send you your password.<br/><br/>
	<form action="mail/action.do" method="post" name="mainform">
		<input type="hidden" name="action" value="reminder"/>
		<table class="comp-block-login" width="38%">
			<tr>
				<td>Login ID</td>
				<td>
					<input class="text2" type="text2" name="login" />
				</td>
			</tr>
			<tr>
				<td height="20" valign="center" colspan="2"><b>OR</b>
				</td>
			</tr>
			<tr>
				<td>Email address</td>
				<td>
					<input class="text2" type="text" name="email" />
				</td>
			</tr>
			<tr>
				<td height="30" valign="bottom" colspan="2"><br/>
					<a class="button-link" onclick="document.forms[0].submit(); return false;" href="#">Submit</a>
				</td>
			</tr>
		</table>
	</form>	
	</div>
	</xsl:template>
</xsl:stylesheet>