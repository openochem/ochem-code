<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="2.0">
	<xsl:include href="../helper.xslt" />
	<xsl:template name="content">
		<title>Action not allowed for anonymous user</title>
		<link rel="stylesheet" type="text/css" href="css/batch.css" />
		<link rel="stylesheet" type="text/css" href="css/smoothness/jquery-ui-1.9.1.custom.css" />
		<style type="text/css">
			.summary {width: 800px; font-size: 120%; margin: 10px;}
			.red {background-color:#FEE; border-bottom: 1px solid black;}
			.summary TD {padding: 10px;}
		</style>
		<center>
		<table>
		<tr><td class="big-padding ui-widget">
				<table class="summary">
					<tr class="red"><td>Action not allowed</td></tr>
					<tr><td>Dear guest, <br/><br/>certain actions (like creating new records) are restricted for guest users.<br/><br/>If you would like to take full advantage of OCHEM, please
						<a href="/user/newuser.do"> register here</a> (it will only take a minute).<br/><br/>
						If you already have an account, please <a href="/login/show.do">log in here</a>.
						Accounts registered with corporate e-mails will be automatically validated. Otherwise contact info@ochem.eu to validate your acount.
						
					</td></tr>
				</table>
		</td></tr>
		</table>
		</center>
	</xsl:template>
</xsl:stylesheet>