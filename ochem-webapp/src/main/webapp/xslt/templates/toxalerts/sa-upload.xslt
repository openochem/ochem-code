<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />
	
	<xsl:template name="content">
		<title>Structural alerts upload</title>
		<table height="100%" width="100%">
			<tr><td class="itunes-up">
				<img src="img/icons/batch_upload.gif"/>
				<h1><span class="toxalerts">ToxAlerts</span>: Upload structural alerts</h1>
				
				</td></tr>
			<tr>
				<td class="itunes-right">
					Here you can upload structural alerts from an Excel file.<br/>
					<br/>
					Obligatory columns: SMARTS, ARTICLEID, PROPERTY, NAME.<br/>
					Optional columns: DESCRIPTION, SMARTS_DESCRIPTION, COMMENT. You can provide additional consa-screen.xsltditions for each alert.<br/><br/>
					
					You may simplify complex SMARTS patterns by using the <a tab="SMARTS substitution variables" href="alerts/variables.do">SMARTS substitution variables</a><br/> 
					If you have any doubts, please download an <a href="documents/toxalert-upload-sample.xls">exemplary Excel file</a>.<br/><br/>
					
					<form action="alerts/uploadSubmit.do" enctype="multipart/form-data" method="post">
					Select an Excel file: <input type="file" name="file"/><br/>
					<input type="checkbox" name="private-upload"/> Upload privately <a class="infolink" title="If checked, the uploaded alerts will be uploaded privately and will be visible only to you"></a>
					<br/><br/>
					<input type="submit" value="Submit"/>
					</form>
					
				</td>
			</tr>
		</table>
		</xsl:template>
</xsl:stylesheet>
