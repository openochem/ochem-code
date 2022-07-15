<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />

	<xsl:template name="content">
		<style type="text/css">
			.comment
			{
				margin-left: 20px;
			}
			
			.user 
			{
				margin-top: 10px;
				color: #666;
			}
			
			.user A {color: #66A;}
			
			.day {
				border-top: 1px solid #AAA; 
				padding-top: 20px;
				color: #AAA;
				font-family: Arial;
				margin-top: 20px;
				margin-bottom: 10px;
			}
			
			.time {
				color: gray;
				margin-right: 10px;
			}
		</style>
		<title>User actions feed</title>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/browsers/useractions-browser.js?ver=1.7.6"></script>
		<table height="100%" width="100%">
			<tr>
				<td class="itunes-up" colspan="2">
					<h1><img src="img/icons/user-32.png"/>User actions feed</h1>
					What are the OCHEM users doing? View the news feed here.
				</td>
			</tr>
			<tr>
				<td class="itunes-right">
					<div id="info">Only the last week's events are shown here.</div>
					<div id="Browser">
					</div>
				</td>
			</tr>
		</table>
		
		<script language="javascript">
			loadUserEvents();
		</script>
		
	</xsl:template>
</xsl:stylesheet>