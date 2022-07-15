<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />
	<xsl:template name="content">
		<style type="text/css">
			DIV.comment {width: 400px;}
		</style>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript"> 
			include.plugins('view'); 
			var sampleBrowser = new Browser(); 
			sampleBrowser.itemElement = "action"; 
			sampleBrowser.itemTemplate = "js/templates/userhistory-item.ejs";
			sampleBrowser.url = "history/listuser.do"; 
			$(document).ready(function() {sampleBrowser.initialize();}); 
			</script>
		<title>User-related changes</title>
		<table width="100%">
			<tr>
				<td class="itunes-up">
				<img src="img/icons/userchanges.png"/>
				<h1>Changes history</h1>
				Tracking user-related changes made to the data
				</td>
			</tr>
			<tr>
				<td class="itunes-right">
				<input type="text" name="query" filter="1"/>
					<a href="javascript:sampleBrowser.request(true)">
						[refresh]
					</a>
				<div class="pager-strip">
					<b class="showed">none</b> of <b class="total">none</b>
				</div>
				<div id="pager">
				</div>
				<div id="Browser">
				</div>
				<div class="pager-strip">
					<b class="showed">none</b> of <b class="total">none</b>
				</div>
			</td>
			</tr>
		</table>
	</xsl:template>
</xsl:stylesheet>