<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />
	
	<xsl:template name="content">
		<title>Message browser</title>
		<style type="text/css">
			.conditions {color: green; font-size: 80%; float: right; clear: right; width: 500px; text-align: right;}
			.article-data {font-size: 80%; margin-top: 7px; margin-left: 0px; margin-bottom: 7px;}
			img {vertical-align: middle;}
		</style>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script language="javascript" src="js/commons/dynamic-select.js"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript" src="js/browsers/message-browser.js"></script>
		<script language="javascript">
			include.plugins('view');
			var messageBrowser = new MessageBrowser();
			$(document).ready(function() {messageBrowser.initialize();});
		</script>
		<style type="text/css">
			.fontweight {font-weight: bold;}
			.fontsize {font-size: 20px;}
			.message-block {border: 1px solid black; width: 95%; background-color: #EFF5FB; margin: 10px 10px 10px 10px; padding: 20px; clear: both;}
			.message-block I {display: block; float: right; width: 200px; text-align: right;}
			.message-body {border: 0px solid black; margin: 15px 0px 15px 20px; text-align: justify}
		</style>
		<table width="100%" height="100%" cellspacing="0">
		<tr>
			<td class="itunes-up" colspan="2">
				<img src="img/icons/e-mail_icon.jpg"/><h1 id="page-title">Message Box</h1>
			</td>
		</tr>
		<tr>
			<td class="itunes-left">
				<br/>
				<div class="fontweight fontsize" id="inbox">
					<a action="inbox" name="inbox" filter="1">INBOX</a>
				</div>
				<br/>
				<div class="fontsize" id="sent">
					<a action="sent" name="sent" filter="1">SENT</a>
				</div>
			</td>
			<td class="itunes-right">
				<div class="upper-command-panel">
					[<a action="unread" name="unread message">show unread</a>]
					[<a action="showall" name="show all message">show all</a>]
				</div>
				<div class="pager-strip">
					<span><b class="showed">none</b> of <b class="total">none</b></span>
					<div id="pager" class="pgr">
					</div>
				</div>
				
				<div id="Browser">
				</div>
				
				<div class="pager-strip">
					<span><b class="showed">none</b> of <b class="total">none</b></span>
					<div id="pager" class="pgr">
					</div>
				</div>
			</td>
		</tr>
	</table>
	</xsl:template>
</xsl:stylesheet>