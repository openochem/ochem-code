<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />
	
	<xsl:template name="content">
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript" src="js/browsers/journal-browser.js"></script>
		<script language="javascript">
			include.plugins('view');
			var sampleBrowser = new JournalBrowser();
			$(document).ready(function() {sampleBrowser.initialize();});
		</script>
		<title>Journal Browser</title>
		<table height="100%" width="100%">
			<tr><td class="itunes-up"><h1><doc term="Browser+of+scientific+publications#Browserofscientificpublications-journalBrowserBrowserofjournals">Journals browser</doc></h1>
			Please search existing journals before adding new ones</td></tr>
			<tr>
				<td class="itunes-right">
				Type part of journal name to filter:
					<input type="text" name="query" filter="1"/>
				
					<a href="javascript:sampleBrowser.request(true)">
						[refresh]
					</a>
					<a action="edit" title="Add new Journal" >
						[Create new <img src="img/icons/new.gif"/>]
					</a>
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