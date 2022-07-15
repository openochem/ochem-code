<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />

	<xsl:template name="content">
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript" src="js/browsers/authors-browser.js"></script>
		<script language="javascript">
			include.plugins('view');
			var sampleBrowser = new AuthorBrowser("author", "author");
			$(document).ready(function() {sampleBrowser.initialize();});
		</script>
		<table height="100%" width="100%">
			<tr><td class="itunes-up"><h1>Authors browser</h1>
			Please search existing authors before adding new ones</td></tr>
			<tr>
				<td class="itunes-right">
				Type part of author last name to filter:
					<input type="text" name="query" filter="1"/>
					<a href="javascript:sampleBrowser.request(true)">
						[refresh]
					</a>
					<a action="edit">
						[new author]
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