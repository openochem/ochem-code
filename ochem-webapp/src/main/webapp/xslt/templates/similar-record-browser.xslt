<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />

	<xsl:template name="content">
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript" src="js/browsers/similar-record-browser.js"></script>
		<script language="javascript">
			include.plugins('view');
			var sampleBrowser = new SimilarRecordBrowser();
			$(document).ready(function() {sampleBrowser.initialize();});
		</script>
		<title>Similar record browser</title>
		<style type="text/css">
			.conditions {color: green; font-size: 80%; float: right; clear: right; width: 500px; text-align: right;}
			.article-data {font-size: 80%; margin-top: 7px; margin-left: 0px; margin-bottom: 7px;}
			.rightfloat {float: right; clear: right; text=align: right}
			img {vertical-align: middle;}
		</style>
		<title>Similar record browser</title>
		<table height="100%" width="100%">
			<tr><td class="itunes-up">
				<xsl:call-template name="area-of-interest"/>
						<img src="img/icons/property.png"/>
						<h1>Similar records browser</h1>
						Records have similar molecular structure and name within the article
				</td></tr>
			<tr>
				<td class="itunes-right" >
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
			<tr>
				<td align="right"><a href="javascript:window.closeTab();">close</a></td>
			</tr>
		</table>
	</xsl:template>
</xsl:stylesheet>