<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />
	
	<xsl:template name="content">
		<title>Batch edit browser</title>
		<style type="text/css">
			.conditions {color: green; font-size: 80%; float: right; clear: right; width: 500px; text-align: right;}
			.article-data {font-size: 80%; margin-top: 7px; margin-left: 0px; margin-bottom: 7px;}
			TABLE.cfilter {margin-top: 5px;}
			img {vertical-align: middle;}
			.time {font-size: 80%;}
			.errorcomment {font-size: 80%; color: #600;}
			.padded {margin-left: 20px !important;}
			.textual-condition {width: 230px;}
			.yui-skin-sam .yui-dialog .ft span.default button {
				color:#000000; !important;
				font-weight:bold; !important;
			}
		</style>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript" src="js/commons/dynamic-select.js"></script>
		
		<script language="javascript" src="js/browsers/batchedit-browser.js"></script>
		<script language="javascript">
			include.plugins('view');
			var batchEdit = new BatchEditBrowser();
			$(document).ready(
				function() {
					batchEdit.initialize();
				}
			);
			
		</script>
		<table width="100%" height="100%" cellspacing="0">
		<tr>
			<td class="itunes-up" colspan="2">
				<h1>Batch edit browser</h1>
				Edit selected molecules in the batch
			</td>
		</tr>
		<tr>
			<td class="itunes-right">
				<a target="_parent" href="epbrowser/show.do?render-mode=full">[Compounds properties browser]</a>
				<div class="pager-strip">
					<b class="showed">none</b> of <b class="total">none</b>
				</div>
				<div id="pager">
				</div>
				<div id="Browser">
				</div>
			</td>
		</tr>
	</table>
	</xsl:template>
</xsl:stylesheet>