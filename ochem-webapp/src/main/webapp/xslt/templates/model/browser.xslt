<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />

	<xsl:template name="content">
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript" src="js/browsers/model-browser.js"></script>
		<script language="javascript">
			include.plugins('view');
			var sampleBrowser = new ModelBrowser("model", "model");
			$(document).ready(function() {sampleBrowser.initialize();});
		</script>
		<title>Model applier</title>
		<table width="100%">
			<tr><td class="itunes-up">
				<h1><doc term="Using+model+browser">Models browser</doc></h1>
			</td></tr>
			<tr><td class="itunes-right">
				<h1>Step 1. Select a model from the list</h1>
				Filter by model name:
				<input type="text" name="query" filter="1"/>
				and property name:
				<input type="text" name="proquery" filter="1"/>
				 Order by:
				 <select name="order" filter="1">
				 	<option value="creation">creation time</option>
					<option value="access">last access time</option>
				 	<option value="modification">last modification time</option>
				 </select> 
				<a href="javascript:sampleBrowser.request(true)">
					[refresh]
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
			</td></tr>
		</table>
	</xsl:template>
</xsl:stylesheet>