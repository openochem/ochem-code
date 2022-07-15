<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />

	<xsl:template name="content">
		<title>Model point profile</title>
		<style type="text/css">
			.conditions {color: green; font-size: 80%; float: right; clear: right; width: 500px; text-align: right;}
			.article-data {font-size: 8pt; margin-top: 7px; margin-left: 0px; margin-bottom: 7px;}
			.selected-point {border: 2px solid #900 !important;}
			img {vertical-align: middle;}
			
			.block-image DIV {font-size: 120%; color: #333; white-space: nowrap; margin-top: 3px;}
			.itunes-right H1 {font-weight: bold; border-bottom: 1px solid gray; font-family: Helvetica;}
		</style>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script language="javascript" src="js/commons/dynamic-select.js"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript" src="js/blocks/neighbours.js"></script>
		<script language="javascript">	
			var mm_id = getParams["mm_id"];
			var ep_id = getParams["ep_id"];
			var model_id = getParams["model_id"];
			var row_num = getParams["row_num"];
			var task_num = getParams["task_num"];
			$(function()
			{
				if (!getParams["predictionVector"])
				{
					$("#neighbours option[value='correlation']").remove();
					$("#neighbours option[value='rank-correlation']").remove();
					$("#neighbours option[value='euclidean']").remove();
				}
			});
			
			
		</script>
		<table width="100%" >
			<tr><td class="itunes-up">
				<h1>Prediction neighbors explorer<a class="infolink" href="https://docs.ochem.eu/display/MAN/Prediction+neighbors" target="_blank"></a></h1>
				The training set compounds nearest to the selected prediction
			</td></tr>
			<tr><td class="itunes-right">
				<h1>The predicted compound</h1>
				<div id="main-profile"></div>
				<div id="neighbours">
					<h1>Nearest training set neighbours</h1>
					Similarity measure: <select name="type" filter="1">
					<option value="structural-similarity">Structural similarity</option>
					<option value="correlation">Correlation (prediction space)</option>
					<option value="rank-correlation">Rank correlation (prediction space)</option>
					<option value="euclidean">Euclidian distance (prediction space)</option>
					</select>
					<div id="nearest-neighbours-browser"></div>		
				</div>	
			</td></tr>		
			<tr><td align="right">
					<a href="javascript:window.closeTab()">close</a>
			</td></tr>
		</table>
	</xsl:template>
	
</xsl:stylesheet>