<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />

	<xsl:template name="content">
		<xsl:variable name="entity">unit</xsl:variable>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript">
			include.plugins('view');
			function ExportHistoryBrowser(controller, item)
			{
				this.controller = controller;
				Browser.call(this);
				this.itemElement = item;
				this.useTableRows = true;
				this.itemTemplate = 'js/templates/exporthistory.ejs';
				this.pager.selectors.pager = ".pgr";
				
				this.onItemDrawn = function()
				{
					this.currentBlock.find('i[title]').tooltip({showURL: false});
				}
				
				this.doDownload = function()
				{
						var win = openTab("Download file", webRoot+this.controller+"/download.do?id="+ this.currentEntity.id);
				}
			}			
			
			var sampleBrowser = new ExportHistoryBrowser("export", "export-action");
			$(document).ready(function() {sampleBrowser.initialize();});
		</script>
		<style type="text/css">
			.maintable TD
			{
				border: 1px solid #DDDDDD;
				vertical-align: top;
				padding: 5px;
			}
			
			.maintable TH
			{
				border: 1px solid #AAAAAA;
				vertical-align: top;
				text-align: center;
				padding: 5px;
				font-weight: bold;
			}
		</style>
		<title>Unit browser</title>
		<table height="100%" width="100%">
			<tr><td class="itunes-up">
				<xsl:call-template name="area-of-interest"/>
					<h1>Export history</h1>
					The history of your exports from OCHEM
				</td></tr>
			<tr>
				<td class="itunes-right" a="{$entity}">
					<div class="pager-strip">
						<span><b class="showed">none</b> of <b class="total">none</b></span>
						<div id="pager" class="pgr">
						</div>
					</div>
					
					<table class="maintable">
					<thead><th>Date</th><th>File name</th><th>Format</th><th>Rows exported</th><th>Bonus cost <a class="infolink" href="https://docs.ochem.eu/display/MAN/Bonus+points+system" target="_blank" title="More information on the bonus ponts on OCHEM"></a></th><th>Free bonuses used</th><th>Downloaded</th><th>Status</th></thead>
					<tbody id="Browser">
					</tbody>
					</table>
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