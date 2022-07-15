<?xml version="1.0" encoding="UTF-8"?>


<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	
	<xsl:include href="../../helper.xslt" />
	<xsl:include href="../model/inc/data-preprocessing-commons.xslt" />
	
	<xsl:template name="content">
	
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script language="javascript" src="js/commons/dynamic-select.js"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript">
			include.plugins('view');
			var browser = new Browser("descriptors");
			browser.url = "descriptorscalculator/list.do";
			browser.itemElement = "descriptorsRow";
			browser.itemTemplate = 'js/templates/descriptorsrow.ejs';
			browser.pager.selectors.pager = ".pgr";
			browser.useTableRows = true;
			browser.options.highlight_row = "no";
			browser.onItemDrawn = function()
			{
				if(this.currentEntity.error != "" &amp;&amp; this.currentEntity.error != undefined)
					this.currentBlock.addClass("highlighterror");
			}
			
			function doShowChange()
			{
				if ($("#show_all").is(":checked"))
					$(".hideable").removeClass("invisible");
				else
					$(".hideable").addClass("invisible");
			}
			
			$(document).ready(function() {
				browser.initialize();
				$("#show_all").change(doShowChange);
			});
		</script>
		<style>
			.maintable 
			{
				width: 100%;
			}
			
			.sumtable
			{
				padding: 50px;
			}
			
			.maintable TD, .sumtable TD
			{
				border: 1px solid #DDDDDD;
				padding: 5px;
			}
			.maintable TH
			{
				border: 1px solid #DDDDDD;
				padding: 5px;
				font-weight: bold;
			}
			
			.filterpanel
			{
				dbackground-color: #EEEEFF;
				padding: 10px;
				margin-bottom: 2px;
			}
			
			.filterpanel LABEL {margin-right: 10px;}
			
			.right
			{
				align:right;
			}
			
			.hideable {text-align: right;}
		</style>
		<table height="100%" width="100%">
			<tr><td class="itunes-up silver">
				<img src="img/icons/mol-descriptors.png"/>
				<h1>DescriptorsCalculator: Calculation results</h1>
				</td></tr>
			<tr>
				<td class="itunes-right">
					<table class="sumtable">
						<tr><td>Total molecules:</td><td><xsl:value-of select="//descriptorsSummary/@total"/></td></tr>
						<tr><td>Correctly calculated:</td><td><xsl:value-of select="//descriptorsSummary/@valid"/></td></tr>
						<tr><td>Errors:</td><td><xsl:value-of select="//descriptorsSummary/@error"/></td></tr>
					</table><br/><br/>
					<a href="descriptorscalculator/download.do" class="fancy-button">Proceed to the descriptors download page</a><br/><br/><br/>
						
					<div class="filterpanel">
						<xsl:if test="//descriptorsSummary/@error &gt; 0">
							<input type="checkbox" name="errors_only" filter="1" id="errors_only"/><label for="errors_only">Display errors only</label>
						</xsl:if>  
						<input type="checkbox" name="show_all" id="show_all"/><label for="show_all">Show descriptor values (first 150 columns)</label>
					</div>
					<div class="pager-strip">
						<span><b class="showed">none</b> of <b class="total">none</b></span>
						<div class="pgr">
						</div>
					</div>
					<table class="maintable"> 
						<thead>
						<tr>
							<th></th>
							<th></th>
							<xsl:for-each select="//descriptorsSummary/columns">
								<th class="hideable invisible"><xsl:value-of select="."/></th>
							</xsl:for-each>
						</tr>
						</thead>
						<tbody id="Browser">
						</tbody>
					</table>
					<div class="pager-strip">
						<span><b class="showed">none</b> of <b class="total">none</b></span>
						<div class="pgr">
						</div>
					</div>
					<br/>
					<a href="descriptorscalculator/download.do" class="fancy-button">Proceed to the descriptors download page</a>
					<br/><br/>
				</td>
			</tr>
		</table>
		
		</xsl:template>
</xsl:stylesheet>