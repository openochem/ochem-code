<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	
	<xsl:include href="../../helper.xslt" />
	<xsl:include href="../model/inc/select-compounds.xslt" />
	
	<xsl:template name="content">
	
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script language="javascript" src="js/commons/dynamic-select.js"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<style type="text/css">
			.comparison-table TH { border: 3px solid white; border-bottom: 1px solid black; font-weight: bold; text-align: center;}
			.comparison-table TH, .comparison-table TD {padding: 6px;}
			.comparison-table TR {border-bottom: 1px solid #AAA !important;}
			.green, .green A {color: #090 !important;}
			.gray, .gray A {color: #555 !important;}
			#MolBrowser DIV {float: left; }
			#MolBrowser DIV IMG {border: 1px solid gray;}
		</style>
		<table height="100%" width="100%">
			<tr><td class="itunes-up">
				<img src="img/icons/compare.png"/>
				<h1><span class="setcompare">SetCompare</span>: Comparison results</h1>
				The comparison summary of the two selected sets
				</td></tr>
			<tr>
				<td class="itunes-right">
					The following table shows the features (molecular descriptors) that were significantly over-represented in one of the two sets.
					<br/>
					It includes appearance counts of the features in each set and the p-Value of such a distribution.<br/>
					<a href="setcomparison/exportResults.do">Export results as a CSV file</a><br/>
					<br/>
					<xsl:if test="//param[@key='errors1'] &gt; 0">
						There are <a href="setcomparison/molecules.do?errors=1&amp;set=1" tab="Failed molecules"><xsl:value-of select="//param[@key='errors1']"/> errors</a> (duplicates, failures) in the smaller set .<br/>
					</xsl:if>
					<xsl:if test="//param[@key='errors2'] &gt; 0">
						There are <a href="setcomparison/molecules.do?errors=1&amp;set=2" tab="Failed molecules"><xsl:value-of select="//param[@key='errors2']"/> errors</a> (duplicates, failures) in the larger set.<br/>
					</xsl:if>
					<div class="pager-strip">
						<span><b class="showed">none</b> of <b class="total">none</b></span>
						<div class="pgr">
						</div>
					</div>
					
					<table class="comparison-table">
						<tr>
							<th>Descriptor</th>
							<th>In set 1<br/>(<xsl:value-of select="//param[@key='set1']"/>)</th>
							<th>In set 2<br/>(<xsl:value-of select="//param[@key='set2']"/>)</th>
							<th>Enrichment factor</th>
							<th>p-Value</th>
						</tr>
						<tbody id="Browser">
						</tbody>
					</table>
					<div class="pager-strip">
						<span><b class="showed">none</b> of <b class="total">none</b></span>
						<div class="pgr">
						</div>
					</div>	
				</td>
			</tr>
		</table>

	<div id="mol-browser" class="invisible">
		<div class="hd">Select your molecule set</div>
		<div class="bd" id="MolBrowser">
			

		</div>
	</div>
		
		<script language="javascript">
			include.plugins('view');
			var resBrowser = new Browser();
			resBrowser.url = "setcomparison/list.do";
			resBrowser.useTableRows = true;
			resBrowser.itemElement="setcompare-result";
			resBrowser.itemTemplate = 'js/templates/setcompare-result.ejs';
			resBrowser.pager.selectors.pager = ".pgr";
			
			$(document).ready(function(){
				resBrowser.initialize();
			});
			
			var inSet1Total = 1*'<xsl:value-of select="//param[@key='set1']"/>'.split(" ")[0];
			var inSet2Total = 1*'<xsl:value-of select="//param[@key='set2']"/>'.split(" ")[0];
			
			
			
		</script>
		
		</xsl:template>
</xsl:stylesheet>
