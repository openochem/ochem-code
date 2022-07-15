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
		</style>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script language="javascript" src="js/commons/dynamic-select.js"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript" src="js/browsers/dotprofile-browser.js"></script>
		<script language="javascript">		
			include.plugins('view');
			var dotBrowser = new DotprofileBrowser();
			var descriptorsBrowser = new DescriptorsBrowser();
			
			$(document).ready(function() 
			{
				dotBrowser.filters.setValue("model_id", getParams["model_id"]);
				dotBrowser.initialize();
				descriptorsBrowser.filters.setFromUrl();
				descriptorsBrowser.initialize();
				tabView = new YAHOO.widget.TabView("demo");
			});		
		</script>
		<table width="100%" >
			<tr><td class="itunes-up">
				<h1>Model point profile</h1>
				Detailed information about one model point
			</td></tr>
			<tr>
				<td>			
					<div id="demo" class="yui-navset"> 
					    <ul class="yui-nav"> 
					    	<li class="selected">
					    		<a href="#tab1"><em>Record</em></a>
					    	</li>
					    	<li>
						       	<a href="#tab2"><em>Model data</em></a>
						    </li>
					    </ul> 
					    <div class="yui-content">
					    	<div class="dotbrowser">
					    		<div class="pager-strip">
								<b class="showed">none</b> of <b class="total">none</b>
								</div>
								<div id="pager">
								</div>		
								<div id="DotBrowser">
								</div>	
							</div>
							<div class="descriptorsbrowser">
								<div class="pager-strip">
									<b class="showed">none</b> of <b class="total">none</b>
								</div>
								<div id="pager">
								</div>
								<div id="DescriptorsBrowser">
								</div>	
							</div>
						</div>	
					</div>							
				</td>
			</tr>		
			<tr>
				<td align="right">
					<a href="javascript:window.closeTab()">close</a>
				</td>
			</tr>
		</table>
	</xsl:template>
	
</xsl:stylesheet>