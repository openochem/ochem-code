<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	
	<xsl:include href="../../helper.xslt" />
	
	<xsl:template name="content">
	
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script language="javascript" src="js/commons/dynamic-select.js"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<style type="text/css">
			#Browser IMG {border: 1px solid gray; margin: 4px;}	
			#Browser DIV {padding: 2px;}
			A IMG {margin-right: 5px;}
		</style>
		<table height="100%" width="100%">
			<tr><td class="itunes-up">
				<img src="img/icons/compare.png"/>
				<h1><span class="setcompare">SetCompare</span>: Selected molecules</h1>
				</td></tr>
			<tr>
				<td class="itunes-right">
					<div id="errors" class="invisible">
						This dialog displays the failed molecules.
					</div>
					<input type="hidden" filter="1" name="name"/>
					<input type="hidden" filter="1" name="set"/>
					<input type="hidden" filter="1" name="errors"/>
					
					
					<a action="openEPBrowser" class="rounded-button"><img src="img/icons/compounds16.png"/>Browse the experimental data for these molecules</a>
					<br/><br/>
					
					<div class="pager-strip">
						<span><b class="showed">none</b> of <b class="total">none</b></span>
						<div class="pgr">
						</div>
					</div>
					<div id="Browser">
					</div>
					<div class="pager-strip">
						<span><b class="showed">none</b> of <b class="total">none</b></span>
						<div class="pgr">
						</div>
					</div>	
				</td>
			</tr>
		</table>
		
		<script language="javascript">
			include.plugins('view');
			
			var molBrowser = new Browser();
			molBrowser.url = "setcomparison/listMolecules.do";
			molBrowser.actionURL = "setcomparison/*.do";
			molBrowser.itemElement="molecule";
			molBrowser.pager.selectors.pager = ".pgr";
			molBrowser.itemTemplate = 'js/templates/setcompare-molecule.ejs';
			
			$(document).ready(function(){
				molBrowser.initialize();
			});
			
			if (getParams["errors"])
				$("#errors").removeClass("invisible");
			
			molBrowser.listenEvent("items_loaded", function(){
				$("#Browser div").css({width: "auto", float: "left"});
			});
			
			molBrowser.onOpenepbrowserSuccess = function(entity)
			{
				window.openTab("Experimental data browser", entity.message.message);
			}
		</script>
		
		</xsl:template>
</xsl:stylesheet>
