<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />
	
	<xsl:template name="content">
		<title>Molecule browser</title>
		<style type="text/css">
			.property {font-size: 13px; width: 50%; float: right; clear: right; text-align: left; align: left;}
			.molecule {font-size: 13px; width: 50%; float: left;}
			.fragmented {color: green;}
			img {vertical-align: middle;}
		</style>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script language="javascript" src="js/commons/dynamic-select.js"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript" src="js/browsers/molecule-browser.js"></script>
		<script language="javascript">
			include.plugins('view');
			var moleculeBrowser = new MoleculeBrowser();
			$(document).ready(function() {
				moleculeBrowser.options.highlight_row = "no";
				moleculeBrowser.initialize();
			});
		</script>
		<table width="100%" height="100%" cellspacing="0">
		<tr>
			<td class="itunes-up" colspan="2">
			<img src="img/icons/compounds.png" id="page-icon"/>
				<xsl:call-template name="area-of-interest"/>
				<h1>Molecules browser</h1>
				search compounds
			</td>
		</tr>
		<tr>
			<td class="itunes-left browserScope" id="EPFilters">
				<div id="MolFiltersDiv">
					<div class="openable opened">
						<h1>Molecule</h1>
						<br/>
						<div class="openable-content">
						Molecule ID <input type="text" name="id" prompt="" value="" filter="1"/><br/><br/>
						Name / InchiKey<input type="text" name="name" prompt="" value="" filter="1"/><br/><br/>
						
						<br/>
						<input type="checkbox" name="with-records" value="true" filter="1" class="checkbox"/>With records
						</div>
						<br/>
						<a href="javascript:void()" class="fb-button" onclick="moleculeBrowser.request(true); return false;">refresh</a>
						<a action="reset" class="fb-button">reset</a>
					</div>
				</div>
			</td>
			<td class="itunes-right">
				<div class="browserScope">
				<div class="upper-command-panel">
						<label>Tags</label>
						<a action="addtagselect" title="Tag selected molecules" ><img src="img/icons/tag.png"/></a>
						<a action="removetagselect" title="Remove tag from selected molecules"><img src="img/icons/tag-orange.png"/></a>
						<label>Records</label>
						<a action="addselect" title="Select all molecules" class="invisible"><img src="img/icons/select_all.gif"/></a>
						<a action="selectpage" title="Select molecules on current page"><img src="img/icons/select_page.gif"/></a>
						<a action="removeselect" title="Unselect all molecules"><img src="img/icons/unselect_all.gif"/></a>
					</div>
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
				</div>
			</td>
		</tr>
	</table>
	
	<div id="waitingDialog"> 
	    <div class="hd">Please wait</div> 
	    <div class="bd" style="text-align: center;"> 
	        Please wait until action is completed.<br/>
	        It may take a while.<br/>
	        <img src="img/roller_small.gif"/> 
	    </div> 
	</div>
	</xsl:template>
</xsl:stylesheet>