<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />	
	<xsl:template name="content">
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript" src="js/commons/mol-editors.js"></script>
		<script type="text/javascript" src="js/blocks/molecule.js?ver=2"></script>
	<script language="javascript">
		include.plugins('view'); 
		var editor = new MoleculeEditor();
		editor.actionURL="molecule/action.do";
		$(document).ready(
		function() {
		    tabView = new YAHOO.widget.TabView("demo");
			editor.initialize();
			editor.tabView = tabView;
			editor.callAction('Get');
			$(".tabs").height($("#bottom").height() - $(".ui-nav").height() - 20);
			$(".tabs").css({overflow: "auto"});
		}
		);
	</script>

	<style type="text/css">
		.conditions {color: green; font-size: 8pt; float: right; clear: right; width: 400px; margin-left: -50px; text-align: right;}
		.article-data {font-size: 8pt; margin-top: 7px; margin-left: 0px; margin-bottom: 7px;}
		img {vertical-align: middle;}
		BODY IFRAME {border: 0px;}
		#bottom {border: 1px solid #999; background-color: #F0F0F0; clear: both; padding: 10px 10px 10px 10px; font-size: 100%;}
		#bottom H1 {font-size: 100%; margin-bottom: 3px;}
		.depiction {overflow:auto;}
		.floating {float: left; clear:none; border: 1px solid black; padding: 2px 2px 2px 2px; margin: 1px 1px 1px 1px;}
		INPUT {border: 1px solid black; padding: 2px 2px 2px 2px; margin: 1px 1px 1px 1px;}
		textarea#styled { width: 290px;	height: 60px; border: 1px solid black; overflow: auto;}
		.fb-button {position: relative; top: 2px;}
		
	</style>
	<title>Molecule editor</title>
	<h1 class="popup-header"><doc term="Molecule+editor">Molecule editor</doc>
		<i>Draw a molecular structure or import it from a variety of formats (SMILES, SDF, MOL)</i></h1>
		<iframe name="iframe" id="iframe" width="1" height="1" onload="editor.iframeLoaded();">
		</iframe>

	<table style="float: left;">
	<tr>
		<td valign="top" width="1px" height="1px">
			<iframe src="../editor.html" id="sketch" class="sketcher-frame" width="380px" height="340px"></iframe>
		</td>
		<td id="bottom" style="vertical-align: top" rowspan="2">
		<form enctype="multipart/form-data" name="molViewer" action="molecule/action.do?action=Submit&amp;out=json" method="post" target="iframe">
		<div class="mainactions">
		<input type="hidden" name="format" value="file"/>
		
			<p>Paste your molecule (SDF, MOL2, SMILES, name)<br/>
			Examples: "Aspirine", "CC=CCCC", "LSD":<br/>
				<textarea name="smiles" id="styled" style="width: 300px;"/>
				<br/>
				<a action="Smiles" class="fb-button">Load pasted molecule</a>
			</p>
			<br/><br/>
			<p>Upload from an SD-file<br/>
				<input type="file" name="file" style="width: 300px;" />
				<br/>
				<a href="#"  class="fb-button" onclick="$(this).closest('form').submit(); return false;">Upload file</a>
			</p>
			<br/><br/>
			<p>Fetch by name from PubChem database<br/>
				<input name="namesearch" type="text" style="width: 300px;" />
				<br/>
				<a action="Searcher"  class="fb-button" search="pubchem-search">Search in PubChem</a>
			</p>
		</div>
		</form>	
		<br/><br/>
		</td>		
		
	</tr>
	<tr>
		<td style="text-align: right; padding: 5px;">
			<div id="bottom2" class="mainactions">
				<input type="button" action="Submit" value="Submit the molecule" class="fancy-button"/>
				<input type="button" onclick="window.closeTab();" value="cancel"/>
			</div>
		</td>
	</tr>
	<tr>
		<td valign="top" colspan="2">
			
		</td>	
	</tr>
	</table>
	
	<div id="demo" class="yui-navset" style="float: left; width: 450px;"> 
	    		<ul class="yui-nav"> 
	        		<li class="selected"><a href="#tab1"><em>Names</em></a></li> 
	        		<li><a href="#tab2"><em>Alternate Depictions</em></a></li> 
	        		<li><a href="#tab3"><em>Properties</em></a></li> 
	    		</ul> 
	    		<div id="add-inf" class="yui-content depiction">
	    			<div id="MolNames" class="tabs">
	    				<div id="MolNameBrowser">
	    				</div>
	    			</div>
	    			<div id="MolDepictions" class="tabs">
	    			</div>
	    			<div id="MolProperties" class="tabs">
						<div class="pager-strip">
							<b class="showed">none</b>
							of
							<b class="total">none</b>
						</div>
						<div id="pager"></div>
						<div id="MolPropertiesBrowser">
	    				</div>
	    			</div>    			    		
	    			    			    		
				</div>
			</div>
	
	</xsl:template>
</xsl:stylesheet>