<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />
	
	<xsl:template name="content">
		<title>Molecule search browser</title>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript" src="js/commons/midnighter.js"></script>
		<script language="javascript" src="js/browsers/search-browser.js"></script>
		<script language="javascript"> 
			include.plugins('view');
			var sampleBrowser = new SearchBrowser(); 
			sampleBrowser.url = "ncbisearch/list.do";
			sampleBrowser.actionURL = "ncbisearch/action.do";
			sampleBrowser.itemTemplate = "js/templates/searchmol.ejs";
			sampleBrowser.itemElement = "molecule"; 
			
			sampleBrowser.onSearchSuccess = function(xml) 
			{
				this.request(); 
			}
			
			sampleBrowser.getActionQuery = function(action)
			{
				if (action == "Search")
				{
					$("#Browser").html('<img src="img/long_green.gif"/>');
					return  "search-string=" + URLEncode(($("input[name='namesearch']").val()).trim()) + 
							"&amp;search-number=" + $("input[name='searchnumber']").val() + 
							"&amp;database=" + $("select[name='db']").val()
							;
				}
			}			
			
			sampleBrowser.onAddSuccess = sampleBrowser.onDeleteSuccess = sampleBrowser.onRenameSuccess = function() {this.request()};
			
			$(document).ready( function()
			{
				sampleBrowser.initialize();
			 	var uri = URLDecode(getParams['search-string']);
			 	$("input[name='namesearchqspr']").attr("value", uri);
			 	$("input[name='searchnumberqspr']").attr("value", "1");
			 	$("input[name='namesearch']").attr("value", uri);
			 	$("input[name='searchnumber']").attr("value", "1");
			});
		</script>
		<table height="100%" width="100%">
			<tr>
				<td class="itunes-up">
					<h1>Molecules search browser</h1>
					Search for molecules by name in Pubchem database
				</td>
			</tr>
			<tr>
				<td class="itunes-right">
					<div>
						Fetch a molecule structure by name from 
						<select name="db" send="1">
							<option value="PubChem">PubChem</option>
							<option value="QSPR">QSPR</option>
						</select>
					 	database
						<br/>
						Name: <input name="namesearch" type="text" size="30" maxlength="75" value=""/>
						No. results <input name="searchnumber" type="text" size="2" maxlength="3" value=""/>
						<a action="Search" search="pubchem">[search]</a>
						<img src="1.gif" width="1" height="1" type="button"/>
						<input type="image" class="invisible" value="submitname" src="submit-button.gif" width="1" height="1" border="0" alt="SUBMIT!" name="image"/>
					</div>
					<div class="pager-strip">
						<b class="showed">none</b> of <b class="total">none</b>
					</div>
					<div id="pager"></div>
					<div id="Browser"></div>
					
					<!--div>Fetch a molecule structure by name from PubChem database<br/>
						Name: <input name="namesearch" type="text" size="20" maxlength="50" value=""/>
						No. results <input name="searchnumber" type="text" size="2" maxlength="3" value=""/>
						<a action="Search" >[search]</a>
						<img src="1.gif" width="1" height="1" type="button"/>
						<input type="image" class="invisible" value="submitname" src="submit-button.gif" width="1" height="1" border="0" alt="SUBMIT!" name="image"/>
					</div>
					<div class="pager-strip">
						<b class="showed">none</b> of <b class="total">none</b>
					</div>
					<div id="pager"></div>
					<div id="Browser"></div-->
				</td>
			</tr>
		</table>
	</xsl:template>
</xsl:stylesheet>