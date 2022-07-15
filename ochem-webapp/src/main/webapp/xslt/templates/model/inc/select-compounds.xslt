<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:template name="select-compounds">
		<xsl:param name="scope" />
		<xsl:param name="brief" />
		<script language="javascript">
			include.plugins('view');
			
			function getDepictionUrl(id)
			{
				return "depiction.jsp?w=250&amp;h=250&amp;id="+id;
			}
			
			function MoleculeSearchBrowser()
			{
				this.controller = "moleculesearch";
				Browser.call(this);
				var self = this;
				this.itemElement = "molecule";
				this.itemTemplate = 'js/templates/moleculesearch.ejs';
				this.itemClass = "molsearchitem";
				this.oldFilters = "";
				this.conditionalRequest = function()
				{
					if (!self.filters.getValue("name"))
					{
						window.alert("Please, enter a molecule name");
						return;
					}
					var newFilters = self.filters.getQueryString();
					if (newFilters != self.oldFilters)
					{
						self.oldFilters = newFilters;
						self.request();
					}
				}
				
				this.onItemsLoaded = function()
				{
					if (self.pager.totalNum == 0)
						self.mainDiv.html(&apos;<div class="search-message">No molecules found.<br/>Maybe there is a typo in the name?</div>&apos;);
					if (self.pager.totalNum == 1)
					{
						self.currentEntity = {id: self.mainDiv.find("[rec-id]").attr("rec-id")};
						self.doItemselect();
					}
				}
				
				this.listenEvent("items_loaded", this.onItemsLoaded);
								
			}
		
			var cpdsForm; 
			(function()
			{
				var form = new AjaxForm();
				var self = form;
				cpdsForm = form;
				var scope = "<xsl:value-of select="$scope"/>";
				if (scope != "")
					form.scope = "#select-compounds-" + scope;
				
				var getByName = function(name)
				{
					return $("[name=" + scope + name + "]");
				}
				
				var checkOption = function(value)
				{
					getByName("data").filter("[value=" + value + "]").attr("checked", "checked");
				}
				
				var getOption = function() {
					return getByName("data").filter(":checked").val();
				}
				
				form.checkOption = function(value) {
					checkOption(value);
				}
	
				form.doSelectset = function()
				{
					var propWin = openTab("Select a set", webRoot+"basket/show.do?render-mode=popup");
					propWin.callback = function(basket)
					{
						self.setValue(scope + "validationsetname", basket.id, basket.name);
						checkOption("validationset");
						propWin.closeTab();
					}
				}
				
				form.doSelecttag = function()
				{
					var propWin = openTab("Select a tag", webRoot+"tags/show.do?render-mode=popup&amp;type=molecule");
					propWin.callback = function(tag)
					{
						self.setValue(scope + "tag", tag.id, tag.name);
						checkOption("tag");
						propWin.closeTab();
					}
				}
				
				form.doSubmit = function()
				{
					var vs = self.getValue(scope + "validationsetname");
					var nm = self.getValue(scope + "n-molecule");
					if (vs == undefined &amp;&amp; vs == "" &amp;&amp; nm == "") 
						return window.alert("Please provide valid set for the model");
			
					document.modeldata.submit();
				}
				
				form.validate = function() {
					var option = getOption();
					console.log(option);
					if (option == "validationset" &amp;&amp; !getByName("validationsetname").val())
						window.alert("You did not select the basket");
					else
						return true;
				}
				
				form.doEditmolecule = function()
				{
					var molWin = openTab("Edit molecule", webRoot+"molecule/show.do?id="+self.getValue(scope + 'n-molecule'));
					molWin.callback = function(newId)
					{
						checkOption("molecule");
						self.setValue(scope + "n-molecule", newId, newId);
						getByName("depiction").attr('src', getDepictionUrl(newId));
						molWin.closeTab();
					}
				}
				
				function fileSelected()
				{
					checkOption("externalfile");
				}
				
				var molBrowser = new MoleculeSearchBrowser();
				molBrowser.scope = molBrowser.filters.scope = 	"#<xsl:value-of select="$scope"/>name";
				molBrowser.container = "<xsl:value-of select="$scope"/>molbrowser";
				molBrowser.doItemselect = function()
				{
					var newId = this.currentEntity.id; 
					self.setValue(scope + "n-molecule", newId, newId);
					getByName("depiction").attr('src', getDepictionUrl(newId));
					getByName("profile-link").removeClass("invisible").attr("href", "molecule/profile.do?depiction=" + newId);
					$("#<xsl:value-of select="$scope"/>molbrowser").html("");
				}
				
				$(document).ready(function(){
					$("#"+scope+"name").keyup(function(e){
						getByName("data").filter("[value=molecule]").attr("checked", "checked");
					});
					
					$("#"+scope+"name").keydown(function(e){
						if (e.which == 13)
						{
							e.stopPropagation();
							e.preventDefault();
							molBrowser.conditionalRequest();
							return false;
						}	
					});
					
					$("#"+scope+"name").parent().find("a").click(function() {
						molBrowser.conditionalRequest();
						return false;
					});
					
										
					$('input[type=file]').on("change", function(){ 
						fileSelected(); 
					});
					
					molBrowser.initialize(true);
				});
			})();
			
			
		</script>
		<style type="text/css">
			INPUT {border: 1px solid black; padding: 2px 2px 2px 2px;}
			.options TD {padding-right: 5px; padding-bottom: 20px;}
			.vtop {vertical-align: top;}
			.molsearchitem
			{
				padding: 5px;
				width: 152px;
				float:left;
			}
			.search-message
			{
				padding-top: 10px;
				color: #555;
				font-style: italic;
				padding-left: 5px;
				padding-right: 5px;
			}
			
			.profile-link {font-size: 9pt; color: #444 !important; display: block;}
			
			.select-compounds SMALL {color: #444;}
		</style>
		<div id="select-compounds-{$scope}" class="select-compounds">
			<input type="hidden" name="compounds-selection-form" value="1"/>
			<input type="hidden" name="{$scope}validationsetname" send="1"/>
			<input type="hidden" name="{$scope}tag" send="1"/>
			<input type="hidden" name="{$scope}tag-title"/>
			<input type="hidden" name="{$scope}n-molecule" send="1" value="23427"/>
			<table class="options">
				<xsl:if test="$brief = ''">
				<tr>
					<td><input type="radio" name="{$scope}data" value="externalfile"/></td>
					<td>Upload compounds from a file<br/> <small>SDF, MOL2, SMILES or an Excel sheet</small></td>
					<td><input type="file" name="{$scope}externalfile" onchange="fileSelected()"/></td>
				</tr>
				<tr class="option-name">			
					<td><input type="radio" name="{$scope}data" value="molecule"/></td>
					<td>Draw Molecule<br/><small>click on depiction to the right to draw</small></td>
					<td>
						<a action="editmolecule"><img name="{$scope}depiction" style="border: 1px solid gray;" src="depiction.jsp?id=23427&amp;w=100&amp;h=100" width="100" height="100" title="Click to draw another molecular structure"/>
						<a href="" class="profile-link invisible" name="{$scope}profile-link" tab="Molecule profile" title="Open the molecule profile for the molecule">[molecule profile]</a>
					</a>
					
					</td>
				</tr>
				<tr>
					<td></td>
					<td>Name/CASRN/SMILES:<br/>
					<small>e.g., "CC=CCC" or "Aspirine"</small></td>
					<td>
						<div class="vtop" id="{$scope}name">
						 <input type="text" id="{$scope}name" name="name" value="" filter="1"/><a href="#">load structure</a>
						<div id="{$scope}molbrowser"></div>
						</div>
					</td>
				</tr>
				</xsl:if>
				<tr>
					<td><input type="radio" name="{$scope}data" value="validationset"/></td>
					<td>Choose a previously prepared set:</td>
					<td><a action="selectset" bindto="{$scope}validationsetname" storein="{$scope}vset-title">[...]</a></td>	
				</tr>
				<tr>
					<td><input type="radio" name="{$scope}data" value="tag"/></td>
					<td>Select molecules by a tag:</td>
					<td><a action="selecttag" bindto="{$scope}tag" storein="{$scope}tag-title">[...]</a></td>	
				</tr>
			</table>
		</div>
	</xsl:template>
	
</xsl:stylesheet>