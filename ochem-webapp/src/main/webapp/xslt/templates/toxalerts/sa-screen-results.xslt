<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	
	<xsl:include href="../../helper.xslt" />
	<xsl:template name="content">
	
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<style type="text/css">
			.itunes-left table TD {padding-bottom: 5px; padding-right: 5px; font-size: 12pt;}
			#EPFilters {position: relative;}
			#top-commands A {margin-right: 20px;}
			#top-commands IMG {height: 24px; width: 24px;}
		</style>
		
		<table height="100%" width="100%">
			<tr><td class="itunes-up" colspan="2">
				<img src="img/icons/sa.png"/>
				<h1><span class="toxalerts">ToxAlerts</span>: Screening results</h1>
				The compounds that matched any alerts grouped by endpoints, publications and by alerts themselves
				</td></tr>
			<tr>
				<td class="itunes-left" id="EPFilters">
					<div id="EPFiltersDiv">
					<a action="hidefilters" class="close-button" title="Hide filters panel">x</a>
					<table>
						<tr><td colspan="3" id="endpoints-span"><h1>Endpoints</h1><br/></td></tr>
						<xsl:for-each select="//others/property">
							<tr><td colspan="2"><nobr><input type="radio" name="alert" filter="1" value="e{@id}"/>&#160;&#160;<xsl:value-of select="@name"/></nobr></td><td align="right"><nobr><xsl:value-of select="@timesUsed"/> compounds</nobr></td></tr>
						</xsl:for-each>
						<tr><td colspan="3" id="publications-span"><h1>Publications</h1><br/></td></tr>
						<xsl:for-each select="//others/article">
							<tr><td colspan="2"><input type="radio" name="alert" filter="1" value="a{@id}"/>&#160;&#160;<a href="article/profile.do?id={@id}" tab="Article profile" title="{@title}"><xsl:value-of select="publication-date/year"/>&#160;<xsl:value-of select="authors/author[1]/LastName"/></a></td><td align="right"><nobr><xsl:value-of select="@basketCount"/> compounds</nobr></td></tr>
						</xsl:for-each>
						<tr><td colspan="3" id="alerts-span"><h1>Detected alerts</h1><br/></td></tr>
						<xsl:for-each select="//others/substructure-alert">
							<tr>
								<td><nobr><input type="radio" name="alert" filter="1" value="{@id}"/>&#160;&#160;<a href="alerts/show.do?id={@id}" tab="Alert details"><xsl:value-of select="name"/></a></nobr></td>
								<td style="text-align: left; padding-right: 45px;"><xsl:value-of select="article/publication-date/year"/>&#160;<xsl:value-of select="article/authors/author[1]/LastName"/></td>
								<td align="right"><nobr><xsl:value-of select="countOccurences"/> compounds</nobr></td>
								</tr>
						</xsl:for-each>
					</table>
					</div>
				</td>
				<td class="itunes-right">
					<div style="position: relative" id="top-commands">
					<a action="showfilters" class="show-button invisible" title="Show filters panel">&lt;</a>
					<a action="openbrowser">View records for the filtered compounds</a>
					<a action="addtag"><img src="img/icons/tag.png"/> Tag the <span class="total"></span> filtered molecules</a>
					<a tab="Export results" href="alerts/exportResults.do"><img src="img/icons/xls.gif"/> Export the screening results</a>
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
				</td>
			</tr>
		</table>
		
		<div id="molViewDialog"> 
		    <div class="hd">Molecule structure</div> 
		    <div class="bd">
		    	<a href="javascript:void()" width="300" height="300" onclick="compoundBrowser.molViewDialogClose(); return false;"><img src=""/></a> 
		    </div> 
		</div>
		
		<script language="javascript">
			include.plugins('view');
			var ScreenResultsBrowser = function()
			{
				var self = this;
				this.url = "alerts/screenlist.do";
				Browser.call(this);
				this.itemTemplate = "js/templates/sa-screen-result.ejs";
				this.itemElement = "row";
				this.pager.selectors.pager = ".pgr";
				
				this.doOpenbrowser = function()
				{
					openTab("Properties of compounds with alerts", "alerts/openPropertyBrowser.do?" + self.filters.getQueryString());
				}
				
				this.doZoom = function(link)
				{
					var size = Math.round(Math.min(window.innerHeight, window.innerWidth) * 0.9);
					var id = this.currentBlock.find('.block-image IMG').attr("mp2");
					var img = $(this.molViewDialog.body).find("IMG");
					img.attr("src","");
					img.attr("src","depiction.jsp?w="+size+"&amp;h="+size+"&amp;mp2="+id);
					img.attr("width",size);
					img.attr("height",size);
				    this.molViewDialog.render();
					this.molViewDialog.show();
				}
				
				this.getRecordIdentifier = function(obj)
				{
					// Will return the compound ID
					return obj.int;
				}
				
				this.doShowfilters = function()
				{
					$("#EPFilters").removeClass("invisible");
					$(".show-button").addClass("invisible");
				}
				
				this.doHidefilters = function()
				{
					$("#EPFilters").addClass("invisible");
					$(".show-button").removeClass("invisible");
				}
				
				this.doAddtag = function()
				{
					var selectWin = openTab("Select a tag", webRoot + "tags/show.do?render-mode=popup");
					selectWin.callback = function(selectedItem)
					{  
						self.ajax.send({
							url:  webRoot + "alerts/addTag.do",
							data: "tag=" + selectedItem.id + "&amp;" + self.filters.getQueryString(),
							success: function()
							{
								window.alert("Tag added successfully!");
								selectWin.closeTab();
							}
						});
					}
				}
				
				this.molViewDialog = new YAHOO.widget.Dialog("molViewDialog", { 
					fixedcenter:true, 
					modal:true, 
					visible:false,
					postmethod: "none" 
			    });
			}
			
			resBrowser = new ScreenResultsBrowser();
			$(document).ready(function() {resBrowser.initialize();});
		</script>
		</xsl:template>
</xsl:stylesheet>
