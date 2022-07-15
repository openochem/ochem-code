<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />

	<xsl:template name="content">
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript">
			 <![CDATA[
			include.plugins('view');
			var sampleBrowser = new Browser();
			sampleBrowser.controller = "tags";
			sampleBrowser.itemElement = "tag";
			sampleBrowser.itemTemplate = "js/templates/tag.ejs";
			sampleBrowser.url = webRoot + "tags/list.do";
			sampleBrowser.actionURL = webRoot + "tags/action.do";
			
			sampleBrowser.onDeleteSuccess = function()
			{
				this.deleteRecord();
			}
			
			sampleBrowser.doProfile = function()
			{
				openTab("Tag profile", webRoot + "tags/profile.do?render-mode=popup&id="+this.currentEntity.id+"&type="+this.currentEntity.type);
			}
			
			sampleBrowser.doShow_mol = function()
			{
				openTab("Tag Molecule",webRoot+"molbrowser/show.do?render-mode=popup&type="+this.currentEntity.type+"&tag="+this.currentEntity.id);
			}
			
			sampleBrowser.doExport = function()
			{
				var win = openTab("Export tagged molecules",webRoot+"tags/export.do?id="+this.currentEntity.id);
			}
			
			sampleBrowser.beforeEdit = function()
			{
				$("a#create-new-tag").attr("query", "type=" + $("select[name=type]").val());
				return true;
			}
			
			$(document).ready(function() {
				sampleBrowser.initialize();
				
				if (getParams["force-type"])
					$("select[name=type]").attr("disabled", "disabled");	
			});
			
			
				]]>
		</script>
		<title>Tags browser</title>
		<table height="100%" width="100%">
			<tr><td class="itunes-up" valign="middle">
				<img src="img/icons/tag.png"/>
				<h1>Tags browser</h1>
				Tags are convenient way to provide meaning to your data. You can apply multiple tags to articles and properties.
			</td></tr>
			<tr>
				<td class="itunes-right">
				Show tags for 
				<select name="type" filter="1">
					<option value="molecule">Molecules</option>
					<option value="property">Properties</option>
				</select>
				
				Type part of tag name to filter:
					<input type="text" name="query" filter="1" search="tag-search"/>
					<a href="javascript:sampleBrowser.request(true)">
						[search]
					</a>
					<a action="edit" title="edit" query="type=property" id="create-new-tag">[create new tag]</a>
					<div class="pager-strip">
						<b class="showed">none</b> of <b class="total">none</b>
					</div>
					<div id="pager">
					</div>
					<div id="Browser">
					</div>
					<div class="pager-strip">
						<b class="showed">none</b> of <b class="total">none</b>
					</div>
				</td>
			</tr>
		</table>
	</xsl:template>
</xsl:stylesheet>