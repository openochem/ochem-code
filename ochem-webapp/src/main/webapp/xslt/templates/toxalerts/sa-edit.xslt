<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />
	
	<xsl:template name="content">
		<style type="text/css">
			.form-table TD {padding: 6px 10px 6px 3px; vertical-align: top;}
			.form-table TEXTAREA {width: 400px; height: 200px;}
			.form-table INPUT[type=text] {width: 400px; }
		</style>
		<link rel="stylesheet" type="text/css" href="css/smoothness/jquery-ui-1.9.1.custom.css" />
		<link rel="stylesheet" type="text/css" href="css/tagit-simple-blue.css" />
		<script language="javascript" src="js/lib/jquery-ui-1.9.1.custom.min.js"></script>
		<script language="javascript" src="js/lib/tagit.js"></script>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<table height="100%" width="100%">
			<tr><td class="itunes-up">
				<img src="img/icons/sa.png"/>
				<h1><span class="toxalerts">ToxAlerts</span>: Modify alert details</h1>
				
				</td></tr>
			<tr>
				<td class="itunes-right">
					<form action="alerts/save.do" method="post" class="formscope">
					<input type="hidden" name="id" value="{substructure-alert/@id}" send="1"/>
					<table class="form-table">
						<tr>
							<td>SMARTS<br/><small>can contain <a href="alerts/variables.do" tab="SMARTS substitution variables">substitution variables</a></small></td>
							<td><input type="text" name="smarts" value="{substructure-alert/smart}" send="1"/><br/>
							<small>you may debug your SMARTS using our simple <a href="alerts/match.do" tab="SMARTS matcher">SMARTS matcher</a></small><br/><br/></td>
						</tr>
						<tr>
							<td>Name</td><td><input type="text" name="name" value="{substructure-alert/name}" send="1"/></td>
						</tr>
						<tr>
							<td>Description (as provided by the author)</td><td><input type="text" name="description" value="{substructure-alert/description}" send="1"/></td>
						</tr>
						<tr>
							<td>SMARTS description (additional clarification)</td><td><input type="text" name="smarts-description" value="{substructure-alert/smartsDescription}" send="1"/></td>
						</tr>
						
						<tr>
							<td>Publication
								<br/><small>a publication where this structural alert has been introduced</small>
							</td>
							<td>
								<a action="selectarticle" bindto="article"><xsl:value-of select="substructure-alert/article/@title"/></a>
								<input type="hidden" name="article" send="1" value="{substructure-alert/article/@id}"/>
							</td>
						</tr>
						<tr>
							<td>Endpoint/Class:<br/>
								<small>e.g., aqueous toxicity, carcinogenicity, etc.</small>
							</td>
							<td>
								<a action="selectproperty" bindto="property"><xsl:value-of select="substructure-alert/property/@name"/></a>
								<input type="hidden" name="property" send="1" value="{substructure-alert/property/@id}"/>
							</td>
						</tr>
						
						<tr>
							<td>Author's comment</td><td>
								<textarea name="comment" value="{substructure-alert/comment}" send="1"><xsl:value-of select="substructure-alert/comment"/></textarea>
							</td>
						</tr>
						
						<tr>
							<td>Parent patterns or groups:<br/>
							</td>
							<td>
								<ul type="text" id="parents" style="width: 300px; height: 100px;" class="tagit">
								</ul>
							</td>
						</tr>
					</table>
					<div class="popup-footer">
						<a action="edit">save</a><a href="javascript:window.closeTab();">cancel</a>
					</div>
					</form>
				</td>
			</tr>
		</table>
		<script language="javascript">
			function AlertForm()
			{
				var self = this;
				this.scope = ".formscope";
				EditForm.call(this);
				this.actionURL = "alerts/save.do";
				
				this.beforeEdit = function()
				{
					var tags = $('#parents').tagit("tags");
					var vals = new Array();
					for (var n in tags)
						vals.push(tags[n].value);
					self.setValue("parents", vals.join(","));
					
					return true;
				}
				
				this.doSelectarticle = function()
				{
					var win = openTab("Select a publication", webRoot+"article/show.do?render-mode=popup");
					win.callback = function(obj)
					{
						self.setValue("article", obj.id, obj.title);
						win.closeTab();
					}	
				}
				
				this.doSelectproperty = function()
				{
					var win = openTab("Select an endpoint or a property", webRoot+"properties/show.do?render-mode=popup");
					win.callback = function(obj)
					{
						self.setValue("property", obj.id, obj.name);
						win.closeTab();
					}	
				}
			}
			
			var initialTags = [];
			
			<xsl:for-each select="/model/substructure-alert/parents/parent">
				initialTags.push({label: '<xsl:value-of select="name"/>', value: '<xsl:value-of select="@id"/>'});
			</xsl:for-each>
			
			
			
			$(function(){
				$('#parents').tagit({tagSource: "/alerts/autocomplete.do?out=json&amp;list-only=1", 
					select: true, 
					sortable: true, 
					allowNewTags: false, 
					initialTags: initialTags,
					triggerKeys: ['comma', 'enter', 'tab'],
					placeholderText: "Start typing an alert or group name to add it to the list"
				});
			});
			
			var form = new AlertForm();
		</script>
		</xsl:template>
</xsl:stylesheet>
