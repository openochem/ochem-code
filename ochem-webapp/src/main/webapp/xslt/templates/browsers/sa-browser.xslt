<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />

	<xsl:template name="content">
		<style type="text/css">
			TABLE.torefactor TD, TABLE.torefactor TH {background-color: #FFD; padding: 10px; text-align: left; border-bottom: 1px solid white;}
			TABLE.torefactor TH {font-style: bold; border-bottom: 1px solid black;}
			TR.ready TD {background-color: #c0e793;}
			.smart {font-size: 130%; }
			.comment {font-size: 80%;margin-bottom: 4px;}
			#top-commands A {margin-right: 20px;}
			
			.compact-alert {border: 1px solid #444; float: left; width: 220px; text-align: center; height: 220px; overflow: hidden; margin: 20px;}
			
			
			a#move_up {
			   position: fixed;
			   top: 10px;
			   display: none;
			
			   width: 90px;
			   height: 21px; 
			   text-align: center;
			   font: 12px Verdana;
			   text-decoration: none;
			   color: #2b587a;
			   background: #e1e7ed;
			   padding-top:5px;
			   opacity:0.9;
			   filter: alpha(opacity=90);
			}
			a#move_up:hover {
			   color: #fff;
			   background: #597da3;
			}
			
			.right-label A {color: #888;}
			
			.alert-name {font-size: 100%;}
		</style>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<title>Model templates</title>
		<a id="move_up" href="#">Go up</a>
		<table height="100%" width="100%">
			<tr><td class="itunes-up" colspan="2">
				<img src="img/icons/sa.png"/><h1><span class="toxalerts">ToxAlerts</span>: Structural alerts browser</h1>
				Here you can browser structural alerts for various toxicological endpoints
				
				</td></tr>
			<tr>
				<td class="itunes-left">
				<h1>FILTERS</h1>
				<br/>
				Article:<br/>
				<select name="article" filter="1">
					<option value="">All articles</option>
					<xsl:for-each select="//available-alerts/articles">
						<option value="{@id}"><xsl:value-of select="publication-date/year"/>&#160;<xsl:value-of select="authors/author[1]/LastName"/></option>
					</xsl:for-each>
				</select>
				<br/><br/>
				Endpoint / Filter type:<br/>
				<select name="property" filter="1">
					<option value="">All endpoints</option>
					<xsl:for-each select="//available-alerts/endpoints">
						<option value="{@id}"><xsl:value-of select="@name"/></option>
					</xsl:for-each>
				</select><br/><br/>
				Name / Alert ID:<br/>
				<input type="text" name="alert-name" filter="1" size="40"/>
				<br/>
				<br/>
				<input type="checkbox" filter="1" name="approved-only"/> Show only approved alerts<br/>
				<div id="filter-selected-only" class="invisible">
					<input type="checkbox" filter="1" name="selected-only"/> Show only <b><span id="selection-size"></span> selected alerts</b>
				</div>
				</td>
				<td class="itunes-right">
					<div style="float: right">
						<a action="compact" title="Switch to compact view"><img src="img/icons/compact-view.png"/></a>
						<xsl:if test="//param[@key='moderator'] = 'true'">
							<a tab="Export alerts" href="alerts/export.do" title="Export all alerts as a CSV file"><img src="img/icons/excel32.gif"/></a>
							<a action="approve" title="Approve the selected records matching current filters"><img src="img/icons/approve32.png"/></a>
							<a action="disapprove" title="Disapprove (unpublish) the selected records matching current filters"><img src="img/icons/disapprove32.png"/></a>
							<input type="checkbox" filter="1" name="awaiting-approval"/>&#160;Awaiting approval&#160;&#160;
						</xsl:if>
					</div>
					<div id="top-commands">
						<a href="alerts/upload.do" tab="Structural alerts upload"><img src="img/icons/batch_upload.gif"/> Upload new alerts</a>
						<a href="alerts/screen.do" tab="Screen compounds against alerts"><img src="img/icons/screen.png"/> Screen compounds</a>
						<a action="addselect" title="Select all records matching current filters"><img src="img/icons/select_all.gif"/></a>
						<a action="removeselect" title="Unselect all records matching current filters"><img src="img/icons/unselect_all.gif"/></a>
					</div>
					<div class="pager-strip">
						<span><b class="showed">none</b> of <b class="total">none</b></span>
						<div class="pgr">
						</div>
					</div>
					<div id="Browser"></div>
					<div class="pager-strip">
						<span><b class="showed">none</b> of <b class="total">none</b></span>
						<div class="pgr">
						</div>
					</div>
				</td>
			</tr>
		</table>
		
		<div id="imageDialog"> 
		    <div class="hd">Upload image</div> 
		    <div class="bd"> 
		        <form method="post" enctype="multipart/form-data" action="alerts/uploadImage.do" target="upload-frame">
		        	<input type="hidden" name="id" id="alert-id" value=""/>
		            <input type="file" name="image"/>
		            <input type="submit" value="Upload image"/>
		        </form> 
		    </div> 
		</div>
		
		<iframe width="1" height="1" name="upload-frame" id="upload-frame"/>
		
		<div id="basicmenu" class="yuimenu browserScope">
		    <div class="bd">
		        <ul class="first-of-type">
			        
			            <li class="yuimenuitem no-trash">
			                <a action="screen" class="yuimenuitemlabel" href="#screen">
			                    Screen a compound set for this alert
			                </a>
			            </li>
			            <li class="yuimenuitem no-trash">
			                <a action="edit" class="yuimenuitemlabel" href="#edit">
			                    Edit this alert
			                </a>
			            </li>
			            <xsl:if test="//param[@key='moderator'] = 'true'">
			            	<li class="yuimenuitem no-trash">
				                <a action="approve" class="yuimenuitemlabel" href="#screen">
				                    Approve this alert
				                </a>
				            </li>	
				            <li class="yuimenuitem no-trash">
				                <a action="disapprove" class="yuimenuitemlabel" href="#screen">
				                    Unapprove this alert
				                </a>
				            </li> 
			            </xsl:if>
			        	
			      </ul>         
			    </div>
			</div>
		
		<script language="javascript">
			
			include.plugins('view');
			var AlertBrowser = function()
			{
				var self = this;
				this.controller = "alerts";
				Browser.call(this);
				this.compact = getParams["compact"] ? true : false;
				this.itemTemplate = this.compact ? "js/templates/sa-record-compact.ejs" : "js/templates/sa-record.ejs";
				this.itemElement = "substructure-alert";
				this.pager.selectors.pager = ".pgr";
				
				this.cleanBasketDialog = new YAHOO.widget.Dialog("imageDialog", { 
			    	width:"325px", 
					fixedcenter:true, 
					modal:true, 
					visible:false
			    });
			    
			    this.doCompact = function()
			    {
			    	self.compact = !self.compact;
			    	self.itemTemplate = this.compact ? "js/templates/sa-record-compact.ejs" : "js/templates/sa-record.ejs";
			    	self.view = false;
			    	self.request(false);	
			    }
			    
			    this.onApproveSuccess = this.onDisapproveSuccess = function()
			    {
			    	self.request(false);
			    }
			    
			    this.onDeleteSuccess = function()
			    {
			    	self.deleteRecord();
			    }
			    
			    this.cleanBasketDialog.render();
			}
			
			var browser = new AlertBrowser();
			var self = browser;
			browser.listenEvent("items_loaded", function()
			{
				// Image upload action
				$("#Browser .block-image IMG").click(function(){
					browser.img = $(this).get(0);
					$("#alert-id").val($(this).parents("[rec-id]").attr("rec-id"));
					browser.cleanBasketDialog.show();
				});
			
				// Expand long SMARTS action	
				$("span.smarts").each(function(){
					var val = $(this).html();
					if (val.length >= 200)
					{
						$(this).html(val.substring(0, 80) + "   ");
						this.fullSMARTS = val;
						var extendLink = $('<a href="#"/>');
						extendLink.html("show full SMARTS");
						$(this).append(extendLink);
						extendLink.click(function(){
							$(this).parent().html($(this).parent().get(0).fullSMARTS);	
							return false;
						});
					}
				});
				
				//if (self.compact)
				//	$("div[rec-id]").each(function(){
				//		$(this).css({position: "relative", left: "" + Math.floor((Math.random()*200)-100) + "px", top: "" + Math.floor((Math.random()*200)-100) + "px"})
				//	});
				
				$("#filter-selected-only").setClass("invisible", browser.selectionSize == 0);
				$("#selection-size").html("" + browser.selectionSize);
			});
			
			browser.onAddselectSuccess = 
			browser.onRemoveselectSuccess = function()
			{
				// Unfreeze interface and reload page
				//this.waitingDialog.cancel();
				self.request();
			}
			
			browser.doScreen = function(link)
			{
				openTab("Screen for an alert", "alerts/screen.do?alert-id=" + this.currentEntity.id);
			}
			
			browser.doRecordmenu = function(link)
			{
				this.recordMenu.cfg.setProperty("context", [$(link).attr('id'), "tl", "tr"]); 
				this.recordMenu.show();
			}
			
			browser.parentSetPosition = browser.setPosition;
			browser.setPosition = function(link)
			{
				// Context menu is outside of the record context
				// So dont change position while clicking menu item link
				// Otherwise you lose focus
				if (!$(link).hasClass("yuimenuitemlabel"))
					this.parentSetPosition(link);
			}
			
			browser.beforeToggleselect = function()
			{
				// Changing checkbox status even before sending request to server
				var newValue = (this.currentEntity.selected == "true") ? "false" : "true";
				this.currentEntity.selected = newValue;
				this.currentBlock.find('img[name="checked"]').attr('src', 
					(newValue == "true") ? "img/icons/checked.gif" : "img/icons/unchecked.gif");
				return true;
			}
			
			browser.onToggleselectSuccess = function(xml)
			{
			}
			
			$(document).ready(function(){
				browser.initialize();
				browser.recordMenu = new YAHOO.widget.Menu("basicmenu");
				browser.recordMenu.render();
			});
			
			$("#upload-frame").load(function(){
			
				if (browser.img)
					reloadImg(browser.img);	
			});
			
			function reloadImg(obj) 
			{
			   var src = obj.src;
			   var pos = src.indexOf('?');
			   if (pos >= 0) {
			      src = src.substr(0, pos);
			   }
			   var date = new Date();
			   obj.src = src + '?v=' + date.getTime();
			   return false;
			}
			
			// "Go up" button. Might be used in the other browsers too
			$(function () 
			{ 
			    $(window).scroll(function () {
			        if ($(this).scrollTop() > 400) $('a#move_up').fadeIn(); 
			        else                           $('a#move_up').fadeOut(400); 
			    });
			    $('a#move_up').click(function () 
			    {
			        $('body,html').animate({ 
			            scrollTop: 0
			        }, 800); 
			        return false;
			    });
			}); 
			
			
			
		</script>
	</xsl:template>
</xsl:stylesheet>