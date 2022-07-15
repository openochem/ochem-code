<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0"
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

	<xsl:include href="../helper.xslt" />

	<xsl:template name="content">
		<link rel="stylesheet" type="text/css" href="css/buttons.css" />
		<link rel="stylesheet" type="text/css" href="css/smoothness/jquery-ui-1.9.1.custom.css" />
		<script type="text/javascript" src="js/commons/actionable.js?ver=1.7.5" />
		<script type="text/javascript" src="js/commons/browser.js" />
		<script type="text/javascript" src="js/commons/tooltip-menu.js" />
		<script type="text/javascript" src="js/lib/json2.js" />
		<script type="text/javascript" src="js/browsers/rgroup-browser.js?ver=1.0"></script>
		<script language="javascript" src="js/commons/dynamic-select.js"></script>
		<script type="text/javascript" src="js/lib/jquery-ui-1.9.1.custom.min.js"></script>
		
		<style type="text/css">
			BODY {
				background-color: #E5E5E5;
			}			
			
			#sketch {border: 1px solid #CCC;}	
			
			.itunes-right TABLE.main TD {vertical-align: top; padding-right: 10px;}

			
			.rgroupitem
			{
				padding: 5px;
				width: 152px;
				float:left;
			}
			
			.rgroupitem IMG
			{
				border: 1px;
				border-style:solid;
				border-color:#ddd;	
			}
			
			.resizable 
			{
				padding: 5px;
				width: 49%;
				float:left;
			}
			
			.ui-resizable-e { 
			    cursor: e-resize; 
			    width: 4px; 
			    right: -2px; 
			    top: 0; 
			    bottom: 0; 
			    background-color: #ddd;
			}
			
			.layouttable
			{
				height: 100%;
				width: 100%;
				table-layout:fixed; 
				min-width: 1100px;
			}
			
			.filters
			{
				background-color: #eeeeee;
				padding: 10px;
				border: 1px;
				border-style:solid;
				border-color:#ddd;
				margin-bottom: 10px;
				min-height: 28px;
			}
			
			.blank
			{
				min-height: 30px;
				margin-bottom: 10px;
				padding: 10px;
			}
			
			.content-message {
				text-align: center;
				font-size: 100%;
				color: #999;
				padding: 50px 0px;
				border: 1px solid #ddd;
			}
			
			#ttmenu-f, #ttmenu-b
			{
				position:absolute;
				opacity: 0.7;
				padding: 5px;
				font-size: 70%;
			}
			
			.groupname
			{
				text-align: center;
				background-color: #eeeeee;
				border-left: 1px solid #ddd;
				border-right: 1px solid #ddd;
				border-bottom: 1px solid #ddd;
			}
			
		</style>
		<script language="javascript">		
		include.plugins('view');
		
		var brFull = new RGroupBrowser();
		var brBasket = new RGroupBasketBrowser();
		brFull.basketBrowser = brBasket;
		brBasket.fullBrowser = brFull;
		
		$(document).ready(function(){
			brFull.initialize();
			brBasket.initialize();
		    $(".panelFull").resizable(
		    {
		        autoHide: true,
		        handles: 'e',
		        resize: function(e, ui) 
		        {
		        	var minWidth = 300;
		        	
		        	if (ui.element.width() &lt; minWidth)
		        		ui.element.css({width: minWidth});
		            
		            var parent = ui.element.parent();
		            var maxWidth = parent.width() - minWidth - (ui.element.outerWidth() - ui.element.width());
		            
		            if (ui.element.width() &gt; maxWidth)
		        		ui.element.css({width: maxWidth});
		        	
		        	var remainingSpace = parent.width() - ui.element.outerWidth();
		              
	                divTwo = ui.element.next();
	                divTwoWidth = ((remainingSpace - (divTwo.outerWidth() - divTwo.width()))/parent.width()*100 - 1)+"%";	                
	                divTwo.width(divTwoWidth);
		        },
		        stop: function(e, ui) 
		        {
		            var parent = ui.element.parent();
		            ui.element.css(
		            {
		                width: ui.element.width()/parent.width()*100+"%",
		            });
		        }
		    });
		});
		</script>
		
		<title>RGroup Browser</title>
		
		<table width="100%">
			<tr><td class="itunes-up silver"><h1>RGroup Browser</h1>A browser where you can handle and organize your RGroup libraries to be used in MolOptimizer screenings.</td></tr>
			<tr><td class="itunes-right">
				<table class="main layouttable">
					<tr><td>
					
					<div>
					<div class="resizable panelFull">
							<div class="filters">
								Click on the available groups on the left panel to add them to a set selected on the right
								</div>	
							<div class="pager-strip">
								<span><b class="showed">none</b> of <b class="total">none</b></span>
								<div class="pager">
								</div>
							</div>
							
							<div id="BrowserFull">
							</div>
							
							<div class="pager-strip">
								<span><b class="showed">none</b> of <b class="total">none</b></span>
								<div class="pager">
								</div>
							</div>

							<div id="ttmenu-f" class="invisible">
								<a action="bsktadd"><img src="img/icons/add_to_set.gif"/> Add to basket</a>
							</div>
					</div>
					<div class="resizable panelBasket">
							<div class="filters">
								You are currently working with: <select id="basket-menu" name="basketid" filter="1" ><option value="-1">Loading...</option></select>
								<a action="createbasketdialog">[New]</a>
								<a action="deletebasketdialog">[Delete]</a>
								<a action="renamebasketdialog">[Rename]</a>
							</div>					
							<div class="pager-strip">
								<span><b class="showed">none</b> of <b class="total">none</b></span>
								<div class="pager">
								</div>
							</div>
							
							<div id="BrowserBasket">
							</div>
							
							<div class="pager-strip">
								<span><b class="showed">none</b> of <b class="total">none</b></span>
								<div class="pager">
								</div>
							</div>
							
							<div id="ttmenu-b" class="invisible">
								<a action="bsktremove"><img src="img/icons/remove_from_set.gif"/> Remove from basket</a>
							</div>							
					</div>
					</div>
					</td></tr>
				</table>
			</td></tr>
			<div id="createbasket" title="New Basket">
				Basket name:
				<input type="text"/>
			</div>
			<div id="renamebasket" title="Rename Basket">
				Basket name:
				<input type="text"/>
			</div>
			<div id="renamegroup" title="Rename Group">
				Group name:
				<input type="text"/>
			</div>
		</table>
	</xsl:template>
</xsl:stylesheet>