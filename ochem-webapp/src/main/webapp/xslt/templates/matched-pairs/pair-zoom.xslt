<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">

	<xsl:template name="pair-zoom">
		<style type="text/css">
			#propTable TD {text-align: center; padding: 2px;}
			#propTable .value {padding: 10px 40px; font-family: Verdana;  font-size: 14pt; white-space: nowrap;}
			#propTable .value.clickable {cursor: pointer;}
			#propTable .value.clickable:hover {background-color: #8D8;}
			.tableSplitter {background-color: gray; height: 2px;}
		</style>
		<script language="javascript">
			$(function() {
				molViewDialog = new YAHOO.widget.Dialog("molViewDialog", { 
					fixedcenter:true, 
					modal:true, 
					visible:false,
					postmethod: "none" ,
					constraintoviewport: true
		    	});
		    	molViewDialog.render();
			});
			
			function attachZoomingClickHandlers(container) {
				container = container || "#Browser";
				console.log(container);
				console.log("Attaching zooming handlers");
				$(container).find(".pair").each(function(){
					$(this).find("A.mol").click(function(){
						var pair = $(this).closest("DIV.pair");
						var molIds = new Array(pair.attr("mol1"), pair.attr("mol2"));
						loadPair(molIds);
						return false;
					});
				});
			}
			
			function loadPair(molIds) {
				$("#molViewDialog").find("IMG").eq(0).attr("src", "depiction.jsp?w=400&amp;h=400&amp;mp2=" + molIds[0] + "&amp;mp2templ=" + molIds[1]);
				$("#molViewDialog").find("IMG").eq(1).attr("src", "depiction.jsp?w=400&amp;h=400&amp;mp2=" + molIds[1] + "&amp;mp2templ=" + molIds[0]);
				molViewDialog.render();
				molViewDialog.show();
				
				loadPairData(molIds);
			}
			
			function loadPairData(molId) {
				$("#propTable TBODY").html("");
				$("IMG.progress").removeClass("invisible");
				var view = new View({element: "propTableLine"});
				var ajax = new QSPR.Ajax();
				ajax.send({
					url: "epbrowser/list.do?name=M" + molId[0] + ",M" + molId[1] + "&amp;pagesize=1000"+ "&amp;hidedummy=1",
					success: function(response) {
						var eps = array(response.list["exp-property"]);
						var propTable = new Array(new Array(), new Array());
						for (var i = 0; i &lt; eps.length; i++)
						{
							var idx = eps[i].molecule.mapping.id == molId[0] ? 0 : 1;
							propTable[idx][eps[i].property.name] = eps[i];
						}
						
						for (var prop in propTable[0])
							if (propTable[1][prop] &amp;&amp; propTable[1][prop].property)
								$("#propTable TBODY").append(view.render([propTable[0][prop], propTable[1][prop]]));
							
						$("#propTable TR.value.clickable").click(function(){
							window.openTab("MMP: data records", "epbrowser/show.do?approval-status=all&amp;property=" + $(this).attr("property") + "&amp;name=" + $(this).attr("mols"));
						});
							
						$(document).trigger("DOM_updated", $("#propTable TBODY"));
					},
					after: function() {
						$("IMG.progress").addClass("invisible");	
					}
				});
			}
			
		</script>
		
		<div id="molViewDialog"> 
		    <div class="hd">Matched pair profile</div> 
		    <div class="bd">
		    	<table id="propTable" width="100%">
		    		<thead>
		    			<tr>
		    				<td width="10"></td>
		    				<td><a href="javascript:void()" onclick="molViewDialog.cancel(); return false;">
					    		<img width="330" height="330" src=""/>
					    		</a>
					    	</td>
		    				<td>
		    					<a href="javascript:void()" onclick="molViewDialog.cancel(); return false;">
						    		<img width="330" height="330" src=""/>
						    	</a> 
		    				</td>
		    			</tr>
		    		</thead>
		    		<tbody>
		    		</tbody>
		    	</table>
		    	<img src="img/long_green.gif" class="invisible progress"/>
		    </div> 
		</div>
		
		<script type="text/template" id="propTableLine">
			<tr class="value clickable" property="[%=data[0].property.id %]" mols="M[%=data[0].molecule.mapping.id %],M[%=data[1].molecule.mapping.id %]">
				<td><nobr>[%=data[0].property.name %]</nobr></td>
				<td>[%=data[0].printableValue %]</td>
				<td>[%=data[1].printableValue %]</td>
			</tr>
		</script>
		
		</xsl:template>
		
		
		
	
</xsl:stylesheet>
