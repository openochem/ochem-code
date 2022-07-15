<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	
	<xsl:include href="../../helper.xslt" />
	
	<xsl:template name="content">
		<script type="text/javascript" src="js/commons/actionable.js" />
		<script type="text/javascript" src="js/commons/browser.js" />
		<script src="js/lib/springy/springy.js"></script>
		<script src="js/lib/springy/springyui.js?ver=2.3.6"></script>
		<style type="text/css">
			.itunes-right {padding: 10px 10px !important;}
			CANVAS {border: 1px solid gray;}
		</style>
		<table height="100%" width="100%">
			<tr><td class="itunes-up silver">
			
				<h1><img src="img/icons/mmp-48.png"/><span class="setcompare">MatchedPairs</span>: Fragment graph (experimental)</h1>
				A graph of molecular fragments based on significant transformations for the selected property. This is an experimental feature.
				</td></tr>
			<tr>
				<td class="itunes-right">
					
					An arrow means a transformation, the arrow's direction shows its effect. 
					Experimental, just trying around. <a class="download">Download the graph as image</a><br/><br/>
					<a class="fb-button invisible" href="#" id="btnShowAll" onclick="showAll(); return false;">Back to all fragments</a><br/>
					<div class="graph"></div>
				</td>
			</tr>
		</table>
		
		<script language="javascript">
			var multiplier = getParams["multiplier"] || 1;
			var imgMultiplier = getParams["imgMultiplier"] || 1;
			var arrMultiplier = getParams["arrMultiplier"] || 1;
			var h = getParams["h"] || 0;
			var filters = new Filters();
			filters.setFromUrl();
			
			function draw() {
				new QSPR.Ajax().send({
					url: "matchedpairs/getFragmentGraph.do?" + filters.getQueryString(),
					success: function(response) {
						var graph = new Springy.Graph();
	
						var nodes = new Object();
	
						var i = 0;
						var w = response.length &gt; 1000 ? 30 : response.length &gt; 100 ? 50 : response.length &gt; 30 ? 70 : 100;
						w = w / 2;
						if (response.length &lt; 10)
							w = 100;
							
						w *= multiplier;
						w *= imgMultiplier;
						var getNode = function(id) {
							if (typeof nodes[id] != 'undefined')
							{
								return nodes[id];
							}
							else
								return nodes[id] = graph.newNode({label: '' + id, image: {src: 'depiction.jsp?alpha=0.9&amp;w=' + w + '&amp;h=' + w + '&amp;mmp_frag=' + id, width: w, height: w}});
							
						}					
						
						for (var i = 0; i &lt; response.length; i++) {
							graph.newEdge(getNode(response[i][0]), getNode(response[i][1]), {color: '#00A0B0'});		
						}
						
						var layout = new Springy.Layout.ForceDirected(
						  graph,
						  400.0, // Spring stiffness
						  400.0, // Node repulsion
						  0.5 // Damping
						);
						
						var renderer = new Springy.Renderer(
						  layout,
						  function clear() {
						    // code to clear screen
						  },
						  function drawEdge(edge, p1, p2) {
						    // draw an edge
						  },
						  function drawNode(node, p) {
						    // draw a node
						  }
						);
								
					 $("DIV.graph > canvas").remove();
					 var canvas = $("canvas.graph").clone();
					 //canvas.height(800);
					 //canvas.width(1024);
					 if (h == 0) {
					 	h = $(document).height() - $(".itunes-up").height() - 100;
					 	h *= multiplier;
					 }
					 
					 $("DIV.graph").append('<canvas class='graph invisible' height="'+h+'" width="'+h * 3 / 2+'"/>');
					 $("DIV.graph > canvas").removeClass("invisible");
					 
					 $("a.download").attr("href", $("DIV.graph > canvas").get(0).toDataURL());
					 $("a.download").get(0).download = "fragment-graph.png";
					 
					  var springy = window.springy = jQuery("DIV.graph > canvas").springy({
					  	renderer: renderer,
					    graph: graph,
					    transparentNodes: true,
					    multiplier: multiplier,
					    bgColor: "#FFFFFF",
					    nodeSelected: function(node){
					    	if (!getParams["noclick"])
					    	{
					    		$("#btnShowAll").removeClass("invisible");
					      		console.log('Node selected: ' + JSON.stringify(node.data));
					      		filters.setValue("fragment", node.data.label);
					      		draw();
					      	}
					    }
					  });
					}
				});
			}
			
			$("a.download").click(function(link) {
				 $("a.download").attr("href", $("DIV.graph > canvas").get(0).toDataURL());
				$("a.download").get(0).download = "fragment-graph.png";
			});
			
			function showAll()
			{
				$("#btnShowAll").addClass("invisible");
				filters.setValue("fragment", "");
				draw();
			}
			
			$(draw);
		</script>
		
		<canvas class='graph invisible' height="800" width="1024"/>
		<div style="position: absolute; height: 1px; bottom: 1px;"/>
		
	</xsl:template>
</xsl:stylesheet>
