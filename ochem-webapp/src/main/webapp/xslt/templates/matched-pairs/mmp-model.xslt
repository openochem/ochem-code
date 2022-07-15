<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	
	<xsl:include href="../../helper.xslt" />
	<xsl:include href="pair-zoom.xslt" />
	
	<xsl:template name="content">
		<script type="text/javascript" src="js/commons/actionable.js" />
		<script type="text/javascript" src="js/blocks/mmp.js" />
		<script type="text/javascript" src="js/commons/browser.js" />
			<script language="javascript" src="js/commons/dynamic-select.js"></script>
		<script language="javascript" type="text/javascript" src="js/lib/excanvas.min.js"></script>
	  	<script language="javascript" type="text/javascript" src="js/lib/jquery.flot-0.6.min.js"></script>
	  	<script language="javascript" src="js/commons/plotting.js"></script>
		
		<style type="text/css">
			#mmp-plot {
				width: 500px;
				height: 500px;
				float: left;
			}
			
			#cluster-map 
			{
				padding-top: 15px;
				padding-bottom: 15px;
			}
			
			.yaxis {
				position: absolute;
				top: 560px;
			    left: 2px;
			    transform: rotate(-90deg);
			    -o-transform: rotate(-90deg);
			    -ms-transform: rotate(-90deg);
			    -moz-transform: rotate(-90deg);
			    -webkit-transform:  rotate(-90deg);
			    transform-origin: 0 0;
			    -o-transform-origin: 0 0;
			    -ms-transform-origin: 0 0;
			    -moz-transform-origin: 0 0;
			    -webkit-transform-origin: 0 0;
			    margin-left: 50px;
			}
			
			.transformation { min-height: 100px;}
			.transformation TD {padding-right: 20px; vertical-align: top;}
			.pair IMG, .transformation IMG {
				width: 100px;
				height: 100px;
				border: 1px solid gray;
			}
			.compact-item {height: 200px;}
			.pair {float: left; margin: 15px 30px 15px 0px; font-size: 10pt;}
			.pair A {color: #808080;}
			.pair {text-align: center;}
			
			TABLE.plot-container {margin-left: 50px;}
			TABLE.plot-container TD {text-align: center; vertical-align: top;}
			
			.roller {position: absolute; left: 250px; top: 430px; font-size: 200%; color: #777;}
			
			.transformations-browser {
				width: 600px;
			}
			
			#TransformationsBrowser > DIV  {padding: 5px 0px; border-bottom: 1px solid gray;}
			#TransformationsBrowser > DIV:hover {cursor: pointer;}
			
			.transformations-browser SELECT {font-size: 10pt;}
			
			.pair .values {font-size: 14pt; color: #222;}
			
			.increased {background-color: #FDD;}
			.decreased {background-color: #DFD;}
			
			#histogram {
				width: 300px;
				height: 150px;
			}
			
			
		</style>
	
		<title>MMP analysis of a model</title>
		<table height="100%" width="100%">
			<tr><td class="itunes-up">
				<h1><img src="img/icons/mmp-48.png"/><span class="setcompare">MatchedPairs</span>: MMP-based model interpretation (experimental) <a class="infolink" target="_blank" href="https://docs.ochem.eu/display/MAN/Molecular+Matched+Pairs" title="Click to read more about molecular matched pairs (aka MMPs) and their use in OCHEM"></a></h1>
				</td></tr>
			<tr>
				<td class="itunes-right">
					This is a very early preview feature. The chart shows the MMP deltas for experimental and predicted values. This should help to identify activity cliffs.<br/>
					Each point is a matched molecular pair (MMP). Click on a point to see the pair details.
					<br/><br/>
					The info is based on <xsl:value-of select="/model/model/training-set/@size"/> records for <b><xsl:value-of select="/model/model/modelMappings[1]/property/@name"/></b>
					<xsl:if test="//others/mmpStatus/unindexedMolecules != 0">
						(excluding <xsl:value-of select="//others/mmpStatus/unindexedMolecules"/> not indexed molecules)
					</xsl:if>
					<br/>
					<div>
						Minimal pair similarity:
						<select name="similarity" filter="1" onchange="changeSimilarity(); return false;">
							<option value="0" selected="1">Any</option>
							<option value="25">25</option>
							<option value="50">50</option>
							<option value="75">75</option>
						</select> 
						<a class="infolink" title="Minimal Tanimoto similarity between two molecules in a pair. Used to hide matched pairs that with too dissimilar molecules."></a>
					</div>
					<br/><br/>
					<div class="yaxis">MMP Delta (predicted)</div>
					<table class="plot-container" width="97%">
						<tr>
							<td height="500">
								<span style="position: absolute" class="invisible roller"><img src="img/roller-big.gif"/><br/>Loading...<br/></span>
								&#916;<sub>pair</sub> chart <a class="infolink" help="delta-chart-help" href="https://docs.ochem.eu/display/MAN/Molecular+Matched+Pairs#MolecularMatchedPairs-ModelsInterpretation" target="_blank"></a>
								<div class="invisible" id="delta-chart-help">
									<b>What is &#916;<sub>pair</sub> chart?</b><br/>
									<i>&#916;<sub>pair</sub> chart</i> shows the deltas (differences in values) for each molecular pair in the training set of your model.<br/><br/>
									Each point is a molecular matched pair, where X is the actual effect of the molecular transformation and Y is the predicted effect.<br/>
									<br/><br/>
									<b>Why?</b><br/>
									The chart allows to identify potential activity cliffs and to distinguish well- and poorly-predicted cliffs. 
									Such activity cliffs allow for interpretation, since they are bounded to a matched pair and, therefore, to a small localized change in the molecular structure.	
								</div>
								<div id="mmp-plot">
								</div>
							</td>
							<td rowspan="3" width="100%" align="left" style="text-align: left; padding-left: 15px;">
								<div class="transformations-browser">
									<nobr>Transformations: minimum 
										<select name="minPairs" filter="1">
											<option value="1">1</option>
											<option value="2">2</option>
											<option value="4" selected="1">4</option>
											<option value="10">10</option>
											<option value="20">20</option>
											<option value="50">50</option>
										</select>
										 pairs, p-value <select name="pValue" filter="1">
											<option value="1">Any</option>
											<option value="0.05">0.05</option>
											<option value="0.01">0.01</option>
											<option value="0.001">0.001</option>
										</select>
										</nobr><br/>
									<div class="pager-strip">
										<span><b class="showed">none</b> of <b class="total">none</b></span>
										<div id="pager" class="pgr">
										</div>
									</div>
									<div id="TransformationsBrowser">
									</div>
									<div>
										<div id="query-status">
										</div>
									</div>
									<div class="pager-strip">
										<span><b class="showed">none</b> of <b class="total">none</b></span>
										<div id="pager" class="pgr">
										</div>
									</div>
									<br/>
									<a action="sendmols" class="fb-button" title="re-calculate MMPs for all molecules">Resubmit all molecules</a>
									<a action="openFragments" class="fb-button" title="Show the graph of molecular fragments derived from MMP transformations">Show fragment graph</a>
<!-- 									<a action="annotate" class="fb-button" title="Save the identified 'significant' transformations as a set for further reuse">Save transformations</a>
 -->									<a action="export" class="fb-button" title="Export pairs in CSV format">Export pairs</a>
								</div>
								<div class="pairs-browser invisible">
									
									<table class="transformation invisible">
										<tr>
											<td style="text-align: left;">
												Transformation details:<br/>
												<img src=""/> &#8594;
												<img src=""/>
												<br/><br/>
												<a href="#" onclick="showAll(); return false;" class="fb-button">Back to all transformations</a>
											</td>
											<td style="text-align: left;">
												<br/>
												Transformation ID: <a tab="Transformation profile" href="matchedpairs/transformationProfile.do?">TR<span id="transformationId"/></a><br/>SMIRKS: <span id="transformationSmirks"/><br/>
											</td>
											<td style="padding-left: 40px;">
												<table class="histogram invisible">
													<tr><td>&#916;<sub>pair</sub> histogram</td></tr>
													<tr><td><div id="histogram"></div></td></tr>
													<tr><td>&#916;<sub>pair</sub> (measured)</td></tr>
												</table>
											</td>
										</tr>
									</table>
									
									<br/>
									
									<xsl:if test="session/user/@login = 'novserj' ">
										<div id="clusterbtn"><a href="#" onclick="enableClusters(true); return false;" class="fb-button">Show structural clusters</a></div>
									</xsl:if>
									
									<div id="cluster-map">
									</div>
									<div style="clear:both"></div>
									
									Pairs having the same transformation as the selected one:
									<div class="pager-strip">
										<span><b class="showed">none</b> of <b class="total">none</b></span>
										<div id="pager" class="pgr">
										</div>
									</div>
									<div id="PairsBrowser">
									</div>
									<div>
										<div id="query-status">
										</div>
									</div>
									<div class="pager-strip">
										<span><b class="showed">none</b> of <b class="total">none</b></span>
										<div id="pager" class="pgr">
										</div>
									</div>
								</div>	
							</td>
						</tr>
						<tr>
							<td>MMP deltas (measured)</td>
						</tr>
						<tr>
							<td height="100%">&#160;</td>
						</tr>
					</table>
					
				</td>
			</tr>
		</table>
		
		<div id="annotatedialog" title="Select annotation set" style="display: none;">
			Save the significant transformations for further reuse:<br/><br/>
			<table>
				<tr><td><input type="radio" name="atype" id="newannotation" value="new" checked="checked"/></td><td><label for="newannotation">New annotation set</label></td></tr>
				<tr><td></td><td><input type="text" name="setnametxt"/></td></tr>
				<tr><td><input type="radio" name="atype" id="eannotation" value="existing"/></td><td><label for="eannotation">Existing annotation set</label></td></tr>
				<tr><td></td><td><select name="setnamesel"/></td></tr>
			</table>
		</div>
		
		<script language="javascript">
			var colors = ["#800000", "#008000", "#0000bd", "#827800", "#262626", "#8f00c7", "#0086fe", "#00fefe", "#fe68fe", "#fe8420", "#70fe00", "#fefe00", "#fed38b", "#a0d681"];
			include.plugins('view');
			var ajax = new QSPR.Ajax();
			var subset = "";
			var currentPair = "";
			var currentTransformation = "";
			var currentSimilarity = 0;
			var clusters = false;
			var plot;
			
			
			function changeSimilarity()
			{
				var newSimilarity = $("[name=similarity]").val();
				if (newSimilarity != currentSimilarity)
				{
					currentSimilarity = newSimilarity;
					
					showAll();
					plot = null;
					loadChart();
					
					transBrowser.filters.setValue("similarity", newSimilarity);					
					transBrowser.request();
				}
			}
			
			function pointClicked(setNum, pointNum)
			{
				var view = new View({element: "mmpValue"});
				var splitter = new View({element: "tableSplitter"}).render({});
				var point = this.dataArray[setNum].data[pointNum];
				var molIds = [point[3], point[4]];
				
				if (setNum == 1)
				{
					loadPair(molIds);
				
					$("#propTable TBODY").append(view.render({property: "Measured", values: point[6]}));
					$("#propTable TBODY").append(view.render({property: "Predicted", values: point[5]}));
					$("#propTable TBODY").append(splitter);
				}
				
				if (setNum == 0)
				{
					console.log(point);
					currentPair = point[2];
					currentTransformation = "";
					selectTransformation();
				}
			}
			
			function selectTransformation() 
			{
				$(".pairs-browser").removeClass("invisible");
				$(".transformations-browser").addClass("invisible");
				pairsBrowser.filters.setValue("subset", subset);
				pairsBrowser.filters.setValue("transformationPair", currentPair);
				pairsBrowser.filters.setValue("transformation", currentTransformation);
				pairsBrowser.request();
				
				$(".histogram").addClass("invisible");
				
				// Display the transformation data
				new QSPR.Ajax().send({
					url: "matchedpairs/getTransformation.do?transformationPair=" + currentPair + "&amp;transformation=" + currentTransformation,
					success: function(r) {
						var tr = r["mmp-transformation"];
						var tb = $(".pairs-browser .transformation");
						tb.removeClass("invisible");
						tb.find("img").eq(0).attr("src", "depiction.jsp?mmp_frag=" + tr.frag1Id);
						tb.find("img").eq(1).attr("src", "depiction.jsp?mmp_frag=" + tr.frag2Id);
						$("#transformationId").html(tr.id);
						$("#transformationSmirks").html(tr.smirks);
						tb.find("A").attr("href", "matchedpairs/transformationProfile.do?id=" + tr.id);
						
						$(".histogram").removeClass("invisible");
						MMP.drawHistogram("subset=" + subset + "&amp;transformation=" + tr.id);
					}
				});
				
				//If applicable, display cluster data
				loadClusters();
				
				loadChart();
			}
			
			function loadClusters()
			{
				$("#cluster-map").html("");
				if (clusters)
				{
					new QSPR.Ajax().send({
						url: "matchedpairs/getClusters.do?transformationPair=" + currentPair + "&amp;transformation=" + currentTransformation + "&amp;subset=" + subset,
						success: function(r) {
							var view = new View({element: "cluster-template"});
							for (var i=0; i&lt;r.length; i++)
								if (r[i][0] != null)
									$("#cluster-map").append(view.render({color:colors[i+1], mp2:r[i][0]}));
								else
									$("#cluster-map").append(view.render({color:colors[i+1], mp2:1}));
						}
					});			
				} 
			}
			
			function enableClusters(enable)
			{
				clusters = enable;
				if (clusters)
					$("#clusterbtn").html('<a href="#" onclick="enableClusters(false); return false;" class="fb-button">Hide structural clusters</a>');
				else
					$("#clusterbtn").html('<a href="#" onclick="enableClusters(true); return false;" class="fb-button">Show structural clusters</a>');
					
				selectTransformation();
			}
			
			function loadChartTransformations() {
				$(".roller").removeClass("invisible");
				var url = "mmpqsar/getDeltaChartTransformations.do?id=" + getParams["id"];
				
				ajax.send({
					url: url,
					success: function(response) {
						$(".roller").addClass("invisible");
						plot = new Plot();
						
						var data = [];
						var transformations = response.transformations;
						subset = response.subset;
						for (var i = 0; i &lt; transformations.length; i++)
						{
							var t = transformations[i];
							data.push([t[3], Math.abs(t[1]), t[0]]);
						}
						
						plot.addData(data);
						
						plot.drawLine(0, plot.minY, 0, plot.maxY);
						plot.drawLine(plot.minX, 0, plot.maxX, 0);
						//savedDimensions = [plot.minX, plot.maxX, plot.minY, plot.maxY];
						
						$("#mmp-plot").off();
						$("#mmp-plot").html("");
						plot.render("#mmp-plot");
						plot.pointClicked = function(setNum, pointNum) {
							var point = this.dataArray[setNum].data[pointNum];
							currentPair = "";
							currentTransformation = point[2];
							console.log(currentTransformation);
							selectTransformation();
						};
					}
				});
			}
			
			
			
			function loadChartPairs() {
			
				$(".roller").removeClass("invisible");
				var url = "mmpqsar/getChart.do?id=" + getParams["id"];
				
				var c = "";				
				if (clusters)
					c = "&amp;clusters=true";
				
					
				if (currentPair)
					url += "&amp;transformationPair=" + currentPair + c;
				if (currentTransformation)
					url += "&amp;transformation=" + currentTransformation + c;
				if (currentSimilarity)
					url += "&amp;similarity=" + currentSimilarity;	
				
				ajax.send({
					url: url,
					success: function(response) {
						$(".roller").addClass("invisible");
						if (!plot)
							plot = new Plot();
						else
						{
							//plot.adaptSize = false;
							plot.dataArray[0].color = "#FCC";
							plot.defaultChartColors = colors;
							while (plot.dataArray.length &gt; 1)
								plot.dataArray.pop();
						}
						
						var data = {};
						var pairs = response.pairs;
						subset = response.subset;
						for (var i = 0; i &lt; pairs.length; i++)
						{
							var p = pairs[i];
							if (data[p[5]] == undefined)
								data[p[5]] = [];
							data[p[5]].push([p[4][1] - p[4][0], 1*p[3][1]-1*p[3][0], p[0], p[1], p[2], p[3], p[4]]);
						}
						
						$.each(data, function(key, value)
						{
							plot.addData(value);
						});
						
						//plot.addData(data);
						plot.drawLine(0, plot.minY, 0, plot.maxY);
						plot.drawLine(plot.minX, 0, plot.maxX, 0);
						savedDimensions = [plot.minX, plot.maxX, plot.minY, plot.maxY];
						
						$("#mmp-plot").off();
						$("#mmp-plot").html("");
						plot.render("#mmp-plot");
						plot.pointClicked = pointClicked;
					}
				});
			}
			
			loadChart = loadChartPairs;//Transformations; 
			//loadChart();
			
			function showAll() 
			{
				$("#clustercheckbox").attr(':checked', false);
				$(".pairs-browser").addClass("invisible");
				$(".transformations-browser").removeClass("invisible");
				currentPair = "";
				currentTransformation = "";
				plot.adaptSize = false;
				plot.dataArray[0].color = "#C00";
				while (plot.dataArray.length &gt; 1)
					plot.dataArray.pop();
					
				plot.drawLine(0, plot.minY, 0, plot.maxY);
				plot.drawLine(plot.minX, 0, plot.maxX, 0);
				
				$("#mmp-plot").off();
				$("#mmp-plot").html("");
				plot.render("#mmp-plot");
				plot.pointClicked = pointClicked;
			}
			
			
		</script>
		
		<script type="text/template" id="mmpValue">
			<tr class="value">
				<td><nobr>[%=data["property"] %]</nobr></td>
				<td>[%=data["values"][0] %]</td>
				<td>[%=data["values"][1] %]</td>
			</tr>
		</script>
		
		<script type="text/template" id="tableSplitter">
			<tr class="tableSplitter">
				<td></td>
				<td></td>
				<td></td>
			</tr>
		</script>
		
		<script language="javascript">
			include.plugins('view');
			function MatchedPairsBrowser() {
				this.scope = ".pairs-browser";
				Browser.call(this);
				this.filters.scope = this.scope;
				this.container = "PairsBrowser";
				this.url = "matchedpairs/getPairs.do";
				//this.filters.setValue("transformation", getParams["id"]);
				this.itemElement = "mmpair";
				this.view = new View({element: "template"});
				this.filters.setFromUrl();
				
				this.listenEvent("items_loaded", function(){
					attachZoomingClickHandlers("#PairsBrowser");
					colorPairs();
				});
			}
			
			function TransformationsBrowser() {
				var self = this;
				this.scope = ".transformations-browser";
				Browser.call(this);
				this.filters.scope = this.scope;
				this.container = "TransformationsBrowser";
				this.url = "mmpqsar/getTransformationsStats.do";
				this.filters.setValue("id", getParams["id"]);
				this.itemElement = "mmp-transformation";
				this.view = new View({element: "transformation-template"});
				this.filters.setFromUrl();
				
				this.annotationSelect = new DynamicSelect('setnamesel', 'mmpqsar/listAnnotationSets.do', 'mmpAnnotationSet');
				
				this.listenEvent("items_loaded", function(){
					self.mainDiv.find("div[rec-id]").click(function(){
						currentTransformation = $(this).attr("rec-id");
						currentPair = "";
						selectTransformation();
					});	
				});
			
				this.doSendmols = function(link) {
						openTab("Sending molecules for indexing", "mmpqsar/indexMolecules.do?" + transBrowser.filters.getQueryString());
					}
				
				this.doExport = function(link) {
						openTab("Export pairs", "mmpqsar/exportPairs.do?" + transBrowser.filters.getQueryString());
					}
				
				this.doOpenfragments = function() {
					openTab("Fragment graph", "matchedpairs/fragmentGraph.do?" + self.filters.getQueryString());
				}
				
				this.doAnnotate = function() 
				{
					this.annotationSelect.update();
					$("#annotatedialog").dialog("open");
				}
				
				this.doAnnotatesubmit = function()
				{
					//Check who is selected
					var name = "";
					
					if ($("[name='atype']:checked").val() == "new")
						name = $("[name='setnametxt']").val();
					else
						name = $("[name='setnamesel'] option:selected").html()
					
					new QSPR.Ajax().send({
						url: "mmpqsar/annotate.do?" + self.filters.getQueryString() + "&amp;name=" + name ,
						success: function(r) 
						{
							$("#annotatedialog").dialog("close");
						}
					});
				}
				
				this.init = function()
				{
					$("#annotatedialog").dialog({
					autoOpen: false,
					height: 300,
					width: 550,
					modal: true,
					buttons: {
						"OK": function() {
							self.doAnnotatesubmit();								
						},
						Cancel: function() {
							$(this).dialog("close");
							}
						}
					});
				}
			}
			
			$(function(){
				pairsBrowser = new MatchedPairsBrowser();
				pairsBrowser.initialize(true);
				
				transBrowser = new TransformationsBrowser();
				transBrowser.init();
				transBrowser.initialize(true)
				transBrowser.request(true, function(){
					loadChart();
				});
			});
		</script>
		
		<script type="text/template" id="cluster-template">
			<div style="padding: 5px; float:left;"><div style="padding: 5px; background-color:[%=color %];"><img src="depiction.jsp?mp2=[%=mp2 %]&amp;w=100&amp;h=100"/></div></div>
		</script>
		
		<script type="text/template" id="transformation-template">
			
			<table class="transformation">
				<tr>
					<td>
						<nobr>
						<img src="depiction.jsp?mmp_frag=[%=frag1Id %]"/> &#8594;
						<img src="depiction.jsp?mmp_frag=[%=frag2Id %]"/>
						</nobr>
					</td>
					<td align="left" style="text-align: left;">
						<span style="font-size: 120%">[%=statistics.pairsCount %] matched pairs</span><br/>
						<nobr>
						[% if (typeof statistics.nPP != 'undefined') { %]
							<span title="[%=statistics.nNP %] molecules changed from 'negative' to 'positive'">&#8593; [%=statistics.nNP %] pairs</span>, 
							<span title="[%=statistics.nPN %] molecules changed from 'positive' to 'negative'">&#8595; [%=statistics.nPN %] pairs</span><br/>
							<span title="[%=statistics.nNN %] molecules remaived 'negative', [%=statistics.nPP %] molecules remained 'positive'">[%=statistics.nNN %] negative and [%=statistics.nPP %] positive molecules unaffected</span>
							<br/>
						[% } else { %]
							&#916;<sub>mean</sub> = [%=statistics.deltaMean %] &#177; [%=statistics.deltaStd %]
						[% }  %]
						</nobr>
						
						 
						<br/>SMIRKS: [%=smirks %]<br/>
					
					</td>
				</tr>
			</table>
								
		</script>
		
		<script type="text/template" id="template">
			<div class="pair" mol1="[%=molecule1.id %]" mol2="[%=molecule2.id %]">
				<table>
					[% if (data.values) { %]
					<tr class="values">
						<td>[%=data.values[0] %]</td>
						<td>[%=data.values[1] %]</td>
					</tr>
					[% } %]
					<tr>
						<td><a title="Click to enlarge" class="mol" href="javascript:void()"><img src="depiction.jsp?mp2=[%=molecule1.id %]"/></a></td>
						<td><a title="Click to enlarge" class="mol" href="javascript:void()"><img src="depiction.jsp?mp2=[%=molecule2.id %]"/></a></td>
					</tr>
					<tr>
						<td><a href="molecule/profile.do?id=[%=molecule1.id %]" tab="Molecule profile">M[%=molecule1.id %]</a></td>
						<td><a href="molecule/profile.do?id=[%=molecule2.id %]" tab="Molecule profile">M[%=molecule2.id %]</a></td>
					</tr>
				</table>
			</div>
		</script>
		
		<xsl:call-template name="pair-zoom"/>
		
	</xsl:template>
</xsl:stylesheet>
