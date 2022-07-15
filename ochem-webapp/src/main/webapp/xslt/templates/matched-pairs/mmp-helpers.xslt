<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
	version="1.0">

	<xsl:template name="transformation-stats-layout">
		<script type="text/javascript" src="js/blocks/mmp.js" />
		<div class="transformations-browser">
			<nobr>
				MMP transformations: minimum
				<select name="minPairs" filter="1">
					<option value="1">1</option>
					<option value="2">2</option>
					<option value="4" selected="1">4</option>
					<option value="10">10</option>
					<option value="20">20</option>
					<option value="50">50</option>
					<option value="100">100</option>
					<option value="200">200</option>
				</select>
				pairs, p-value
				<select name="pValue" filter="1">
					<option value="1">Any</option>
					<option value="0.05" selected="1">0.05</option>
					<option value="0.01">0.01</option>
					<option value="0.001">0.001</option>
					<option value="0.0001">0.0001</option>
				</select>
				,
				max. atoms
				<select name="maxAtoms" filter="1">
					<option value="4">4</option>
					<option value="5">5</option>
					<option value="6">6</option>
					<option value="7">7</option>
					<option value="8">8</option>
					<option value="9">9</option>
					<option value="10" selected="1">Any</option>
				</select>
			</nobr>
			<br />
			<br />
			<div class="pager-strip">
				<span>
					<b class="showed">none</b>
					of
					<b class="total">none</b>
				</span>
				<div id="pager" class="pgr" />

			</div>
			<div id="TransformationsBrowser" />

			<div>
				<div id="query-status" />

			</div>
			<div class="pager-strip">
				<span>
					<b class="showed">none</b>
					of
					<b class="total">none</b>
				</span>
				<div id="pager" class="pgr" />

			</div>
			<br />
			<a action="sendmols" class="fb-button" title="re-calculate MMPs for all molecules">Resubmit all</a>
			<a action="openFragments" class="fb-button"
				title="Show the graph of molecular fragments derived from MMP transformations">Show fragment graph</a>
<!-- 			<a action="annotate" class="fb-button"
				title="Save the identified 'significant' transformations as a set for further reuse">Save transformations</a>
 -->			<a action="export" class="fb-button" title="Export pairs in CSV format">Export pairs</a>
			<a action="cleancache" class="fb-button" title="Clean cache">Clean cache</a>
		</div>
		<div class="pairs-browser invisible">
			<table class="transformation invisible">
				<tr>
					<td style="text-align: left;">
						Transformation details:
						<br />

						<span class="frag-container">
							<img src="" />
							<a action="openFragments" title="Show in the fragment graph">FG</a>
						</span>
						&#8594;
						<span class="frag-container">
							<img src="" />
							<a action="openFragments" title="Show in the fragment graph">FG</a>
						</span>
						<br />
						<br />
						<a href="#" onclick="transformationStats.showAll(); return false;"
							class="fb-button">Back to all transformations</a>
					</td>
					<td style="text-align: left;">
						<br />
						Transformation ID:
						<a class="trLink" tab="Transformation profile" href="matchedpairs/transformationProfile.do?">
							TR
							<span id="transformationId" />
						</a>
						<br />
						SMIRKS:
						<span id="transformationSmirks" />
						<br />
					</td>
					<td style="padding-left: 40px;">
						<table class="histogram invisible">
							<tr>
								<td>
									&#916;
									<sub>pair</sub>
									histogram
								</td>
							</tr>
							<tr>
								<td>
									<div id="histogram" />
								</td>
							</tr>
							<tr>
								<td>
									&#916;
									<sub>pair</sub>
									(measured)
								</td>
							</tr>
						</table>
					</td>
				</tr>
			</table>
			<br />
			<br />
			Pairs having the same transformation as the selected one:
			<div class="pager-strip">
				<span>
					<b class="showed">none</b>
					of
					<b class="total">none</b>
				</span>
				<div id="pager" class="pgr" />

			</div>
			<div id="PairsBrowser" />

			<div>
				<div id="query-status" />

			</div>
			<div class="pager-strip">
				<span>
					<b class="showed">none</b>
					of
					<b class="total">none</b>
				</span>
				<div id="pager" class="pgr" />

			</div>
		</div>

		<div id="annotatedialog" title="Select annotation set" style="display: none;">
			Save the significant transformations for further reuse:
			<br />
			<br />
			<table>
				<tr>
					<td>
						<input type="radio" name="atype" id="newannotation" value="new"
							checked="checked" />
					</td>
					<td>
						<label for="newannotation">New annotation set</label>
					</td>
				</tr>
				<tr>
					<td />
					<td>
						<input type="text" name="setnametxt" />
					</td>
				</tr>
				<tr>
					<td>
						<input type="radio" name="atype" id="eannotation" value="existing" />
					</td>
					<td>
						<label for="eannotation">Existing annotation set</label>
					</td>
				</tr>
				<tr>
					<td />
					<td>
						<select name="setnamesel" />
					</td>
				</tr>
			</table>
		</div>
	</xsl:template>

	<xsl:template name="transformation-stats-scripts">
		<script type="text/javascript" src="js/commons/actionable.js" />
		<script type="text/javascript" src="js/commons/browser.js?ver=2.3.2" />
		<script language="javascript" src="js/commons/dynamic-select.js" />
		<script language="javascript" type="text/javascript" src="js/lib/excanvas.min.js" />
		<script language="javascript" type="text/javascript"
			src="js/lib/jquery.flot-0.6.min.js" />
		<script language="javascript" src="js/commons/plotting.js" />
		<style type="text/css">
			#mmp-plot {
			width: 500px;
			height: 400px;
			float: left;
			}

			.yaxis {
			top: 450px;
			}

			.yaxis {
			position: absolute;
			top: 500px;
			left: 2px;
			transform: rotate(-90deg);
			-o-transform: rotate(-90deg);
			-ms-transform: rotate(-90deg);
			-moz-transform: rotate(-90deg);
			-webkit-transform: rotate(-90deg);
			transform-origin: 0 0;
			-o-transform-origin: 0 0;
			-ms-transform-origin: 0 0;
			-moz-transform-origin: 0 0;
			-webkit-transform-origin: 0 0;
			margin-left: 50px;
			}

			.transformation {text-align: center; min-height: 100px;}
			.transformation TD
			{padding-right: 20px; vertical-align: top;}
			.pair IMG, .transformation
			IMG {
			width: 100px;
			height: 100px;
			border: 1px solid gray;
			}
			.compact-item {height: 200px;}
			.pair {float: left;
			margin: 15px 30px 15px 0px; font-size: 10pt;}
			.pair A {color:
			#808080;}
			.pair {text-align: center;}

			TABLE.plot-container {margin-left: 50px;}
			TABLE.plot-container TD {text-align: center;
			vertical-align: top;}

			.roller {position: absolute; left: 280px; top: 350px; font-size: 200%;
			color: #777;}

			.transformations-browser {
			width: 600px;
			}

			.classification .transformations-browser {
			width: 900px !important;
			}

			TABLE.plot-container.classification {margin-left: -2px;}

			#TransformationsBrowser > DIV {padding: 5px 0px; border-bottom: 1px solid gray;}
			#TransformationsBrowser > DIV:hover {cursor: pointer;}

			.transformations-browser SELECT {font-size: 10pt;}

			.pair .values {font-size: 14pt; color: #222;}

			.increased {background-color: #FDD;}
			.decreased {background-color: #DFD;}

			#histogram {
			width: 300px;
			height: 150px;
			}

			SPAN.frag-container {position: relative;}
			.frag-container A {display: none; visibility:
			hidden; position: absolute; right: 10px; top: -30px;
			background-color: rgba(90, 90, 255, 0.7); text-decoration: none;
			padding: 4px; color: white; border-radius: 5px;}
			.frag-container:hover A {display: block; visibility: visible;}


		</style>

		<script language="javascript">
			include.plugins('view');

			var transformationStats = new function() {
			var plot = "";
			var selfTS = this;
			var ajax = new QSPR.Ajax();

			this.subset = "";
			this.currentPair = "";
			this.currentTransformation = "";
			this.currentSimilarity = 0;

			this.changeSimilarity = function()
			{
			var newSimilarity = $("[name=similarity]").val();
			if (newSimilarity != selfTS.currentSimilarity)
			{
			selfTS.currentSimilarity = newSimilarity;

			selfTS.showAll();
			transBrowser.filters.setValue("similarity", newSimilarity);
			transBrowser.request(true, function(){
			selfTS.plot = null;
			selfTS.loadChart();
			});
			}
			}

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

			this.doOpenfragments = function(link) {
			openTab("Fragment graph", "matchedpairs/fragmentGraph.do?" +
			transBrowser.filters.getQueryString() + "&amp;fragment=" +
			$(link).parent(".frag-container").find("img").attr("frag"));
			return false;
			}
			}

			function TransformationsBrowser()
			{
			var self = this;
			this.scope = ".transformations-browser";
			Browser.call(this);

			this.resetFiltersOnLoad = false;
			this.filters.scope = this.scope;
			this.container = "TransformationsBrowser";
			this.url = "mmpqsar/getTransformationsStats.do";
			this.filters.setValue("id", getParams["id"]);
			this.itemElement = "mmp-transformation";
			this.view = new View({element: "transformation-template"});
			this.filters.setFromUrl();

			this.annotationSelect = new DynamicSelect('setnamesel', 'mmpqsar/listAnnotationSets.do',
			'mmpAnnotationSet');

			this.listenEvent("items_loaded", function(){
			selfTS.subset = self.filters.getValue("subset");
			self.mainDiv.find("div[rec-id]").click(function(){
			selfTS.currentTransformation = $(this).attr("rec-id");
			selfTS.currentPair = "";
			selfTS.selectTransformation();
			});
			});

			this.doCleancache = function(link) {
			openTab("MMPs for a basket", "mmpqsar/cleanCache.do?" +
			transBrowser.filters.getQueryString());
			}

			this.doSendmols = function(link) {
			openTab("Sending molecules for indexing", "mmpqsar/indexMolecules.do?" +
			transBrowser.filters.getQueryString());
			}

			this.doOpenfragments = function(link) {
			openTab("Fragment graph", "matchedpairs/fragmentGraph.do?" +
			transBrowser.filters.getQueryString() + "&amp;fragment=" +
			$(link).parent(".frag-container").find("img").attr("frag"));
			}

			this.doExport = function(link) {
			openTab("Export pairs", "mmpqsar/exportPairs.do?" +
			transBrowser.filters.getQueryString());
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
			url: "mmpqsar/annotate.do?" + self.filters.getQueryString() +
			"&amp;name=" + name ,
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

			this.selectTransformation = function() {
			$(".pairs-browser").removeClass("invisible");
			$(".transformations-browser").addClass("invisible");
			pairsBrowser.filters.setValue("subset", selfTS.subset);
			pairsBrowser.filters.setValue("transformationPair", selfTS.currentPair);
			pairsBrowser.filters.setValue("transformation", selfTS.currentTransformation);
			pairsBrowser.request();

			$(".histogram").addClass("invisible");

			// Display the transformation data
			new QSPR.Ajax().send({
			url: "matchedpairs/getTransformation.do?transformationPair=" +
			selfTS.currentPair + "&amp;transformation=" +
			selfTS.currentTransformation,
			success: function(r) {
			var tr = r["mmp-transformation"];
			var tb = $(".pairs-browser .transformation");
			tb.removeClass("invisible");
			tb.find("img").eq(0).attr("src", "depiction.jsp?mmp_frag=" + tr.frag1Id).attr("frag", tr.frag1Id);
			tb.find("img").eq(1).attr("src", "depiction.jsp?mmp_frag=" +
			tr.frag2Id).attr("frag", tr.frag2Id);
			$("#transformationId").html(tr.id);
			$("#transformationSmirks").html(tr.smirks);
			tb.find("A.trLink").attr("href",
			"matchedpairs/transformationProfile.do?id=" + tr.id);

			$(".histogram").removeClass("invisible");
			MMP.drawHistogram("subset=" + selfTS.subset + "&amp;transformation=" + tr.id);
			}
			});
			selfTS.loadChart();
			}

			var pointClicked = function(setNum, pointNum)
			{
			var view = new View({element: "mmpValue"});
			var splitter = new View({element: "tableSplitter"}).render({});
			var point = this.dataArray[setNum].data[pointNum];
			var molIds = [point[3], point[4]];

			if (setNum == 1)
			{
			loadPair(molIds);

			$("#propTable TBODY").append(view.render({property: "Measured", values:
			point[6]}));
			$("#propTable TBODY").append(view.render({property: "Predicted", values:
			point[5]}));
			$("#propTable TBODY").append(splitter);
			}

			if (setNum == 0)
			{
			//loadPair(molIds);

			selfTS.currentPair = point[2];
			selfTS.currentTransformation = "";
			selfTS.selectTransformation();
			//selfTS.loadChart();
			}
			}

			this.loadChart = function() {

			if (MMP.isClassification)
			{
			$(".plot-container").addClass("classification");
			return;
			}

			$("#mmp-plot").parent().removeClass("invisible");
			$(".roller").removeClass("invisible");
			var url = "mmpqsar/" + selfTS.chartURL + ".do?" +
			window.location.search.substring(1);
			if (selfTS.currentPair)
			url += "&amp;transformationPair=" + selfTS.currentPair;
			if (selfTS.currentTransformation)
			url += "&amp;transformation=" + selfTS.currentTransformation;
			if ($("[name=similarity]").val())
			url += "&amp;similarity=" + $("[name=similarity]").val();
			ajax.send({
			url: url,
			success: function(response) {
			$(".roller").addClass("invisible");
			var plot = selfTS.plot;
			if (!plot)
			plot = new Plot();
			else
			{
			//plot.adaptSize = false;
			plot.dataArray[0].color = "#FCC";
			plot.defaultChartColors[1] = "#080";
			while (plot.dataArray.length &gt; 1)
			plot.dataArray.pop();
			}

			var data = [];
			var pairs = response.pairs;
			selfTS.subset = response.subset;
			for (var i = 0; i &lt; pairs.length; i++)
			{
			var p = pairs[i];
			data.push(selfTS.getPointArray(p));
			}

			plot.addData(data);
			plot.drawLine(0, plot.minY, 0, plot.maxY);
			plot.drawLine(plot.minX, 0, plot.maxX, 0);
			savedDimensions = [plot.minX, plot.maxX, plot.minY, plot.maxY];

			$("#mmp-plot").off();
			$("#mmp-plot").html("");
			plot.render("#mmp-plot");
			plot.pointClicked = pointClicked;

			selfTS.plot = plot;
			}
			});
			}

			this.showAll = function() {
			$(".pairs-browser").addClass("invisible");
			$(".transformations-browser").removeClass("invisible");
			}

			$(function(){
			pairsBrowser = new MatchedPairsBrowser();
			pairsBrowser.initialize(true);

			transBrowser = new TransformationsBrowser();

			transBrowser.filters.setValue("similarity",
			$("[name=similarity]").val());
			transBrowser.init();
			transBrowser.initialize(true);
			transBrowser.request(true, function(){
			selfTS.loadChart();
			});
			});
			}();
		</script>

		<script type="text/template" id="transformation-template">
			<table class="transformation">
				<tr>
					<td>
						<nobr>

							<span class="frag-container">
								<img src="depiction.jsp?mmp_frag=[%=frag1Id %]" frag="[%=frag1Id %]" />
								<a action="openFragments" title="Show in the fragment graph">FG</a>
							</span>
							&#8594;
							<span class="frag-container">
								<img src="depiction.jsp?mmp_frag=[%=frag2Id %]" frag="[%=frag2Id %]" />
								<a action="openFragments" title="Show in the fragment graph">FG</a>
							</span>
						</nobr>
					</td>
					<td align="left" style="text-align: left;">
						<span style="font-size: 120%">[%=statistics.pairsCount %] matched pairs</span>
						<br />
						<nobr>
							[% if (typeof statistics.nPP != 'undefined') {
							MMP.isClassification = true; %]
							<span
								title="[%=statistics.nNP %] molecules changed from 'negative' to 'positive'">&#8593; [%=statistics.nNP %] pairs</span>
							,
							<span
								title="[%=statistics.nPN %] molecules changed from 'positive' to 'negative'">&#8595; [%=statistics.nPN %] pairs</span>
							<br />
							<span
								title="[%=statistics.nNN %] molecules remaived 'negative', [%=statistics.nPP %] molecules remained 'positive'">[%=statistics.nNN %] negative and [%=statistics.nPP %]
								positive molecules unaffected</span>
							<br />
							[% } else { %]
							&#916;
							<sub>mean</sub>
							= [%=statistics.deltaMean %] &#177; [%=statistics.deltaStd %]
							[% } %]
						</nobr>
						<br />
						SMIRKS: [%=smirks %]
						<br />

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
						<td>
							<a title="Click to enlarge" class="mol" href="javascript:void()">
								<img src="depiction.jsp?mp2=[%=molecule1.id %]" />
							</a>
						</td>
						<td>
							<a title="Click to enlarge" class="mol" href="javascript:void()">
								<img src="depiction.jsp?mp2=[%=molecule2.id %]" />
							</a>
						</td>
					</tr>
					<tr>
						<td>
							<a href="molecule/profile.do?id=[%=molecule1.id %]" tab="Molecule profile">M[%=molecule1.id
								%]</a>
						</td>
						<td>
							<a href="molecule/profile.do?id=[%=molecule2.id %]" tab="Molecule profile">M[%=molecule2.id
								%]</a>
						</td>
					</tr>
				</table>
			</div>
		</script>

		<script type="text/template" id="mmpValue">
			<tr class="value">
				<td>
					<nobr>[%=data["property"] %]</nobr>
				</td>
				<td>[%=data["values"][0] %]</td>
				<td>[%=data["values"][1] %]</td>
			</tr>
		</script>

		<script type="text/template" id="tableSplitter">
			<tr class="tableSplitter">
				<td />
				<td />
				<td />
			</tr>
		</script>

	</xsl:template>




</xsl:stylesheet>
