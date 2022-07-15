<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	
	<xsl:include href="../../helper.xslt" />
	<xsl:include href="pair-zoom.xslt" />
	
	<xsl:template name="content">
		<script type="text/javascript" src="js/commons/actionable.js" />
		<script type="text/javascript" src="js/blocks/mmp.js" />
		<script type="text/javascript" src="js/commons/browser.js" />
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
			
			.transformation { min-height: 100px; margin-bottom: 10px;}
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
	
		<table height="100%" width="100%">
			<tr><td class="itunes-up">
				<h1><img src="img/icons/mmp-48.png"/><span class="setcompare">MatchedPairs</span>: MMP-based model interpretation (experimental) <a class="infolink" target="_blank" href="https://docs.ochem.eu/display/MAN/Molecular+Matched+Pairs" title="Click to read more about molecular matched pairs (aka MMPs) and their use in OCHEM"></a></h1>
				</td></tr>
			<tr>
				<td class="itunes-right">
					
					
					<table>
						<tr>
							<td>
								<div id="mmp-plot" style="width:1200px;height:800px">
								</div>
							</td>
							<td>
								<div id="transformation">
								</div>
							</td>
						</tr>
					</table>
					
				</td>
				
			</tr>
		</table>
		
		<script language="javascript">
			include.plugins('view');
			$(function(){
				var ajax = new QSPR.Ajax();
				
				var map = new Object();
				
				ajax.send({
					url: "mmpqsar/getSignificanceDataClass.do?out=json",
					success: function(response) {
						var plot = new Plot();
						
						for (var i = 0; i &lt; response.length; i++)
						{
							var dt = [];
							var s = response[i];
							for (var t = 0; t &lt; s.length; t++)
							{
								var tData = s[t];
								dt.push([tData[2], tData[1], tData]);
								map["" + tData[0] + "," + i] = tData;
							}
							plot.addData(dt);
							plot.pointClicked = function(seriesNum, pointNum) {
								var point = plot.dataArray[seriesNum].data[pointNum];
								
								new QSPR.Ajax().send({
									url: "matchedpairs/getTransformation.do?transformation=" + point[2][0],
									success: function(r) {
										var tr = r["mmp-transformation"];
										$("#transformation").html("");
										for (var k = 0; k &lt; 2; k++)
										{
											var row = map["" + tr.id + "," + k];
											tr.statistics = row[3];
											tr.type = k;
											tr.sigLevel = row[2].toFixed(1);
											var html = new View({element: "transformation-template"}).render(tr);
											$("#transformation").append(html);
										}
									}
								});
								
								
							}
						}
						
						plot.render("#mmp-plot");
					}
				});
			});
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
						<span style="font-size: 120%">[%=statistics.pairsCount %] pairs [% if (type == 0) { %](measured only)[% } else { %](measured and predicted)[% } %]</span><br/>
						<nobr>
						[% if (typeof statistics.nPP != 'undefined') { %]
							<span title="[%=statistics.nNP %] 'negative' molecules out of [%=statistics.nNP + statistics.nNN%] changed to 'positive'">&#8593; [%=statistics.nNP %] pairs out of [%=statistics.nNP + statistics.nNN%] negatives</span>, 
							<span title="[%=statistics.nPN %] 'positive' molecules out of [%=statistics.nPN + statistics.nPP%] changed to 'negative'">&#8595; [%=statistics.nPN %] pairs out of [%=statistics.nPN + statistics.nPP%] positives</span>,
							<span>(sig. level [%=sigLevel %])</span><br/>
							<br/>
						[% } else { %]
							&#916;<sub>mean</sub> = [%=statistics.deltaMean %] &#177; [%=statistics.deltaStd %] (sig. level [%=sigLevel %])
						[% }  %]
						</nobr>
						
						 
						<!-- <br/>SMIRKS: [%=smirks %]<br/>  -->
					
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
