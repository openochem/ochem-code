<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	
	<xsl:include href="../../helper.xslt" />
	<xsl:include href="pair-zoom.xslt" />
	
	<xsl:template name="content">
		<script type="text/javascript" src="js/commons/actionable.js" />
		<script type="text/javascript" src="js/commons/browser.js" />
		<script type="text/javascript" src="js/blocks/mmp.js" />
		<script language="javascript" type="text/javascript" src="js/lib/excanvas.min.js"></script>
	  	<script language="javascript" type="text/javascript" src="js/lib/jquery.flot-0.6.min.js"></script>
	  	<script language="javascript" src="js/commons/plotting.js"></script>
		<style type="text/css">
			.transformation {text-align: center; min-height: 100px;}
			.theader { font-size: 100%; border-bottom: 1px solid #666;}
			.theader TD {padding-top: 5px;}
			.theader:not(:first) TD {padding-top: 5px;}
			.transformation TD {padding-right: 20px; vertical-align: top; white-space: nowrap;}
			.pair IMG, .transformation .depiction IMG {
				width: 100px;
				height: 100px;
				border: 1px solid gray;
			}
			.compact-item {height: 200px;}
			.pair {float: left; margin: 15px 30px 15px 0px;; font-size: 10pt;}
			.pair A {color: #808080;}
			.pair {text-align: center;}
			
			TR.annotation TD {white-space: nowrap;}
			TR.annotation.selected {background-color: #EEFFEE; border: 1px solid black;}
			TR.annotation:hover {cursor: pointer; background-color: #EEFFEE; border: 1px solid gray;}
			
			#histogram {
				width: 300px;
				height: 150px;
			}
			
			.pair .values {font-size: 14pt; color: #222;}
			
			.increased {background-color: #FDD;}
			.decreased {background-color: #DFD;}
			
		</style>
	
		<table height="100%" width="100%">
			<tr><td class="itunes-up">
				<h1><img src="img/icons/mmp-48.png"/><span class="setcompare">MatchedPairs</span>: Molecular transformation profile<a class="infolink" target="_blank" href="https://docs.ochem.eu/display/MAN/Molecular+Matched+Pairs" title="Click to read more about molecular matched pairs (aka MMPs) and their use in OCHEM"></a></h1>
				A molecular transformation with all the corresponding molecular matched pairs 
				</td></tr>
			<tr>
				<td class="itunes-right">
				
					<table class="transformation">
						<tr>						
							<td class="depiction">
								<xsl:if test="/model/mmp-transformation-profile/mol1">
									<b>Evidence for molecules:</b><br/>
									<img src="depiction.jsp?mp2={/model/mmp-transformation-profile/mol1/@id}"/> &#8594;
									<img src="depiction.jsp?mp2={/model/mmp-transformation-profile/mol2/@id}"/>
									<br/><b>Transformation:</b><br/>
								</xsl:if>
								<img src="depiction.jsp?mmp_frag={//mmp-transformation/frag1Id}"/> &#8594;
								<img src="depiction.jsp?mmp_frag={//mmp-transformation/frag2Id}"/>
							</td>
							<td align="left">
								<xsl:if test="/model/mmp-transformation-profile/mol1">
									<br/><img src="/img/blank.gif" height="100px"></img><br/><br/>
								</xsl:if>
								Transformation ID: TR<xsl:value-of select="//mmp-transformation/id"/><br/>
								SMIRKS: <xsl:value-of select="//mmp-transformation/smirks"/><br/>
							</td>
							<td align="left" style="padding-left: 50px;">
								<xsl:if test="//mmp-transformation/annotations">
									<b>Affected properties</b><a class="infolink" help="affected-properties-help"></a>
									<div id="affected-properties-help" class="invisible">
										The list of properties/activities affected by this transformation. OCHEM automatically calculates statistics based on experimental data and identifies transformations that tend to change molecular properties.
										<br/><br/>
										Click on the property to see only the relevant molecular pairs.
									</div>
									<table>
									
									<xsl:for-each select="//mmp-transformation/annotations">
										<xsl:if test="not(preceding-sibling::annotations[1]/annotationSet/@id = annotationSet/@id)">
											<tr class="theader"><td colspan="2"><i><xsl:value-of select="annotationSet/name"/></i></td></tr>
										</xsl:if>
										
										<tr property="{property/@id}" class="annotation">
										<xsl:if test="not(annotationSet/@id = '0')">
											<xsl:attribute name="nodata">true</xsl:attribute>
										</xsl:if>
										<td>
											<xsl:if test="increasing = 'true'">
												<img src="img/icons/arrowup-green-16.png" /> 
											</xsl:if>
											<xsl:if test="increasing != 'true'">
												<img src="img/icons/arrowdown-red-16.png" /> 
											</xsl:if>
											<span><xsl:value-of select="property/@name"/></span>
										</td>
										<td>
											<xsl:choose>
												<xsl:when test="property/@type = 0">
													&#916;<sub>pair</sub> = <xsl:value-of select="subsetStats/deltaMean"/> &#177; <xsl:value-of select="subsetStats/deltaStd"/>
												</xsl:when>
												<xsl:otherwise>
													<span title="{subsetStats/@nNP} molecules changed from 'negative' to 'positive'">&#8593; <xsl:value-of select="subsetStats/@nNP"/> pairs</span>, 
													<span title="{subsetStats/@nPN} molecules changed from 'positive' to 'negative'">&#8595; <xsl:value-of select="subsetStats/@nPN"/> pairs</span>,
													<span title="{subsetStats/@nNN} molecules remaived 'negative', {subsetStats/@nPP} molecules remained 'positive'"><xsl:value-of select="subsetStats/@nNN"/> negative and <xsl:value-of select="subsetStats/@nPP"/> positive molecules unaffected</span>
												</xsl:otherwise>
											</xsl:choose>
											
										</td>
										</tr>
									</xsl:for-each>
									</table>
								</xsl:if>
							</td>
							<td style="padding-left: 20px; width: 310px; height: 160px;">
								<table class="histogram invisible">
									<xsl:if test="/model/mmp-transformation-profile/property">
										<tr><td>Property <b><xsl:value-of select="/model/mmp-transformation-profile/property/@name"/></b></td></tr>
									</xsl:if>
									<tr><td>&#916;<sub>pair</sub> histogram</td></tr>
									<tr><td><div id="histogram"></div></td></tr>
									<tr><td>&#916;<sub>pair</sub> (<span id="xAxis"></span>)</td></tr>
								</table>
							</td>
						</tr>
					</table>
					<br/>
					<div id="BrowserContainer">
						<b>Molecular pairs matching this transformation.</b> <br/>
						Minimal similarity:
							<select name="similarity" filter="1">
								<option value="0">Any</option>
								<option value="25">25</option>
								<option value="50">50</option>
								<option value="75">75</option>
							</select>
						
						<xsl:if test="/model/mmp-transformation-profile/propertyChangeDirection = '1'">
							<input type="checkbox" name="onlyIncreasing" checked="checked" filter="1" /> Only relevant pairs
						</xsl:if>
						<xsl:if test="/model/mmp-transformation-profile/propertyChangeDirection = '-1'">
							<input type="checkbox" name="onlyDecreasing" checked="checked" filter="1"/> Only relevant pairs
						</xsl:if>
						
						<div class="pager-strip">
							<span><b class="showed">none</b> of <b class="total">none</b></span>
							<div id="pager" class="pgr">
							</div>
						</div>
						<div id="Browser">
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
					<div id="NoBrowserContainer" class="invisible">
						<i>No matched pairs information for this annotation</i>
					</div>
				</td>
			</tr>
		</table>
		
		<script language="javascript">
			include.plugins('view');
			function MatchedPairsBrowser() {
				Browser.call(this);
				this.url = "matchedpairs/getPairs.do";
				this.filters.setValue("transformation", getParams["id"]);
				this.itemElement = "mmpair";
				this.view = new View({element: "template"});
				this.filters.setFromUrl();
				
				this.listenEvent("items_loaded", function(){attachZoomingClickHandlers();
					colorPairs();
				});
			}
			
			$(function(){
				if ($(".theader").length == 1)
					$(".theader").remove();
			
				var browser = new MatchedPairsBrowser();
				browser.initialize();
				
				$("[property]").click(function(item){
					if (!$(item.currentTarget).attr("nodata"))
					{
						$("#NoBrowserContainer").addClass("invisible");
						$("#BrowserContainer").removeClass("invisible");
					
						browser.filters.setValue("property", $(this).attr("property"));
						$("[property]").removeClass("selected");
						$(".histogram").addClass("invisible");
						$(this).addClass("selected");
						browser.request(true, function() {
							MMP.drawHistogram(browser.filters.getQueryString());
						});
					} else
					{
						$("#BrowserContainer").addClass("invisible");
						$(".histogram").addClass("invisible");
						$("#NoBrowserContainer").removeClass("invisible");
					}
				});
				
				<xsl:if test="/model/mmp-transformation-profile/property">
					MMP.drawHistogram(browser.filters.getQueryString());
				</xsl:if>
			});
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
