<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	
	<xsl:include href="../../helper.xslt" />
	<xsl:include href="pair-zoom.xslt" />
	
	<xsl:template name="content">
		<script type="text/javascript" src="js/commons/actionable.js" />
		<script type="text/javascript" src="js/commons/browser.js" />
		<style type="text/css">
			.transformation IMG, .pair IMG {
				width: 100px;
				height: 100px;
				border: 1px solid gray;
			}
			
			.compact-item {height: 200px;}
			
			.pair {float: left; margin-right: 30px;}
			
			.transformation {text-align: center; min-height: 100px; border-right: 1px solid gray;}
			.pair A, .transformation A {color: #808080;}
			.pair {text-align: center;}
			
		</style>
	
		<title>Molecular transformations browser</title>
		<table height="100%" width="100%">
			<tr><td class="itunes-up">
				<h1><img src="img/icons/mmp-48.png"/><span class="setcompare">MatchedPairs</span>: Molecular transformations browser (experimental)<a class="infolink" target="_blank" href="https://docs.ochem.eu/display/MAN/Molecular+Matched+Pairs" title="Click to read more about molecular matched pairs (aka MMPs) and their use in OCHEM"></a></h1>
				Molecular Matched Pairs Analysis (MMPA) is an experimental feature that is currently under active development.
				</td></tr>
			<tr>
				<td class="itunes-right">
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
				</td>
			</tr>
		</table>
		
		<script language="javascript">
			include.plugins('view');
			function MatchedPairsBrowser() {
				var self = this;
				Browser.call(this);
				this.url = "matchedpairs/transformationsList.do";
				this.itemElement = "mmp-transformation";
				this.view = new View({element: "template"});
				this.filters.setFromUrl();
				
				this.listenEvent("items_loaded", function(){attachZoomingClickHandlers();});
				this.listenEvent("items_loaded", function() {
					mmpSubset = self.filters.getValue("subset");
					$(".transformation-link").each(function(){
						$(this).attr("href", $(this).attr("href") + "&amp;" + window.location.search.substring(1));
					});
				});
				
			}
			
			$(function(){
				var browser = new MatchedPairsBrowser();
				browser.initialize();
			});
		</script>
		
		<script type="text/template" id="template">
			<table class="compact-item">
				<tr>
					<td width="300" class="transformation">	
						<img src="depiction.jsp?mmp_frag=[%=frag1Id %]"/> &#8594;
						<img src="depiction.jsp?mmp_frag=[%=frag2Id %]"/><br/><br/>
						<a href="matchedpairs/transformationProfile.do?id=[%=id %]" tab="Transformation profile">[transformation profile]</a>
					</td>
					<td>
						[%=pairsCount %] matched pairs <a href="matchedpairs/transformationProfile.do?id=[%=id %]" tab="Transformation profile" class="transformation-link">[view all]</a><br/>
						[%
							if (data.pair)
							{
								var pairsArr = array(data.pair);
								for (var i = 0; i != pairsArr.length; i++)
								{
									%]
									<div class="pair" mol1="[%=pairsArr[i].molecule1.id %]" mol2="[%=pairsArr[i].molecule2.id %]">
										<table>
											<tr>
												<td title="Click the molecule to enlarge it"><a class="mol" href="javascript:void()"><img src="depiction.jsp?mp2=[%=pairsArr[i].molecule1.id %]"/></a></td>
												<td title="Click the molecule to enlarge it"><a class="mol" href="#"><img src="depiction.jsp?mp2=[%=pairsArr[i].molecule2.id %]"/></a></td>
											</tr>
											<tr>
												<td><a href="molecule/profile.do?id=[%=pairsArr[i].molecule1.id %]" tab="Molecule profile" title="Open molecule profile">M[%=pairsArr[i].molecule1.id %]</a></td>
												<td><a href="molecule/profile.do?id=[%=pairsArr[i].molecule2.id %]" tab="Molecule profile" title="Open molecule profile">M[%=pairsArr[i].molecule2.id %]</a></td>
											</tr>
										</table>
										
										
									</div>
									[%
								}
							} 
						%]
					</td>
				</tr>
			
			</table>
		</script>
		
		<xsl:call-template name="pair-zoom"/>
		
	</xsl:template>
</xsl:stylesheet>
