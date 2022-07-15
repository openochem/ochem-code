<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />
	<xsl:include href="matched-pairs/pair-zoom.xslt" />

	<xsl:template name="content">
		<script type="text/javascript" src="js/commons/actionable.js" />
		<script type="text/javascript" src="js/commons/browser.js" />
		<xsl:variable name="mp2id" select="molecule-profile/molecule/mapping/@id"/>
		<style type="text/css">
			.itunes-right H2 {font-weight: normal; font-family: Arial; margin-top: 20px;}
			.itunes-right H1 { font-family: Georgia; border-bottom: 1px solid; color: #666; }
			.section {margin-left: 20px;}
			.itunes-right IMG {border: 1px solid gray;}
			TABLE.properties TD {padding-right: 10px;}
			#descriptors-info DIV {margin-top: 10px;}
			#descriptors-info TABLE TD {font-size: 80%; padding: 5px; background-color: #F5F5F5; border: 1px solid white;}
			#descriptors-info B {display: block; font-size: 80%;}
			
			.pair {float: left; margin: 10px 5px 10px 5px;}
			.pair A {color: #808080;}
			.pair {text-align: center; font-size: 10pt;}
		</style>
		<table height="100%" width="100%">
			<tr><td class="itunes-up silver">
				<h1><img src="img/icons/molecule.png"></img>Molecule profile</h1>
				A brief overview of the data, depictions and properties related to the molecule
				</td></tr>
			<tr>
				<td class="itunes-right">
					
					<table width="1000">
						<tr>
							<td valign="top" style="padding-right: 25px;">
								<img src="depiction.jsp?id={molecule-profile/molecule/@id}"/>
							</td>
							<td valign="top" width="100%">
								
								Internal identifier: <b>M<xsl:value-of select="$mp2id"/></b><br/>
								<h1>Data availability</h1>
								<xsl:choose>
									<xsl:when test="molecule-profile/recordsCount &gt; 0">
									<h2>This molecule has <a href="epbrowser/show.do?name=M{$mp2id}&amp;approval-status=all" tab="Molecule records"><xsl:value-of select="molecule-profile/recordsCount"/> records</a> for the following properties:</h2>
									<div class="section">
										<table class="properties">
										<xsl:for-each select="molecule-profile/properties">
											<tr><td><xsl:value-of select="@name"/></td><td><a href="epbrowser/show.do?property={@id}&amp;name=M{$mp2id}&amp;approval-status=all" tab="Molecule records"><xsl:value-of select="@timesUsed"/> records</a></td></tr>
										</xsl:for-each>
										</table>
									</div>
									</xsl:when>
									<xsl:otherwise>
										There is no experimental data available for this molecule.
									</xsl:otherwise>
								</xsl:choose>
								<xsl:if test="molecule-profile/baskets">
									<h2>This molecule is referenced from <xsl:value-of select="count(molecule-profile/baskets)"/> basket(s):</h2>
									<div class="section">
									<xsl:for-each select="molecule-profile/baskets"> 
										<a tab="Basket profile" href="basket/edit.do?id={@id}"><xsl:value-of select="@name"/></a><br/>
									</xsl:for-each>
									</div>
									
									<xsl:if test="molecule-profile/modelsCount &gt; 0 or molecule-profile/pendingModelsCount &gt; 0">
										<br/>
										This molecule was used for creation of <xsl:value-of select="molecule-profile/modelsCount"/> models
										<xsl:if test="molecule-profile/pendingModelsCount &gt; 0">
											(+<xsl:value-of select="molecule-profile/pendingModelsCount"/> pending models)	
										</xsl:if>
										
									</xsl:if>
								</xsl:if>
								<h1>General information</h1>
								<table class="properties">
									<tr><td>Name</td><td>=</td><td><xsl:value-of select="molecule-profile/canonicalName"/></td></tr>
									<tr><td>Molecular Weight</td><td>=</td><td><xsl:value-of select="molecule-profile/molecule/molWeight"/></td></tr>
									<tr><td>InChI Key</td><td>=</td><td><xsl:value-of select="molecule-profile/inchi"/></td></tr>
									<tr><td>SMILES</td><td>=</td><td><xsl:value-of select="molecule-profile/smiles"/></td></tr>
									<tr><td>Formula</td><td>=</td><td><xsl:value-of select="molecule-profile/formula"/></td></tr>
								</table>
								
								<xsl:if test="molecule-profile/molecule/mapping/mmpaIndexStatus = 3">
									<h1>Molecular matched pairs<a class="infolink" help="mmp-help"></a></h1>
									
									<div id="mmp-help" class="invisible">
										A <i>matched molecular pair</i> is a pair of molecules that have only a minor single-point difference (a fragment or a functional group of less than 10 atoms)
									</div>
									<a class="fb-button pairs-link" href="#" onclick="showMatchedPairs(); return false;">show matched pairs for this molecule</a>
									
									
									<div id="pairs-browser" class="invisible">
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
								</xsl:if>
								
								<xsl:if test="/model/session/user/ochem-labs = 'true'">
									<h1>Stored molecular descriptors</h1>
									
									<a class="fb-button" href="#" onclick="requestDescriptors(); return false;">make request</a>
									<div id="descriptors-info">
									</div>
								</xsl:if>
								
							</td>
						</tr>
					</table>
				</td>
			</tr>
		</table>
		
		<script type="text/template" id="template">
			[% var mol = molecule1.id == "<xsl:value-of select="$mp2id"/>" ? molecule2 : molecule1 %]
			<div class="pair" mol1="{$mp2id}" mol2="[%=mol.id%]">
				<a class="mol" href="#"><img src="depiction.jsp?mp2=[%=mol.id%]" title="Click to show the transformation details"/></a><br/>
				<a href="molecule/profile.do?id=[%=mol.id %]" tab="Molecule profile" title="Open molecule profile">M[%=mol.id %]</a>
				via 
				<a href="matchedpairs/transformationProfile.do?id=[%=transformation.id %]" tab="Transformation profile" title="Open transformation profile">TR[%=transformation.id %]</a>
			</div>
								
		</script>
		
		<script language="javascript">
		
			include.plugins('view');
			
			function MatchedPairsBrowser() {
				Browser.call(this);
				this.url = "matchedpairs/getPairs.do";
				this.filters.setValue("molecule", <xsl:value-of select="$mp2id"/>);
				this.itemElement = "mmpair";
				this.view = new View({element: "template"});
				this.filters.setFromUrl();
				this.filters.values.pagesize = 10;
				
				this.listenEvent("items_loaded", function(){attachZoomingClickHandlers()});
			}
			
			$(function(){
				pairsBrowser = new MatchedPairsBrowser();
				//browser.initialize();
			});
			
			function showMatchedPairs() {
				$("#pairs-browser").removeClass("invisible");
				$(".pairs-link").addClass("invisible");
				pairsBrowser.initialize();
			}
			
			
		
			function requestDescriptors()
			{
				$("#descriptors-info").html('<img src="img/long_green.gif"/>');
				var ajax = new QSPR.Ajax();
				ajax.send({
					url: "molecule/getDescriptors.do",
					dataType: "json",
					data: "mp2=<xsl:value-of select="$mp2id"/>&amp;out=json",
					success: function(obj) {
						$("#descriptors-info").html("");
						if (!obj.others.cacheEntry)
							$("#descriptors-info").html("There are no stored descriptors for this molecule");
							
						var entries = array(obj.others.cacheEntry);
						for (var i in entries)
						{
							var entry = entries[i];
							console.log(entry);
							
							if (entry.error)
								$("#descriptors-info").append("Stored failure: " + entry.error);
							else
							{
								var names = array(entry.names);
								var values = array(entry.values);
								
								var d = $("<div/>");
								//d.css({overflow: "auto", width: "400px"});
								var t = $("<table/>");
								var header = $("<tr/>");
								var valuesRow = $("<tr/>");
								for (var i = 0; i &lt; names.length; i++)
								{
									header.append("<td>" + names[i] + "</td>");
									valuesRow.append("<td>" + values[i] + "</td>")
								}
								t.append(header);
								t.append(valuesRow);
								d.append('<b></b>');
								d.find("b").html(entry.config.type + (entry.user ? " (uploaded by " + entry.user + ")" : " (public storage)") );	
								d.append(t);
								$("#descriptors-info").append(d);
								$("html, body").animate({ scrollTop: $(document).height() }, 500);
							}
								
						}
					},
					error: function(msg){
						$("#descriptors-info").html("Its a pity, but we could not retrieve the descriptors at the moment.");
						console.log(msg);
					}
				});	
			}
		</script>
		
		<xsl:call-template name="pair-zoom"/>
		
	</xsl:template>
</xsl:stylesheet>