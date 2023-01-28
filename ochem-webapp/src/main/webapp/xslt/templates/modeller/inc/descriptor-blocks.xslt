<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:template name="descriptor-blocks">
		<style type="text/css">
			.params {margin-left: 20px; background-color: #FFC; padding: 10px;}
			.params LABEL {width: 300px; float: left;}
			.params INPUT[type="checkbox"] {margin-right: 4px; float: left; width: 20px;}
			INPUT[type="checkbox"] {margin-right: 2px;}
			.width40 {width: 40px !important;}
			.params A {font-size: 80%; margin-bottom: 5px;}
			#alvaDesc-params,#RDKIT-params  {width: 650px; overflow: auto;}
			#GSFrag-params {width: 340px; overflow: auto;}
			.mappings {border-left: 1px solid black; padding: 5px 20px 5px 20px; margin-left: 5px;}
			#MOPAC-params {width: 650px; overflow: auto;}
			#MOPAC2016-params {width: 650px; overflow: auto;}
			#EState-atomic-params {width: 650px; overflow: auto;}
			#EState-params {width: 650px; overflow: auto;}
			#CDK2-params {width: 650px; overflow: auto;}
			#PaDEL2"-params {width: 650px; overflow: auto;}
			#MORDRED-params {width: 650px; overflow: auto;}
			#Fragmentor-params, #Fragmentor-2011-params, #QNPR-params, #Random-params {margin-left: 20px; background-color: #FFC; padding: 10px; width: 650px; }
			#SilicosItScaffold-params {width: 340px; overflow: auto;}
			.faded {opacity: 0.3;}
			
			.two-columns {
				border-collapse: separate;
    			border-spacing: 3px;
			}
			.two-columns TD {vertical-align: top; padding-right: 35px; min-width: 600px;
				border-radius: 5px;
				border: 1px solid #BBB;
				-moz-border-radius: 5px;
				padding: 20px;
			}
			
			.two-columns TD H3 {
				margin-top: 0px;
			}
			
		</style>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/browser.js?ver=1.8.6"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script language="javascript" src="js/commons/dynamic-select.js"></script>
		<script language="javascript" src="js/blocks/select-model.js"></script>	
		<title>Model builder - Select descriptors</title>
		<h1><doc term="Molecular+descriptors">Select the molecular descriptors</doc></h1>
		
	<xsl:choose>
		<xsl:when test="//param[@key='availability-unknown']">
			<font color="red">Warning: Could not retrieve the list of available descriptor types. Some listed descriptors might be not available.</font> 
		</xsl:when>
		<xsl:otherwise>
			<script language="javascript">
			// Automatically determine the list of supported descriptors
			$(document).ready(function(){
				<xsl:for-each select="//param[@key='available-task']">
					$("input[name='<xsl:value-of select="."/>']").attr("available", "1");
				</xsl:for-each>
				$("input[descriptor]").filter("[available!=1]")
					.attr("disabled", "1")
					.removeAttr("checked").change()
					.parent("div").addClass("faded").attr("title", "This descriptor type is not supported by your installation").append('<i style="color:red">Not supported by your installation</i>');
					;
			});
			
		</script>
		</xsl:otherwise>
	</xsl:choose>
		
	<table class="two-columns">
		<tr>
			<td>
				<h3>Recommended descriptor types (2D)</h3>
		
		<div><doc term="OEState" hide="true"><input type="checkbox" name="OEstate" checked="checked" details="1" descriptor="1"/> OEState</doc></div>
		<div class="params invisible" id="OEstate-params">
			<input type="checkbox" name="oestate_bond" checked="true"/><label>Bonds Indices</label><br/>
			<input type="checkbox" name="oestate_count"/><label>Counts only</label><br/>
		</div>
		
		<div><doc term="ALogPS" hide="true"><input type="checkbox" checked="checked" name="ALogPS" descriptor="1"/> ALogPS (2)</doc></div>

		<div><doc term="MOLD2" hide="true"><input type="checkbox" name="MOLD2" descriptor="1"/> Mold2 (777)</doc></div>

		<div><doc term="CDDD" hide="true"><input type="checkbox" name="CDDD" descriptor="1"/> CDDD</doc></div>

		<div><doc term="JPlogP" hide="true"><input type="checkbox" name="JPlogP" descriptor="1"/> JPlogP</doc></div>

		<div><doc term="SIRMS" hide="true"><input type="checkbox" name="SIRMS" details="1" descriptor="1"/> SIRMS</doc></div>
		<div class="params invisible" id="SIRMS-params">
			<div>
				Atom labeling:<br/>
<!-- Not yet implemented with CDK 
				<input type="checkbox" name="sirms_charge"  checked="true"/><label>Atom charges</label><br/>
				<input type="checkbox" name="sirms_logp"  checked="true"/><label>Lipophilicity</label><br/>
				<input type="checkbox" name="sirms_refractivity"  checked="true"/><label>Refractivity</label><br/>
 -->				<input type="checkbox" name="sirms_hb"  checked="true"/><label>Hydrogen acceptors and donors</label><br/>
				<input type="checkbox" name="sirms_elm"/><label>Atom elements</label><br/>
				<input type="checkbox" name="sirms_none"/><label>Only topogy</label><br/>
				<br/>
			</div>
			<input type="checkbox" name="sirms_noh" checked="true"/><label>Exclude all H-atoms</label><br/>
			<input type="checkbox" name="sirms_mixtures"/><label>Treat pure compounds as mixtures</label><br/>
			Fragments from <input name="sirms-min-frag" type="text" value="1"/> to <input name="sirms-max-frag" type="text" value="4"/><br/>
		</div>
				
		<div><doc term="ISIDA+Fragments" hide="true"><input type="checkbox" name="Fragmentor" details="1" descriptor="1"/> ISIDA fragments</doc></div>
		<div class="invisible" id="Fragmentor-params">	
			Fragments from <input name="min-length" type="text" value="2"/> to <input name="max-length" type="text" value="4"/><br/>
			Type of fragments:
			<select name="fragment-type">
				<option value="1">Sequences of atoms only</option>
				<option value="2">Sequences of bonds only</option>
				<option value="3" selected="true">Sequences of atoms and bonds</option>
			</select>
		</div>

		<div><doc term="MAP4" hide="true"><input type="checkbox" name="MAP4" details="1" descriptor="1"/> The in Hashed Atom Pair fingerprint (MAP4)</doc></div>
		<div class="params invisible" id="MAP4-params">
			<select name="dimensions">
				<option value = "512">512</option>
				<option value = "1024"  selected="1">1024</option>
				<option value = "2048">2048</option>
			</select>
			MAP4 Radius: 
			<select name="radius">
				<option value = "2" selected="1">2</option>
				<option value = "4">4</option>
				<option value = "3">3</option>
			</select><br/>
		</div>
	
		<div><doc term="GSFrag" hide="true"><input type="checkbox" name="GSFrag" details="1" descriptor="1"/> GSFragment (1138)</doc></div>
		<div class="params invisible" id="GSFrag-params">
			<input type="checkbox" name="gsfrag-normal" checked="true"/><label>GSFrag (247)</label>
			<input type="checkbox" name="gsfragl" checked="true"/><label>GSFragL (891)</label><br/>
		</div>

		<div><doc term="QNPR" hide="true"><input type="checkbox" name="QNPR" details="1" descriptor="1"/> QNPR</doc></div>
		<div class="invisible" id="QNPR-params">
			Fragments from <input name="qnpr-min-frag" type="text" value="1"/> to <input name="qnpr-max-frag" type="text" value="3"/><br/>
		</div>

		<div><doc term="ToxAlerts%3A+Database+of+structural+alerts" hide="true"><input type="checkbox" name="StructuralAlerts" details="1" descriptor="1"/> Structural alerts (ToxAlerts and Functional Groups)</doc></div>
		<div class="params invisible" id="StructuralAlerts-params">
			<table cellpadding="3" cellspacing="3">
				<tr><td colspan="2"><br/><b>Select alerts:</b></td></tr>
						<tr><td>
						Publication&#160;</td>
							<td>
								<select name="article" filter="1">
									<option value="">All articles</option>
									<xsl:for-each select="//available-alerts/articles">
										<option value="{@id}"><xsl:value-of select="publication-date/year"/>&#160;<xsl:value-of select="authors/author[1]/LastName"/></option>
									</xsl:for-each>
								</select>
							</td></tr>
						<tr><td>Endpoint&#160;</td>
							<td><select name="property" filter="1">
									<option value="">All endpoints</option>
									<xsl:for-each select="//available-alerts/endpoints">
										<option value="{@id}"><xsl:value-of select="@name"/></option>
									</xsl:for-each>
								</select>
							</td></tr>
						<tr class="invisible" id="alert-filter"><td>A single alert selected:&#160;</td>
							<td>
								<b><span id="alert-name"></span></b>
								<input type="hidden" name="alert-id"/>
							</td>
						</tr>
						<tr><td colspan="2"><input type="checkbox" name="approved-only"/>&#160;Only approved alerts</td></tr>
						<xsl:if test="//available-alerts/selectionSize &gt; 0">
						<tr id="selected-filter"><td colspan="2"><input type="checkbox" name="selected-only"/>&#160;Only <b><xsl:value-of select="//available-alerts/selectionSize"/> selected alerts</b></td>
						</tr>
						</xsl:if>
					</table>
		</div>

		<h3><br/>Recommended descriptor types (3D)</h3>

		<xsl:if test="/model/@inhouse = 'false'">
		<div><doc term="alvaDesc" hide="true"><input type="checkbox" name="alvaDesc" details="1" descriptor="1"/> alvaDesc v.2.0.4 (5666/3D)</doc></div>
		<div class="params invisible" id="alvaDesc-params">
			<a href="javascript:selectBlocks(true, 'alvaDesc-params')">[select all]</a><a href="javascript:selectBlocks(false, 'alvaDesc-params')">[select none]</a>
			<a href="javascript:selectBlocksThree(true, 'alvaDesc-params')">[select 3D]</a><a href="javascript:selectBlocksThree(false, 'alvaDesc-params')">[unselect 3D]</a>
			<br/>
			<input type="checkbox" name="alva1" checked="true"/><label>Constitutional descriptors (50)</label>
			<input type="checkbox" name="alva2" checked="true"/><label>Ring descriptors (35)</label>
			<input type="checkbox" name="alva3" checked="true"/><label>Topological indices (79)</label>
			<input type="checkbox" name="alva4" checked="true"/><label>Walk and path counts (46)</label>
			<input type="checkbox" name="alva5" checked="true"/><label>Connectivity indices (37)</label>
			<input type="checkbox" name="alva6" checked="true"/><label>Information indices (51)</label>
			<input type="checkbox" name="alva7" checked="true"/><label>2D matrix-based descriptors (608)</label>
			<input type="checkbox" name="alva8" checked="true"/><label>2D autocorrelations (213)</label>
			<input type="checkbox" name="alva9" checked="true"/><label>Burden eigenvalues (96)</label>
			<input type="checkbox" name="alva10" checked="true"/><label>P_VSA-like descriptors (69)</label>
			<input type="checkbox" name="alva11" checked="true"/><label>ETA indices (40)</label>
			<input type="checkbox" name="alva12" checked="true"/><label>Edge adjacency indices (324)</label>
			<input type="checkbox" name="alva13" checked="true" id="3D"/><label>Geometrical descriptors (3D, 38)</label>
			<input type="checkbox" name="alva14" checked="true" id="3D"/><label>3D matrix-based descriptors (3D, 132)</label>
			<input type="checkbox" name="alva15" checked="true" id="3D"/><label>3D autocorrelations (3D, 80)</label>
			<input type="checkbox" name="alva16" checked="true" id="3D"/><label>RDF descriptors (3D, 210)</label>
			<input type="checkbox" name="alva17" checked="true" id="3D"/><label>3D-MoRSE descriptors (3D, 224)</label>
			<input type="checkbox" name="alva18" checked="true" id="3D"/><label>WHIM descriptors (3D, 114)</label>
			<input type="checkbox" name="alva19" checked="true" id="3D"/><label>GETAWAY descriptors (3D, 273)</label>
			<input type="checkbox" name="alva20" checked="true" id="3D"/><label>Randic molecular profiles (3D, 41)</label>
			<input type="checkbox" name="alva21" checked="true" id="3D"/><label>Functional group counts (3D, 154)</label>
			<input type="checkbox" name="alva22" checked="true"/><label>Atom-centred fragments (115)</label>
			<input type="checkbox" name="alva23" checked="true"/><label>Atom-type E-state indices (346)</label>
			<input type="checkbox" name="alva24" checked="true"/><label>Pharmacophore descriptors (165)</label>
			<input type="checkbox" name="alva25" checked="true"/><label>2D Atom Pairs (1596)</label>
			<input type="checkbox" name="alva26" checked="true" id="3D"/><label>3D Atom Pairs (3D, 36)</label>
			<input type="checkbox" name="alva27" checked="true" id="3D"/><label>Charge descriptors (3D, 15)</label>
			<input type="checkbox" name="alva28" checked="true" id="3D"/><label>Molecular properties (3D, 27)</label>
			<input type="checkbox" name="alva29" checked="true"/><label>Drug-like indices (30)</label>
			<input type="checkbox" name="alva30" checked="true" id="3D"/><label>CATS 3D (3D, 300)</label>
			<input type="checkbox" name="alva31" checked="true" id="3D"/><label>WHALES (3D, 33)</label>
			<input type="checkbox" name="alva32" checked="true"/><label>MDE (19)</label>
			<input type="checkbox" name="alva33" checked="true"/><label>Chirality (70)</label>
		</div>
		</xsl:if>


		<div><doc term="CDK2" hide="true"><input type="checkbox" name="CDK2" details="1" descriptor="1"/> CDK 2.7.1 descriptors (256/3D)</doc></div>
		<div class="params invisible" id="CDK2-params">
			<a href="javascript:selectBlocks(true, 'CDK2-params')">[select all]</a><a href="javascript:selectBlocks(false, 'CDK2-params')">[select none]</a><br/>
			<input type="checkbox" name="cdk21" checked="true"/><label>constitutional descriptors (38)*</label>
			<input type="checkbox" name="cdk22" checked="true"/><label>topological descriptors (195)</label><br/>
			<input type="checkbox" name="cdk23" checked="true"/><label>geometric descriptors (49)</label>
			<input type="checkbox" name="cdk24" checked="true"/><label>electronic descriptor (23)</label>
			<input type="checkbox" name="cdk25" checked="true"/><label>hybrid descriptor (34)</label>
			<br/><br/>*Some descriptors belong to two or more packages.
			<br/>Time to calculate descriptors per molecule (max: 99 minutes)<br/>
			<label>Time (in minutes):</label>	<input class="width40" type="text" name="cdk2Timeout" value="10"/>
		</div>


		<div><doc term="RDKIT" hide="true"><input type="checkbox" name="RDKIT" details="1" descriptor="1"/> RDKit descriptors <i>(3D)</i></doc></div>
		<div class="params invisible" id="RDKIT-params">
			<a href="javascript:selectBlocks(true, 'RDKIT-params')">[select all]</a><a href="javascript:selectBlocks(false, 'RDKIT-params')">[select none]</a>
			<a href="javascript:selectBlocksThree(true, 'RDKIT-params')">[select 3D]</a><a href="javascript:selectBlocksThree(false, 'RDKIT-params')">[unselect 3D]</a>
			<br/>
			<input type="checkbox" name="rdkitb1" checked="true"/><label>Scalars (53)</label>
			<input type="checkbox" name="rdkitb2" checked="true"/><label>Scalars secondary (61)</label>
			<input type="checkbox" name="rdkitb3" checked="true"/><label>2D auto-correlations (192)</label>
			<input type="checkbox" name="rdkitb4" checked="true" id="3D"/><label>3D auto-correlations  <i>(3D)</i> (80)</label>
			<input type="checkbox" name="rdkitb5" checked="true"/><label>Topological (see Topological below)</label>
			<input type="checkbox" name="rdkitb6" checked="true" id="3D"/><label>GETAWAY <i>(3D)</i> (272)</label>
			<input type="checkbox" name="rdkitb7" checked="true" id="3D"/><label>Morse <i>(3D)</i> (224)</label>
			<input type="checkbox" name="rdkitb8" checked="true" id="3D"/><label>RDF <i>(3D)</i> (210)</label>
			<input type="checkbox" name="rdkitb9" checked="true" id="3D"/><label>WHIM <i>(3D)</i> (114)</label>
			<input type="checkbox" name="rdkitb10" checked="true"/><label>Morgan (ECFP-like) (see Morgan below)</label>
			<input type="checkbox" name="rdkitb11" checked="true"/><label>MACCS keys (166)</label>
			<input type="checkbox" name="rdkitb12"/><label>Atom pairs</label>
			<input type="checkbox" name="rdkitb13"/><label>Sheridan BT pairs</label>
			<input type="checkbox" name="rdkitb14"/><label>Sheridan BP pairs</label>
			<input type="checkbox" name="rdkitb15" checked="true"/><label>Topological Torsions</label>
			<input type="checkbox" name="rdkitb16" checked="true"/><label>Synthesability score (1)</label>
			<input type="checkbox" name="rdkitb17"/><label>Avalon descriptors (see Avalon below)</label>
			<input type="checkbox" name="rdkitb18"/><label>Default RDKIT fingerprint (2048)</label>
			<br/>&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;--------------------------  Additional parameters --------------------------  <br/>

			WHIM threshold:	<input class="width40" type="text" name="whim_thresh"  value="0.1"/>
			&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;
			Topological bits: 
			<select name="topological_nbits">
				<option value = "512">512</option>
				<option value = "1024"  selected="1">1024</option>
				<option value = "2048">2048</option>
				<option value = "4096">4096</option>
				<option value = "8192">8192</option>
			</select>
			<br/>&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;---------------------  Parameters of Morgan descriptors ---------------------  <br/>
			ECFP Bits: 
			<select name="morgan_nbits">
				<option value = "512">512</option>
				<option value = "1024"  selected="1">1024</option>
				<option value = "2048">2048</option>
				<option value = "4096">4096</option>
				<option value = "8192">8192</option>
			</select>
			&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;
			ECFP Radius: 
			<select name="morgan_radius">
				<option value = "2" selected="1">2</option>
				<option value = "3">3</option>
				<option value = "4">4</option>
			</select><br/>	
			<input type="checkbox" name="morgan_fcfp"/><label>Calculate pharmacophore features</label>
			<input type="checkbox" name="morgan_counts"/><label>Use counts </label>	
			<br/>&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;---------------------  Parameters of Avalon descriptors ---------------------  <br/>
			Avalon Bits: 
			<select name="avalon_nbits">
				<option value = "512">512</option>
				<option value = "1024"  selected="1">1024</option>
				<option value = "2048">2048</option>
				<option value = "4096">4096</option>
				<option value = "8192">8192</option>
			</select>
			<input type="checkbox" name="avalon_counts"/><label>Use counts </label><br/>	
	
		</div>
	
		<div><doc term="MORDRED" hide="true"><input type="checkbox" name="MORDRED" details="1" descriptor="1"/> MORDRED descriptors (1826/3D)</doc></div>
		<div class="params invisible" id="MORDRED-params">
			<select name="mordred">
				<option value = "D2" selected="1">2D only (1613)</option>
				<option value = "D3">All (1826)</option>
			</select><br/>
		</div>

		<div><doc term="MOPAC2016-derived+descriptors" hide="true"><input type="checkbox" name="MOPAC2016" descriptor="1"/> MOPAC2016 descriptors <i>(35/3D)</i></doc></div>

		<div><doc term="KrakenX" hide="true"><input type="checkbox" name="KrakenX" descriptor="1"/> KrakenX descriptors (MOPAC2016 derived)<i>(124/3D)</i></doc></div>

		<div><doc term="PyDescriptor" hide="true"><input type="checkbox" name="PyDescriptor" details="1" descriptor="1"/> PyDescriptor descriptors <i>(16251/3D)</i></doc></div>
			<div class="params invisible" id="PyDescriptor-params">
			<label>Time (in minutes):</label>	<input class="width40" type="text" name="pydescriptorTimeout"  value="10"/>
		</div>
	
			<div><doc term="MERA+descriptors" hide="true"><input type="checkbox" name="Mera" descriptor="1"/> MERA descriptors <i>(529/3D)</i></doc></div>
		
		<div><doc term="MERSY+descriptors" hide="true"><input type="checkbox" name="Mersy" descriptor="1"/> MERSY descriptors <i>(42/3D)</i></doc></div>

		<div><doc term="Inductive Descriptors" hide="true"><input type="checkbox" name="InductiveDescriptors" details="1" descriptor="1"/> 'Inductive' descriptors <i>(54/3D)</i></doc></div>
		<div class="params invisible" id="InductiveDescriptors-params">
			<a href="javascript:selectBlocks(true, 'InductiveDescriptors-params')">[select all]</a><a href="javascript:selectBlocks(false, 'InductiveDescriptors-params')">[select none]</a><br/>
			<input type="checkbox" name="inductiveMElectoneg" checked="true"/><label>electronegativity-based</label><br/>
			<input type="checkbox" name="inductiveMHardness"  checked="true"/><label>hardness-based</label><br/>
			<input type="checkbox" name="inductiveMSoftness"  checked="true"/><label>softness-based</label><br/>
			<input type="checkbox" name="inductiveMCharge"    checked="true"/><label>charge-based</label><br/>
			<input type="checkbox" name="inductiveMSigma"     checked="true"/><label>inductive sigma-based</label><br/>
			<input type="checkbox" name="inductiveMRs"        checked="true"/><label>steric Rs-based</label><br/>
		</div>
		
				
		<div><doc term="Spectrophores" hide="true"><input type="checkbox" name="Spectrophores" details="1" descriptor="1"/> Spectrophores (144/3D)</doc></div>
		<div class="params invisible" id="Spectrophores-params">
			Required accuracy: 
			<select name="spectrophores-accuracy">
				<option value = "1">1</option>
				<option value = "2">2</option>
				<option value = "5">5</option>
				<option value = "10">10</option>
				<option value = "15">15</option>
				<option value = "20" selected="1">20</option>
				<option value = "30">30</option>
				<option value = "36">36</option>
				<option value = "45">45</option>
				<option value = "60">60</option>
			</select><br/>
			Stereospecificity: 
			<select name="spectrophores-cage">
				<option value = "0" selected="1">No (48)</option>
				<option value = "1">Unique (72)</option>
				<option value = "2">Mirror (72)</option>
				<option value = "3">Unique+Mirror (144)</option>
			</select><br/>
			Resolution in Angstrom: 
			<select name="spectrophores-resolution">
				<option value = "0.1"  selected="1">0.1</option>
				<option value = "0.2">0.2</option>
				<option value = "0.5">0.5</option>
				<option value = "1">1</option>
				<option value = "2">2</option>
				<option value = "3">3</option>
				<option value = "5">5</option>
				<option value = "10">10</option>
				<option value = "20">20</option>
			</select><br/>
		</div>
		</td>
						
		<td>
			<h3>Predictions by OCHEM's featured models<a help="models-help" class="infolink"></a></h3>
			<xsl:for-each select="//others/model">
				<input type="checkbox" name="model-id" value="{@id}"/><label><xsl:value-of select="@name"/></label><br/>
			</xsl:for-each>
			<br/>
			<input type="checkbox" name="other-models" details="1"/>
			<label>Outputs of other OCHEM models</label>
				<div class="params invisible" id="other-models-params">
				<div><a action="add">[Add models]</a></div>
				<div id="ModelsBrowser"></div>
				</div>
					
			<div id="models-help" class="invisible">
				OCHEM allows to to use model predictions as molecular descriptors.<br/>
				The approach of using outputs of models as an input for a new model is sometimes referred to as <i>feature nets</i>.<br/><br/>
				You can use the featured models from the list below as well your own OCHEM models.
			</div>
				
			<br/><br/>
		<h3>Obsolete/Additional descriptor types</h3>
		
		<div><doc term="MOPAC-derived+descriptors" hide="true"><input type="checkbox" name="MOPAC" descriptor="1"/> MOPAC 7.1 descriptors <i>(25/3D)</i></doc></div>
		
			</td>
		</tr>
	</table>
				
		<xsl:if test="/model/others/desc-config-entry">
			<h3>Stored descriptors (<a href="descriptorsstorage/show.do" tab="Descriptors storage">explore the descriptors storage</a>)</h3>
			<xsl:for-each select="/model/others/desc-config-entry">
				<div><input type="checkbox" name="Stored-{objectID}" details="1"/><xsl:value-of select="type"/></div>
			</xsl:for-each>
			<br/>
		</xsl:if>
	<h3>Special descriptors (scaffolds, fingerprints):</h3>
		
		<div><doc term="Silicos-It+scaffolds" hide="true"><input type="checkbox" name="SilicosItScaffold" details="1" descriptor="1"/> Silicos-It Scaffolds</doc></div>
			<div class="params invisible" id="SilicosItScaffold-params">
			<input type="checkbox" name="silicos_scaffolds" checked="true"/><label>Bemis-Murcko scaffolds</label>
			<input type="checkbox" name="silicos_frameworks" checked="true"/><label>Bemis-Murcko frameworks</label><br/>
			<input type="checkbox" name="silicos_oprea" checked="true"/><label>Oprea frameworks</label><br/>
			<input type="checkbox" name="silicos_schuffenhauer" checked="true"/><label>Schuffenhauer scaffolds</label><br/>
		</div>

		<div><doc term="ECFP+fingerprints" hide="true"><input type="checkbox" name="ECFP" details="1" descriptor="1"/> ECFP Fingerprints</doc></div>
			<div class="params invisible" id="ECFP-params">
			ECFP Bits: 
			<select name="ecfp_nbits">
				<option value = "512">512</option>
				<option value = "1024"  selected="1">1024</option>
				<option value = "2048">2048</option>
			</select>
			ECFP Radius: 
			<select name="ecfp_radius">
				<option value = "2" selected="1">2</option>
				<option value = "3">3</option>
				<option value = "4">4</option>
			</select><br/>	
			<input type="checkbox" name="ecfp_fcfp"/><label>Calculate functional groups</label>
			<input type="checkbox" name="ecfp_counts"/><label>Use counts </label><br/>
		</div>		
		
		<div><doc term="MolPrint" hide="true"><input type="checkbox" name="MolPrint" details="1" descriptor="1"/> MolPrint Fingerprints </doc></div>
			<div class="params invisible" id="MolPrint-params">
			Depth of the fingerprint vector: <input class="width40" type="text" name="molprint-depth"  value="2"/><br/>
		</div>
				
		<xsl:if test="//param[@key='use-mixtures-descriptors'] = 'true'">
		<hr/>Mixture processing options: 
			<select name="mixtures">
				<option value="fraction">Descriptors calculated for the whole mixture plus molar fractions are added as descriptors</option>
				<option value="concatenate">Descriptors are concatenated (only binary mixtures)</option>			
				<option value="wsum">Descriptors are averaged (weighted by molar fraction; any mixtures)</option>			
				<option value="wsumdiff">Sum and absolute differences of weighted by molar fraction descriptors (binary mixtures)</option>
			</select>
		</xsl:if>
		
		<hr/>

		<div id="Conditions" class="invisible">
			<br/>Conditions of experiments<br/>
			<div class="params" id="Conditions-Browser">
			</div>
		</div>
		
		<h3>Under development: can change anytime and backward compatibility is not guaranteed. Use at your own risk!</h3>
	
		<div><doc term="PaDEL2" hide="true"><input type="checkbox" name="PaDEL2" details="1" descriptor="1"/> PaDEL2 descriptors (3D)</doc></div>
		<div class="params invisible" id="PaDEL2-params">
			<a href="javascript:selectBlocks(true, 'PaDEL2-params')">[select all]</a><a href="javascript:selectBlocks(false, 'CDK2-params')">[select none]</a><br/>
			<input type="checkbox" name="padel1" checked="true"/><label>2D descriptors</label>
			<input type="checkbox" name="padel2" checked="true"/><label>3D descriptors (3D)</label><br/>
			<input type="checkbox" name="padel3" checked="true"/><label>fingerprints</label>
		<br/><br/>
		</div>
	
		<div><doc term="SIGMA profiles" hide="true"><input type="checkbox" name="SIGMA"  details="1" descriptor="1"/> COSMO-RS SIGMA profiles (MOPAC2016 derived, <i>(3D)</i>)</doc></div>
		<div class="params invisible" id="SIGMA-params">
		Select SIGMA profile width:
			<select name="sigmawidth">
				<option value = "0.005">0.005 (11 bins)</option>
				<option value = "0.001"  selected="1">0.001 (51 bins)</option>
				<option value = "0.0005">0.0005 (101 bins)</option>
			</select>
		</div>
				
		<div><doc term="EPA" hide="true"><input type="checkbox" name="EPA" details="1"  descriptor="1"/> EPA (T.E.S.T) descriptors (797)</doc></div>
		<div class="params invisible" id="EPA-params">
			<input type="checkbox" name="epa3D"/><label>Use also 3D descriptors (total 1141)</label><br/>
		</div>
	
		<br/>
		<script language="javascript">
		

			function selectBlocksThree(flag, param)
			{
				var boxes = $("#"+param).find("input[id='3D']");
				if (flag)
					boxes.attr("checked", "checked");
				else
					boxes.removeAttr("checked");
			}
						
			function selectBlocks(flag, param)
			{
				var boxes = $("#"+param).find("input[type='checkbox']");
				if (flag)
					boxes.attr("checked", "checked");
				else
					boxes.removeAttr("checked");
			}

		var condBrowser = null;
		<xsl:if test="/model/model/training-set/@id">
			condBrowser = new Browser();
			condBrowser.url = "basket/listConditions.do";
			condBrowser.filters.setValue("id", <xsl:value-of select="/model/model/training-set/@id"/>);
			condBrowser.container = "Conditions-Browser";
			condBrowser.scope = "#Conditions-Browser";
			condBrowser.itemElement = "property";
			condBrowser.itemTemplate = "js/templates/models/conscriptor.ejs";
			condBrowser.listenEvent("items_loaded", function(){
				if (condBrowser.pager.totalNum &gt; 0)
					$("#Conditions").removeClass("invisible");
			});
			
			condBrowser.onItemDrawn = function()
			{
				this.currentBlock.find(".regression").setClass("invisible", this.currentEntity.type != 0);
				this.currentBlock.find(".classification").setClass("invisible", this.currentEntity.type != 1);
				var dom = this.currentBlock.get(0);
				
				if (this.currentEntity.type == 1)
				{
					var optionsSelect = new DynamicSelect('option-' + this.currentEntity.id, 'properties/listoptions.do', 'option', this.currentBlock);
					optionsSelect.update("id=" + this.currentEntity.id + "&amp;basket=<xsl:value-of select="/model/model/training-set/@id"/>");
				}
				else
				{
					var unitsSelect = new UnitSelect('unit-' + this.currentEntity.id, 'unit/list.do', 'unit', this.currentBlock);
					unitsSelect.update("category=" + this.currentEntity.unitCategory.id, this.currentEntity.defaultUnit.id);
				}
				
				
				this.currentBlock.find("input[name=condition]").change(function(){
					$(this).siblings(".details").setClass("invisible", !$(this).is(":checked"));
				});
			}
			
			condBrowser.doDetails = function()
			{
				this.currentBlock.find(".mappings").removeClass("invisible");
			}
			
			condBrowser.doAddmapping = function()
			{
				var div = $("<div> = </div>");
				var selectRule = "id=" + this.currentEntity.id + "&amp;basket=<xsl:value-of select="/model/model/training-set/@id"/>";
				div.prepend(sel1 = $("<select></select>"));
				div.append(sel2 = $("<select></select>"));
				sel1.attr("name", "mapping-" + this.currentEntity.id + "-1");
				sel2.attr("name", "mapping-" + this.currentEntity.id + "-2");
				this.currentBlock.find(".mappings").append(div);
				var dSel1 = new DynamicSelect(sel1, 'properties/listoptions.do', 'option', this.currentBlock);
				dSel1.update(selectRule);
				var dSel2 = new DynamicSelect(sel2, 'properties/listoptions.do', 'option', this.currentBlock);
				dSel2.update(selectRule);
			}
			
			condBrowser.listenEvent("items_loaded", function(){
				condBrowser.mainDiv.find("input[name=condition]").change();
			});
			</xsl:if>
			
			$(document).ready(function()
			{
				var updateVisibility = function(){
					var details = $("#"+$(this).attr('name')+"-params");
				if (details.length > 0)
					if ($(this).is(':checked'))
						details.removeClass("invisible");
					else
						details.addClass("invisible");
						
				};
				$("input[details]").click(updateVisibility);
				$("input[details]").each(function() {updateVisibility.call(this);});
				
				$("input[type='checkbox']").click(function() {
					if ($("input:checkbox:checked:visible").length == -1)
						$("input[type='submit']").attr("disabled", "disabled");
					else
						$("input[type='submit']").removeAttr("disabled");
				});

				if (condBrowser)
					condBrowser.initialize();
			});
			
			function selectExpValuesBasket()
			{
				var win = openTab("Select a basket with experimental values", "basket/show.do");
				win.callback = function(basket)
				{
					$("#exp-values-basket").html(basket.name);
					$("[name='exp-values-basket-id']").val(basket.id);
					
					$("#exp-values-properties").html("Available properties:");
					$("#exp-values-properties").append("<br/>");
					var props = basket.properties;
					if (!props.length)
						props = new Array(props);
					for (var i = 0; i &lt; props.length; i++)
					{
						var prop = props[i];
						var div = $("<div/>");
						div.html('<input type="checkbox" value="' + prop.id + '" name="exp-values-property-id"/>');
						div.append(prop.name);
						$("#exp-values-properties").append(div);
					}
					$("#exp-values-properties input[type=checkbox]").setChecked(true);
					win.closeTab();
				}		
			}
			
			function deleteEtmTemplate(link)
			{
				$(link).parents("span").remove();
			}
			
			$(function(){
				$(document).on("change", "input[name='model-id']", function(){
					$("#scenario").setClass("invisible", $("input[name='model-id']:checked").length == 0);
				});
				$("input[name='model-id']").change();
			});
			
		</script>
		
		<xsl:if test="//attachment/configuration/descriptors">
			<script language="javascript">
				var setDragonBlocks = function(scope, prefix, value, totalBlocks)
				{
					if (value == '')
						return;
					value = 1 * value;
					var block = $(scope);
					for (var i = 0; i &lt; totalBlocks; i++)
					{
						console.log(value);
						block.find("input[name="+prefix+(i + 1) + "]").setChecked(value % 2 == 1);
						value &gt;&gt;= 1;
					}
				}
				
				<xsl:if test="//attachment/configuration/descriptors/mixtures">
					var mixtures = '<xsl:value-of select="//attachment/configuration/descriptors/mixtures"/>';
					mixtures = mixtures.toLowerCase();
					setValue("mixtures", mixtures);
				</xsl:if>

				<xsl:if test="//attachment/configuration/descriptors/allowMerge">
				setCheckbox("allow_merge", '<xsl:value-of select="//attachment/configuration/descriptors/allowMerge"/>');
				</xsl:if>				

				
				<xsl:if test="//attachment/configuration/descriptors/types[type='OEstate']">
				setCheckbox("oestate_bond", '<xsl:value-of select="//attachment/configuration/descriptors/types/configuration/bond"/>');
				setCheckbox("oestate_count", '<xsl:value-of select="//attachment/configuration/descriptors/types/configuration/count"/>');
				</xsl:if>
	
		
				<xsl:if test="//attachment/configuration/descriptors/types[type='SilicosItScaffold']">
				setCheckbox("silicos_scaffolds", '<xsl:value-of select="//attachment/configuration/descriptors/types/configuration/scaffolds"/>');
				setCheckbox("silicos_frameworks", '<xsl:value-of select="//attachment/configuration/descriptors/types/configuration/frameworks"/>');
				setCheckbox("silicos_oprea", '<xsl:value-of select="//attachment/configuration/descriptors/types/configuration/oprea"/>');
				setCheckbox("silicos_schuffenhauer", '<xsl:value-of select="//attachment/configuration/descriptors/types/configuration/schuffenhauer"/>');
				</xsl:if>
				
				<xsl:if test="//attachment/configuration/descriptors/types[type='GSFrag']">
				setCheckbox("gsfrag-normal", '<xsl:value-of select="//attachment/configuration/descriptors/types/configuration/gsfrag"/>');
				setCheckbox("gsfragl", '<xsl:value-of select="//attachment/configuration/descriptors/types/configuration/gsfragl"/>');
				</xsl:if>
				
				setDragonBlocks("#RDKIT-params", "rdkitb", '<xsl:value-of select="//attachment/configuration/descriptors/types[type='RDKIT']/configuration/dragonBlocks"/>', 18);
				setDragonBlocks("#alvaDesc-params", "alva", '<xsl:value-of select="//attachment/configuration/descriptors/types[type='alvaDesc']/configuration/dragonBlocks"/>', 30);
				setDragonBlocks("#PaDEL2", "padel", '<xsl:value-of select="//attachment/configuration/descriptors/types[type='PaDEL2']/configuration/dragonBlocks"/>', 3);
	
				<xsl:if test="//attachment/configuration/descriptors/types[type='ECFP']">
					setValue('ecfp_nbits', '<xsl:value-of select="//attachment/configuration/descriptors/types[type='ECFP']/configuration/MORGAN_NBITS"/>');
					setValue('ecfp_radius', '<xsl:value-of select="//attachment/configuration/descriptors/types[type='ECFP']/configuration/MORGAN_RADIUS"/>');
					setCheckbox('ecfp_fcfp', '<xsl:value-of select="//attachment/configuration/descriptors/types[type='ECFP']/configuration/MORGAN_FCFP"/>');
					setCheckbox('ecfp_counts', '<xsl:value-of select="//attachment/configuration/descriptors/types[type='ECFP']/configuration/MORGAN_COUNTS"/>');
				</xsl:if>
				
				<xsl:if test="//attachment/configuration/descriptors/types[type='RDKIT']">
					setValue('whim_thresh', '<xsl:value-of select="//attachment/configuration/descriptors/types[type='RDKIT']/configuration/WHIM_THRESHOLD"/>');
					setValue('topological_nbits', '<xsl:value-of select="//attachment/configuration/descriptors/types[type='RDKIT']/configuration/TOPOLOGICAL_NBITS"/>');
					setValue('morgan_nbits', '<xsl:value-of select="//attachment/configuration/descriptors/types[type='RDKIT']/configuration/MORGAN_NBITS"/>');
					setValue('morgan_radius', '<xsl:value-of select="//attachment/configuration/descriptors/types[type='RDKIT']/configuration/MORGAN_RADIUS"/>');
					setCheckbox('morgan_fcfp', '<xsl:value-of select="//attachment/configuration/descriptors/types[type='RDKIT']/configuration/MORGAN_FCFP"/>');
					setCheckbox('morgan_counts', '<xsl:value-of select="//attachment/configuration/descriptors/types[type='RDKIT']/configuration/MORGAN_COUNTS"/>');
					setCheckbox('avalon_nbits', '<xsl:value-of select="//attachment/configuration/descriptors/types[type='RDKIT']/configuration/AVALON_NBITS"/>');
					setValue('avalon_counts', '<xsl:value-of select="//attachment/configuration/descriptors/types[type='RDKIT']/configuration/AVALON_COUNTS"/>');
				</xsl:if>
				
				<xsl:if test="//attachment/configuration/descriptors/types[type='Fragmentor']">
					setValue('min-length', '<xsl:value-of select="//attachment/configuration/descriptors/types[type='Fragmentor']/configuration/minFragmentLength"/>');
					setValue('max-length', '<xsl:value-of select="//attachment/configuration/descriptors/types[type='Fragmentor']/configuration/maxFragmentLength"/>');
					setValue('fragment-type', '<xsl:value-of select="//attachment/configuration/descriptors/types[type='Fragmentor']/configuration/fragmentType"/>');
				</xsl:if>
				
				<xsl:if test="//attachment/configuration/descriptors/types[type='SIGMA']">
					setCheckbox("sigmawidth", '<xsl:value-of select="//attachment/configuration/descriptors/types[type='SIGMA']/configuration/sigmawidth"/>');
					setCheckbox("normalise", '<xsl:value-of select="//attachment/configuration/descriptors/types[type='SIGMA']/configuration/normalise"/>');
				</xsl:if>

				<xsl:if test="//attachment/configuration/descriptors/types[type='EPA']">
					setCheckbox("epa3D", '<xsl:value-of select="//attachment/configuration/descriptors/types[type='EPA']/configuration/epa3D"/>');
				</xsl:if>

				<xsl:if test="//attachment/configuration/descriptors/types[type='MAP4']">
					setCheckbox("dimensions", '<xsl:value-of select="//attachment/configuration/descriptors/types[type='MAP4']/configuration/dimensions"/>');
					setCheckbox("radius",     '<xsl:value-of select="//attachment/configuration/descriptors/types[type='MAP4']/configuration/radius"/>');
				</xsl:if>
				
				<xsl:if test="//attachment/configuration/descriptors/types[type='CDK2']">
					setCheckbox("cdk21", '<xsl:value-of select="//attachment/configuration/descriptors/types[type='CDK2']/configuration/descriptorTypes = 'constitutionalDescriptor'"/>');
					setCheckbox("cdk22", '<xsl:value-of select="//attachment/configuration/descriptors/types[type='CDK2']/configuration/descriptorTypes = 'topologicalDescriptor'"/>');
					setCheckbox("cdk23", '<xsl:value-of select="//attachment/configuration/descriptors/types[type='CDK2']/configuration/descriptorTypes = 'geometricalDescriptor'"/>');
					setCheckbox("cdk24", '<xsl:value-of select="//attachment/configuration/descriptors/types[type='CDK2']/configuration/descriptorTypes = 'electronicDescriptor'"/>');
					setCheckbox("cdk25", '<xsl:value-of select="//attachment/configuration/descriptors/types[type='CDK2']/configuration/descriptorTypes = 'hybridDescriptor'"/>');
					setCheckbox("cdk26", '<xsl:value-of select="//attachment/configuration/descriptors/types[type='CDK2']/configuration/descriptorTypes = 'proteinDescriptor'"/>');
					setValue("cdk2Timeout", '<xsl:value-of select="//attachment/configuration/descriptors/types[type='CDK2']/configuration/moleculeTimeout"/>');
				</xsl:if>
							
				<xsl:if test="//attachment/configuration/descriptors/types[type='InductiveDescriptors']">
					setCheckbox("inductiveMElectoneg", '<xsl:value-of select="//attachment/configuration/descriptors/types[type='InductiveDescriptors']/configuration/calcMElectronegativity"/>');
					setCheckbox("inductiveMHardness", '<xsl:value-of select="//attachment/configuration/descriptors/types[type='InductiveDescriptors']/configuration/calcMHardness"/>');
					setCheckbox("inductiveMSoftness", '<xsl:value-of select="//attachment/configuration/descriptors/types[type='InductiveDescriptors']/configuration/calcMSoftness"/>');
					setCheckbox("inductiveMCharge", '<xsl:value-of select="//attachment/configuration/descriptors/types[type='InductiveDescriptors']/configuration/calcMPartialCharges"/>');
					setCheckbox("inductiveMSigma", '<xsl:value-of select="//attachment/configuration/descriptors/types[type='InductiveDescriptors']/configuration/calcMInductiveParam"/>');
					setCheckbox("inductiveMRs", '<xsl:value-of select="//attachment/configuration/descriptors/types[type='InductiveDescriptors']/configuration/calcMStericParam"/>');
				</xsl:if>
				
				<xsl:if test="//attachment/configuration/descriptors/types[type='Spectrophores']">
					setValue('spectrophores-accuracy', '<xsl:value-of select="//attachment/configuration/descriptors/types[type='Spectrophores']/configuration/accuracy"/>');
					setValue('spectrophores-cage', '<xsl:value-of select="//attachment/configuration/descriptors/types[type='Spectrophores']/configuration/cage"/>');
					setValue('spectrophores-resolution', '<xsl:value-of select="//attachment/configuration/descriptors/types[type='Spectrophores']/configuration/resolution"/>');
				</xsl:if>

				<xsl:if test="//attachment/configuration/descriptors/types[type='MORDRED']">
					setValue('mordred', '<xsl:value-of select="//attachment/configuration/descriptors/types[type='MORDRED']/configuration/mordred"/>');
				</xsl:if>
				
				<xsl:if test="//attachment/configuration/descriptors/types[type='SIRMS']">
					setCheckbox('sirms_none', '<xsl:value-of select="//attachment/configuration/descriptors/types[type='SIRMS']/configuration/descriptorTypes = 'none'"/>');
					setCheckbox('sirms_elm', '<xsl:value-of select="//attachment/configuration/descriptors/types[type='SIRMS']/configuration/descriptorTypes = 'elm'"/>');
					setCheckbox('sirms_charge', '<xsl:value-of select="//attachment/configuration/descriptors/types[type='SIRMS']/configuration/descriptorTypes = 'CHARGE'"/>');
					setCheckbox('sirms_logp', '<xsl:value-of select="//attachment/configuration/descriptors/types[type='SIRMS']/configuration/descriptorTypes = 'LOGP'"/>');
					setCheckbox('sirms_hb', '<xsl:value-of select="//attachment/configuration/descriptors/types[type='SIRMS']/configuration/descriptorTypes = 'HB'"/>');
					setCheckbox('sirms_refractivity', '<xsl:value-of select="//attachment/configuration/descriptors/types[type='SIRMS']/configuration/descriptorTypes = 'REFRACTIVITY'"/>');
					setCheckbox('sirms_noh', '<xsl:value-of select="//attachment/configuration/descriptors/types[type='SIRMS']/configuration/noH"/>');
					setCheckbox('sirms_mixtures', '<xsl:value-of select="//attachment/configuration/descriptors/types[type='SIRMS']/configuration/mixtures"/>');
					setValue('sirms-min-frag', '<xsl:value-of select="//attachment/configuration/descriptors/types[type='SIRMS']/configuration/minFragment"/>');
					setValue('sirms-max-frag', '<xsl:value-of select="//attachment/configuration/descriptors/types[type='SIRMS']/configuration/maxFragment"/>');
				</xsl:if>				
				
				setValue('qnpr-min-frag', '<xsl:value-of select="//attachment/configuration/descriptors/types[type='QNPR']/configuration/minFragmentLength"/>');
				setValue('qnpr-max-frag', '<xsl:value-of select="//attachment/configuration/descriptors/types[type='QNPR']/configuration/maxFragmentLength"/>');
			</script>
		</xsl:if>
		
	</xsl:template>
</xsl:stylesheet>
