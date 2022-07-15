<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />
	<xsl:template name="content">
	<style>
	h1 {
		color:#7777FF;
		display:block;
		font-family:Georgia,Helvetica;
		font-size:24px;
	}
	
	.acknowledgements td {
		padding: 10px;
		width: 700px;
	}
	
	body.popup {
		padding: 20px !important;
	}
	
	.acknowledgements TD.icon	{
		border: 0px solid;
		color: black;
		text-align: center;
		vertical-align: middle;
		width: 130px; height: 130px;
		padding: 0px;
	}
	
	</style>
	<h1 id="page-title">Our acknowledgements</h1><br/>
	<table cellspacing="5" class="acknowledgements">
	<tr>
	<td class="icon"><a target="_blank" href="http://cdk.sourceforge.net/"><img src="img/logo-cdk.png"/></a></td>
	<td>The <a target="_blank" href="http://cdk.sourceforge.net/">Chemistry Development Kit</a> (CDK) is a Java library for structural chemo- and bioinformatics. It is now developed by more than 50 developers all over the world and used in more than 10 different academic as well as industrial projects world wide.</td>
	</tr>
	<tr>
	<td class="icon"><a target="_blank" href="http://www.mn-am.com"><img src="img/mn_am.png"/></a></td>
	<td><a target="_blank" href="http://www.mn-am.com">Molecular Networks GmbH and Altamira, LLC</a> provides multifaceted, innovative software to the chemical, biotechnology and pharmaceutical industry. The company's suite of chemoinformatics applications covers many different areas including: handling of chemical information, design of new chemical entities and prediction of physicochemical and biological properties of chemical compounds. </td>
	</tr>
	<tr>
	<td class="icon"><a target="_blank" href="http://chm.kode-solutions.net/"><img src="img/logo-kode.png"/></a></td>
	<td><a target="_blank" href="http://chm.kode-solutions.net/">Kode Chemoinformatics</a> provides services and products in the field of chemoinformatics, with particular expertise on chemical data analysis using chemometric instruments and machine learning techniques (in particular, QSAR/QSPR modelling for drug design and eco-toxicological screening). Kode develop and distribute chemoinformatic software applications, and is the official distributor of Dragon 7.
	</td>
	</tr>
	<tr>
	<td class="icon"><a target="_blank" href="http://www.alvascience.com"><img src="img/alvascience_logo.jpg" width="200" /></a></td><td><a target="_blank" href="http://www.alvascience.com">Alvascience</a> 
	Alvascience provides software products and targeted consulting to reduce the gap between science and IT. Alvascience expertise is in chemoinformatics, QSAR, molecular modelling and data science. 
	</td>
	</tr>
	<tr>
	<td class="icon"><a target="_blank" href="http://peter-ertl.com/jsme/"><img src="img/jsmelogo.png"/></a></td>
	<td><a target="_blank" href="http://peter-ertl.com/jsme">The JSME Molecule Editor</a> is a free, feature rich molecule editor written in JavaScript. JSME supports drawing and editing of molecules and reactions on desktop computers, as well as on handheld devices including iPhone, iPad and Android smartphones and tablets.</td>
	</tr>
	<tr>
	<td class="icon"><a target="_blank" href="http://www.inchi-trust.org/"><img src="img/InChI.png"/></a></td>
	<td>The <a target="_blank" href="http://www.inchi-trust.org/">InChI Trust</a> develops and supports the non-proprietary IUPAC InChI standard and promotes its uses to the scientific community.</td>
	</tr>	
	</table>
	</xsl:template>
	
</xsl:stylesheet>
