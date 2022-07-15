<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />
	
	<xsl:template name="content">
	
	<style type="text/css">
		.comp-block-image img
		{
			border: 1px solid black;
			width: 150px;
			height: 150px;
		}
		
		.comp-rblock-image img
		{
			border: 1px solid black;
			width: 100px;
			height: 100px;
		}
	</style>
	
	<h1 align="center"><b>Composer</b></h1>		
	<center>
	<table width="95%" class="comp-block-image">
      <tr>
        <td align="center"><a onclick="MComposer.editMolecule(this,0);" href="#"><img id="baseMol" src="img/indicator.white.gif"/></a></td>
        <td align="center" valign="middle">	<a href="#" class="comp-button-link" onclick="MComposer.update(); return false;">Create molecule &gt;</a></td>
        <td align="center"><a onclick="MComposer.editMolecule(this,0);" href="#"><img id="new-mol" src="img/indicator.white.gif"/></a></td>
      </tr>
      <tr>
        <td align="center">Base Mol</td>
        <td></td>
        <td align="center">Created Mol</td>
      </tr>      
      <tr>
        <td align="center"></td>
        <td></td>
        <td align="center"><a href="#" class="comp-button-link" onclick="MComposer.submit(); return false;">Save Molecule</a></td>
      </tr>
    </table>
	</center>
	<blockquote>
		<a href="#" class="comp-button-link" onclick="MComposer.addRblock(); return false;">Add R group</a>
	</blockquote>	
	<table>
	<tr id="r-groups" valign="top"></tr>
	</table>	
	<script language="javascript" src="js/blocks/composer.js"/>
	</xsl:template>
</xsl:stylesheet>