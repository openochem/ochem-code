<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />

	<xsl:template name="content">
		<title>Basket editor</title>
		<style type="text/css">
			#form-table TD {padding: 5px 10px 5px 10px;}
			#form-table TD {border: 1px solid #999;}
			.stats TH {font-size: 70%;}
			.stats TD {border: 1px solid #BBB !important;}
			.compact-item TD {border: none;}
			INPUT, TEXTAREA {border: 1px solid black; padding: 2px 2px 2px 2px; margin: 1px 1px 1px 1px;}
			.right {align:right;}
			.hidden {display: none;}
			.message {font-size: 11px; font-style: italic;}
			
			.actions {padding: 10px; background-color: #F9F9F9; float: left; margin-top: 10px; border: 1px solid #EEE; margin-bottom: 10px; width: 95%;}
			.actions B {margin-bottom: 5px; display: block;}
			
			.plaintable TD {padding: 5px 15px !important; border: none !important; vertical-align: top;}
			.plaintable A {margin-bottom: 3px; display: block;}
			.plaintable A IMG {margin-right: 3px; white-space: nowrap;}
			
			.non-editable {border: 0px none;}
		</style>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script language="javascript" src="js/commons/dynamic-select.js"></script>
		<script language="javascript" src="js/browsers/tags-inline-browser.js"></script>
		<script language="javascript" src="js/blocks/basket-edit.js"></script>
		<h1 class="popup-header">Basket editor
		<i>Add new basket or edit exiting basket</i></h1>
		<table width="900" id="form-table">
		<tr class="EditForm"><td>
		<xsl:choose>
			<xsl:when test="//basket/ownerName = /model/session/user/@login or not(//basket)">
				<input type="hidden" name="id" value="{basket/@id}" send="1" filter="1"/>
				Name: <br/><input type="text" onkeyup="return checkBasketName(this)" name="name" value="{basket/@name}" send="1" style="width: 90%;"/>
				<div id="basketName" class="message edit-elem">(min. 2 characters)</div><br/>
				<br/>
				Description (optional):<br/>
				<textarea name="description" send="1" style="width: 90%; height: 50px;"><xsl:value-of select="basket/description"/></textarea>
				<br/>	
				Excluded implicit records (under development):<br/>
				<textarea name="excludedrecords" send="1" style="width: 90%; height: 50px;"><xsl:value-of select="basket/excludedrecords"/></textarea>
				<br/>	
			</xsl:when>
			<xsl:otherwise>
				Basket name: <b><xsl:value-of select="basket/@name"/></b><br/>
				This basket belongs to user: <xsl:value-of select="basket/ownerName"/>
				<xsl:if test="basket/description">
				Description:<br/>
				<xsl:value-of select="basket/description"/>
				</xsl:if>
			</xsl:otherwise>
		</xsl:choose>
		
		
		<xsl:if test="basket/@id">
			<div class="actions">
			<b>Actions</b>
			<table class="plaintable">
				<tr>
					<td><a href="basket/clone.do?basket={basket/@id}" tab="Copy of {basket/@name}"><img src="img/icons/clone.gif"/>Create a copy of this basket</a>
						<a href="basket/findPrimaryRecords.do?basket={basket/@id}" title="Create a basket with primary records of this basket" tab="Primary records basket" class="edit-elem">Create a primary records basket</a>
						<a href="basket/addRecords.do?id={basket/@id}" tab="Add/delete records from a basket" class="edit-elem">Add or delete particular records</a>
						<xsl:if test="basket/propertyUsed/property/@type = 0 and count(basket/propertyUsed/property) = 1">
							<a href="basket/discretize.do?id={basket/@id}" tab="Discretize the basket"  class="edit-elem"><img src="img/icons/discretize.gif"/>Discretize the numerical values</a>
						</xsl:if>
						<xsl:if test="/model/session/user/ochem-labs = 'true'">
							<a href="basket/solventExractSubmit.do?id={basket/@id}" title="Convert basket to use Solvent as Condition" tab="Extract solvents" class="edit-elem">Extract solvents</a>
						</xsl:if>
					</td>
					<td>
						<xsl:if test="basket/@models + basket/@pending-models &gt; 0">
								<a tab="Model summary" href="multiplemodels/show.do?set={basket/@id}" waitMsg="It may take a while to generate a report for a large number of models. Please, wait patiently."><img src="img/icons/table.png"/>&#160;Models summary for <xsl:value-of select="basket/@models + basket/@pending-models"/> models</a> 
						</xsl:if>
					<a href="basket/split.do?id={basket/@id}"  class="edit-elem" tab="Split the basket into two sets" hint="Split the basket into two sets. This is usefull to, for example, automatically create a randomized training and validation sets">Split the basket into two sets</a>
					<a href="baskettransformer/show.do?id={basket/@id}"  class="edit-elem" tab="Transform a basket using OScript">Transform the basket using OScript</a>
					<a href="basket/exportBasket.do?id={basket/@id}" tab="Basket export"><img src="img/icons/xls.gif"/>Export this basket into Excel, CSV or SDF</a>
					<!-- <a href="matchedpairs/transformations.do?basket={basket/@id}" tab="Basket MMPs"><img src="img/icons/mmp-16.png"/>Matched molecular pairs for this basket</a> -->
					</td>
				</tr>
			</table>
			
			
			</div>
			
		</xsl:if>
		
		</td></tr>
		<tr><td id="Properties">
				<input type="hidden" name="id" value="{basket/@id}" filter="1"/>
					<b>Statistics of the basket</b><br/><br/>
					<table class="stats">
						<tr>
						<th>Properties </th>
						<th>Records</th>
						<xsl:if test="basket/excludedCount &gt; 0">
							<th>Thereof excluded</th>
						</xsl:if>
						<th>Unique compounds</th>
						<th></th>
						</tr>
						<xsl:for-each select="basket/propertyUsed/property">
							<tr>
								<td><a href="properties/show.do?id={@id}" tab="property"><xsl:value-of select="@name"/></a></td>
								<td><a href="epbrowser/show.do?basket-select={//basket/@id}&amp;property={@id}" tab="records"><xsl:value-of select="@timesUsed"/> records</a></td>
								<xsl:if test="//basket/excludedCount &gt; 0">
									<td><a href="epbrowser/show.do?basket-select={//basket/@id}&amp;property={@id}&amp;basket-excluded=1" tab="Excluded records"><xsl:value-of select="countExcluded"/> records</a></td>
								</xsl:if>
								<td><xsl:value-of select="countUniqueCompounds"/> compounds</td>
								<td><a tab="MMPs for a basket" href="mmpqsar/basket.do?basket={//basket/@id}&amp;property={@id}"><img src="img/icons/mmp-16.png"/> Show MMPs</a></td>
							</tr>
							<xsl:for-each select="option">
								<xsl:if test="countOfRecords &gt; 0 or countUniqueCompounds &gt; 0">
								<tr>
									<td style="text-align: right;"><xsl:value-of select="@name"/></td>
									<td><a tab="{../@name} - {@name}" href="epbrowser/show.do?basket-select={//basket/@id}&amp;property={../@id}&amp;option={@id}"><xsl:value-of select="countOfRecords"/> records</a></td>
									<xsl:if test="//basket/excludedCount &gt; 0">
										<td><a tab="Excluded {../@name} - {@name}" href="epbrowser/show.do?basket-select={//basket/@id}&amp;property={../@id}&amp;option={@id}&amp;basket-excluded=1"><xsl:value-of select="countExcluded"/> records</a></td>
									</xsl:if>
									<td><xsl:value-of select="countUniqueCompounds"/> compounds</td>
									<td></td>
								</tr>
								</xsl:if>
							</xsl:for-each>
						</xsl:for-each>
						<tr>
							<td>Total Compounds (ignoring stereo-chemistry)</td>
							<td>
								<xsl:attribute name="colspan">
									<xsl:choose>
										<xsl:when test="//basket/excludedCount &gt; 0">4</xsl:when>
										<xsl:otherwise>3</xsl:otherwise>
									</xsl:choose>
								</xsl:attribute>
								<xsl:value-of select="basket/totalUniqueCompounds_SC"/> (<xsl:value-of select="basket/totalUniqueCompounds"/>) compounds
							</td>					
						</tr>
					</table>
					<br/>
					<table class="stats">
						<tr>
						<th>Articles </th>
						<th> Count</th>
						</tr>
						<xsl:for-each select="basket/articleUsed/article">
							<tr>
								<td><a href="article/profile.do?id={@id}" tab="article"><xsl:value-of select="@title"/></a></td>
								<td><a href="epbrowser/list.do?basket-select={//basket/@id}&amp;article={@id}" tab="records"><xsl:value-of select="@basketCount"/></a></td>
							</tr>
						</xsl:for-each>
					</table>
					<br/>
					<xsl:if test="basket/tagUsed/tag">
					<table class="stats">
						<tr>
						<th>Tags </th>
						<th> Count</th>
						</tr>
						<xsl:for-each select="basket/tagUsed/tag">
							<tr>
								<td><a href="tags/show.do?id={@id}&amp;type=molecule" tab="tag"><xsl:value-of select="@name"/></a></td>
								<td><a href="epbrowser/list.do?basket-select={//basket/@id}&amp;tag={@id}" tab="records"><xsl:value-of select="@basketCount"/></a></td>
							</tr>
						</xsl:for-each>
					</table>
					</xsl:if>
		</td></tr>
		</table>
		<div class="EditForm popup-footer edit-elem">
			<a action="edit" restrict="form">save</a>
			<a href="javascript:window.closeTab();">cancel</a>
		</div>
		<script language="javascript">
			$(document).ready(function(){
			 if (currentUser != "<xsl:value-of select="//basket/ownerName"/>" &amp;&amp; "<xsl:value-of select="//basket/@id"/>" &gt; 0)
			{
			 	$(".edit-elem").addClass("invisible");
			 	
				$("input[type=text]").addClass("non-editable");
				$("textarea").addClass("non-editable");
			 }
			<xsl:for-each select="//others/property">
				conditionsBrowser.drawFromJSON({id:"<xsl:value-of select="@id"/>", name:"<xsl:value-of select="@name"/>"});
			
			</xsl:for-each>
			});
		</script>
	</xsl:template>
	
</xsl:stylesheet>