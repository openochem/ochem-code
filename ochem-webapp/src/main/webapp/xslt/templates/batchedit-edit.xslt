<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />
	
	<xsl:template name="content">
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script language="javascript" src="js/commons/dynamic-select.js"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript" src="js/blocks/batchedit-edit.js"></script>
		<script language="javascript" src="js/blocks/skimming.js"></script>
		<script language="javascript">				
			var images = new Array();
			<xsl:for-each select="list/compressedEP/imageId">
				images.push(<xsl:value-of select="."/>);
			</xsl:for-each>
			var selectedUnitId  = '<xsl:value-of select="list/compressedEP/ep/unit/@id"/>';
			$(document).ready(
				function(){
					var finalString = skim(images);
				}
			);
			
		</script>
		
		<title>Batch record editor</title>
		
		<style type="text/css">
			.outer {width: 100%;}
			.outer TD {}
			.inner TD {padding: 2px 2px 2px 2px;}
			.centered {border: 1px solid #999;}
			TD.inner {border: 1px solid #999; padding: 10px;}
			INPUT {border: 1px solid black; padding: 2px 2px 2px 2px;}
			.narrow {width:100px;}
			.outer IMG {border: 1px solid black;}
			.inner TR TD:first-child {width: 100px;}
			TEXTAREA {width: 500px; height: 50px;}
			#Property-info IMG {border: none !important;}
			#Property-info A {text-decoration: none;}
			.cond-text-value {width: 90%;} 
			
		</style>
		<title>Batch editor</title>
		<h1 class="popup-header">Batch compounds properties editor
		<i>Add and update information about the compounds properties</i></h1>
		
		<table class="outer">
			<tr>
				<td rowspan="5" class="formscope centered">
					<input type="hidden" name="id" value="{list/compressedEP/ep/@id}" send="1"/>
					<input type="hidden" name="key" value="{list/compressedEP/@compressionKey}" send="1"/>
					
					<!-- these three are needed to reconstruct the compression key -->
					<input type="hidden" name="property" send="1" value="{list/compressedEP/ep/property/@id}"/>
					<input type="hidden" name="article" send="1" value="{list/compressedEP/ep/article/@id}"/>
					<input type="hidden" name="unit" send="1" value="{list/compressedEP/ep/unit/@id}"/>
					
					<input type="hidden" name="unitcategory" value="{list/compressedEP/ep/property/unitCategory/@id}"/>
					<input type="hidden" name="unitname" send="1" value="{list/compressedEP/ep/unit/@name}"/>
					
					<input type="hidden" name="newproperty" send="1" value=""/>
					<!-- newunit and newoption are in the down -->
					
					<input type="hidden" name="newarticle" send="1" value=""/>
					
					
					<div id="skim"></div>
				</td>
				<td class="inner formscope">
					<table>
						<tr>
						<td><a title="apply this property to all records">
						<input type="checkbox" name="apply_prop" value="true" send="1" /></a> Property:</td>
						<td width="100">
							<a action="editproperty" name="property-link">
								<xsl:choose>
									<xsl:when test="list/compressedEP/ep/property/@name"><xsl:value-of select="list/compressedEP/ep/property/@name"/></xsl:when>
									<xsl:otherwise>[...]</xsl:otherwise>
								</xsl:choose>
							</a>
						</td>
						
						<td>
							<div id="Quantitive">
								<xsl:if test="list/compressedEP/ep/property/@qualitive='true'">
									<xsl:attribute name="class">invisible</xsl:attribute>
								</xsl:if>
								<select name="predicate" selected="{list/compressedEP/ep/@predicate}" disabled=""></select>
								
								<input type="text" name="value" value="{list/compressedEP/ep/value}" class="w50 right" disabled=""/>
								
<!-- 								<select name="n-newunit" send="1" width="60"></select> -->
								
								<select name="newunit" send="1" width="60">
								<xsl:for-each select="list/compressedEP/ep/property/unitCategory/unit">
									<option value="{@id}">
										<xsl:if test="@id=../../../unit/@id">
											<xsl:attribute name="selected">selected</xsl:attribute>
										</xsl:if>
										<xsl:value-of select="@name"/>
									</option>
								</xsl:for-each>
							</select>								
							</div>
							<div id="Qualitive">
								<xsl:if test="count(list/compressedEP/ep/property/option) = '0'">
									<xsl:attribute name="class">invisible</xsl:attribute>
								</xsl:if>
								<select name="newoption" send="1">
									<xsl:if test="list/compressedEP/ep/property/option/@multi='true'">
										<option value="-1" selected="selected">
											do not change -- (multiple options selected)
										</option>
									</xsl:if>	
									<xsl:for-each select="list/compressedEP/ep/property/option">
										<option value="{@id}">
<!-- 											<xsl:if test="@id=../../option/@id"> -->
<!-- 												<xsl:attribute name="selected">selected</xsl:attribute> -->
<!-- 											</xsl:if> -->
											<xsl:value-of select="@name"/>
										</option>
									</xsl:for-each>
								</select>
							</div>
						</td>
						</tr>
					</table>
				</td>
			</tr>
			<tr>
			<td class="inner formscope">
				<table>
					<tr>
						<td><a title="apply this value to all records">
						<input type="checkbox" name="apply_article" value="true" send="1" /></a> Article:</td>
						<td>
						<a action="editarticle" name="article-link" edit="select-article">
							<xsl:choose>
								<xsl:when test="list/compressedEP/ep/article/@title"><xsl:value-of select="list/compressedEP/ep/article/@title"/></xsl:when>
								<xsl:otherwise>[...]</xsl:otherwise>
							</xsl:choose>
						</a>
					</td></tr>
					<tr>
					<td><a title="apply this value to all records">
					<input type="checkbox" name="apply_page" value="true" send="1" /></a> Page:</td>
					<td><input type="text" name="newpage" value="{list/compressedEP/ep/art-page-num}" class="w50 right" send="1" title="Put a value here to change the page for all records. Leave it blank to keep values untouched." /></td></tr>
					<tr>
					<td><a title="apply this value to all records">
					<input type="checkbox" name="apply_line" value="true" send="1" /></a> Line:</td>
					<td><input type="text" name="newline" value="{list/compressedEP/ep/art-line-num}" class="w50 right" send="1" title="Put a value here to change the line for all records. Leave it blank to keep values untouched." /></td></tr>
					<tr>
					<td><a title="apply this value to all records">
					<input type="checkbox" name="apply_table" value="true" send="1" /></a> Table:</td>
					<td><input type="text" name="newtable" value="{list/compressedEP/ep/art-table-num}" class="w50 right" send="1" title="Put a value here to change the table for all records. Leave it blank to keep values untouched." /></td>
					</tr>
				</table>
			</td>
			</tr>
			<tr>
			<td class="inner formscope">
				<table>
					<tr><td><a title="apply this evidence to all records">
					<input type="checkbox" name="apply_evi" value="true" send="1" /></a> 
					Evidence:</td>
						<td>
						<select name="newevidence" send="1">
							<option value="-1">
<!-- 								<xsl:if test="list/compressedEP/@epEvidence = -1"> -->
									<xsl:attribute name="selected">selected</xsl:attribute>
<!-- 								</xsl:if> -->
								do not change <!-- keep: "keep <xsl:value-of select="list/compressedEP/@epEvidence"/>" -->
							</option>
							<option value="0">
<!-- 								<xsl:if test="list/compressedEP/@epEvidence = 0"> -->
<!-- 									<xsl:attribute name="selected">selected</xsl:attribute> -->
<!-- 								</xsl:if> -->
								No evidence specified
							</option>
							<option value="1">
<!-- 								<xsl:if test="list/compressedEP/@epEvidence = 1"> -->
<!-- 									<xsl:attribute name="selected">selected</xsl:attribute> -->
<!-- 								</xsl:if> -->
								Measured in this article
							</option>
							<option value="2">
<!-- 								<xsl:if test="list/compressedEP/@epEvidence = 2"> -->
<!-- 									<xsl:attribute name="selected">selected</xsl:attribute> -->
<!-- 								</xsl:if> -->
								Measured in this article (to be verified)
							</option>
							<option value="3">
<!-- 								<xsl:if test="list/compressedEP/@epEvidence = 3"> -->
<!-- 									<xsl:attribute name="selected">selected</xsl:attribute> -->
<!-- 								</xsl:if> -->
								Error in the record (see discussion)
							</option>
							<option value="4">
<!-- 								<xsl:if test="list/compressedEP/@epEvidence = 4"> -->
<!-- 									<xsl:attribute name="selected">selected</xsl:attribute> -->
<!-- 								</xsl:if> -->
								Invalid record
							</option>
						</select>
						</td>
					</tr>
				</table>
			</td>
			</tr>
			<tr><td class="inner conditionsscope">
				<!-- CONDITIONS BLOCK -->
				<!--input type="hidden" name="n-article" filter="1" value="{list/compressedEP/ep/article/@id}"/>
				<input type="hidden" name="n-property" value="{list/compressedEP/ep/property/@id}" filter="1"/>
				<input type="hidden" name="n-unit" value="{list/compressedEP/ep/unit/@id}" filter="1"/-->
				<!--input type="hidden" name="id" value="{list/compressedEP/@id}" filter="1"/-->
				<input type="hidden" name="set-id" value="{list/compressedEP/ep/conditions/@id}" filter="1"/>
				<table><tr>
					<td class="formscope" width="80">
						<a title="apply these conditions to all records">
						<input type="checkbox" name="apply_cond" value="true" send="1"/></a> Conditions:</td>
					<td><a action="new">[add]</a></td>
				</tr></table>
				<div id="ConditionsBrowser" class="formscope"></div>
			</td></tr>
			<tr><td class="inner formscope">
					<xsl:choose>
					<xsl:when test="list/compressedEP/ep/@rights=0">					
						<input type="checkbox" name="apply_public" value="true" send="1" /> Check to make these records public (i.e., visible to other users). 
					</xsl:when>
					<xsl:otherwise>						
						<input type="checkbox" name="apply_public" value="true" send="1" /> Some of these records might be hidden. Check to make these records also public. 
					</xsl:otherwise>
					</xsl:choose>
			</td></tr>
		</table>
		
		<div class="formscope popup-footer"><a action="edit">save</a><a href="javascript:window.closeTab();">cancel</a></div>
		
		<div id="waitingDialog"> 
	    	<div class="hd">Please wait</div> 
	    	<div class="bd" style="text-align: center;"> 
	        	Please wait until action is completed.<br/>
	        	It may take a while.<br/>
	        	<img src="img/roller_small.gif"/> 
	    	</div> 
		</div>
		
	</xsl:template>
</xsl:stylesheet>