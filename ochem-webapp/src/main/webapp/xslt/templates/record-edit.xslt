<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />
	
	<xsl:template name="content">
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script language="javascript" src="js/commons/dynamic-select.js"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript" src="js/blocks/record-edit.js"></script>
		
		<title>Record editor</title>
		
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
		<script language="javascript">
			var dt = ([
				<xsl:for-each select="//moleculename">
					{"id":<xsl:value-of select="@id"/>,"name":&quot;<xsl:value-of select="@name"/>&quot;,"validation":<xsl:value-of select="@validation"/>,"user":&quot;<xsl:value-of select="@user"/>&quot;}
					<xsl:if test="position() != last()">,</xsl:if>
				</xsl:for-each>
				]);
			nameBrowser.items = dt;	
		</script>
		
		<h1 class="popup-header">Compound property editor<i>Add and update information about the compounds properties</i></h1>
		<table class="outer">
			<tr>
				<td rowspan="5" class="centered" valign="top" width="250" style="padding: 3px;">
					<div class="formscope">
						<input type="hidden" name="id" value="{exp-property/@id}" send="1"/>
						<input type="hidden" name="n-molecule" send="1" value="{exp-property/molecule/@id}"/>
						<input type="hidden" name="n-molecule-mapping1" send="1" value="{exp-property/molecule/mapping/mapping1/@id}"/>
						<input type="hidden" name="n-molecule-mapping2" send="1" value="{exp-property/molecule/mapping/@id}"/>
						<input type="hidden" name="n-article" send="1" value="{exp-property/article/@id}"/>
						<input type="hidden" name="n-property" value="{exp-property/property/@id}" send="1"/>
						<input type="hidden" name="options-count" value="{exp-property/property/options-count}" send="1"/>
						<input type="hidden" name="original-option" value="{exp-property/option/@id}"/>
						<input type="hidden" name="original-option-name" value="{exp-property/option/@name}"/>
						<input type="hidden" name="count-similar-ep" value="{exp-property/@count}"/>
					</div>
					<table>
						<tr><td class="formscope" width="250">
							<xsl:choose>					        	 	
				        	 	<xsl:when test="exp-property/molecule">
				        	 		<a action="editmolecule"><img id="depiction" src="depiction.jsp?id={exp-property/molecule/@id}" width="150" height="150"/></a>					        					        	 				
								</xsl:when>
								<xsl:otherwise>
									<a action="editmolecule"><img id="depiction" src="img/click-to-draw.png" alt="" width="150" height="150"/></a>					
								</xsl:otherwise>
							</xsl:choose>
						</td></tr>
						<tr><td id="modify" class="inner formscope invisible highlighted">
							<!-- Modify Similar Records BLOCK -->
							<input type="checkbox" name="modify-similar-ep" send="1" value="true"/>
							You have <xsl:value-of select="exp-property/@count - 1"/> more records in this article with identical molecule structure. 
							Please check to modify other records as well
							<br/><a action="viewrecord">[View Record]</a>
						</td></tr>
						<tr><td class="namescope">
							<!-- NAME BLOCK -->
							<input type="hidden" name="id" value="{exp-property/@id}" filter="1"/>	
							<input type="hidden" name="mol-id" filter="1" value="{exp-property/molecule/@id}"/>					
							<table><tr>
								<td width="60"><b>Names:</b></td>
								<td><a action="new">[add]</a></td>
								<td><a action="checknames" class="search-name">[check]</a></td>	
							</tr></table>
							<div id="NameBrowser"></div>
						</td></tr>						
						<tr><td class="synonymscope">
							<!-- SYNONYM BLOCK -->
							<input type="hidden" name="mol-id" value="{exp-property/molecule/@id}" filter="1"/>
							<input type="hidden" name="id" value="{exp-property/@id}" filter="1"/>
							<input type="hidden" name="initial-category" value="{exp-property/property/unitCategory/@id}"/>
							<input type="hidden" name="initial-unit" value="{exp-property/unit/@id}"/>							
							<div width="80"><b>Synonyms:</b></div>
							<div class="pager-strip">
								<span><b class="showed">none</b> of <b class="total">none</b></span>
								<div id="pager" class="pgr">
								</div>
							</div>
							<div id="SynonymBrowser"></div>
							<div class="pager-strip">
								<span><b class="showed">none</b> of <b class="total">none</b></span>
								<div id="pager" class="pgr">
								</div>
							</div>
						</td></tr>						
					</table>
				</td>
				
				<td class="inner formscope">
					<table>
						<tr>
							<td>Property:</td>
							<td width="140">
								<a action="editproperty" name="property-link" edit="select-property">
									<xsl:choose>
										<xsl:when test="exp-property/property/@name"><xsl:value-of select="exp-property/property/@name"/></xsl:when>
										<xsl:otherwise>[...]</xsl:otherwise>
									</xsl:choose>
								</a>
							</td>	
							<td>
								<div id="Quantitive">
									<xsl:if test="exp-property/property/@qualitive='true'">
										<xsl:attribute name="class">invisible</xsl:attribute>
									</xsl:if>

									<select name="n-predicate" send="1" width="60">
										<xsl:for-each select="others/predicate">
											<option value="{@id}">
												<xsl:if test="@id=../../exp-property/predicate/@id">
													<xsl:attribute name="selected">selected</xsl:attribute>
												</xsl:if>
												<xsl:attribute name="shortname"><xsl:value-of select="@shortName"/></xsl:attribute>
												<xsl:attribute name="name"><xsl:value-of select="@name"/></xsl:attribute>
												<xsl:value-of select="@name"/>
											</option>
										</xsl:for-each>
									</select>

									<input type="text" name="n-value" value="{exp-property/value}" class="w50 right" send="1"/>
									<span id="second-value">
										<span id="predicate">&#x00B1;</span>
										<input type="text" name="n-second-value" value="{exp-property/secondValue}" class="w50 right" send="1"/>
									</span>
									<select name="n-unit" send="1" width="60">
										<xsl:for-each select="exp-property/property/unitCategory/unit">
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
									<xsl:if test="exp-property/property/options-count = '0'">
										<xsl:attribute name="class">invisible</xsl:attribute>
									</xsl:if>
									<select name="n-option" send="1">
										<xsl:for-each select="exp-property/property/option">
											<option value="{@id}">
												<xsl:if test="@id=../../option/@id">
													<xsl:attribute name="selected">selected</xsl:attribute>
												</xsl:if>
												<xsl:value-of select="@name"/>
											</option>
										</xsl:for-each>
									</select>
									<span></span>
								</div>
							</td>
							<td id="Property-info">
								<xsl:if test="not(exp-property/property)">
										<xsl:attribute name="class">invisible</xsl:attribute>
									</xsl:if>
								<a action="describeproperty" title="Property detailed info"><img src="img/icons/edit.gif"/></a>
								<!-- wiki_info --> <!-- <a action="wikiproperty" title="Wiki page for this property"><img src="img/icons/wiki.gif"/></a> -->
							</td>
						</tr>
					</table>
				</td>
			</tr>
			<tr>
				<td class="inner formscope">
					<table>
						<tr><td>Article:</td>
							<td>
							<a action="editarticle" name="article-link" edit="select-article">
								<xsl:choose>
									<xsl:when test="exp-property/article/@title"><xsl:value-of select="exp-property/article/@title"/></xsl:when>
									<xsl:otherwise>[...]</xsl:otherwise>
								</xsl:choose>
							</a>
						</td></tr>
						<tr><td>Page:</td><td><input type="text" name="n-page" value="{exp-property/art-page-num}" class="w50 right" send="1"/></td></tr>
						<tr><td>Line:</td><td><input type="text" name="n-line" value="{exp-property/art-line-num}" class="w50 right" send="1"/></td></tr>
						<tr><td>Table:</td><td><input type="text" name="n-table" value="{exp-property/art-table-num}" class="w50 right" send="1"/></td></tr>
						<tr><td>N (Mol ID):</td><td><input type="text" name="n-art-mol-id" value="{exp-property/art-mol-id}" class="w50 right" send="1"/></td></tr>
					</table>
				</td>
			</tr>
			<tr>
				<td class="inner formscope">
					<table>
						<tr><td>Evidence:</td>
							<td>
							<select name="evidence" send="1">
							<option value="0" style="background-color: #FEE;">No evidence specified</option>
							<option value="1">
								<xsl:if test="(exp-property/@id=exp-property/@connected_id) and not(exp-property/ep_status)">
									<xsl:attribute name="selected">selected</xsl:attribute>
								</xsl:if>
								Measured in this article
							</option>
							<option value="2">
								<xsl:if test="not(exp-property/@id=exp-property/@connected_id) and (exp-property/@connected_id)">
									<xsl:attribute name="selected">selected</xsl:attribute>
								</xsl:if>
								Measured in another article
							</option>
							<option value="3">
								<xsl:if test="exp-property/ep_status = 1">
									<xsl:attribute name="selected">selected</xsl:attribute>
								</xsl:if>
								Measured in this article (to be verified)
							</option>
							<option value="4">
								<xsl:if test="exp-property/ep_status = 0">
									<xsl:attribute name="selected">selected</xsl:attribute>
								</xsl:if>
								Error in the record (see discussion)
							</option>
							<option value="5">
								<xsl:if test="exp-property/ep_status = 3">
									<xsl:attribute name="selected">selected</xsl:attribute>
								</xsl:if>
								Invalid record
							</option>
							</select>
							</td>
						</tr>
						<tr><td>Comment:</td>
							<td>
								<textarea name="n-comment" send="1"><xsl:value-of select="exp-property/other"/></textarea>
							</td>
						</tr>
						<tr name="reference">
							<input type="hidden" name="n-reference" send="1" value="{exp-property/@connected_id}"/>
							<td>Reference:</td>
							<td>
							<a action="editreference" name="reference-link">
								[Record <xsl:value-of select="exp-property/@connected_id"/>]
							</a>
							</td>
						</tr>
					</table>
				</td>
			</tr>
			<tr><td class="inner conditionsscope">
				<!-- CONDITIONS BLOCK -->
				<input type="hidden" name="id" value="{exp-property/@id}" filter="1"/>
				<input type="hidden" name="set-id" value="{exp-property/conditions/@id}" filter="1"/>
				<table><tr>
					<td width="80">Conditions:</td>
					<td><a action="new">[add]</a></td>
				</tr></table>
				<div id="ConditionsBrowser" class="formscope"></div>
			</td></tr>
			<tr><td class="inner formscope">
				<xsl:choose>
				<xsl:when test="exp-property/@id=-1">					
					<input type="checkbox" name="hide" value="true" send="1"/> Check to make the record HIDDEN (i.e., non visible to other users). 
				</xsl:when>
				<xsl:otherwise>
					<xsl:choose>
					<xsl:when test="exp-property/@rights=0">					
						<input type="checkbox" name="hide" value="true" send="1" checked="checked"/> Check to make the record HIDDEN (i.e., non visible to other users). 
					</xsl:when>
					<xsl:otherwise>This record is public</xsl:otherwise>
					</xsl:choose>
				</xsl:otherwise>
				</xsl:choose>
				<br/>
			</td></tr>
		</table>
		<div class="formscope popup-footer"><a action="edit">save</a><a href="javascript:window.closeTab();">cancel</a></div>
		
		<div id="recordDialog"> 
		    <div class="hd">Message</div> 
	   		 <div class="bd"> 
		        <form> 
		        	<input type="hidden" name="rec-id" value=""></input> 
		            <div id="dupRecords">duplicate records</div>
		        </form> 
		    </div> 
		</div>
		
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