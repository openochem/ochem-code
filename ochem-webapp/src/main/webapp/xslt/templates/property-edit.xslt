<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />

	<xsl:template name="content">
		<title>Property editor</title>
		<style type="text/css">
			TD {padding: 4px 4px 4px 4px;}
			TD {border: 1px solid #999; padding: 15px;}
			.compact-item TD {border: none;}
			INPUT, TEXTAREA {border: 1px solid black; padding: 2px 2px 2px 2px; margin: 1px 1px 1px 1px;}
			.right {align:right;}
			.hidden {display: none;}
			.message {font-size: 11px; font-style: italic;}
			A IMG {margin-right: 5px;}
		</style>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script language="javascript" src="js/commons/dynamic-select.js"></script>
		<script language="javascript" src="js/browsers/tags-inline-browser.js"></script>
		<script language="javascript" src="js/blocks/property-edit.js"></script>
		
		<h1 class="popup-header"><doc term="Property+browser#Propertybrowser-Propertyprofile">Property profile</doc>
		<i>View or/and modify an OCHEM property </i></h1>
		
		<table width="780" class="EditTable">
		<tr class="EditForm"><td>
		<input type="hidden" name="id" value="{property/@id}" send="1" filter="1"/>
		<input type="hidden" name="options-count" value="{property/options-count}" send="1" filter="1"/>
		<input type="hidden" name="initial-unit" value="{//property/defaultUnit/@id}"/>
		Name:<input type="text" onkeyup="return checkPropName(this)" name="name" value="{property/@name}" send="1"/>
		<xsl:if test="property/@id">
			<!-- docs_info --> <!-- <a tab="Wiki" href="wikipage/action.do?entities=property&amp;id={property/@id}"><img src="img/icons/wiki.gif"/></a>  -->
		</xsl:if>&#160;&#160;&#160;&#160;
		Type: 
		<select name="type" selected-value="{property/@type}" send="1">
			<xsl:if test="property/@id">
				<xsl:attribute name="disabled">disabled</xsl:attribute>
			</xsl:if>
			<option value="0">Numeric</option>
			<option value="1">Classification</option>
			<xsl:if test="/model/session/user/superuser = 'true'">
				<option value="2">Textual</option>
			</xsl:if>
		</select>&#160;&#160;&#160;&#160;
		<xsl:if test="property/isCondition = 'false'">
			<xsl:if test="property/@rights = 2">
				
				<xsl:choose>
					<xsl:when test="property/approved = 'true'">
						<nobr>
						<img src="img/icons/public.png"/>&#160;Public property, 
						<xsl:choose>
							<xsl:when test="property/moderator">
								moderated by <a href="user/profile.do?login={property/moderator/@login}" tab="User profile"><xsl:value-of select="property/moderator/@login"/></a>
							</xsl:when>
							<xsl:otherwise>
								unmoderated
							</xsl:otherwise>
						</xsl:choose>
						</nobr>
					</xsl:when>
					<xsl:otherwise><img src="img/icons/awaiting-approval.png"/>&#160;Public property
						<xsl:choose>
							<xsl:when test="/model/session/@moderator = 'true'">&#160;&#160;<input type="checkbox" name="approve" send="1"/>&#160;Approve</xsl:when>
							<xsl:otherwise> (awaiting approval)</xsl:otherwise>
						</xsl:choose>
					</xsl:otherwise>
				</xsl:choose>
			</xsl:if>
			
			<xsl:if test="property/@rights = 0">
				<img src="img/icons/ghost.png"/>&#160;Hidden property&#160;&#160;<input type="checkbox" name="public" send="1"/>&#160;Make public&#160;<a help="public-hint" class="question-hint">[?]</a>
				<div class="invisible" id="public-hint">
					Currently, this property is private and visible to you only.<br/>
					If you publish this property, it will become available to all users after approval by OCHEM administrators.<br/><br/>
					<b>In case this property gets approved, you will become a moderator of the data for this property</b>
				</div>
			</xsl:if>
		</xsl:if>
		
		<div id="propName" class="message">(min. 2 characters and max. 40 characters)</div>
		
		<xsl:if test="property/parent">
			<br/>This property is a part of the group <a tab="Group of properties" href="properties/edit.do?id={property/parent/@id}"><xsl:value-of select="property/parent/@name"/></a><br/>
		</xsl:if>
		<xsl:if test="property/isDirectory = 'true' and property/@id">
			<br/><img src="img/icons/folder.gif"/> This is a group that contains <a tab="Children of {property/@name}" href="properties/show.do?parent={property/@id}"><xsl:value-of select="property/@children-count"/> properties</a>
		</xsl:if>
		<xsl:if test="property/@property-record != 0">
			<br/>There are <a href="epbrowser/show.do?property={property/@id}" tab="The records of a property"><xsl:value-of select="property/@property-record"/> records</a> for this property
			<br/><a href="editor/filteredDataReport.do?property={property/@id}" tab="Data report for {property/@name}"><img src="img/icons/analytics-16.png"/>Summary of data records for this property</a>
			<xsl:if test="/model/session/user/ochem-labs = 'true'">
			<br/><a href="matchedpairs/transformations.do?property={property/@id}" tab="MMPs for {property/@name}"><img src="img/icons/mmp-16.png"/>View molecular matched pairs for this property</a>
			</xsl:if>
		</xsl:if>
		<xsl:if test="(property/@property-record = 0) and property/@id">
			<br/>There are no records for this property yet
		</xsl:if>
		</td></tr>
		
		<xsl:if test="/model/session/user/superuser = 'true'">
		<tr class="EditForm" id="Units">
			<td>
				Property weight for bonus points: <input type="text" name="bonusPointsWeight" send="1" value="{property/bonusPointsWeight}" class="right"/>
				<br/>
				<div class="message">As OCHEM editor, you are authorized to specify weight for each property to represent the value of the data.
				Higher weight will mean more bonus points for the users who upload such data.</div>	
			</td>
		</tr>
		</xsl:if>
		
		<tr class="EditForm" id="Units"><td>
			<div>
				System of units
					<select name="category" send="1">
						<xsl:for-each select="others/unitcategory">
							<option value="{@id}">
								<xsl:if test="/model/property/unitCategory/@id = @id">
									<xsl:attribute name="selected">selected</xsl:attribute>
								</xsl:if>
								<xsl:value-of select="@name"/>
							</option>
						</xsl:for-each>
					</select>
				<label>Default unit:</label> 
					<select name="unit" send="1">
						<xsl:for-each select="property/unitCategory/unit">
							<option value="{@id}">
							<xsl:if test="@id=//property/defaultUnit/@id">
								<xsl:attribute name="selected">selected</xsl:attribute>
							</xsl:if>
							<xsl:value-of select="@name"/>
							</option>
						</xsl:for-each>
					</select>
				<xsl:if test="property/used-unit">
					<br/>Following units have been used for this property: 
					<xsl:for-each select="property/used-unit">
						<a href="unit/edit.do?id={@id}" tab="Unit profile">
							<xsl:if test="@name = ''">
								empty unit
							</xsl:if>
							<xsl:value-of select="@name"/></a> (<a href="epbrowser/show.do?property={../@id}&amp;unit={@id}" tab="Records of property+unit"><xsl:value-of select="countInProperty"/></a>), 
					</xsl:for-each>
				</xsl:if>
				
			</div>
		</td></tr>
		<tr><td id="Tags">
				<input type="hidden" name="property-id" value="{property/@id}" filter="1"/>
				<div style="float: left; width: 100px;">Tags <a action="add">[+]</a></div>
					<div id="TagsBrowser" style="float: left;">
					</div>
		</td></tr>
		<tr id="Conditions"><td>
			<xsl:if test="property/isCondition = 'true'">
				<xsl:attribute name="class">invisible</xsl:attribute>
			</xsl:if>
			<div style="float: left; width: 200px;">
				Obligatory Conditions<a action="add">[+]</a></div>
			<div id="ConditionsBrowser" style="float: left;" class="EditForm">
			</div><br/>
		</td></tr>
		<tr id="MiscConditions" class="invisible"><td>
				<div style="float: left;">
					Conditions used:  
				</div>
				<div class="list">
						
				</div>	
		</td></tr>		
		<tr class="EditForm"><td>
			Aliases: <input type="text" style="width: 100%" send="1" name="aliases" value="{property/aliases}"/><br/>
			<div id="propDesc2" class="message">(comma separated list of the other names, i.e. synonyms, for this property)</div><br/>
			Description:<br/>
			<textarea name="description" send="1" onkeyup="return checkPropDesc(this)" style="width: 100%; height: 80px;">
				<xsl:value-of select="property/description" />
			</textarea>
			<div id="propDesc" class="message">(min. 50 characters)</div>
		</td></tr>
		<tr id="Options">
		<td>
			<xsl:choose>
				<xsl:when test="property/options-count &lt; 100">
					<a action="add" class="rounded-button">add new option</a>
				</xsl:when>
				<xsl:otherwise>
					This property has <xsl:value-of select="property/options-count"/> options
				</xsl:otherwise>
			</xsl:choose>
			
			<xsl:if test="property/@id &gt; 0">
				<a href="propertyoptions/show.do?property={property/@id}" tab="Property options" class="rounded-button">Open the options browser</a>
			</xsl:if>
			<br/><br/>
			<div id="Browser">
			</div>
		</td>
		</tr>
		</table>
		<div class="EditForm popup-footer">
			<a action="edit" restrict="form">save</a>
			<a href="javascript:window.closeTab();">cancel</a>
		</div>
		
		
		<div id="confirm-dialog"> 
		    <div class="hd">Confirmation request</div> 
		    <div class="bd">
		    	<div id="confirm-message" style="height: 300px; overflow: auto;">
		    	</div> 
		    </div> 
		</div>
		
		<div id="waiting-dialog"> 
		    <div class="hd">Please wait</div> 
		    <div class="bd" style="text-align: center;"> 
		        Please wait until action is completed.<br/>
		        It may take a while.<br/>
		        <img src="img/roller_small.gif"/> 
		    </div>
		</div>
		
		<script language="javascript">
			$(document).ready(function(){
			<xsl:for-each select="//others/property">
				conditionsBrowser.drawFromJSON({id:"<xsl:value-of select="@id"/>", name:"<xsl:value-of select="@name"/>"});
			</xsl:for-each>
			
			<xsl:if test="property/@id &gt; 0">
				var ajax = new QSPR.Ajax();
				
				// Load the slow condition statistics dynamically for a better user experience
				ajax.send({
					url: "properties/getUsedConditions.do?id=<xsl:value-of select='property/@id'/>",
					success: function(response){
						if (!response.property.conditionsUsed)
							return;
						var conds = array(response.property.conditionsUsed.condition);
						if (conds &amp;&amp; conds.length > 0)
						{
							$("#MiscConditions").removeClass("invisible");
							for (var i = 0; i &lt; conds.length; i++)
								$("#MiscConditions .list").append(getTemplate("condition-template").render(conds[i]));
							
							$(document).trigger("DOM_updated", $("#MiscConditions"));
						}
					},
					after: function(){}
				});
				</xsl:if>
			});
			
		</script>
		
		<script type="text/template" id="condition-template">
			<nobr><a tab="Property profile: [%=name %]" href="properties/edit.do?id=[%= id %]" id="[%= id %]">[%= name %]</a> (<a tab="Filtered records" id="[%= id %]" href="epbrowser/show.do?property={property/@id}&amp;cond-id-1=[%= id %]">[%= timesUsed %]</a>),</nobr>			
		</script>
	</xsl:template>
	
	
</xsl:stylesheet>