<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"  version="2.0">

	<xsl:output method="html" 
		omit-xml-declaration="yes" 
		indent="yes"
		 doctype-public="-//W3C//DTD XHTML 1.0 Strict//EN"
        doctype-system="http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd"
		/>

	<xsl:template match="main-menu">
	<div id="yui-main">
		<div class="yui-b">
			<div id="{@id}" class="yuimenubar yuimenubarnav">
				<div class="bd">
					<div class="font-bar">
						<a href="javascript:changeFontSize(5);" title="Increase font size">A+</a>
						<a href="javascript:changeFontSize(-5);" title="Decrease font size">a-</a>
						<a href="Privacy_Policy.htm" title="Read our privacy statement">Privacy statement</a>
					</div>
					<ul class="first-of-type">
						<xsl:for-each select="item">
							<li class="yuimenubaritem first-of-type">
								<a class="yuimenubaritemlabel" href="{@href}" target="{@target}">
									<xsl:value-of select="@title" />
								</a>
								<xsl:apply-templates />
							</li>
						</xsl:for-each>
					</ul>
				</div>
			</div>
		</div>
	</div>
	<script language="javascript">
			YAHOO.util.Event.onDOMReady(
			function()
			{
		    	var oMenuBar = new YAHOO.widget.MenuBar("<xsl:value-of select='@id'/>", { autosubmenudisplay: true, hidedelay: 750, lazyload: true });
        		oMenuBar.render();
        	}
        	);
		</script>
	</xsl:template>


	<xsl:template match="wiki">
		<a tab="Wiki" href="wikipage/action.do?name={@term}" title="Wikipedia article"><img src="img/icons/wiki2.gif"/></a>
	</xsl:template>
	
	<xsl:template match="doc">
			<xsl:variable name="finaltitle">
				<xsl:choose>
					<xsl:when test="@title"><xsl:value-of select="@title"/> (Click to read more)</xsl:when>
					<xsl:otherwise>Click to read</xsl:otherwise>
				</xsl:choose>
			</xsl:variable>
			
			<xsl:choose>
				<xsl:when test="@hide='true'">
					<div class="doc"><xsl:apply-templates/><a class="infolink blank" href="https://docs.ochem.eu/display/MAN/{@term}" title="{$finaltitle}" target="_blank"></a></div>
				</xsl:when>
				<xsl:otherwise>
					<div class="doc"><xsl:apply-templates/><a class="infolink" href="https://docs.ochem.eu/display/MAN/{@term}" title="{$finaltitle}" target="_blank"></a></div>
				</xsl:otherwise>
			</xsl:choose>
	</xsl:template>

	<xsl:template match="sub-menu">
		<div id="{@id}" class="yuimenu">
			<div class="bd">
				<xsl:if test="item">
					<ul>
						<xsl:apply-templates select="item"/>	
					</ul>
				</xsl:if>
				<xsl:apply-templates select="items"/>
			</div>
		</div>
	</xsl:template>
	
	<xsl:template match="items">
		<ul>
			<xsl:apply-templates select="item"/>	
		</ul>
	</xsl:template>
	
	<xsl:template match="item">
		<li class="yuimenuitem">
			<a class="yuimenuitemlabel" href="{@href}" target="{@target}"><xsl:if test="@img"><img src="{@img}"/></xsl:if><xsl:value-of select="@title"/></a>
			<xsl:apply-templates />
		</li>
	</xsl:template>
	
	<xsl:template match="item[@type='bold']">
		<li class="yuimenuitem">
			<a class="yuimenuitemlabel" href="{@href}" target="{@target}"><b><xsl:if test="@img"><img src="{@img}"/></xsl:if><xsl:value-of select="@title"/></b></a>
			<xsl:apply-templates />
		</li>
	</xsl:template>

	<xsl:template match="popup-menu">
		<script type="text/javascript">
			YAHOO.util.Event.onDOMReady(function () 
			{
				var <xsl:value-of select="@id"/>ItemData = 
				[
                	<xsl:apply-templates mode="popup"/>
                ];
                
                var <xsl:value-of select="@id"/> = 
                new YAHOO.widget.ContextMenu
                (
                      "<xsl:value-of select="@id"/>",
                	{
              	   		trigger: "clones",
                	 itemdata: oFieldContextMenuItemData,
                	 lazyload: true,
                 	effect: 
                 	{ 
                 		effect: YAHOO.widget.ContainerEffect.FADE,
                 		duration: 0.1
                	}                                                     
               }
               );

			});

		</script>
	</xsl:template>

	<xsl:template match="item" mode="popup">
		{text: '<xsl:value-of select="." />'}
		<xsl:if test="position() != last()">,</xsl:if>
	</xsl:template>
	

	<xsl:template match="@*|node()" priority="-2"><xsl:copy><xsl:apply-templates select="@*|node()" /></xsl:copy></xsl:template>
	
	<xsl:template match="script|iframe" priority="-1">
		<xsl:copy>	
			<xsl:copy-of select="@*"/>
			<xsl:comment/>
			<xsl:apply-templates/>
		</xsl:copy>
	</xsl:template>


</xsl:stylesheet>