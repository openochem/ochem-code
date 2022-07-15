<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />

	<xsl:template name="content">
		<title>Molecule sets</title>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript" src="js/browsers/basket-browser.js"></script>
		<table height="100%" width="100%">
			<tr><td class="itunes-up" colspan="2">
				<h1><doc term="Basket+browser">Basket browser</doc></h1>
				Browse, Compare or Join molecule sets
			</td></tr>
			<tr>
				<td class="itunes-right">
					<div class="upper-command-panel">
						Filter by name: <input filter="1" type="text" name="name"/>
						<a action="edit" title="Create new molecule set">
							[Create new <img src="img/icons/new.gif"/>]
						</a>
					<a class="invisible" action="compare" title="Compare selected molecule sets to find duplicates">
						<img src="img/icons/compare.gif"/>
					</a>
					<a class="invisible" action="combine" title="Join selected molecules set ">
						<img src="img/icons/union.gif"/>
					</a>
					<input type="checkbox" name="public" value="true" filter="1"/>Show public sets
					<xsl:if test="/model/session/user/group">
						<input type="checkbox" name="group" value="true" filter="1"/>
						Sets of group <xsl:value-of select="/model/session/user/group/name"/>
					</xsl:if> 
					</div>
					<input type="hidden" name="name" value="" filter="1"/>
					<input type="hidden" name="showsystem" value="1" filter="1"/>
					<div class="pager-strip">
						<span><b class="showed">none</b> of <b class="total">none</b></span>
						<div id="pager" class="pgr">
						</div>
					</div>
					<div id="aaa">&#160;</div>
					<div id="Browser">
					</div>
					<div class="pager-strip">
						<b class="showed">none</b> of <b class="total">none</b>
					</div>
				</td>
			</tr>
		</table>
		<div id="basicxlsmenu" class="yuimenu browserScope">
    	<div class="bd">
        <ul class="first-of-type">
        	<li class="yuimenuitem">
                <a action="basket" class="yuimenuitemlabel">
                    Export this basket to a file (Excel, CSV or SDF)
                </a>
            </li>
            <li class="yuimenuitem">
                <a action="indices" class="yuimenuitemlabel">
                   Calculate descriptors for this basket
                </a>
            </li>
          </ul>
          <ul>
            <li class="yuimenuitem">
                <a action="upload" class="yuimenuitemlabel">
                    Upload records into this basket
                </a>
            </li>
            </ul>
    	</div>
		</div>
		<div id="basicbasketmenu" class="yuimenu browserScope">
    	<div class="bd">
        <ul class="first-of-type">
        	<li class="yuimenuitem">
                <a action="nameclick" class="yuimenuitemlabel">
                    Rename the basket
                </a>
            </li>
            <li class="yuimenuitem">
                <a action="statistics" class="yuimenuitemlabel">
                   Show statistics of basket
                </a>
            </li>
          </ul>
    	</div>
		</div>
	</xsl:template>
</xsl:stylesheet>