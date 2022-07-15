<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:template name="structure-optimisation">
		<h1>Select a tool to optimize molecule structures</h1>
		<input type="radio" name="optimisation" value = "0" id="none">
			<xsl:choose>
				<xsl:when test="//param[@key='recommendation'] = 0">
					<xsl:attribute name="checked">checked</xsl:attribute>
				</xsl:when>
				<xsl:otherwise>
					<xsl:attribute name="deprecated">1</xsl:attribute>
				</xsl:otherwise>
			</xsl:choose>
		</input>
		<label for="none"> No optimisation</label><br/>

		<xsl:if test="//method != 'DIMENET'">
			<div id="deprecated" class="warning">
				The descriptors you selected require structure optimisation. Selecting no 3D optimisation  may result into errors.
			</div>
		</xsl:if>

		<xsl:choose>
			<xsl:when test="/model/@inhouse = 'true'">
				<input type="radio" name="optimisation" value = "1" id="corina"/>
				<label for="corina"> Optimise with Corina</label><br/>
		
				<input type="radio" name="optimisation" value = "6" id="balloon">
						<xsl:attribute name="checked">checked</xsl:attribute>
				</input>
				<label for="balloon"> Optimise with BALLOON</label><br/>
			</xsl:when>
			<xsl:otherwise>
				<input type="radio" name="optimisation" value = "1" id="corina">
						<xsl:attribute name="checked">checked</xsl:attribute>
				</input>
				<label for="corina"> Optimise with Corina</label><br/>
		
				<input type="radio" name="optimisation" value = "6" id="balloon"/>
				<label for="balloon"> Optimise with BALLOON</label><br/>
			</xsl:otherwise>
		</xsl:choose>

		<input type="radio" name="optimisation" value = "4" id="babel"/>
		<label for="babel"> Optimise with OpenBabel</label><br/>

		<input type="radio" name="optimisation" value = "5" id="obgen"/>
		<label for="obgen"> Optimise with OBGEN (in OpenBabel)</label><br/>

		<xsl:if test="(//user/ochem-labs = 'false')">
			<input type="radio" name="optimisation" value = "7" id="CCCD"/>
			<label for="cccd">Cambridge Crystallographic Data Centre</label><br/>
		</xsl:if>

		<script language="javascript">
			<xsl:if test="//attachment">
				$("input[type=radio]").removeAttr("checked");
				$("input[value=0]").attr("checked", "checked");
			</xsl:if>

			$(document).ready(function()
			{
				$("[name=optimisation]").change(function(){
					$("#deprecated").setClass("invisible", $("[name=optimisation]").filter(":checked").filter("[deprecated]").length == 0);
				});

				$("[name=optimisation]").change();
			});	
			
		</script>		
	</xsl:template>
	
</xsl:stylesheet>