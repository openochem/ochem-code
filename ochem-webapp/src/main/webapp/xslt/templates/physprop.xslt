<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0"
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

	<xsl:include href="../helper.xslt" />

	<xsl:template name="content">
		
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script language="javascript" src="js/commons/mol-editors.js"></script>
		<script language="javascript" src="js/blocks/physprop.js?ver=2.0.8"></script>
		<script language="javascript">document.domain = 'ochem.eu';</script>
		
		<style type="text/css">
			#outer TD {vertical-align: top; padding-right: 10px; }
			.outer-cell {background-color: #FAFAFA; border: 3px solid white; padding: 5px 10px; }
			#outer H1 {font-weight: bold; font-family: Arial; font-size: 14pt; margin-bottom: 15px; padding-left: 3px; color: #111;}
			#results #inner div {padding-left: 3px;}
			#results #inner TD {padding: 5px 15px; background-color: #EEE; border: 3px solid white;}
			TD.value {text-align: right; white-space: nowrap; vertical-align: middle !important;}
			.normalButton, #submitButton {border: 1px solid gray; font-size: 14pt; padding: 5px; background-color: #CCC;}
			#submitButton:hover {background-color: #9E9;}
			#submitButton {font-weight: bold; float: right;background-color: #ACA;}
			.buttons {margin-top: 10px;}
			.foot {font-size: 8pt;}
			.sketcher-frame {
			}
			
			.running {
				font-size: 90%;
				color: #666;
			}
			
			#status {color: #666; margin-top: 5px;}
			
			.model-name INPUT {margin-right: 5px;}
			.model-name {white-space: nowrap;}
			
			#comparison IMG {border-right: 1px solid #DDD; margin-right: 5px;}
			#comparison DIV {float: left; font-size: 90%; border: 1px solid gray; width: 400px; margin-right: 15px; margin-bottom: 10px; position: relative;}
			#comparison TABLE {margin-top: 2px;}
			#comparison TABLE TD {white-space: nowrap;}
			H1 {margin-top: 0px;}
			
			#comparison .delete {position: absolute; top: 0px; right: 3px; font-size: 150%; color: #500; text-decoration: none !important;}
			
		</style>
		
		<div>
		<table id="outer">
			<tr><td class="outer-cell">
					<!-- initialize with phenol data, as shown in the depiction -->
					<!-- phenol -->
					<!-- <input type="hidden" name="moldata" send="1" value="Oc1ccccc1%0A123%0A%20%0A%20%207%20%207%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200999%20V2000%0A%20%20%20%201.2124%20%20%20%204.2000%20%20%20%200.0000%20O%20%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%0A%20%20%20%201.2124%20%20%20%200.0000%20%20%20%200.0000%20C%20%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%0A%20%20%20%200.0000%20%20%20%200.7000%20%20%20%200.0000%20C%20%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%0A%20%20%20%202.4249%20%20%20%200.7000%20%20%20%200.0000%20C%20%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%0A%20%20%20%200.0000%20%20%20%202.1000%20%20%20%200.0000%20C%20%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%0A%20%20%20%202.4249%20%20%20%202.1000%20%20%20%200.0000%20C%20%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%0A%20%20%20%201.2124%20%20%20%202.8000%20%20%20%200.0000%20C%20%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%0A%20%201%20%207%20%201%20%200%20%200%20%200%20%200%0A%20%202%20%203%20%202%20%200%20%200%20%200%20%200%0A%20%202%20%204%20%201%20%200%20%200%20%200%20%200%0A%20%203%20%205%20%201%20%200%20%200%20%200%20%200%0A%20%204%20%206%20%202%20%200%20%200%20%200%20%200%0A%20%205%20%207%20%202%20%200%20%200%20%200%20%200%0A%20%206%20%207%20%201%20%200%20%200%20%200%20%200%0AM%20%20END%0A"/> -->
					<!-- caffein -->
					<input type="hidden" name="moldata" send="1" value="Cn1cnc2c1c(%3DO)n(C)c(%3DO)n2C%0A123%0A%20%0A%2014%2015%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200999%20V2000%0A%20%20%20%200.3903%20%20%20%205.2641%20%20%20%200.0000%20C%20%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%0A%20%20%20%205.7917%20%20%20%204.2000%20%20%20%200.0000%20C%20%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%0A%20%20%20%203.3668%20%20%20%200.0000%20%20%20%200.0000%20C%20%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%0A%20%20%20%203.3668%20%20%20%205.6000%20%20%20%200.0000%20O%20%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%0A%20%20%20%205.7917%20%20%20%201.4000%20%20%20%200.0000%20O%20%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%0A%20%20%20%200.0000%20%20%20%202.8000%20%20%20%200.0000%20C%20%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%0A%20%20%20%200.8229%20%20%20%201.6674%20%20%20%200.0000%20N%20%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%0A%20%20%20%203.3668%20%20%20%204.2000%20%20%20%200.0000%20C%20%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%0A%20%20%20%204.5792%20%20%20%202.1000%20%20%20%200.0000%20C%20%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%0A%20%20%20%202.1544%20%20%20%202.1000%20%20%20%200.0000%20C%20%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%0A%20%20%20%202.1544%20%20%20%203.5000%20%20%20%200.0000%20C%20%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%0A%20%20%20%200.8229%20%20%20%203.9326%20%20%20%200.0000%20N%20%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%0A%20%20%20%204.5792%20%20%20%203.5000%20%20%20%200.0000%20N%20%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%0A%20%20%20%203.3668%20%20%20%201.4000%20%20%20%200.0000%20N%20%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%0A%20%201%2012%20%201%20%200%20%200%20%200%20%200%0A%20%202%2013%20%201%20%200%20%200%20%200%20%200%0A%20%203%2014%20%201%20%200%20%200%20%200%20%200%0A%20%204%20%208%20%202%20%200%20%200%20%200%20%200%0A%20%205%20%209%20%202%20%200%20%200%20%200%20%200%0A%20%206%20%207%20%202%20%200%20%200%20%200%20%200%0A%20%206%2012%20%201%20%200%20%200%20%200%20%200%0A%20%207%2010%20%201%20%200%20%200%20%200%20%200%0A%20%208%2011%20%201%20%200%20%200%20%200%20%200%0A%20%208%2013%20%201%20%200%20%200%20%200%20%200%0A%20%209%2013%20%201%20%200%20%200%20%200%20%200%0A%20%209%2014%20%201%20%200%20%200%20%200%20%200%0A%2010%2011%20%202%20%200%20%200%20%200%20%200%0A%2010%2014%20%201%20%200%20%200%20%200%20%200%0A%2011%2012%20%201%20%200%20%200%20%200%20%200%0AM%20%20END%0A"/>
						<h1>Select a molecule for prediction</h1>
				
						<div id="dep">
							<div><a action="editmolecule"><img id="depiction" src="" style="border-width:4px;border-style:groove;color:#b0c4de;"/></a></div>
						</div>
						<div id="mol-editor" valign="top">
							<iframe src="../editor.html" id="sketch" class="sketcher-frame"></iframe>
							<div class="buttons">
								<input class="normalButton" id="selectButton" type="button" action="select" value="Select" />
								<input class="normalButton" id="cancelButton" type="button" action="cancel" value="Cancel" />
							</div>
						</div>
						<div class="buttons">
							<input class="normalButton" id="drawButton" type="button" action="editmolecule" title="Draw a molecule in a structure editor" value="Draw a molecule" />
							<input id="submitButton" type="button" action="submit" value="Predict!" />
						</div>
			</td>		
			<td class="outer-cell">		
			<h1>Prediction results</h1>
				<div id="results">
					<table>
						<tr>
							<td id="inner">
								<table>
									<xsl:for-each select="//others/model">
											<xsl:for-each select="modelMappings">
												<tr>
													<td class="model-name">
														
														<input type="checkbox" model-id="{../@id}"/>
														<span>
															<xsl:attribute name="title">Used model: <xsl:value-of select="../@name"/></xsl:attribute>
															<xsl:choose>
															<xsl:when test="count(../modelMappings) > 1">
																<xsl:value-of select="property/@name"/>
															</xsl:when>
															<xsl:otherwise>
																<xsl:value-of select="../@featuredName"/>
															</xsl:otherwise>
														</xsl:choose>
														</span>
														
														
													</td>
													<td class="value">-</td>
												</tr>
											</xsl:for-each>
									</xsl:for-each>
								</table>
								<div class="invisible foot">&#42; within a confidence interval of 66%</div>
							</td>
						</tr>
					</table>
					<div class="running invisible">
					<br/>We are currently calculating your predictions. It will take a few more moments.
					<br/> Thank you for your patience!
					</div>
					<div id="status"></div>
					
					
				</div>
				<div>
					<input id="stopButton" type="button" action="stop" value="Stop" class="invisible"/>
				</div>
			</td>
			<td id="comparison" class="invisible outer-cell">
				<h1>Predictions comparison</h1>
			</td>
			</tr>
		</table>		
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