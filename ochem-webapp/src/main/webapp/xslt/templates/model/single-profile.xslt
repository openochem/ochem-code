<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE xsl:stylesheet [
  <!ENTITY return "&#13;">
]>

<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt"/>
	<xsl:include href="profile-common.xslt"/>
	<xsl:template name="content">
		
		
		<style type="text/css">
	  		table.tiny {margin: 10px 10px 10px 10px;}
	  		table.tiny TD, table.tiny TH {padding: 5px; border: 1px solid black;}
	  		table.tiny TD {text-align: center;}
	  		.warning {color:#F00; font-size: 90%;}
	  		.right-note {float: right; text-align: right; width: 300px; font-size: 90%;}
	  		.ad B {display: block; margin-top: 10px;}
	  		.classificationSummary TD {padding: 15px; text-align: center;}
	  		.confusion-matrix{background-color: #C1E1FB}
			.real{background-color: #BBEEBB}
			.predicted{background-color: #EEAAAA}
			.stath{background-color: #EEEEAA}
			.statv{background-color: #EEEEEE}
			.dot {margin-right: 3px;}
			TABLE.scatter-plots TD {padding-right: 20px;}
			.actions {margin: 10px;}
			.actions SPAN {margin-right: 20px;}
			
			#size-summary {padding: 5px; background-color: #dbe3ec; border: 1px solid gray; margin-top: 10px;}
	  	</style>
	  	<script language="javascript" type="text/javascript" src="js/lib/excanvas.min.js"></script>
	  	<script language="javascript" type="text/javascript" src="js/lib/jquery.flot-0.6.min.js"></script>
	  	<script language="javascript" src="js/commons/plotting.js"></script>
	  	<script language="javascript" src="js/blocks/model-profile.js"></script>
	  	<title>Model profile</title>
	  	<table width="100%">
	  	<xsl:if test="not(model/taskId)">
					<tr><td class="itunes-up">
						<h1><doc term="Statistical+parameters">Model profile</doc></h1>
						Statistical parameters, tables, charts - all the information related to the model.
					</td></tr>
		</xsl:if>
		<tr >
			<td class="itunes-right">
			
			
			<div id="demo" class="yui-navset">
				    <ul class="yui-nav"> 
				    	<li class="selected">
				    		<a href="#tab1"><em>Overview</em></a>
				    	</li>
				        <xsl:if test="others/ad">
					        <li>
					        	<a href="#tab33"><em>Applicability domain</em></a>
					        </li>
				        </xsl:if>
				        <xsl:if test="model/qmrfUrl">
				        	<li>
					        	<a href="#tab34"><em>QMRF data</em></a>
					        </li>	
				        </xsl:if>
				    </ul> 
					<div class="yui-content">
						<div>
							<div class="right-note">
								<xsl:call-template name="lf2br">
       								<xsl:with-param name="StringToTransform" select="model/description"/>
   								</xsl:call-template>
   								<br/><br/>
   								<xsl:if test="model/timeToComplete">
   									<i>Calculated in <xsl:value-of select="model/timeToComplete"/> seconds</i><br/>
   								</xsl:if>
   								<i><a href="javascript:showSizes()" title="Click to see the detailed size summary">Size: <xsl:value-of select="round(model/size div 1024)"/> Kb</a></i>
   								<div id="size-summary" class="invisible">
   								</div>
							</div>
							Model name: 
							<span id="model-name" publicid="{model/publicId}"><xsl:value-of select="model/@name"/></span>&#160; 
							<xsl:if test="model/owner = 'true'">
								<a href="javascript:rename()">[rename]</a>
							</xsl:if>
								
							<xsl:if test="model/published='true'">
								, published in <a tab="Model's article" href="article/profile.do?id={model/article/@id}"><xsl:value-of select="model/article/@title"/></a>
								<br/>Public ID is <a href="model/{model/publicId}" target="_blank" name="link_public_id"><xsl:value-of select="model/publicId"/></a>
							</xsl:if>
							<xsl:if test="model/publicId and model/published != 'true'">
								<br/>Temporal Public ID: <a href="model/{model/publicId}"><xsl:value-of select="model/publicId"/></a> - use this link to share this temporal model (will be deleted after 3 months unless published) 
							</xsl:if>
							<br/><br/>
							<span id="property-name">Predicted property: <b><xsl:value-of select="model/modelMappings/property/@name"/></b></span> <!-- docs_info --> <!--  <a tab="Property description" href="wikipage/action.do?entities=property&amp;id={model/modelMappings/property/@id}"><xsl:value-of select="model/modelMappings/property/@name"/></a><br/> -->
							<xsl:if test="model/modelMappings/unit/@id != '9'">
								modeled in <xsl:value-of select="model/modelMappings/unit/@name"/>
							</xsl:if>
							<br/>
							<span id="method-name">Training method: <xsl:value-of select="model/template/@name"/></span><br/>
							<table id="rmse" class="tiny">
								<tr>
									<th rowspan="2">Data Set</th>
									<xsl:if test="//others/statistics[2]">
										<th colspan="5" style="text-align: center;">Original</th>
										<th colspan="5" style="text-align: center;">Recalculated</th>
									</xsl:if>
								</tr>
								<tr>
									<td>#</td>
									<xsl:choose>
										<xsl:when test="model/modelMappings/property/@qualitive = 'true'">
											<td>Accuracy</td>
											<td>Balanced Accuracy</td>
											<td>MCC</td>
											<td>AUC</td>
											<xsl:if test="//others/statistics[2]">
												<td>#</td>
												<td>Accuracy</td>
												<td>Balanced Accuracy</td>
												<td>MCC</td>
												<td>AUC</td>
											</xsl:if>
										</xsl:when>
										<xsl:otherwise>
											<td>R<sup>2</sup></td>
											<td title="Coefficient of determination">q<sup>2</sup></td>
											<td>RMSE</td>
											<td>MAE</td>
											<xsl:if test="//others/statistics[2]">
												<td>#</td>
												<td>R<sup>2</sup></td>
												<td title="Coefficient of determination">q<sup>2</sup></td>
												<td>RMSE</td>
												<td>MAE</td>
											</xsl:if>
										</xsl:otherwise>
									</xsl:choose>
								</tr>
								<tr>
									<th><img src="img/icons/red-dot.gif" class="dot"/>Training set: 
										<a tab="Training set" href="basket/edit.do?id={model/training-set/@id}">
											<xsl:value-of select="model/training-set/@name"/>
										</a> 
									</th>
									<td><a tab="Training set" href="model/browseRecords.do?mm_id={model/modelMappings/@id}&amp;set=training"><xsl:value-of select="//others/statistics[1]/set[1]/@size"/> records</a>
										<xsl:call-template name="numRecords"><xsl:with-param name="statisticsSet" select="//others/statistics[1]/set[1]"/></xsl:call-template>
									</td>
									<xsl:choose>
										<xsl:when test="model/modelMappings/property/@qualitive = 'true'">
											<xsl:call-template name="classificationColumns"><xsl:with-param name="statisticsSet" select="//others/statistics[1]/set[1]"/></xsl:call-template>
											<xsl:if test="//others/statistics[2]">
												<td>
													<a tab="Training set" href="epbrowser/show.do?basket-select={model/training-set/@id}&amp;property={model/modelMappings/property/@id}&amp;excluded=-{model/@id}"><xsl:value-of select="//others/statistics[2]/set[1]/@size"/> records</a>
													<xsl:call-template name="numRecords"><xsl:with-param name="statisticsSet" select="//others/statistics[2]/set[1]"/></xsl:call-template>
												</td>
												<xsl:call-template name="classificationColumns"><xsl:with-param name="statisticsSet" select="//others/statistics[2]/set[1]"/></xsl:call-template>
											</xsl:if>
										</xsl:when>
										<xsl:otherwise>
											<xsl:call-template name="regressionColumns"><xsl:with-param name="statisticsSet" select="//others/statistics[1]/set[1]"/></xsl:call-template>
											<xsl:if test="//others/statistics[2]">
												<td><a tab="Training set" href="model/browseRecords.do?mm_id={model/modelMappings/@id}&amp;set=training&amp;recalculated=1"><xsl:value-of select="//others/statistics[2]/set[1]/@size"/> records</a>
												<xsl:call-template name="numRecords"><xsl:with-param name="statisticsSet" select="//others/statistics[2]/set[1]"/></xsl:call-template>
												</td>
												<xsl:call-template name="regressionColumns"><xsl:with-param name="statisticsSet" select="//others/statistics[2]/set[1]"/></xsl:call-template>
											</xsl:if>
										</xsl:otherwise>
									</xsl:choose>
								</tr>
								<xsl:if test="//others/statistics/set[@setId='validation']">
									<tr>
										<th>
											<img src="img/icons/green-dot.gif" class="dot"/>Test set: 
											<xsl:choose>
												<xsl:when test="count(model/validation-sets/validation-set) &gt; 1">
													<select name="validation" filter="1">
														<xsl:for-each select="model/validation-sets/validation-set">
															<option value="{position()}"><xsl:value-of select="@name"/></option>
														</xsl:for-each>
													</select>
												</xsl:when>
												<xsl:otherwise>
													<a tab="Test set" href="basket/edit.do?id={model/selectedValidationSet/@id}"><xsl:value-of select="model/selectedValidationSet/@name"/></a>
												</xsl:otherwise>
											</xsl:choose> 
											<a href="javascript:deleteValidationSet();" class="delete-link" title="Delete this validation set from the model">[x]</a>
										</th>
										<td>
											<a tab="Test set" href="model/browseRecords.do?mm_id={model/modelMappings/@id}&amp;set=validation&amp;validationSetNum={//others/statistics[1]/validationSetId}"><xsl:value-of select="//others/statistics[1]/set[@setId='validation']/@size"/> records</a>
											<xsl:call-template name="numRecords"><xsl:with-param name="statisticsSet" select="//others/statistics[1]/set[@setId='validation']"/></xsl:call-template>
										</td>
										<xsl:choose>
										<xsl:when test="model/modelMappings/property/@qualitive = 'true'">
											<xsl:call-template name="classificationColumns"><xsl:with-param name="statisticsSet" select="//others/statistics[1]/set[@setId='validation']"/></xsl:call-template>
											<xsl:if test="//others/statistics[2]">
												<td><a tab="Test set" href="epbrowser/show.do?basket-select={model/selectedValidationSet/@id}&amp;property={model/modelMappings/property/@id}&amp;excluded=-{model/@id}"><xsl:value-of select="//others/statistics[2]/set[@setId='validation']/@size"/> records</a>
												<xsl:call-template name="numRecords"><xsl:with-param name="statisticsSet" select="//others/statistics[2]/set[@setId='validation']"/></xsl:call-template>
												</td>
												<xsl:call-template name="classificationColumns"><xsl:with-param name="statisticsSet" select="//others/statistics[2]/set[@setId='validation']"/></xsl:call-template>
											</xsl:if>
										</xsl:when>
										<xsl:otherwise>
											<xsl:call-template name="regressionColumns"><xsl:with-param name="statisticsSet" select="//others/statistics[1]/set[@setId='validation']"/></xsl:call-template>
											<xsl:if test="//others/statistics[2]">
												<td><a tab="Test set" href="model/browseRecords.do?mm_id={model/modelMappings/@id}&amp;set=validation&amp;validationSetNum={//others/statistics[2]/validationSetId}&amp;recalculated=1"><xsl:value-of select="//others/statistics[2]/set[@setId='validation']/@size"/> records</a>
												<xsl:call-template name="numRecords"><xsl:with-param name="statisticsSet" select="//others/statistics[2]/set[@setId='validation']"/></xsl:call-template>
												</td>
												<xsl:call-template name="regressionColumns"><xsl:with-param name="statisticsSet" select="//others/statistics[2]/set[@setId='validation']"/></xsl:call-template>
											</xsl:if>
										</xsl:otherwise>
									</xsl:choose>
									</tr>
								</xsl:if>
								<xsl:if test="//others/statistics/set[@setId='excluded'] and //others/statistics/set[@setId='excluded']/@size &gt; 0">
									<tr>
										<th><img src="img/icons/blue-dot.gif" class="dot"/>Excluded from training set 
										</th>
										<td><a tab="Excluded from training set" href="model/browseRecords.do?mm_id={model/modelMappings/@id}&amp;set=excluded"><xsl:value-of select="//others/statistics[1]/set[@setId='excluded']/@size"/> records</a>
										<xsl:call-template name="numRecords"><xsl:with-param name="statisticsSet" select="//others/statistics[1]/set[@setId='excluded']"/></xsl:call-template>
										</td>
										<xsl:choose>
										<xsl:when test="model/modelMappings/property/@qualitive = 'true'">
											<xsl:call-template name="classificationColumns"><xsl:with-param name="statisticsSet" select="//others/statistics[1]/set[@setId='excluded']"/></xsl:call-template>
											<xsl:if test="//others/statistics[2]">
												<td><a tab="Excluded from training set" href="epbrowser/show.do?basket-select={model/training-set/@id}&amp;property={model/modelMappings/property/@id}&amp;excluded={model/@id}"><xsl:value-of select="model/modelMappings/@excludedSize"/><xsl:value-of select="//others/statistics[2]/set[@setId='excluded']/@size"/> records</a>
												<xsl:call-template name="numRecords"><xsl:with-param name="statisticsSet" select="//others/statistics[2]/set[@setId='excluded']"/></xsl:call-template>
												</td>
												<xsl:call-template name="classificationColumns"><xsl:with-param name="statisticsSet" select="//others/statistics[2]/set[@setId='excluded']"/></xsl:call-template>
											</xsl:if>
										</xsl:when>
										<xsl:otherwise>
											<xsl:call-template name="regressionColumns"><xsl:with-param name="statisticsSet" select="//others/statistics[1]/set[@setId='excluded']"/></xsl:call-template>
											<xsl:if test="//others/statistics[2]">
												<td><a tab="Excluded from training set" href="model/browseRecords.do?mm_id={model/modelMappings/@id}&amp;set=excluded&amp;recalculated=1"><xsl:value-of select="//others/statistics[2]/set[@setId='excluded']/@size"/> records</a>
												<xsl:call-template name="numRecords"><xsl:with-param name="statisticsSet" select="//others/statistics[2]/set[@setId='excluded']"/></xsl:call-template>
												</td>
												<xsl:call-template name="regressionColumns"><xsl:with-param name="statisticsSet" select="//others/statistics[2]/set[@setId='excluded']"/></xsl:call-template>
											</xsl:if>
										</xsl:otherwise>
									</xsl:choose>
									</tr>
								</xsl:if>
								<xsl:if test="//others/statistics/set[@setId='modified-on-the-fly'] and //others/statistics/set[@setId='modified-on-the-fly']/@size &gt; 0">
									<tr>
										<th>
											<img src="img/icons/yellow-dot.gif" class="dot"/>Modified (need recalculation) <a help="modified-points-help">[?]</a><br/>
											<div id="modified-points-help" class="invisible">
												You can dynamically exclude the points from the training set.<br/>
												Until the model is recalculated, such modified points will be marked yellow.<br/><br/>
												Once the model is recalculated, the points will be colored as blue (excluded points) or as red (training set points)
											</div>
											<small><i>These compounds were excluded or included after the model was calculated.<br/>The model needs recalculation to update these points</i></small>
										</th>
										<td>
											<xsl:if test="//others/statistics[1]/set[@setId='modified-on-the-fly']/@size &gt; 0">
												<a tab="Excluded from training set" href="model/browseRecords.do?mm_id={model/modelMappings/@id}&amp;set=modified-on-the-fly"><xsl:value-of select="//others/statistics[1]/set[@setId='modified-on-the-fly']/@size"/> records</a>
												<br/><a href="#" onclick="restorePoints(false); return false;" title="Restore these modified points to make them consistent with the model">[restore]</a>
											</xsl:if>
										</td>
										<xsl:choose>
										<xsl:when test="model/modelMappings/property/@qualitive = 'true'">
											<xsl:call-template name="classificationColumns"><xsl:with-param name="statisticsSet" select="//others/statistics[1]/set[@setId='modified-on-the-fly']"/></xsl:call-template>
											<xsl:if test="//others/statistics[2]">
												<xsl:call-template name="classificationColumns"><xsl:with-param name="statisticsSet" select="//others/statistics[2]/set[@setId='modified-on-the-fly']"/></xsl:call-template>
											</xsl:if>
										</xsl:when>
										<xsl:otherwise>
											<xsl:call-template name="regressionColumns"><xsl:with-param name="statisticsSet" select="//others/statistics[1]/set[@setId='modified-on-the-fly']"/></xsl:call-template>
											<xsl:if test="//others/statistics[2]">
												<td>
													<xsl:if test="//others/statistics[2]/set[@setId='modified-on-the-fly']/@size &gt; 0">
														<a tab="Excluded from training set" href="model/browseRecords.do?mm_id={model/modelMappings/@id}&amp;set=modified-on-the-fly&amp;recalculated=1"><xsl:value-of select="//others/statistics[2]/set[@setId='modified-on-the-fly']/@size"/> records</a>
														<br/><a href="#" onclick="restorePoints(true); return false;" title="Restore these modified points to make them consistent with the model">[restore]</a>
													</xsl:if>
												</td>
												<xsl:call-template name="regressionColumns"><xsl:with-param name="statisticsSet" select="//others/statistics[2]/set[@setId='modified-on-the-fly']"/></xsl:call-template>
											</xsl:if>
										</xsl:otherwise>
									</xsl:choose>
									</tr>
								</xsl:if>
							</table>
							<table class="scatter-plots">
								<tr>
									<td></td>
									<td>
										<xsl:if test="//others/statistics[2]">
											<img src="img/icons/left_arrow.png"/><a href="#" onclick="replaceOriginalModel(); return false;">Store the recalculated model as reference (original model)</a><br/>
											<img src="img/icons/red_cross.gif"/><a href="#" onclick="discardRecalculatedModel(); return false;">Discard the recalculated model</a><br/>
										</xsl:if>
									</td>
								</tr>
								<tr>
									<xsl:choose>
       								<xsl:when test="model/modelMappings/property/@qualitive != 'true'" >
										<td>
											<div style="position: relative;">
												<div id="placeholder1" style="width:650px;height:500px;"></div>
												<a href="javascript:plot1.toggleFullScreen();" title="Toggle full screen mode" class="full-screen-link"><img src="img/icons/full-screen.gif"/></a>
											</div>
											<div style="width:650px; text-align: center;">Measured value</div>
											<div id="placeholder1-selection"></div>
										</td>
										<xsl:if test="//others/statistics[2]">
											<td>
												<div style="position: relative;">
													<div id="placeholder2" style="width:650px;height:500px;"></div>
													<a href="javascript:plot2.toggleFullScreen();" title="Toggle full screen mode" class="full-screen-link"><img src="img/icons/full-screen.gif"/></a>
												</div>
												<div style="width:650px; text-align: center;">Measured value</div>
												<div id="placeholder2-selection"></div>
											</td>
										</xsl:if>
				       				</xsl:when>
				       				<xsl:otherwise>
										<td><div id="placeholder1"></div></td>
										<xsl:if test="//others/statistics[2]">
											<td><div id="placeholder2"></div></td>
										</xsl:if>
										<br/>
										<a href="#" id="roccurve">Show ROC curves</a>
										<div id="rocarea"></div>
									</xsl:otherwise>
									</xsl:choose>
								</tr>
							</table>
							<xsl:if test="others/statistics/set/point[@error]">
							<div class="warning">
								Number of compounds ignored because of <a href="javascript:errorsShow(false)">errors in original model</a> = 
								<xsl:apply-templates select="//others/statistics[1]" mode="errorNum"/>
								<xsl:if test="//others/statistics[2]">
								, 
								<a href="javascript:errorsShow(true)"> errors in recalculated model</a> = 
								<xsl:apply-templates select="//others/statistics[2]" mode="errorNum"/>.
								</xsl:if>
							</div>
							</xsl:if>
							<xsl:if test="others/statistics/deletedRecordsCount &gt; 0">
							<div class="warning">
								There were <a href="epbrowser/show.do?recordsDeletedFromModel={model/modelMappings/@id}" tab="Records, deleted after the model was created"><xsl:value-of select="others/statistics/deletedRecordsCount"/> records</a> that were deleted from the training set after the model was calculated
							</div>
							</xsl:if>
							
							<xsl:if test="count(others/statistics/set[hasPredicates='true']) &gt; 0">
								<input type="checkbox" name="consider-predicates" filter="1"/><label>Account for predicates ("&lt;", "&gt;" or intervals)</label>
							</xsl:if>
							
							<xsl:if test="others/statistics/hasDuplicates &gt; 0">
								<span><a href="model/profile.do?id={model/@id}&amp;exclude-duplicates=1&amp;render-mode=popup">
                                   [Exclude duplicated records]</a></span>
							</xsl:if>							 
								
							<xsl:if test="count(others/statistics/set[classificationSummary/tp]) &gt; 0">
								<input type="checkbox" name="consider-predicates" filter="1"/><label>Use optimal threshold for balanced accuracy</label>
							</xsl:if>
							
							<div class="actions">
								<span><a href="model/exportModel.do?id={model/@id}&amp;render-mode=popup" tab="Download descriptors and model statistics"><img src="img/icons/xls.gif"/>  Export this model</a></span>
								<xsl:if test="/model/model/owner = 'true'">
									<span><a href="model/createModelCopy.do?id={model/@id}" tab="Model copy"><img src="img/icons/clone.gif"/>  Create a copy of this model</a></span>
								</xsl:if>
								<span><a href="javascript:xmlDialog.show();" title="Show the full configuration of the model in an XML format"><img src="img/icons/xml.jpg"/>  View configuration XML</a></span>
								<span><a href="model/exportModelXml.do?id={model/@id}" title="Export the full configuration of the model in an XML format"><img src="img/icons/xml.jpg"/>  Export configuration XML</a></span>
								<span><a href="mmpqsar/model.do?id={model/@id}&amp;property={model/modelMappings/property/@id}" tab="MMP-based model analysis" title="Matched molecular pairs analysis of this model"><img src="img/icons/mmp-16.png"/>   MMP-based analysis (experimental)</a></span>
								
								<xsl:if test="//user/superuser = 'true'">								
									<span><a href="model/exportStandaloneModel.do?public_id={model/publicId}"><img src="img/icons/xml.jpg"/>  Export standalone model</a></span>
									<span><a href="model/exportModelFile.do?public_id={model/publicId}"><img src="img/icons/xml.jpg"/>  Export model file</a></span>
								</xsl:if>
								
								<xsl:if test="model/published='true' and model/approved='false'">
								<span>
									<a onclick="modelApprove.dialog(); return false;" href="javascript:void()"><img src="img/icons/approve16.png"/>  Approve this model</a>
								</span>
								</xsl:if>
								
								<xsl:if test="//session/developer = 'true'">
									<a href="model/attachment.do?id={model/@id}" target="_blank">Export model XML</a>
								</xsl:if>
							</div>
						</div>
						<xsl:if test="others/ad">
							<div>
								Williams plot with <i><xsl:value-of select="others/ad/@name"/></i> used as a distance to model.
								<xsl:if test="//others/ad/useStandardizedResiduals = 'true'">
									This plot uses standartized (studentized) residuals. <a tab="Wiki" class="wikilink"  href="http://en.wikipedia.org/wiki/Studentized_residual"><img src="img/icons/wiki2.gif"/></a> 
								</xsl:if><br/>
								<table class="ad">
									<tr>
									<td valign="top">
									<b>Distance to model<a class="infolink" href="https://docs.ochem.eu/display/MAN/Applicability+domain+assessment" title="Distance ot model is a measure of prediction uncertainty. Click to read more." target="_blank"></a></b>
									<select name="dmname" filter="1">
										<xsl:for-each select="others/statistics[last()]/set[1]/dm">
											<option value="{.}">
												<xsl:if test=".=//others/ad/@name">
													<xsl:attribute name="selected">selected</xsl:attribute>
												</xsl:if>
												<xsl:value-of select="."/>
											</option>
										</xsl:for-each>
										<xsl:if test="count(others/statistics[1]/set[1]/dm) > 10">
											<option value="all">Show all
												<xsl:if test="count(//others/ad) > 1">
													<xsl:attribute name="selected">selected</xsl:attribute>
												</xsl:if>
											</option>
										</xsl:if>
									</select>
									<b>Averaging type</b>
									<select name="averaging" filter="1">
										<option value="">Default</option>
										<option value="mgd">Bin-based (MGD)</option>
										<option value="ma">Sliding window</option>
										<option value="cumulative">Cumulative</option>
									</select>
									<xsl:if test="//others/ad/ad-configuration/averagingType = 'ma'">
										<b>Window size:</b>
										<select id="wsize" name="wsize" filter="1">
											<option value="3">3% of the set</option>
											<option value="5">5% of the set</option>
											<option value="10" selected="selected">10% of the set</option>
											<option value="20">20% of the set</option>
										</select>
										</xsl:if>
									
									<b>X axis</b>
									<select name="xaxis" filter="1">
										<option value="dm"><xsl:value-of select="//others/ad/@name"/></option>
										<option value="percents">Percentage of compounds</option>
									</select>
									<br/>
									<xsl:if test="//others/ad/ad-configuration/averagingType = 'mgd'">
										<input type="checkbox" name="show-negative" filter="1"/><label>Show negative residuals</label><br/>
									</xsl:if>
									<br/>
									
									</td>
									<td>
										<div id="leverage" style="width:600px; height:250px; margin-left: 5px;"></div>
										<div style="width:600px; margin-left: 5px; text-align: center;" id="x-axis-description">Distance to model</div>
										<div id="leverage-selection"></div>
									</td>
									<td rowspan="2" valign="top" style="padding-left: 20px;">
										<xsl:if test="//others/statistics/set[@setId='validation']">
											<a href="#" onclick="showEstimatedStatistics(); return false;">Show the estimated predictive statistics for the validation sets</a> 
											<a class="infolink" help="help-1"></a>
											<div id="help-1" class="invisible">
												This feature allows to compare the <i>estimated</i> prediction quality for the validation sets with the <i>actual</i> one.<br/>
												This information enables you to judge how the fairness of the prediction quality estimates.
											</div>
											<div id="estimated-statistics">
											</div>
										</xsl:if> 
									 </td>
								</tr>
								<tr>
									<td colspan="2">
										<xsl:if test="//others/ad/ad-configuration/averagingType = 'mgd'">
											Select outliers with p-value less or equal than 
											<select name="p-value">
												<option value="">---</option>
												<option value="0.2">0.2</option>
												<option value="0.1">0.1</option>
												<option value="0.05">0.05</option>
												<option value="0.01">0.01</option>
												<option value="0.005">0.005</option>
												<option value="0.001">0.001</option>
												<option value="0.0001">0.0001</option>
												<option value="0.00001">0.00001</option>
												<option value="0.000001">0.000001</option>
												<option value="0.0000001">0.0000001</option>
											</select>
										</xsl:if>
									</td>
								</tr>
								</table>
							</div>
						</xsl:if>
					</div>
				</div>
				
			<xsl:if test="not(model/taskId)"><xsl:call-template name="recalculate"/></xsl:if>
			</td>
		</tr>
		</table>		
		
		<div id="xmlDialog">
			<div class="hd">XML Configuration of the model</div> 
   		 	<div class="bd">
   		 		<div style="overflow: auto; height: 400px;">
					<pre>
						<xsl:value-of select="model/configurationXml"/>
					</pre>
			</div>
			</div>
		</div>
		
		<script language="javascript">
			mid = "";
			qualitativeProperty = true;
			valSetId = '<xsl:value-of select="//others/statistics/set[@setId='validation']/basketId"/>';
			
			<xsl:if test="model/@id">
				mid = "&amp;id=<xsl:value-of select="model/@id"/>&amp;mm_id=<xsl:value-of select="model/modelMappings/@id"/>&amp;mapping_id=<xsl:value-of select="model/modelMappings/@id"/>&amp;model_id=<xsl:value-of select="model/@id"/>";
				modelId = <xsl:value-of select="model/@id"/>;
			</xsl:if>
			
			
			
			<xsl:if test="model/modelMappings/property/@qualitive != 'true'">
				qualitativeProperty = false;
				var dt;
				var dtDeleted;
				var plot1 = new Plot();
				plot1.selector = "#placeholder1";
				SelectionHandler.call(plot1);
				
				<xsl:for-each select="//others/statistics[1]/set">
					<xsl:call-template name="rpData"/>
					plot1.addData(dt, setColors["<xsl:value-of select="@setId"/>"]);
					plot1.currentData.setId = <xsl:value-of select="position()-1"/>;
					if (dtDeleted.length > 0)
					{
						plot1.addData(dtDeleted, "#666");
						plot1.currentData.setId = <xsl:value-of select="position()-1"/>;
					}
				</xsl:for-each>
				
				
				var plot2;
				<xsl:if test="//others/statistics[2]">
					var plot2 = new Plot();
					plot2.selector = "#placeholder2";
					SelectionHandler.call(plot2);
					<xsl:for-each select="//others/statistics[2]/set">
						<xsl:call-template name="rpData"/>
						plot2.addData(dt, setColors["<xsl:value-of select="@setId"/>"]);
						plot2.currentData.setId = <xsl:value-of select="position()-1"/>;
						if (dtDeleted.length > 0)
						{
							plot2.addData(dtDeleted, setColors.deleted);
							plot2.currentData.setId = <xsl:value-of select="position()-1"/>;
						}
					</xsl:for-each>
					
					PlotSynchronizer.call(plot1);
					PlotSynchronizer.call(plot2);
					plot1.createPointHash();
					plot2.createPointHash();
					plot1.friendlyPlots.push(plot2);
					plot2.friendlyPlots.push(plot1);
					
				</xsl:if>
			</xsl:if>
			
			// Draw Wilkinson-plot	
			var leveragePlot;
			<xsl:if test="others/statistics[last()]/set/dm">
	       		leveragePlot = new Plot();
	       		leveragePlot.selector = "#leverage";
	       		leveragePlot.options.legend.position = "ne";
	       		SelectionHandler.call(leveragePlot);
	       		
	       		var AD = undefined;
	       		var dtDeleted = ([]);
	       		<xsl:for-each select="//others/ad">
	       			var useStandartizedResiduals = <xsl:value-of select="useStandardizedResiduals"/>;
	       		{
	       			<xsl:if test="ad-configuration/averagingType = 'mgd'">
	       				// MGD
	       				
	       				<xsl:for-each select="//others/statistics[last()]/set">
	       					var setIdentifier = "<xsl:value-of select="@setId"/>";
	       					<xsl:if test="adConfiguration">
	       						<xsl:apply-templates select="." mode="WilkinsonData"/>
		       					
			       				leveragePlot.addData(dt, setColors[setIdentifier]).setLabel(setIdentifier);
			       				leveragePlot.currentData.setId = 0;
			       			</xsl:if>
		       			</xsl:for-each>
		       			
		       			<xsl:for-each select="//others/statistics[last()]/set">
	       					var setIdentifier = "<xsl:value-of select="@setId"/>";
	       					<xsl:if test="adConfiguration">
		       					AD = undefined;
		       					<xsl:apply-templates select="adConfiguration" mode="adData"/>
			       				leveragePlot.defaultColor = "#000";
			       				// Draw MGD lines
			       				
			       				if (setIdentifier == "training")
			       				{
			       					<xsl:if test="//others/ad/dmThreshold">
				       					var dmThreshold = <xsl:value-of select="//others/ad/dmThreshold"/>;
				       					<xsl:if test="adConfiguration/percents">
				       						dmThreshold = 95;
				       					</xsl:if>
				       					leveragePlot.drawLine(dmThreshold, leveragePlot.minY, dmThreshold, leveragePlot.maxY);
				       					leveragePlot.currentData.lines = {lineWidth: 1};
					       				leveragePlot.currentData.shadowSize = 0;
					       				leveragePlot.fillArea(-1, -1, dmThreshold, leveragePlot.maxY);
				       				</xsl:if>
			       				
			       					if (!useStandartizedResiduals)
			       					{
					       				// Draw MGD "steps" for AD of a regression model
					       				var line = function(min, max, error)
					       				{
					       					leveragePlot.drawLine(min, error, max, error);
					       					<xsl:if test="//others/ad/showNegativeValues = 'true'">
					       					leveragePlot.drawLine(min, -error, max, -error);
					       					</xsl:if>
					       				}
					       				trIntervals = AD.xAxis;
					       				trErrors = AD.errors;
						       			line(leveragePlot.minX, AD.intervals[0], AD.errors[0]);
						       			
						       			for (var i = 0; i &lt; AD.xAxis.length - 1; i++)
						       				line(AD.xAxis[i], AD.xAxis[i+1], AD.errors[i+1]);
						       			line(AD.xAxis[i], AD.xAxis[i+1], AD.errors[i+1]);
					       			}
					       			else
					       			{
					       				// We use standardtized residuals. Draw the specific features of the chart (2.5 std lines, warning leverage, etc.)
					       				var std = 1;
					       				var line = function(level)
					       				{
					       					leveragePlot.drawLine(0, level, leveragePlot.maxX, level);
					       					leveragePlot.currentData.lines = {lineWidth: 1};
					       					leveragePlot.currentData.shadowSize = 0;
					       				}
					       				line(2.5 * std);
					       				line(-2.5 * std);
					       				line(3 * std);
					       				line(-3 * std);
					       			}
					       			
					       			
					       			
					       			<xsl:if test="//others/ad/showNegativeValues = 'true'">
					       			leveragePlot.drawLine(leveragePlot.minX, 0, leveragePlot.maxX, 0);
					       			leveragePlot.currentData.lines = {lineWidth: 1};
					       			</xsl:if>
				       			}
			       			</xsl:if>
		       			</xsl:for-each>
		       			
		       			
		       			
		       			if (dtDeleted.length > 0)
		       				leveragePlot.addData(dtDeleted, setColors.deleted).setLabel("deleted");
		       			
	       			</xsl:if>
	       			<xsl:if test="ad-configuration/averagingType != 'mgd'">
	       				// Draw the reliability-vs-accuracy chart for a classification model
	       				
	       				<xsl:for-each select="(//others/statistics)[last()]/set">
	       					var setIdentifier = "<xsl:value-of select="@setId"/>";
	       					AD = undefined;
	       					<xsl:if test="adConfiguration">
			       				<xsl:apply-templates select="adConfiguration" mode="adData"/>
			       				var dt = [];
			       				for (var i = 0; i &lt; AD.xAxis.length - 1; i++)
					       			dt.push([AD.xAxis[i], AD.errors[i], AD.pointNums[i]]);
					       		leveragePlot.addData(dt, setColors[setIdentifier]).setLabel(setIdentifier + ", <xsl:value-of select="@name"/>");
					       		leveragePlot.setPoints({line: true});
					       		leveragePlot.currentData.lines = {lineWidth: 1};
				       		</xsl:if>
			       		</xsl:for-each>
	       			</xsl:if>
	       		}
	       		</xsl:for-each>
       		</xsl:if>
       		
			      		
       		$(document).ready(function()
       		{
       			$("#roccurve").click(function(){
       				if ($("#roccurve").html() == "Show ROC curves")
       				{
       					$("#roccurve").html("Hide ROC curves");
	       				var valSetId = '<xsl:value-of select="//model/selectedValidationSet/@id"/>';
	       				$("#rocarea").html('<img src="img/roller.gif"/>');
	       				var ajax = new QSPR.Ajax("model/rocCurve.do?validation=" + valSetId + mid);
						ajax.send({
							data: "",
							success: function(reply)
							{
								if ($("#roccurve").html() == "Show ROC curves")
									return;
								$("#rocarea").html('<div id="roc" style="width: 450px; height: 450px;"></div>');
								var rocCurve = new Plot();
		       					rocCurve.selector = "#roc";
		       					rocCurve.options.legend.position = "ne";
		       					rocCurve.drawLine(0, 0, 1, 1);
						       	rocCurve.currentData.lines = {lineWidth: 1};
						       	rocCurve.currentData.shadowSize = 0;
						       	var rocData = array(reply.others.rocCurve);
		       					for (var set = 0; set &lt; rocData.length; set++)
		       					{
			       					var dt = [];
			       					var rocPoints = array(rocData[set].points);
				       				for (var i = 0; i &lt; rocPoints.length; i++)
						       			dt.push([parseFloat(rocPoints[i].x), parseFloat(rocPoints[i].y), parseInt(rocPoints[i].id)]);
				       				rocCurve.addData(dt, setColors[rocData[set].setId]).setLabel(rocData[set].setId);
				       				rocCurve.currentData.setId = set;
				       				rocCurve.setPoints({line: true});
						       		rocCurve.currentData.lines = {lineWidth: 1};
					       		}
			       				rocCurve.render("#roc");
							}
						});
					} else
					{
						$("#roccurve").html("Show ROC curves");
						$("#rocarea").html("");
					}
					return false;
				});
       		
       		
       		
       			xmlDialog = new YAHOO.widget.Dialog("xmlDialog", { 
			    	width:"700px", 
					fixedcenter:true, 
					modal:true, 
					visible:false 
			    });
			    
			    xmlDialog.render();
       			
       			tabView = new YAHOO.widget.TabView("demo");
				
				if (getParams["dmname"] || getParams["wsize"] || getParams["xaxis"])
					tabView.set("activeIndex", 1);
					
       			if (leveragePlot)
				{
					leveragePlot.render("#leverage");
					//selectOutliers(leveragePlot, trIntervals, trErrors);
					leveragePlot.pointClicked = plot2clicked;
					var quantiles = new Array();
					quantiles["0.0000001"] = 5.326;
					quantiles["0.000001"] = 4.891;
					quantiles["0.00001"] = 4.417;
					quantiles["0.0001"] = 3.890;
					quantiles["0.001"] = 3.290;
					quantiles["0.005"] = 2.807;
					quantiles["0.01"] = 2.575;
					quantiles["0.02"] = 2.326;
					quantiles["0.05"] = 1.960;
					quantiles["0.1"] = 1.645;
					quantiles["0.2"] = 1.282;
					$('select[name="p-value"]').change(function()
					{
						if ($(this).val() == "")
							leveragePlot.unselectAll();
						else
							selectOutliers(leveragePlot, trIntervals, trErrors, quantiles[$(this).val()], $(this).val());	
					});
				}
				
				$("[filter]").change(function()
				{
					var input = $(this);
					var name = input.attr("name");
					var append = "";
					if (input.attr("type") == "checkbox")
					{
						if (input.is(":checked"))
							append = name + "=1";
					}
					else
					{
						var value = input.attr("value");
						append = name + "=" + value;
					}
					
					var rg = new RegExp(name + "=[^&amp;]*");
					var newLocation = window.location.href.replace(rg, "") + "&amp;" + append + "&amp;tab=" + tabView.get("activeIndex");
					if (newLocation.indexOf("render-mode") == -1)
						newLocation += "&amp;render-mode=popup";
					window.location.href = newLocation;
				});
				
				$("[filter]").each(function()
				{
					var value = getParams[$(this).attr("name")];
					if (value)
					{
						$(this).val(value);
						if ($(this).attr("type") == "checkbox")
							$(this).attr("checked", "checked");
					}
					if ($(this).attr("name") == "xaxis" &amp;&amp; (value == "percents"))
						$("#x-axis-description").html("Distance to model (in the percentage scale)");
				});
				
				if (getParams["tab"])
					tabView.set("activeIndex", getParams["tab"]);
				
       			<xsl:choose>
       				<xsl:when test="model/modelMappings/property/@qualitive != 'true'" >
       					plot1.addBaseLine().render("#placeholder1");
       					
						plot1.pointClicked = plot1clicked;
						if (plot2 != undefined)
						{
							plot2.addBaseLine().render("#placeholder2");
							plot2.pointClicked = plot2clicked;
						}
						
       				</xsl:when>
       				<xsl:otherwise>
       					var num = 1;
       					var optionsName = new Array();
	       				var mm_id = <xsl:value-of select="//modelMappings/@id"/>
	       				var model_id = <xsl:value-of select="//model/@id"/>
       					<xsl:for-each select="//others/option">
       						optionsName.push('<xsl:value-of select="@name"/>');
       					</xsl:for-each>
       					<xsl:for-each select="//others/statistics">
       					    var statName = "";
       						var placeholder = $("#placeholder"+num);
       						if(num == 1)
       							statName = "(Original)";
       						else
       							statName = "(Recalculated)";
       						var mainTable = $('<table class="classificationSummary"></table>');
							var maintr;
							mainTable.append(maintr = $("<tr></tr>"));
							var set = 0;
							<xsl:for-each select="set">
								var setName = "";
								if(set == 0)							
									setName = "Training "+statName;
								else if(set == 1)
									setName = "Test "+statName;
								else
									setName = "Excluded "+statName;
								// This is a classification model. Draw a confision matrix
								var maintd;
								var table ;
								var tr;
								var td;
								maintr.append(maintd = $("<td></td>"));
								
								var maxOption = 0;
								<xsl:for-each select="classificationSummary/nodes">
									if(<xsl:value-of select="@real"/> &gt; maxOption)
										maxOption = parseInt(<xsl:value-of select="@real"/>);
								</xsl:for-each>
								
								//one additional for table header		
								maxOption = maxOption + 1;
									
								maintd.append(table = $('<table class="tiny"></table>'));
								
								
								var specsens = $('<tr>0</tr>');
								
								for (var i = 0; i &lt;= maxOption+1; i++) {
									table.append(tr = $("<tr></tr>"));
									for (var k = 0; k &lt;= maxOption+1; k++)
										if (i != maxOption + 1 || k != maxOption + 1)
											tr.append('<td>0</td>');
								}
								var col = maxOption +1;
								var inf = $('<tr></tr>');
								inf.append('<td colspan="'+(col+1)+'">'+setName+'</td>');
								table.append(inf);
								
								placeholder.append(mainTable);
								
								//set table row header
								table.find("tr").eq(maxOption + 1).append("<td>&#160;</td>");
								table.find("tr").eq(0).find("td").eq(0).html("Real&#8595;/Predicted&#8594;");
								for (var i = 1; i &lt;= maxOption; i++)
								{
									table.find("tr").eq(0).find("td").eq(i).html(optionsName[(i*1)-1]).addClass('predicted');
									table.find("tr").eq(i).find("td").eq(0).html(optionsName[(i*1)-1]).addClass('real');
									table.find("tr").eq(i).find("td").eq(i).addClass('confusion-matrix');
									
									table.find("tr").eq(i).find("td").eq(maxOption + 1).addClass('statv');
								}
								table.find("tr").eq(maxOption + 1).find("td").filter(":not(:first)").addClass('statv');
								
								table.find("tr").eq(maxOption + 1).find("td").eq(0).html("Precision").addClass('stath');
								table.find("tr").eq(0).find("td").eq(maxOption + 1).html("Hit rate").addClass('stath');
								
								<xsl:for-each select="classificationSummary/nodes">
									var real = <xsl:value-of select="@real"/>;
									var predicted = <xsl:value-of select="@predicted"/>;
									var pointSelectors = [];
									if (getParams['best-predictions'])
										pointSelectors.push("best-predictions=" + getParams['best-predictions']);
									if (getParams['dm-threshold'])
										pointSelectors.push("dm-threshold=" + getParams['dm-threshold']);
									var pointSelectorStr = pointSelectors.join("&amp;");
									var htmlText = '<a title="Click to see the records" tab="Records from a model" href="epbrowser/show.do?mm_id='+mm_id+'&amp;real='+real+'&amp;predicted='+predicted+'&amp;num='+num+'&amp;set=' + set + '&amp;validationSetNum=' + $('select[name=validation]').val()+ '&amp;' + pointSelectorStr + '">'+<xsl:value-of select="@count"/>+'</a>';
									table.find("tr").eq(parseInt(real)+1).find("td").eq(parseInt(predicted)+1).html(htmlText);
								</xsl:for-each>
								
								<xsl:for-each select="classificationSummary/sensitivityByClass">
									table.find("tr").eq(<xsl:value-of select="position()"/>).find("td").eq(maxOption + 1).html(<xsl:value-of select="@formatted-value"/>);
								</xsl:for-each>

								<xsl:for-each select="classificationSummary/valueByClass">
									table.find("tr").eq(maxOption + 1).find("td").eq(<xsl:value-of select="position()"/>).html(<xsl:value-of select="@formatted-value"/>);
								</xsl:for-each>
								
								$(document).trigger("DOM_updated", $(document));
								$("#placeholder" + num + " a[title]").tooltip({showURL: false});
								set++;
							</xsl:for-each>
       						num++;
       					</xsl:for-each>
       				</xsl:otherwise>
       			</xsl:choose>
       			
       			
       		});
       		
       		
			if (typeof String.prototype.startsWith != 'function') 
   			{
				String.prototype.startsWith = function (str)
				{
					return this.indexOf(str) == 0;
				};
			}
       		
			
		</script>
	
	</xsl:template>
	
	<xsl:template name="numRecords">
		<xsl:param name="statisticsSet"/>
		<xsl:if test="$statisticsSet/@size != $statisticsSet/@matchedPoints">
			<br/><small title="Only {$statisticsSet/@matchedPoints} points were included for statistics calulation">(<xsl:value-of select="$statisticsSet/@matchedPoints"/> selected)</small>
		</xsl:if>
	</xsl:template>
	
	<xsl:template name="classificationColumns">
		<xsl:param name="statisticsSet"/>
		<td>
			<xsl:value-of select="$statisticsSet/classificationSummary/accuracyTotal/@formatted-value"/>%
			<xsl:if test="$statisticsSet/classificationSummary/accuracyTotal/@formatted-std != 'NaN'">
				± <xsl:value-of select="$statisticsSet/classificationSummary/accuracyTotal/@formatted-std"/>
			</xsl:if>
		</td>
		<td>
			<xsl:value-of select="$statisticsSet/classificationSummary/accuracyBalanced/@formatted-value"/>%
			<xsl:if test="$statisticsSet/classificationSummary/accuracyBalanced/@formatted-std != 'NaN'">
				± <xsl:value-of select="$statisticsSet/classificationSummary/accuracyBalanced/@formatted-std"/>
			</xsl:if>
		</td>
		<td>
			<xsl:value-of select="$statisticsSet/classificationSummary/mcc/@formatted-value"/>
			<xsl:if test="$statisticsSet/classificationSummary/mcc/@formatted-std != 'NaN'">
				± <xsl:value-of select="$statisticsSet/classificationSummary/mcc/@formatted-std"/>
			</xsl:if>
		</td>
		<td>
			<xsl:value-of select="$statisticsSet/classificationSummary/auc/@formatted-value"/>
			<xsl:if test="$statisticsSet/classificationSummary/auc/@formatted-std != 'NaN'">
				± <xsl:value-of select="$statisticsSet/classificationSummary/auc/@formatted-std"/>
			</xsl:if>
		</td>
	</xsl:template>	
	
	
	<xsl:template name="regressionColumns">
		<xsl:param name="statisticsSet"/>
		<td>
			<xsl:value-of select="$statisticsSet/r2"/>
			<xsl:if test="$statisticsSet/r2std">
				± <xsl:value-of select="$statisticsSet/r2std"/>
			</xsl:if>
		</td>
		<td>
			<xsl:value-of select="$statisticsSet/q2"/>
			<xsl:if test="$statisticsSet/q2std">
				± <xsl:value-of select="$statisticsSet/q2std"/>
			</xsl:if>
		</td>
		<td>
			<xsl:value-of select="$statisticsSet/rmse"/>
			<xsl:if test="$statisticsSet/rmsestd">
				± <xsl:value-of select="$statisticsSet/rmsestd"/>
			</xsl:if>
		</td>
		<td>
			<xsl:value-of select="$statisticsSet/mae"/>
			<xsl:if test="$statisticsSet/maestd">
				± <xsl:value-of select="$statisticsSet/maestd"/>
			</xsl:if>
		</td>	
	</xsl:template>
	
	<xsl:template name="factors">
		<div class="tc">
		<table id="regfactors">
			<tr><th>Descriptor</th><th>Regression factor(original)</th><th>Regression factor(recalculated)</th></tr>
			<xsl:for-each select="//chosenDescriptors/value">
				<tr><td><xsl:value-of select="."/></td><td><xsl:value-of select="format-number(//factors/value[position()], '###,##0.00')"/></td><td><xsl:value-of select="format-number(//recalculatedFactors/value[position()], '###,##0.00')"/></td></tr>
			</xsl:for-each>
		</table>
		</div>
	</xsl:template>
	
	<xsl:template name="rpData">
			dt = ([
				<xsl:for-each select="point">
					<xsl:if test="not(@error)">
					[<xsl:value-of select="real"/>, <xsl:value-of select="predicted"/>, <xsl:value-of select="id"/>],
					</xsl:if>
				</xsl:for-each>
				[0, 0, 0]
				]);	
				dt.pop();
				
			dtDeleted = ([
				<xsl:for-each select="point">
					<xsl:if test="not(@error) and @deleted">
					[<xsl:value-of select="real"/>, <xsl:value-of select="predicted"/>, <xsl:value-of select="id"/>],
					</xsl:if>
				</xsl:for-each>
				[0, 0, 0]
				]);	
			dtDeleted.pop();
	</xsl:template>
	
	<xsl:template match="*" mode="adData">
		var AD = new Object();
		AD.intervals = [
			<xsl:for-each select="interval">
				<xsl:value-of select="."/><xsl:if test="position() != last()">, </xsl:if>
			</xsl:for-each>
			];
		AD.errors = [
			<xsl:for-each select="error">
				<xsl:value-of select="."/><xsl:if test="position() != last()">, </xsl:if>
			</xsl:for-each>
			];
		
		AD.percents = [
			<xsl:for-each select="percents">
				<xsl:value-of select="."/><xsl:if test="position() != last()">, </xsl:if>
			</xsl:for-each>
			];
			
		AD.xAxis = AD.percents.length > 0 ? AD.percents : AD.intervals;
		
		AD.pointNums = [
			<xsl:for-each select="epId">
				<xsl:value-of select="."/><xsl:if test="position() != last()">, </xsl:if>
			</xsl:for-each>
			];
	</xsl:template>
	
	<xsl:template match="*" mode="WilkinsonData">
			dt = ([
				<xsl:for-each select="point[ad]">
					[<xsl:value-of select="ad/@dmValue"/>, <xsl:value-of select="ad/@error"/>, <xsl:value-of select="id"/>],
				</xsl:for-each>
				[0, 0, 0]
				]);	
			dt.pop();
			
			dt.sort(function(a,b){
				return a[0]-b[0];
			});
			<xsl:if test="adConfiguration/percents">
				var realCount = 0;
				var previous = -100500;
				for (var i=0; i &lt; dt.length; i++)
				{
					if (Math.abs(dt[i][0] - previous) &gt; 0.0001)
					{
						previous = dt[i][0];
						realCount += 1;
					}
				}
				
				var currentCount = 0;
				previous = -100500;
				for (var i=0; i &lt; dt.length; i++)
				{
					if (Math.abs(dt[i][0] - previous) &gt; 0.0001)
					{
						previous = dt[i][0];
						currentCount += 1;
					}	
					dt[i][0] = currentCount * 100 / realCount;
				}
			</xsl:if>
			
			<xsl:for-each select="point[ad]">
				<xsl:if test="not(@error) and @deleted">
					dtDeleted.push([<xsl:value-of select="ad/@dmValue"/>, <xsl:value-of select="ad/@error"/>, <xsl:value-of select="id"/>]);
				</xsl:if>
			</xsl:for-each>
	</xsl:template>
	
	<xsl:template match="*" mode="errorNum">
				<xsl:value-of select="count(set/point[@error])"/>
	</xsl:template>
	<!--name="lf2br" IS REMOVED -->
	<xsl:template name="lf2br">
	  <xsl:param name="StringToTransform"/>
	  <xsl:choose>
	    <xsl:when test="contains($StringToTransform,'&#10;')">
	      <xsl:value-of select="substring-before($StringToTransform,'&#10;')"/>
	       <br/>
	      <xsl:call-template name="lf2br">
	        <xsl:with-param name="StringToTransform" select="substring-after($StringToTransform,'&#10;')"/>
	      </xsl:call-template>
	    </xsl:when>
	    <xsl:otherwise>
	      <xsl:value-of select="$StringToTransform"/>
	    </xsl:otherwise>
	  </xsl:choose>
	</xsl:template>
</xsl:stylesheet>