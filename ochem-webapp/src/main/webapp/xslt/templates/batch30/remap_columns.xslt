<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
	version="2.0">
	<xsl:include href="../../helper.xslt" />

	<xsl:template name="content">
		<script type="text/javascript" src="js/lib/jquery-ui-1.10.3.custom.min.js" />
		<link rel="stylesheet" type="text/css" href="css/smoothness/jquery-ui-1.10.3.custom.min.css" />
		<link rel="stylesheet" type="text/css" href="css/batch.css" />
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script language="javascript" src="js/blocks/batch-upload-30.js"></script>

		<script language="javascript">
			remapping = new ColumnRemappingActionable();
			$(document).ready(function(){
				$("#tabs li > a").each(function(){ this.href = location.href + this.hash });
				remapping.init();
				$("#tabs").tabs();
				$("#tabs").on( "tabsactivate", function(event, ui){	remapping.scopeBlock = ui.newPanel;	});
				<xsl:if test="//file/@selectedSheet != ''">
					$( "#tabs" ).tabs( "option", "active", <xsl:value-of select="//file/@selectedSheet" /> );
				</xsl:if>
				remapping.callAction("validate", null, null);
			});
		</script>
		<table class="layouttable">
			<tr><td class="itunes-up silver">
					<img src="img/icons/batchupload.png" />
					<h1><doc term="Batch+data+upload">Batch Upload 3.0 - File preview and column remapping</doc></h1>
					Preview your data, select the sheet and the columns you would like to upload
			</td></tr>
			<tr>
				<td style="padding: 10px 10px 10px 10px; vertical-align:top;">
					<div id="tabs">
						<ul>
							<xsl:for-each select="//sheets">
								<xsl:variable name="currentSheet">
									<xsl:value-of select="position()-1" />
								</xsl:variable>
								<li>
									<a href="#tabs-{position()}">
										<em>
											<xsl:value-of select="@name" />
										</em>
									</a>
								</li>
							</xsl:for-each>
						</ul>
						<xsl:for-each select="//sheets">
							<xsl:variable name="currentSheet">
								<xsl:value-of select="position()-1" />
							</xsl:variable>
							<div class="preview" id="tabs-{position()}">
								<form action="batchupload30/remap_columns_submit.do"
									method="post">
									<input type="hidden" send="1" name="sheet" value="{$currentSheet}" />
									<table>
										<tr>
											<xsl:for-each select="columns">
												<xsl:choose>
												<xsl:when test="@hidden = 'false'">
													<td>
														<table>
															<tr>
																<td>
																	<xsl:choose>
																		<xsl:when test="@type = 'undefined'">
																			<xsl:attribute name="class">red</xsl:attribute>
																		</xsl:when>
																		<xsl:when test="@type = 'property_and_value'">
																			<xsl:attribute name="class">dgreen</xsl:attribute>
																		</xsl:when>
																		<xsl:when test="@type = 'condition_and_value'">
																			<xsl:attribute name="class">dgreen</xsl:attribute>
																		</xsl:when>
																		<xsl:otherwise>
																			<xsl:attribute name="class">green</xsl:attribute>
																		</xsl:otherwise>
																	</xsl:choose>
																	<input type="checkbox" send="1" name="select{position()-1}"
																		id="select{position()-1}">
																		<xsl:if test="@ignore = 'false'">
																			<xsl:attribute name="checked">checked</xsl:attribute>
																		</xsl:if>
																	</input>
																	<a action="menu">
																		<xsl:value-of select="@name" />
																	</a>
																	<xsl:if test="@subscriptName">
																		(<xsl:value-of select="@subscriptName" />)
																	</xsl:if>																	
																	<br/>
																	<input type="hidden" send="1" name="name" value="{@name}" />
																</td>
															</tr>
															<tr>
																<td class="small">
																	<xsl:if test="@defaultValue != ''">
																		Default:
																		<xsl:value-of select="@defaultValue" />
																	</xsl:if>
																	&#160;
																</td>
															</tr>
															<xsl:for-each select="sampleValues">
																<tr>
																	<td class="grey">
																		<xsl:choose>
																			<xsl:when test=".=''">
																				&#160;
																			</xsl:when>
																			<xsl:otherwise>
																				<xsl:value-of select="." />
																			</xsl:otherwise>
																		</xsl:choose>
																	</td>
																</tr>
															</xsl:for-each>
														</table>
													</td>
												</xsl:when>
												<xsl:otherwise>
													<input type="hidden" send="1" name="select{position()-1}" value="true" />
													<input type="hidden" send="1" name="name" value="{@name}" />
												</xsl:otherwise>
												</xsl:choose>
											</xsl:for-each>
										</tr>
									</table>
									<br />
									<div name="adddummy" class="invisible">
										<input type="checkbox" send="1" name="adddummy"/> Upload molecules without experimental data
										<a class="infolink" title="Check this if you would upload a set of molecules rather than a set of experimental measurements"></a>
									</div>
									<div name="messages">
										<xsl:for-each select="messages">
											<div>
												<xsl:choose>
													<xsl:when test="type = 'error'">
														<b class="error">
															<xsl:value-of select="message" />
														</b>
													</xsl:when>
													<xsl:when test="type = 'warning'">
														<b class="warn">
															<xsl:value-of select="message" />
														</b>
													</xsl:when>
													<xsl:otherwise>
														<b class="notice">
															<xsl:value-of select="message" />
														</b>
													</xsl:otherwise>
												</xsl:choose>
											</div>
										</xsl:for-each>
									</div>
									<br />
									Green titles indicate recognized columns, red titles indicate
									errors.
									Please click on the red columns and select whether the column
									indicates a property, condition or another column type like
									name, value or molecule,
									then select the matching entity and confirm your selection by
									clicking on the green button on the left.
									<br />
									If you have irrelevant columns in your sheet, you can leave
									them red and they will be ignored in the further process.
									If you need help, feel free to drop us an e-mail at
									<a style="color: #3333BB; text-decoration: none;" href="mailto:info@ochem.eu">info@ochem.eu.</a>
									<br />
									<br />
									<div name="error">
										<xsl:if test="@sheetStatus != 'error'">
											<xsl:attribute name="class">invisible</xsl:attribute>
										</xsl:if>
										<b class="error">Please resolve the critical errors with the sheet
											to proceed with the upload</b>
										<br />
										<input type="submit" name="submit" disabled="true"
											value="Upload this sheet" />
									</div>
									<div name="warning">
										<xsl:if test="@sheetStatus != 'warning'">
											<xsl:attribute name="class">invisible</xsl:attribute>
										</xsl:if>
										<b class="warn">There are warnings regarding this sheet</b>
										<br />
										<input type="submit" name="submit" value="Upload this sheet" />
									</div>
									<div name="notice">
										<xsl:if test="@sheetStatus != 'notice'">
											<xsl:attribute name="class">invisible</xsl:attribute>
										</xsl:if>
										<input type="submit" name="submit" value="Upload this sheet" />
									</div>
								</form>
							</div>
						</xsl:for-each>
					</div>
				</td>
			</tr>
			<tr>
				<td class="bubutton itunes-right big-padding ui-widget right">
					<a name="cancel" href="batchupload30/cancel.do">Cancel Batch Upload</a>
					<a name="downloadexcel" href="batchupload30/report.do">Download Excel file</a>
				</td>
			</tr>
		</table>
		<ul id="menu" style="display: none;">
			<li>
				<a action="mapkc">Known column</a>
			</li>
			<li>
				<a action="mappr">Property</a>
			</li>
			<li>
				<a action="mapco">Condition</a>
			</li>
		</ul>
		<div id="knowncolumns" title="Select column name" style="display: none;">
			<select>
				<xsl:for-each select="others/list[@name='known']/string">
					<option value="{.}">
						<xsl:value-of select="." />
					</option>
				</xsl:for-each>
			</select>
		</div>
		<div id="selenium-batch-upload-page-2" class="invisible"/>
	</xsl:template>
</xsl:stylesheet>