<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="2.0">
	<xsl:include href="../helper.xslt" />
	<xsl:include href="inc/license.xslt" />

	<xsl:template name="content">
		<style type="text/css">
			.style {padding: 10px; width: 900px; margin: 10px 10px 10px 50px;}
			.trStyle {border: 1px none #999;}
			.trStyle > TD {padding-top: 10px; padding-bottom: 20px;}
			.trStyle.header {font-weight: bold; border-bottom: 1px solid black; text-align: left !important; }
			.trStyle.header TD {padding-top: 15px; padding-bottom: 0px;}
			TD {padding: 3px;}
			LABEL {width: 200px}
			INPUT, TEXTAREA {border: 1px solid black; padding: 2px 2px 2px 2px; margin: 1px 1px 1px 1px;}
			.warning {color: red !important;}
			.warn {color: red !important; font-size: 11px; font-style: italic;}
			.redirect {background-color: #FFDDDD; font-size: 90%; padding: 10px; color: red;}
			.message {font-size: 11px; font-style: italic;}
			<xsl:if test="not(@id)">
				.required {background-color: #FFA;}
			</xsl:if>
					
		</style>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script language="javascript" src="js/blocks/user.js"></script>
		<xsl:apply-templates select="/model/user"/>
	</xsl:template>
	<xsl:template match="/model/user">
		<title>Account</title><table height="100%" width="100%">
			<tr>
				<td class="itunes-up" colspan="2">
				<img src="img/icons/article.png"/>
					<h1>User account</h1>
					Details of your personal OCHEM account
				</td>
			</tr>
			<tr>
				<td class="itunes-right">
					<table class="style">
						<tr>
							<td>
								<input type="hidden" name="userId" send="1" value="{@id}"/>
								<div class="warning invisible"><xsl:value-of select="/model/user/message"/></div>
								<xsl:if test="not(Organisation) and @id">
									<div class="redirect">Dear <xsl:value-of select="/model/session/user/Title"/>&#160;<xsl:value-of select="/model/session/user/FirstName"/>&#160;<xsl:value-of select="/model/session/user/LastName"/>,<br/><br/>
									 	 You have been redirected here due to update your profile with the form of organization you are working in. Please select the organization form you fit in best.<br/><br/>
									  	 Thank you,<br/>
									  	 your OCHEM Team.
									</div>
								</xsl:if>	
							</td>
						</tr>
						
						<tr class="trStyle header" style="text-align: center; font-size: 16px"><td>Registration Information</td></tr>
						
						<tr class="trStyle">
							<td>
								<table>
									<tr>
										<td>Login*</td>
										<td><input type="text" n-disable="1" name="login" send="1" value="{@login}" maxlength="20" class="required"/>
											<div class="message">(min. 4 characters and max. 20 characters)</div>
											<div id="check"><a action="check" class="button-link">Check availability </a></div>
										</td>
										<td class="warn" name="warn_login"></td>
									</tr>
									<tr>
										<td width="150px">e-mail*</td>
										<td><input type="text" n-disable="1" name="emailId" send="1" value="{E-mail}" class="required"/></td>
										<td class="warn" name="warn_email"></td>
									</tr>
									<xsl:if test="referral">
									<tr>
										<td>Invited by</td>
										<td><img src="img/icons/user-16.png"/><xsl:value-of select="referral/@login"/>
										</td>
									</tr>
								</xsl:if>
							<xsl:choose>
								<xsl:when test="@id">
									<tr>
										<td colspan="2"><div id="pwd"><a action="change">Change Password</a></div></td>
									</tr>
								</xsl:when>
								<xsl:otherwise>
								
									<tr>
										<td>Password*</td>
										<td><input type="password" name="passwd" send="1" class="required"/></td>
										<td class="warn" name="warn_passwd"></td><td>Password can contain only letters and numbers.</td>
									</tr>
									<tr>
										<td>Confirm password*</td>
										<td><input type="password" name="c_passwd" send="1" class="required"/></td>
										<td class="warn" name="warn_c_passwd"></td>
									</tr>
								
								</xsl:otherwise>
							</xsl:choose>
								
								</table>
							</td>
						</tr>
						
						<tr class="trStyle header" style="text-align: center; font-size: 16px"><td>Personal Information</td></tr>
					
						<tr class="trStyle">
							<td>
								<table>
									<tr>
										<td width="150px">Title*</td>
										<td>
											<xsl:variable name="titel" select="Title" />
											<select name="title" send="1">
												<xsl:for-each
													select="tokenize('-- please select --,Mr.,Miss.,Mrs.,Ms.,Dr.,Prof.,MD', ',')">
													<option>
														<xsl:if test="$titel = .">
															<xsl:attribute name="selected">selected</xsl:attribute>
														</xsl:if>
														<xsl:attribute name="value"><xsl:value-of
															select="." /></xsl:attribute>
														<xsl:value-of select="." />
													</option>
												</xsl:for-each>
											</select>
										</td>
										<td class="warn" name="warn_title"></td>
									</tr>
									<tr>
										<td width="150px">First name*</td>
										<td>
											<input type="text" name="firstname" value="{FirstName}" send="1"
												maxlength="50" class="required" />
										</td>
										<td class="warn" name="warn_firstname"></td>
									</tr>
									<tr>
										<td width="150px">Last name*</td>
										<td>
											<input type="text" name="lastname" value="{LastName}" send="1"
												maxlength="50" class="required" />
										</td>
										<td class="warn" name="warn_lastname"></td>
									</tr>
									<tr>
										<td width="150px">Affiliation*</td>
										<td>
											<input type="text" name="affiliation" value="{Affiliation}"
												send="1" maxlength="200" />
										</td>
									</tr>
									<tr>
										<td>Form of organization*</td>
										<td>
											<xsl:variable name="orga" select="Organisation" />
											<select name="organisation" send="1">
												<xsl:for-each
													select="tokenize('-- please select --,Commercial,Non profit - Governmental,Academic,Self-employed', ',')">
													<option>
														<xsl:if test="$orga = .">
															<xsl:attribute name="selected">selected</xsl:attribute>
														</xsl:if>
														<xsl:attribute name="value"><xsl:value-of
															select="." /></xsl:attribute>
														<xsl:value-of select="." />
													</option>
												</xsl:for-each>
											</select>
										</td>
										<td class="warn" name="warn_orga"></td>
									</tr>
							
								</table>
							</td>
						</tr>
						
						<tr class="trStyle">
							<td>
								<table>
									<tr><td width="150px">City</td><td><input type="text" name="city" value="{City}" send="1" maxlength="25"/></td></tr>
									<tr><td>State</td><td><input type="text" send="1" name="state" value="{State}" maxlength="25"/></td></tr>
									<tr><td>Country</td><td><input type="text" send="1" name="country" value="{Country}" maxlength="35"/></td></tr>
									<tr><td>Zip</td><td><input type="text" name="zip" value="{Zip}" send="1"/></td></tr>
									<tr><td>Phone</td><td><input type="text" name="telephone" value="{TelePhone}" send="1" maxlength="50"/></td></tr>
								</table>
							</td>
						</tr>
						
						<tr class="trStyle">
							<td>
								<table>
									<tr><td width="150px">Position</td><td><input type="text" name="occupation" value="{Occupation}" send="1" maxlength="100"/></td></tr>
									<tr><td>Web</td><td><input type="text" name="website" value="{WebSite}" send="1" maxlength="100"/></td></tr>
								</table>
							</td>
						</tr>		
						
					<xsl:if test="@id">
						<tr class="trStyle header"><td>Session information</td></tr>
						<tr class="trStyle"><td>Session ID: <xsl:value-of select="//session/@id"/></td></tr>
						<tr class="trStyle"><td>Session GUID: <xsl:value-of select="//session/guid"/></td></tr>
						<tr class="trStyle"><td><a class="fancy-button" action="update">Update</a></td></tr>
					</xsl:if>
					<xsl:if test="not(@id)">
					
						<xsl:if test="not(//param[@key='inhouse'])">
						<tr class="trStyle" style="text-align: center; font-size: 16px"><td>Terms of Service</td></tr>
						<tr class="trStyle">
							<td>
								<table>
									<tr>
										<td>
											<textarea rows="20" cols="180" style="width: 95%; height: 200px; text-align: justify;" readonly="readonly" onfocus="this.rows=10">
												<xsl:call-template name="license"/>
											</textarea>
										</td>
									</tr>
									<tr>
										<td>
											<b>By clicking on 'I accept' below I acknowledge that I have read and fully understand the foregoing information and agree to abide by 
											<a href="license_agreement.htm" tab="License agreement">License agreement</a> above and the 
											<a href="Privacy_Policy.htm" tab="Privacy policy">Privacy Policy</a>.</b>
										</td>
									</tr>
								</table>
							</td>
						</tr>
						</xsl:if>
						<tr class="trStyle">
						
							<xsl:choose>
								<xsl:when test="//param[@key='inhouse']">
									<td><a class="fancy-button" action="submit">Create my account</a></td>		
								</xsl:when>
								<xsl:otherwise>
									<td><a class="button-link" action="submit">I accept. Create my account.</a><a href="login/show.do?render-mode=popup" class="button-link">I reject.</a></td>	
								</xsl:otherwise>
							</xsl:choose>
						</tr>	
					
					</xsl:if>
				
					
					</table>	
				</td>
			</tr>
		</table>
		
	</xsl:template>	
</xsl:stylesheet>