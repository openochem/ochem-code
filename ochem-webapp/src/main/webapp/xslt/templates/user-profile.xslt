<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="2.0">
	<xsl:include href="../helper.xslt" />
	<xsl:include href="inc/data-report.xslt" />
	<xsl:template name="content">
	<script language="javascript" src="js/commons/actionable.js"></script>
	<style>
	h1 {
		color: #333;
		display:block;
		font-family: Helvetica;
		font-size:24px;
		margin-bottom: 10px;
	}
	
	.acknowledgements td {
		padding: 10px;
		width: 700px;
	}
	
	body.popup {
		padding: 40px !important;
		
	}
	
	.user {padding-top: 20px; padding-bottom: 20px;}
	.avatar {margin-right: 20px;}
	.user {border-bottom: 1px solid gray; border-top: 1px solid gray;}
	.user TABLE TD {padding-right: 10px; padding-bottom: 3px;}
	
	TD.contributions {vertical-align: top; padding: 20px; width: 50%; border: 10px solid white; background-color: #FFE;}
	
	.avatar {border: 1px solid #333;}
	
	</style>
	<title>User profile</title>
	
	<xsl:if test="/model/userProfile/firstLogin = 'true'">
			<p>
			Welcome to OCHEM! This is your profile. You can review your information and data here.
			</p>
		</xsl:if>
	
	<div class="user">
	
		<h1>General info</h1>
		<img src="img/no-photo.png" align="left" class="avatar"/>
		<table style="font-size: 110%">
			<tr><td>Username</td><td><b><xsl:value-of select="/model/userProfile/user/@login"/></b></td></tr>
			<tr><td>Real name</td><td><b><xsl:value-of select="/model/userProfile/user/Title"/>&#160;<xsl:value-of select="/model/userProfile/user/FirstName"/>&#160;<xsl:value-of select="/model/userProfile/user/LastName"/></b></td></tr>
			<xsl:if test="/model/session/user/superuser = 'true' and /model/userProfile/user/E-mail">
			<tr><td>E-mail</td><td><xsl:value-of select="/model/userProfile/user/E-mail"/></td></tr>
			</xsl:if>
			<tr><td>Affiliation</td>
				<td>
					<b><xsl:value-of select="/model/userProfile/user/Affiliation"/></b>
					<xsl:if test="not(/model/userProfile/user/Affiliation)">
						not provided
					</xsl:if>
				</td>
			</tr>
			<xsl:if test="(/model/session/user/@login = /model/userProfile/user/@login) or (/model/session/user/superuser = 'true')">
				<tr><td>Bonus points</td><td><a title="Click to see the bonus points overview" tab="Bonus points overview" href="bonuses/overview.do?user={/model/userProfile/user/@login}"><xsl:value-of select="/model/userProfile/bonusPoints"/></a><a class="infolink" href="https://docs.ochem.eu/display/MAN/Bonus+points+system" target="_blank" title="More information on the bonus ponts on OCHEM"></a></td></tr>
			</xsl:if>
		</table>
		<br/>Joined OCHEM on <xsl:value-of select="/model/userProfile/user/registrationTime"/><br/>
		<xsl:if test="/model/userProfile/user/rank = 0">
			<font style="color: #900;">This user has not yet been validated by OCHEM administration team.</font>
		</xsl:if>
		<xsl:if test="/model/userProfile/user/rank = 1">
			<font style="color: #090;">This user has been validated by the OCHEM administartion.</font>
		</xsl:if>
		
		<br/><br/>
		<xsl:choose>
			<xsl:when test="/model/session/user/@login = /model/userProfile/user/@login">
				<a class="fb-button" href="user/show.do">Edit my account details</a>
				<a class="fb-button" href="user/contributions.do" tab="User contributions: {/model/userProfile/user/@login}">View my contributions</a>
			</xsl:when>
			<xsl:otherwise>
				<a class="fb-button" href="user/contributions.do?introducer={/model/userProfile/user/@login}" tab="User contributions: {/model/userProfile/user/@login}">View user's contributions</a>
				<a class="fb-button" href="dialogue/dialogue.do?user={/model/userProfile/user/@login}" tab="Dialogue with {/model/userProfile/user/@login}">Message this user</a>
			</xsl:otherwise>	
		</xsl:choose>
		<xsl:if test="/model/session/user/superuser = 'true'">
			<a class="fb-button" href="useractions/show.do?login={/model/userProfile/user/@login}" tab="User activity: {/model/userProfile/user/@login}">View user's activity</a>
			<a class="fb-button" href="oauth/create.do?user={/model/userProfile/user/@login}" tab="Create App: {/model/userProfile/user/@login}">My Applications</a>
		</xsl:if>
		
		<br style="clear: both"/>
		
	</div>
	
	<xsl:if test="(/model/userProfile/pendingTasksToBeDeleted &gt; 0) or (/model/userProfile/modelsToBeDeleted &gt; 0)">
		<div class="user" style="border-top: none">
			<img src="img/icons/exclamation.gif"/>Some of your data was not used for a long time and is marked for automatic deletion:<br/>
			<br/>
			<table>
			<xsl:if test="(/model/userProfile/pendingTasksToBeDeleted &gt; 0)">
				<tr>
					<td><a tab="Tasks marked for deletion" href="pendingtasks/tasks.do?toBeDeleted=1"><xsl:value-of select="/model/userProfile/pendingTasksToBeDeleted"/> pending tasks</a> to be deleted</td>
					<td>
						<a class="fb-button" action="deleteMarkedTasks">Delete all now</a>
						<a class="fb-button" action="keepMarkedTasks">Keep all</a>
						<a class="fb-button" tab="Tasks marked for deletion" href="pendingtasks/tasks.do?toBeDeleted=1">Review individually</a>
					</td>
				</tr>
			</xsl:if>
			<xsl:if test="(/model/userProfile/modelsToBeDeleted &gt; 0)">
				<tr>
					<td>
						<a tab="Models marked for deletion" href="model/select.do?toBeDeleted=1"><xsl:value-of select="/model/userProfile/modelsToBeDeleted"/> models</a> to be deleted
					</td>
					<td>
						<a class="fb-button" action="deleteMarkedModels">Delete all now</a> 
						<a class="fb-button" action="keepMarkedModels">Keep all</a> 
						<a class="fb-button" tab="Models marked for deletion" href="model/select.do?toBeDeleted=1">Review individually</a> 
					</td>
				</tr>
			</xsl:if>
			</table>
		</div>
	</xsl:if>
	
	<div class="user" style="border-top: none">
		<h1>User contributions</h1>
		<table width="100%">
			<tr>
				<td class="contributions">
					<b>Data contributed by <xsl:value-of select="/model/userProfile/user/@login"/></b><img id="progress" src="img/roller_transparent.gif"/><br/>
					
					<table style="margin-left: 20px;" id="user-records">
					</table>
					
				</td>
				<td class="contributions">
					<b>Models contributed by <xsl:value-of select="/model/userProfile/user/@login"/></b><br/>
					<xsl:choose>
						<xsl:when test="/model/userProfile/user/public-models-count &gt; 0">
							This user has developed <b><xsl:value-of select="/model/userProfile/user/public-models-count"/></b> public models:<br/><br/>
								<table style="margin-left: 20px;">
								<xsl:for-each select="/model/userProfile/user/public-models/model">
									<tr>
										<td>
											Model
										</td>
										<td title="{@name}"> for <nobr><b><xsl:value-of select="modelMappings[1]/property/@name"/></b></nobr>
										</td>
										<td>
											based on <nobr><a href="epbrowser/show.do?basket-select={training-set/@id}" tab="Training set of the model"><xsl:value-of select="training-set/@size"/> records</a></nobr>
										</td>
										<td>
											<a href="model/profile.do?public_id={publicId}" tab="Model profile">[view model]</a>
										</td>
									</tr> 
								</xsl:for-each>
								</table>
						</xsl:when>
						<xsl:otherwise>
							This user did not create any public models yet.
						</xsl:otherwise>
					</xsl:choose>
				</td>
				</tr>
		</table>
		
		
	</div>
	
	<script id="template" type="text/template">
		<tr>
			<td>[%=property.name %]</td>
			<td align="right">
				<a href="epbrowser/show.do?property=[%=property.id %]&amp;introducer={/model/userProfile/user/@login}&amp;approval-status=all" tab="[%=property.name %] records by {/model/userProfile/user/@login}">[%= (0 + 1*metrics.awaitingApproval + 1*metrics.approvedRecords + 1*metrics.privateRecords) %] records</a>
			</td>
		</tr>
	</script>
	
	<script language="javascript">
		include.plugins('view');
		var a = new Actionable();
		a.actionURL = "model/markedForDeletionAction.do";
		
		a.beforeDeletemarkedmodels = a.beforeDeletemarkedtasks = function() {
			return window.confirm("Are you sure? Most probably you are, just asking.");
		}
		
		a.onDeletemarkedmodelsSuccess = a.onKeepmarkedmodelsSuccess = a.onDeletemarkedtasksSuccess = a.onKeepmarkedtasksSuccess = function() {
			window.alert("Your request has been successfully processed");
			window.location.reload();
		}
		
		$(a.initialize);
		
		$(function(){
			(new QSPR.Ajax()).send({
				url: "/editor/getReportData.do?groupings=property&amp;introducer=<xsl:value-of select="/model/userProfile/user/@login"/>",
				success: function(response) {
					var view = new View({element: "template"});
					var a = array(response.awaitingDataReport.dataTree.child);
					
					for (var i = 0; i &lt; a.length; i++)
						$("#user-records").append(view.render(a[i]));
						
						if (a.length == 0)
							$("#user-records").parent().append("The user did not contribute any data yet");
						
					$(document).trigger("DOM_updated", $("#user-records"));
				},
				after: function() {
					$("#progress").remove();
				}
			});
		});
		
	</script>
	</xsl:template>
	
</xsl:stylesheet>
