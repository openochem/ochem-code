<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />

	<xsl:template name="content">
		<title>Registered users</title>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript">
			include.plugins('view');
			
			var rankTitles = new Array("unverified", "verified", "superuser");
			var ranks = new Array(0, 1, 10);
			
			var sampleBrowser = new Browser();
			sampleBrowser.itemElement = "user";
			sampleBrowser.itemTemplate = "js/templates/user.ejs";
			sampleBrowser.url = webRoot + "user/list.do";
			sampleBrowser.actionURL = webRoot + "user/action.do";
			sampleBrowser.onChangerankSuccess = function(){
				sampleBrowser.request(false);
			};
			sampleBrowser.onDeleteSuccess = sampleBrowser.onSuspendSuccess = sampleBrowser.onUnsuspendSuccess = function()
			{
				sampleBrowser.request();
			};
			sampleBrowser.listenEvent("items_loaded", function(){
				$(document).find(".rank").change(function(){
					sampleBrowser.setPosition($(this));
					sampleBrowser.callAction("changerank");
				});
			});
			$(document).ready(function() {sampleBrowser.initialize();});
		</script>
		<style type="text/css">
			.rank-0 {color: #666;}
			.rank-10 {background-color: #FEE;}
		</style>
		<table height="100%" width="100%">
			<tr>
				<td class="itunes-up" colspan="2">
					<h1>Registered users</h1>
					Browse and manage registered users.
				</td>
			</tr>
			<tr>
				<td class="itunes-right">
					User Name <input type="text" name="username" filter="1"/>
					Organisation type:
					<select name="organisation" filter="1">
						<option value="">All organisation types</option>
						<option value="Commercial">Commercial</option>
						<option value="Academic">Academic</option>
						<option value="Non profit - Governmental">Non profit - Governmental</option>
						<option value="Self-employed">Self-employed</option>Academic
						<option value="empty">Not specified</option>
					</select>
					Sort by:
					<select name="sort" filter="1">
						<option value="registrationTime">Joined date</option>
						<option value="activitiesCountTotal">Activity</option>
						<option value="activitiesCountMonth">Activity (last month)</option>
						<option value="latestActivityTime">Last activity time</option>
					</select>
					<a href="javascript:sampleBrowser.request(true)">
						[search]
					</a>
					<div class="pager-strip">
						<b class="showed">none</b> of <b class="total">none</b>
					</div>
					<div id="pager">
					</div>
					
					<div id="Browser">
					</div>
					
					<div class="pager-strip">
						<b class="showed">none</b> of <b class="total">none</b>
					</div>
				</td>
			</tr>
		</table>
	</xsl:template>
</xsl:stylesheet>