<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />
	<xsl:template name="content">
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/browser.js?ver=1.8.6"></script>
		<script type="text/javascript" src="js/blocks/datacube.js?ver=1.8.6.1" />
		<script type="text/javascript" src="js/commons/floating-headers.js"></script>
		<style type="text/css">
			TABLE.striped TD, TABLE.striped TH {padding: 5px 15px; border-bottom: 1px solid white; font-family: Verdana; font-size: 10pt;}
			TABLE.striped TH {font-weight; bold;}
			TABLE.striped TR:nth-child(even) {background-color: rgba(0, 0, 0, 0.05);}
			TABLE.striped TH {font-style: bold; border-bottom: 1px solid black;}
			TR.ready TD {background-color: #c0e793;}
			.rounded-button {background-color: #5B5 !important; color: white !important;}
			.rounded-button:hover {background-color: #595 !important;}
			
			TD.highlighted-column {
				background-color: rgba(0, 0, 0, 0.03);
			}
			
			.depth-1, .depth-0 {font-weight: bold;}
			.depth-0 TD {font-size: 105% !important;}
			.depth-1 TD:first-child {padding-left: 30px;}
			.depth-2 TD:first-child {padding-left: 60px;}
			.depth-3 TD:first-child {padding-left: 90px;}
			.depth-4 TD:first-child {padding-left: 120px; font-size: 75%;}
			
			.groupings {
				border-radius: 5px;
				background-color: #FFD;
				padding: 10px;
				border: 1px solid gray;
				width: 500px;
			}
			
			.groupings DIV.content {margin-bottom: 5px; font-size: 80%;}
			A.img {text-decoration: none !important; }
			A.img IMG {opacity: 0.5;}
			A.img:hover IMG {opacity: 1.0;}
			
			.report A {color: #11A; width:1200px;}
			
			.report TR:hover {background-color: #CEC !important;}
			
			.summary {width: 800px; font-size: 120%; margin: 10px;}
			.summary TD {padding: 10px;}
			.approve, .reject {font-size: 90%; }
			.approve IMG, .reject IMG {margin-right: 10px; align: left;}
			
		</style>
		<center>
		<table class="summary">
		<xsl:choose>
			<xsl:when test="count(//entry) &gt; 0">
					<xsl:choose>
					<xsl:when test="//param[@key='approve'] = 'approve'">
						<tr><td class="approve">
						<img src="img/icons/approve32.png"/>You are going to approve the records and assign bonus points to record introducers</td></tr>
					</xsl:when>
					<xsl:when test="//param[@key='approve'] = 'reject'">
						<tr><td class="reject"><img src="img/icons/disapprove32.png"/>You are going to reject the records and suggest the record introducers to review their data.<br/>
						<br/>You are welcome to explain the reasons for rejection to the data introducer by droping him a message.
						The user can resubmit his data for approval at any time.
						</td></tr>
					</xsl:when>
					<xsl:otherwise>
					<tr><td class="approve">
						You are going to unapprove (resend for moderation) the following records</td></tr>
					</xsl:otherwise>
					</xsl:choose>
				<tr><td>
				<table class="striped report">
					<thead>
						<tr>
							<th>Property<br/>weight</th>
							<th>Record count</th>
						</tr>
					</thead>
					<tbody>
					</tbody>
				</table>
				<div class="progress invisible"><img src="/img/long_green.gif"/></div>
				</td></tr>
				<xsl:choose>
					<xsl:when test="//param[@key='approve'] = 'approve'">
						<tr><td class="approve">
							Awarded bonuses: 
								<xsl:for-each select="//userBonuses/entry">
									<xsl:value-of select="value"/> to user <xsl:value-of select="key/@login"/>
									<xsl:if test="position() != count(//userBonuses/entry)">,</xsl:if>
								</xsl:for-each>
							<br/>
							For your moderation efforts, OCHEM grants you <b><xsl:value-of select="//basketBonusReport/moderatorAcceptBonuses"/> bonus points</b>.
						</td></tr>
						<tr><td class="act">
							<a action="approve" class="fancy-button" style="width: 400px; text-align: center;">Approve the records and award bonus points</a>
						</td></tr>					
					</xsl:when>
					<xsl:when test="//param[@key='approve'] = 'reject'">
						<tr><td class="reject">
							For your moderation efforts, OCHEM grants you <b><xsl:value-of select="//basketBonusReport/moderatorRejectBonuses"/> bonus points</b>.
						</td></tr>					
						<tr><td class="act">
							<a action="reject" class="fancy-button-red" style="width: 300px; text-align: center;">Reject the records</a>
						</td></tr>
					</xsl:when>
					<xsl:otherwise>
						<tr><td class="act">
							<a action="unapprove" class="fancy-button-grey" style="width: 300px; text-align: center;">Unapprove the records</a>
						</td></tr>
					</xsl:otherwise>
				</xsl:choose>
			</xsl:when>
			<xsl:otherwise>
				<tr><td>
					<xsl:choose>
					<xsl:when test="//param[@key='approve'] = 'approve'">
						You have selected no unapproved records
					</xsl:when>
					<xsl:when test="//param[@key='approve'] = 'reject'">
						You have selected no unapproved records
					</xsl:when>
					<xsl:otherwise>
					You have selected no approved records
					</xsl:otherwise>
					</xsl:choose>
				</td></tr>
				<tr><td class="act">
					<a action="close" class="fancy-button" style="width: 400px; text-align: center;">Close this tab</a>
				</td></tr>		
			</xsl:otherwise>
		</xsl:choose>
		</table>
		</center>
		<script language="javascript">
			<xsl:choose>
					<xsl:when test="//param[@key='approve'] = 'approve'">
						var report = new DataCubeReport("editor/getApprovePreview.do?out=json&amp;approved=false");
					</xsl:when>
					<xsl:when test="//param[@key='approve'] = 'reject'">
						var report = new DataCubeReport("editor/getApprovePreview.do?out=json&amp;approved=false");
					</xsl:when>
					<xsl:otherwise>
						var report = new DataCubeReport("editor/getApprovePreview.do?out=json&amp;approved=true");
					</xsl:otherwise>
			</xsl:choose>
			
			report.drawNode = function(o) {
				if (o.user)
					return getTemplate("user-template").render(o.user);
				else if (o.article)
					return o.article.title;
				else if (o.property)
					return getTemplate("property-template").render(o.property);
				else 
					return "No";
			}
			
			report.dataDrawn.register(function() {	
				console.log(report.responseOthers.dataReportRequest);
			});
			
			report.drawMetrics = function(tr, metrics, queryString, node) {
				console.log(queryString);
				tr.append(getTemplate("metrics-template").render({metrics: metrics, queryString: queryString, property: node.property, basketId:"<xsl:value-of select="//basketBonusReport/basket/@id"/>"}));		
			}
			
			var a = new Actionable();
			a.scope = $(".act");
			a.doApprove = function(link)
			{
				a.ajax.send({
					url: "editor/approve.do",
					data: "basket_id=<xsl:value-of select="//basket/@id"/>",
					success: function()
					{
						// Sometimes user get confised and think that nothing happened. Show a clarification message.
						window.alert("The records you selected have been marked as approved. Thank you!");
						
						window.callback();
						window.closeActiveTab(); 
					}
				});
			}
			
			a.doReject = function(link)
			{
				a.ajax.send({
					url: "editor/reject.do",
					data: "basket_id=<xsl:value-of select="//basket/@id"/>",
					success: function()
					{
						window.callback();
						window.closeActiveTab(); 
					}
				});
			}
			
			a.doUnapprove = function(link)
			{
				a.ajax.send({
					url: "editor/unapprove.do",
					data: "basket_id=<xsl:value-of select="//basket/@id"/>",
					success: function()
					{
						window.callback();
						window.closeActiveTab(); 
					}
				});
			}
			
			a.doClose = function(link)
			{
				window.callback();
				window.closeActiveTab(); 
			}
		</script>
		
		<script type="text/template" id="metrics-template">
			[%
				function printLink(metrics, metricsName, queryString, basketId) {
					if (metrics[metricsName] > 0) {
						%]<a tab="Details" href="/editor/browseData.do?basket-select=[%=basketId %]&amp;metrics=[%=metricsName %]&amp;[%=queryString %]">[%=metrics[metricsName] %]</a>[%
					}
				}
			%]
			<td class="right">[% if (property != undefined) { %]
					[%=property.bonusPointsWeight %]
				[% } %]
			</td>
			<xsl:choose>
					<xsl:when test="//param[@key='approve'] = 'approve'">
						<td class="right  highlighted-column">[% printLink(metrics, "awaitingApproval", queryString, basketId) %]</td>
					</xsl:when>
					<xsl:when test="//param[@key='approve'] = 'reject'">
						<td class="right  highlighted-column">[% printLink(metrics, "awaitingApproval", queryString, basketId) %]</td>
					</xsl:when>
					<xsl:otherwise>
						<td class="right  highlighted-column">[% printLink(metrics, "approvedRecords", queryString, basketId) %]</td>
						<td class="right  highlighted-column">[% printLink(metrics, "rejected", queryString, basketId) %]</td>
					</xsl:otherwise>
			</xsl:choose>
			
		</script>
		
		<script type="text/template" id="user-template">
			<a class="img" tab="Messages with [%=login%]" title="Message [%=login%]" href="dialogue/dialogue.do?user=[%=login%]"><img src="img/icons/mail_icon.png"/></a><a tab="User profile" href="user/profile.do?login=[%=login%]">[%=login%]</a>
		</script>
		
		<script type="text/template" id="property-template">
			<a href="properties/edit.do?id=[%=id %]" tab="Property profile: [%=name%]">[%=name %]</a>
			[% if (data.moderator &amp;&amp; !getParams["moderator"]) { %]<a class="img" tab="Dialogue with [%=moderator.login %]" href="dialogue/dialogue.do?user=[%=moderator.login %]" title="Message the moderator of this property ([%=moderator.login %])"><img src="img/icons/mail_icon.png"/></a>[% } %]
			<xsl:if test="//param[@key='global-report']">
			[% if (approved != 'true') { %]
			<small><a action="approve" property="[%=id %]"  title="This property has NOT yet approved by the OCHEM editor. Click to mark is as approved."> [approve]</a></small>
			[% } %]
			</xsl:if>
		</script>
	</xsl:template>
</xsl:stylesheet>