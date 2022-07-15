<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:template name="data-report">
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
			.depth-4 TD:first-child {padding-left: 120px;}
			.depth-4 {font-size: 80%;}
			
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
			
			.report A {color: #11A;}
			
			.report TR:hover {background-color: #CEC !important;}
			
		</style>
		<div class="groupings">
			<div class="content">Select the desired data groupings. The groupings can be re-ordered.</div>
			<div class="grouping"><input type="checkbox" name="isPrimary"/> Primary record status</div>
			
			<div class="grouping">
				<input type="checkbox" name="property" checked="checked">
					<xsl:if test="//param[@key='moderator-report']">
						<xsl:attribute name="disabled">1</xsl:attribute>
					</xsl:if>
				</input> Property</div>
			<div class="grouping cond"><input type="checkbox" name="cond-"/> <span>Experimental condition <a class="infolink" title="You can group data by a qualitative condition value (check the option to select the condition)"></a></span></div>
			<div class="grouping"><input type="checkbox" name="introducer" checked="checked"/> Data introducer</div>
			<xsl:if test="not(//param[@key='global-report'])">
				<div class="grouping"><input type="checkbox" name="article"/> Publication</div>
			</xsl:if>
		</div>
		<br/>
		<table class="striped report">
			<thead>
			<tr>
				<xsl:if test="//param[@key='global-report']">
					<th>Moderator</th>
					<th>Property<br/>weight</th>
				</xsl:if>
				<th class="private-records invisible">Private records</th>
				<th>Awaiting approval</th>
				<th>Rejected records</th>
				<th class="approved-records invisible">Approved records</th>
				<th class="private-records invisible">Private models</th>
				<th>Models awaiting approval</th>
				<th>Approved models</th>
			</tr>
			</thead>
			<tbody>
			</tbody>
			
		</table>
		<div class="progress invisible"><img src="/img/long_green.gif"/></div>
		
		<script language="javascript">
			if (getParams["my"])
				$("input[name='introducer']").closest("div").remove();	
			
			var report = new DataCubeReport("editor/getReportData.do?out=json");
			report.drawNode = function(o) {
				if (o.user)
					return getTemplate("user-template").render(o.user);
				else if (o.article)
					return o.article.title;
				else if (o.property)
					return getTemplate("property-template").render(o.property);
				else 
					return getTemplate("generic-template").render(o);
			}
			
			report.dataDrawn.register(function() {	
				console.log(report.responseOthers.dataReportRequest);
				if (report.responseOthers.dataReportRequest.countPrivateRecords == "true")
					$(".private-records").removeClass("invisible");
				if (report.responseOthers.dataReportRequest.countApprovedRecords == "true")
					$(".approved-records").removeClass("invisible");
			});
			
			report.drawMetrics = function(tr, metrics, queryString, node) {
				console.log(queryString);
				tr.append(getTemplate("metrics-template").render({metrics: metrics, queryString: queryString, property: node.property}));		
			}
			
			var a = new Actionable();
			a.scope = $(".report");
			a.doAssign = function(link)
			{
				var win = openTab("Select a user", "user/all.do");
				win.callback = function(entity)
				{
					a.ajax.send({
						url: "editor/assignModerator.do",
						data: "property=" + $(link).attr("property") + "&amp;login=" + entity.login,
						success: function()
						{
							var a = $('<a></a>').attr("href", "user/profile.do?login=" + entity.login);
							a.html(entity.login);
							$(link).parent().html('').append(a);
							win.closeTab();	
						}
					});
					
				}
			}
			
			a.doApprove = function(link)
			{
				a.ajax.send({
						url: "properties/action.do",
						data: "action=approve&amp;id=" + $(link).attr("property"),
						success: function()
						{
							$(link).parent().html('');
						}
					});
			}
			
			$(".cond INPUT").change(function(){
				if ($(this).is(":checked"))
				{
					var win = openTab("Select a qualitative condition for grouping the data", "properties/show.do?condition=1");
					win.callback = function(entity)
					{
						$(".cond INPUT").attr("name", "cond-" + entity.id);
						$(".cond SPAN").html(entity.name);
						report.load();
						win.closeTab();
					}	
				}
			});
			
			$(document).bind('DOM_updated', null, function()
			{
				a.attachActions($(".report"));
				ActivateFloatingHeaders("TABLE.report");
			});
			
		</script>
		
		<script type="text/template" id="metrics-template">
			[%
				function printLink(metrics, metricsName, queryString) {
					if (metrics[metricsName] > 0) {
						%]<a tab="Details" href="/editor/browseData.do?metrics=[%=metricsName %]&amp;[%=queryString %]">[%=metrics[metricsName] %]</a>[%
					}
				}
			%]
			<xsl:if test="//param[@key='global-report']">
				<td>[% if (property != undefined) { %]
					[% if (property.moderator) {%]
						<a tab="Moderator profile" href="user/profile.do?login=[%=property.moderator.login %]" property="[%=property.id %]">[%=property.moderator.login %]</a><a action="assign" title="Reassign a moderator">[~]</a>
					[% } else { %]
						unmoderated <a href="" property="[%=property.id %]" title="Assign a moderator" action="assign">[+]</a>
					[% }} %]
				</td>
				<td class="right">[% if (property != undefined) { %]
						[%=property.bonusPointsWeight %]
					[% } %]
				</td>
			</xsl:if>
			<td class="right private-records invisible  highlighted-column">[% printLink(metrics, "privateRecords", queryString ) %]</td>
			<td class="right  highlighted-column">[% printLink(metrics, "awaitingApproval", queryString) %]</td>
			<td class="right  highlighted-column">[% printLink(metrics, "rejected", queryString) %]</td>
			<xsl:if test="not(//param[@key='global-report'])">
				<td class="right  highlighted-column">[% printLink(metrics, "approvedRecords", queryString) %]</td>
			</xsl:if>
			
			<td class="right invisible private-records">[% printLink(metrics, "privateModels", queryString) %]</td>
			<td class="right">[% printLink(metrics, "awaitingModels", queryString) %]</td>
			<td class="right">[% printLink(metrics, "approvedModels", queryString) %]</td>
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
		
		<script type="text/template" id="generic-template">
			[% if (title == "false") { %]
				Non-primary records 
			[% } else if (title == "true") { %]
				Primary records [% } else { %]
				[%=title %][%
			} %]
		</script>
	</xsl:template>
</xsl:stylesheet>