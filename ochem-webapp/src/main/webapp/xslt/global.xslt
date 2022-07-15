<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"  version="2.0">
	<xsl:output method="xml" omit-xml-declaration="no" indent="yes"/>
	
	<xsl:template match="model">
		<html>
			<head>
				<base href="{@web-root}"/>
				<meta name="Online Chemical Modeling Environment by QSPR Team" http-equiv="content-type" content="text/html; charset=UTF-8" />
				<xsl:call-template name="imports"/>
				<title>Online Chemical Modeling Environment</title>
			</head>
			<body class="yui-skin-sam qspr">
				<iframe id="yui-history-iframe"></iframe>
				<input id="yui-history-field" type="hidden"/>
				<table class="maintable" align="center" border="0" cellspacing="0" cellpadding="0">
					<tr>
						<td valign="center" background="img/bg-top.gif"
							height="70" align="left">
							<img src="img/logo.gif"/>
						</td>
						<td class="top top-right" background="img/bg-top.gif" align="right" valign="top" style="padding-right: 10px; padding-top: 10px;">
							<div id="browser-version"></div>
							<xsl:choose>
								<xsl:when test="session">
									<xsl:choose>
										<xsl:when test="session/user">
											Welcome, Dear <xsl:value-of select="session/user/Title"/> <xsl:value-of select="session/user/LastName"/>! 
											<a href="dialogue/dialogue.do" title="Check messages"><img src="img/icons/mail_icon.png"/>
											
											<b id="message-count"><xsl:if test="session/user/@message &gt; 0">(<xsl:value-of select="session/user/@message"/> new mails!)</xsl:if></b>
											
											</a>
											
											<a href='user/profile.do?login={session/user/@login}'>My account</a>
										</xsl:when>
										<xsl:otherwise>Welcome, Guest!</xsl:otherwise>
									</xsl:choose>
									<a href="login/logout.do?render-mode=redirect">Logout</a>
								</xsl:when>
								<xsl:otherwise>
									<img src="img/icons/login.png"/><a href="login/show.do">log in</a><xsl:if test="//model/allowRegistration = 'true'"><a href="user/newuser.do">create account</a></xsl:if>
								</xsl:otherwise>
							</xsl:choose>
							
							<div style="position: absolute; right: 3px; top: 3px; font-size: 80%; color: #966;" title="The current version of the OCHEM"><xsl:value-of select="version-info"/></div>
						</td>
					</tr>
					<tr>
						<td colspan="2">
						<main-menu id="mm">
							<item href="home/show.do" title="Home">
								<sub-menu id="mm_home">
									<items>
										<item href="home/show.do" title="Home page"/>
									</items>
									<items>
										<item href="https://docs.ochem.eu/display/MAN" title="OCHEM Manual" target="_blank"/>
										<item href="mailto:info@ochem.eu" title="Request support" target="_blank"/>
										<item href="static/tutorials.do" title="Tutorials and Use Cases"/>
									</items>
									<items>
										<item href="static/acknowledgements.do" title="Acknowledgements"/>
										<item href="static/reference.do" title="Citation"/>
									</items>
									<items>
										<xsl:choose>	
										<xsl:when test="session">
											<xsl:choose>
												<xsl:when test="session/user">
													<item href='user/profile.do?login={/model/session/user/@login}' img="img/icons/login.png" title="My account"/>
													<item href="dialogue/dialogue.do" img="img/icons/mail_icon.png" title="Check messages"/>
												</xsl:when>
												<xsl:otherwise></xsl:otherwise>
											</xsl:choose>
											 <item href="login/logout.do" title="Logout"/>
										</xsl:when>
										<xsl:otherwise>
											<item href="login/show.do" img="img/icons/login.png" title="Log in"/>
											<xsl:if test="//model/allowRegistration = 'true'">
												<item href="user/newuser.do" title="Register an account"/>
											</xsl:if>
										</xsl:otherwise>
										</xsl:choose>
									</items>
								</sub-menu>
							</item>
							<item href="#" title="Database">
								<sub-menu id="mm_articles">
									<items>
										<item href="epbrowser/show.do" title="Compound properties"/>
										<item href="basket/show.do" title="Baskets"/>
										<item href="batchupload30/show.do" title="Batch data upload"/>
									</items>
									<items>
										<item href="properties/show.do" title="Properties"/>
										<item href="properties/show.do?condition=true" title="Conditions"/>
										<item href="unit/show.do" title="Units"/>
										<item href="unit/systems.do" title="Systems of units"/>
										<item href="article/show.do" title="Articles/Books"/>
										<item href="journal/show.do" title="Journals"/>
										<item href="alerts/home.do" title="ToxAlerts">
											<sub-menu id="mm_toxalerts">
												<items>
													<item href="alerts/home.do" title="ToxAlerts home"/>
													<item href="alerts/show.do" title="View alerts"/>
													<item href="alerts/screen.do" title="Screen compounds against alerts"/>
													<item href="alerts/upload.do" title="Upload new alerts"/>
												</items>
											</sub-menu>
										</item>
										<item href="matchedpairs/transformations.do" title="MatchedPairs"  img="img/icons/mmp-16.png">
											<sub-menu id="mm_matchedpairs">
												<xsl:if test="/model/session/user/ochem-labs = 'true'">
												<items>
													<item href="matchedpairs/transformations.do" title="All molecular transformations"/>
												</items>
												<items>
													<item href="matchedpairs/status.do" title="Indexing status (technical info)"/>
												</items>
												<items>
													<item href="matchedpairs/clearCache.do" title="Clear caches"/>
												</items>
												</xsl:if>
											</sub-menu>
										</item>	
									</items>
									<items>
										<item href="epbrowser/show.do?trash=1" title="Trash" img="img/icons/trash_small.gif"/>
									</items>
								</sub-menu>
							</item>
							<item href="#" title="Models">
								<sub-menu id="mm_models">
									<items>
										<item href="modelconfigurator/choose.do" title="Create a model" type="bold"/>
										<item href="multiplemodels/create.do" title="Create multiple models"/>
										<item img="img/icons/mol-descriptors-16.png" href="descriptorscalculator/show.do" title="Calculate descriptors"/>
									</items>
									<items>
										<item href="pendingtasks/tasks.do" title="View pending tasks"/>
									</items>
									<items>
										<item href="model/select.do" title="Apply a model"/>
										<item href="predictor/show.do" title="Open predictor"/>
									</items>
									<items>
										<item href="setcomparison/select.do" title="SetCompare utility"/>
										<item href="descriptorsstorage/show.do" title="Descriptors storage"/>
									</items>
								</sub-menu>
							</item>
							
								<item href="#" title="Administration">
									<sub-menu id="administration">
										<items>	
											<item href="systemstatus/monitors.do" title="OCHEM load"/>
											<item href="systemstatus/meta.do" title="Metaserver status"/>
											<item href="systemstatus/show.do" title="Server tests"/>
										</items>
										<items>
											<item href="molbrowser/show.do" title="Molecules"/>
											<item href="epbrowser/show.do?trash=1&amp;globaltrash=1" title="Global trash" img="img/icons/trash16.png"/>
										</items>
									</sub-menu>
								</item>

						</main-menu>
						</td>
					</tr>
						<tr>
							<xsl:if test="not(announcement/value)">
								<xsl:attribute name="class">invisible</xsl:attribute>
							</xsl:if>
							<td colspan="2" class="announcement"><xsl:value-of select="announcement/value"/></td>
						</tr>
					<tr>
						<td height="100%" class="content" valign="top" colspan="2">
							<div style="position: relative;" class="aaaa">
								<div id="main-tab-container" class="yui-navset main-tabs">
									<ul class="yui-nav invisible">
										<li class="selected">
											<a href="#maintab1">
												<em>undefined</em>
											</a>
										</li>
									</ul>
									<div class="yui-content">
										<div id="tab1">
											<div class='frame-status invisible'><span></span><br/><img src='img/roller_transparent.gif'/></div>
											<iframe frameborder="0">
												<xsl:attribute name="src">
													<xsl:choose>
														<xsl:when test="contains(innerUrl, '?')"><xsl:value-of select="innerUrl"/>&amp;render-mode=popup</xsl:when>
														<xsl:otherwise><xsl:value-of select="innerUrl"/>?render-mode=popup</xsl:otherwise>
													</xsl:choose>	
												</xsl:attribute>
												<!-- <xsl:apply-templates select="content" /> -->
											</iframe>
										</div>
									</div>
								</div>
							</div>
						</td>
					</tr>
				</table>
				<div id="bottom-aligner">ABC</div>
				<xsl:if test="/model/@web-root = 'https://ochem.eu/'">
			<!--##IVT DISABLED
					<script language="javascript">
						
						 var _gaq = _gaq || [];
			 			_gaq.push(['_setAccount', 'UA-17599118-1']);
						 _gaq.push(['_setDomainName', '.ochem.eu']);
			 				_gaq.push(['_trackPageview']);
						
						 (function() {
						   var ga = document.createElement('script'); ga.type = 'text/javascript'; ga.async = true;
						   ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js';
						   var s = document.getElementsByTagName('script')[0]; s.parentNode.insertBefore(ga, s);
						 })();
				
					</script>
			IVT DISABLED##-->
				</xsl:if>
			</body>
		</html>
	</xsl:template>
	
	<xsl:template match="model[@render-mode='popup']">
		<html>
			<head>
				<base href="{@web-root}"/>
				<meta http-equiv="content-type" content="text/html; charset=UTF-8" />
				<xsl:call-template name="imports"/>
			</head>

			<body class="yui-skin-sam popup">
				<iframe id="yui-history-iframe"></iframe>
				<input id="yui-history-field" type="hidden"/>
				<xsl:apply-templates select="content" />
				<xsl:if test="/model/@web-root = 'https://ochem.eu/'">
			<!--##IVT DISABLED
					<script language="javascript">
						
						 var _gaq = _gaq || [];
			 			_gaq.push(['_setAccount', 'UA-17599118-1']);
						 _gaq.push(['_setDomainName', '.ochem.eu']);
			 				_gaq.push(['_trackPageview']);
						
						 (function() {
						   var ga = document.createElement('script'); ga.type = 'text/javascript'; ga.async = true;
						   ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js';
						   var s = document.getElementsByTagName('script')[0]; s.parentNode.insertBefore(ga, s);
						 })();
					</script>
			IVT DISABLED##-->
				</xsl:if>
				<div class="invisible" id="developer-info">
					Template name: <xsl:value-of select="/model/templateName"/>
				</div>
			</body>
		</html>
	</xsl:template>
	
	<xsl:template match="model[@render-mode='iprior']">
		<html>
			<head>
				<base href="{@web-root}"/>
				<meta http-equiv="content-type" content="text/html; charset=UTF-8" />
				<xsl:call-template name="imports"/>
				<style type="text/css">
					BODY.iprior {padding: 30px;}
					.logo {border-bottom: 2px solid #666; padding-bottom: 10px;}
					#content {padding: 10px;}
				</style>
			</head>

			<body class="iprior">
				<iframe id="yui-history-iframe"></iframe>
				<input id="yui-history-field" type="hidden"/>
				<div class="logo">
					<img src="img/eADMET-Logo-360dpi_Noborder.jpg" alt="eADMET" style="hight:60px;width:220px"/>
				</div>
				<div id="content">
					<xsl:apply-templates select="content" />
				</div>
				<!--##qspr
					<script language="javascript">
						
						 var _gaq = _gaq || [];
			 			_gaq.push(['_setAccount', 'UA-17599118-1']);
						 _gaq.push(['_setDomainName', '.ochem.eu']);
			 				_gaq.push(['_trackPageview']);
						
						 (function() {
						   var ga = document.createElement('script'); ga.type = 'text/javascript'; ga.async = true;
						   ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js';
						   var s = document.getElementsByTagName('script')[0]; s.parentNode.insertBefore(ga, s);
						 })();
				
					</script>
				qspr##-->
			</body>
		</html>
	</xsl:template>
	
	<xsl:template match="@*|node()" priority="-1">
		<xsl:copy>
			<xsl:apply-templates select="@*|node()" />
		</xsl:copy>
	</xsl:template>
	
	<xsl:template match="content">
		<xsl:apply-templates/>
	</xsl:template>
	
	<xsl:template name="imports">
		
		<link rel="stylesheet" type="text/css" href="css/main.css?ver=2.3.5" />
		<link rel="stylesheet" type="text/css" href="css/forms.css" />
		<link rel="stylesheet" type="text/css" href="css/tablecloth.css" />
		<link rel="stylesheet" type="text/css" href="css/itunes.css" />
		<link rel="stylesheet" type="text/css" href="css/custom.css" />
		<link rel="stylesheet" type="text/css" href="css/buttons.css" />
		<link rel="stylesheet" type="text/css" href="css/smoothness/jquery-ui-1.10.3.custom.css" />
		
		<link rel="stylesheet" type="text/css" href="js/lib/yahoo/reset/reset.css" />
		<link rel="stylesheet" type="text/css" href="js/lib/yahoo/autocomplete/assets/skins/sam/autocomplete-skin.css" />
		<link rel="stylesheet" type="text/css" href="js/lib/yahoo/fonts/fonts-min.css" />
		<link rel="stylesheet" type="text/css" href="js/lib/yahoo/container/assets/skins/sam/container.css" />
		<link rel="stylesheet" type="text/css" href="js/lib/yahoo/tabview/assets/skins/sam/tabview.css" />
		<link rel="stylesheet" type="text/css" href="js/lib/yahoo/menu/assets/skins/sam/menu.css" />
		<link rel="stylesheet" type="text/css" href="js/lib/yahoo/datatable/assets/skins/sam/datatable.css" />
		
		
		

		<!-- <link rel="stylesheet" type="text/css" href="css/yahooskin/menu-skin.css" /> -->
		
		<script type="text/javascript" src="js/lib/yahoo/yahoo.ochem.custom.js"/>
		<script type="text/javascript" src="js/lib/jquery-1.8.2.min.js" />
		<script type="text/javascript" src="js/lib/jquery-ui-1.10.3.custom.min.js" />
		
		<!--   <script type="text/javascript" src="js/lib/jquery.dimensions.pack.js" /> -->
		<!-- <script type="text/javascript" src="js/lib/jquery.tooltip.pack.js" />  -->
  		<script type="text/javascript" src="js/lib/jmvc/include.js" />
		<script type="text/javascript" src="js/lib/xmljson.js" />
		
		<script type="text/javascript" src="js/commons/midnighter.js?ver=1.8.10" />
		<script type="text/javascript" src="js/commons/datasources.js" />
		<script type="text/javascript" src="js/commons/autocomplete.js" />
		<script type="text/javascript" src="js/commons/initialize.js" />
		<script type="text/javascript" src="js/commons/pager.js?ver=2.3.5" />
		<script type="text/javascript" src="js/commons/long-operation.js?ver=1.5.4" />
		<script language="javascript">
			var webRoot = "<xsl:value-of select="@web-root"/>";
			var currentUser = "<xsl:value-of select="session/user/@login"/>";
			var currentUserRank = 1.0 * "<xsl:value-of select="session/user/rank"/>";
			var currentUserLimit = "<xsl:value-of select="session/limit"/>";
			var maximumTaskPriority = "<xsl:value-of select="session/@max-priority"/>";
			var isSuperUser = (currentUserRank >= 10);
			var isModerator = ("true" == "<xsl:value-of select="session/@moderator"/>");
			var isMobile = false;
			var unreadMessagesCount = '<xsl:value-of select="session/user/@message"/>' || 0;
			<xsl:if test="session/isMobile = 'true'">
				isMobile = true;
			</xsl:if>
		</script>
	
	</xsl:template>
</xsl:stylesheet>