<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />
	
	<xsl:template name="content">
		
		<link rel="stylesheet" type="text/css" href="css/smoothness/jquery-ui-1.9.1.custom.css" />
		<link rel="stylesheet" type="text/css" href="css/tagit-simple-blue.css" />
		<script language="javascript" src="js/lib/jquery-ui-1.9.1.custom.min.js"></script>
		<script language="javascript" src="js/lib/tagit.js"></script>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script type="text/javascript" src="js/lib/jstree/jquery.jstree.js"></script>
		<style type="text/css">
			#tree {
				height: 500px; overflow: auto;
			}
			
			.jstree-default li.folder-false > A ins {
			    background-image: url("/img/icons/text_file.jpg") !important;
			    background-position: 0px 0px !important;
			}
		</style>
		<table height="100%" width="100%">
			<tr><td class="itunes-up">
				<img src="img/icons/sa.png"/>
				<h1><span class="toxalerts">ToxAlerts</span>: Table of contents</h1>
				
				</td></tr>
			<tr>
				<td class="itunes-right">
					<a class="rounded-button" action="addchild">Add a folder</a>
					<a class="rounded-button" action="delete">Delete</a>
					<a class="rounded-button" action="rename">Rename</a><br/><br/>
					<div id="tree">
					<ul id="root">
						<xsl:apply-templates select="/model/tree-node"/>
					</ul>
					</div>
				</td>
			</tr>
		</table>
		
		<script language="javascript">
			var tree;
			$(function(){
				tree = $("#tree");
				tree.jstree({
					core: {
						load_open: true
					},
					dnd: {
						"drop_finish": function(data)
						{
							window.alert("" + data.o + " " + data.r);
						}
					},
					crrm : {
				        move : {
				            "check_move" : function (m) {
				            	return true;
				           }
				        }
				    },
					"plugins" : [ "themes", "html_data", "dnd", "crrm", "ui" ]
				});
				
				tree.bind("loaded.jstree", function (event, data) {
        			tree.jstree("open_all");
    			});
    			
    			tree.bind("move_node.jstree", function(event, data){
    				console.log(data.rslt);
    				var params = new Array();
    				params.push("action=move");
    				params.push("id=" + (data.rslt.o.attr("alert-id")));
    				params.push("parent=" + (data.rslt.np.attr("alert-id")));
    				params.push("type=" + data.rslt.p);
    				params.push("replaced-node=" + data.rslt.or.attr("alert-id"));
    				actionable.ajax.send({
    					url: actionable.actionURL,
    					data: params.join("&amp;"),
    					success: function(){}
    				});
				});
				
				tree.bind("dblclick.jstree", function (event) {
				   var node = $(event.target).closest("li");
				   openTab("Alert details", node.children("a").attr("href"));
				});
    			
    			$("li[folder='false']").addClass("alert");
    				actionable.ajax.waitingDialog = createWaitingDialog();
				});
			
			var actionable = new Actionable();
			actionable.actionURL = "alerts/contentsAction.do"
			actionable.beforeAddchild = function()
			{
				var userInput = window.prompt('Please, enter the name of the folder', '');
				if (userInput) {
					actionable.name = userInput;
					return true;
				}
			}
			
			actionable.beforeRename = function()
			{
				var userInput = window.prompt('Please, enter the new name of the folder', '');
				if (userInput) {
					actionable.name = userInput;
					return true;
				}
			}
			
			actionable.getActionQuery = function()
			{
				return "id=" + $('#tree').jstree('get_selected').attr('alert-id') + "&amp;name=" + URLEncode(actionable.name);
			}
			
			actionable.onDeleteSuccess = function()
			{
				$('#tree').jstree('get_selected').remove();
			}
			
			actionable.onRenameSuccess = function()
			{
				$('#tree').jstree('rename_node', $('#tree').jstree('get_selected') , actionable.name);
			}
			
			actionable.onAddchildSuccess = function()
			{
				var newNode = { data: actionable.name };
				tree.jstree("create_node", $('#tree').jstree('get_selected'), 0, newNode, false, false);
			}
			
			
		</script>
		
		<div id="waitingDialog"> 
		    <div class="hd">Please wait</div> 
		    <div class="bd" style="text-align: center;"> 
		        Please wait until action is completed.<br/>
		        It may take a while.<br/>
		        <img src="img/roller_small.gif"/> 
		    </div> 
		</div>
		
		</xsl:template>
		
		<xsl:template match="tree-node">
			<li alert-id="{substructure-alert/@id}" class="folder-{substructure-alert/@folder}">
				<a href="alerts/show.do?id={substructure-alert/@id}"><xsl:value-of select="substructure-alert/name"/></a>
				<xsl:if test="children/tree-node">
					<ul><xsl:apply-templates select="children/tree-node"/></ul> 
				</xsl:if>
			</li>
		</xsl:template>
</xsl:stylesheet>
