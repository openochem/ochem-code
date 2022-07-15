function PropertyRemappingActionable()
{
	
	Actionable.call(this);
	EventDispatcher.call(this);
	
	var self = this;
	
//	this.doKnown = function(link)
//	{
//		this.knownDialog.show();
//	}
	
	this.setContainerValue = function(link, value)
	{
		$(link).parents("div.container:first").find("input[type='hidden']:first").attr('value', value);
		$(link).parents("div.container:first").find(".error:first").html("");
		$(link).parents("div.container:first").find(".warn:first").html("");
		$(link).parents("div.container:first").find(".notice:first").html("");
		$(link).html(value);		
	}
	
	this.doArticle = function(link)
	{	
		var artWin = openTab("Select article", webRoot + "article/show.do?render-mode=popup");
		artWin.callback = function(article)
		{
			self.setContainerValue(link, "A"+article["id"]);
//			$(link).parents("div.container:first").find("input[type='hidden']:first").attr('value', "A"+article["id"]);
//			$(link).parents("div.container:first").find(".error:first").html("");
//			$(link).html("A"+article["id"]);
			//Remap/refresh units and options
			artWin.closeTab();
		}
	}
	
	this.doBasket = function(link)
	{	
		var baskWin = openTab("Select molecule set", webRoot + "baskets/show.do?render-mode=popup");
		baskWin.callback = function(basket)
		{
			self.setContainerValue(link, basket["name"]);
//			$(link).parents("div.container:first").find("input[type='hidden']:first").attr('value', "A"+article["id"]);
//			$(link).parents("div.container:first").find(".error:first").html("");
//			$(link).html(basket["name"]);
			//Remap/refresh units and options
			baskWin.closeTab();
		}
	}
	
	this.doProperty = function(link)
	{	
		var propWin = openTab("Select property", webRoot + "properties/show.do?render-mode=popup");
		propWin.callback = function(property)
		{
			self.setContainerValue(link, property["name"]);
//			$(link).parents("div.container:first").find("input[type='hidden']:first").attr('value', property["name"]);
//			$(link).parents("div.container:first").find(".error:first").html("");
//			$(link).html(property["name"]);
			//Remap/refresh units and options
			propWin.closeTab();
		}
	}

	this.doCondition = function(link)
	{	
		var condWin = openTab("Select condition", webRoot + "properties/show.do?render-mode=popup&condition=true");
		condWin.callback = function(condition)
		{
			self.setContainerValue(link, condition["name"]);
//			$(link).parents("div.container:first").find("input[type='hidden']:first").attr('value', condition["name"]);
//			$(link).parents("div.container:first").find(".error:first").html("");
//			$(link).html(condition["name"]);
			//Remap/refresh units and options
			condWin.closeTab();
		}
	}
	
	this.doUnit = function(link)
	{
		var pname = $($(link).parents("div.container").get(1)).find("input[type='hidden']:first").attr('value');
		var unitWin = openTab("Select unit ("+pname+")", webRoot + "unit/show.do?render-mode=popup&pname="+URLEncode(pname));

		unitWin.callback = function(unit)
		{
			self.setContainerValue(link, unit["name"]);
//			$(link).parents("div.container:first").find("input[type='hidden']:first").attr('value', unit["name"]);
//			$(link).parents("div.container:first").find(".error:first").html("");
//			$(link).html(unit["name"]);
			unitWin.closeTab();
		}
	}
	
	this.doOption = function(link)
	{
	}
	
//	this.prepareLink = function(link)
//	{
//		self.link = link;
//		self.cell = $(self.link).parents("th");
//		self.label = $(self.cell).contents().find("label").get(0);
//		
//		self.name = $(self.cell).contents().find("input[type='hidden']").get(0);
//		self.value = $(self.cell).contents().find("input[type='hidden']").get(1);
//		
//		self.comment = $(self.cell).contents().find("nobr[name='comment']").get(0);
//		
//		self.unitlink = $(self.cell).contents().find("a[action='unitmap']");
//		self.articlelink = $(self.cell).contents().find("a[action='articlemap']");
//	}
//	
//	this.doRecordmenu = function(link)
//	{
//		self.prepareLink(link);
//		
//		self.recordMenu.cfg.setProperty("context", [$(link).attr('id'), "tl", "tr"]); 
//		self.recordMenu.show();
//	}
//	
//	this.doUnitmap = function(link)
//	{
//		self.prepareLink(link);
//		
//		var unitWin = openTab("Select unit", webRoot + "unit/show.do?render-mode=popup&pname="+URLEncode($(self.name).attr('value')));
//
//		unitWin.callback = function(newUnit)
//		{
//			$(self.value).attr('value', newUnit["name"]);
//			
//			var unitName = newUnit.name;
//			if (unitName == '')
//				unitName = 'dimensionless';
//							
//			$(self.unitlink).html('['+unitName+']');
//			$(self.articlelink).html('');
//
//			unitWin.closeTab();
//		}
//
//	}
//	
//	this.doArticlemap = function(link)
//	{
//		self.prepareLink(link);
//		
//		var articleWin = openTab("Select article", webRoot + "article/show.do?render-mode=popup");
//
//		articleWin.callback = function(newArticle)
//		{
//			var articleName = "A"+newArticle["id"];
//			$(self.value).attr('value', articleName);
//			
//			$(self.unitlink).html();
//			$(self.articlelink).html('['+articleName+']');			
//			
//			articleWin.closeTab();
//		}
//
//	}	
//	this.handleDialogSubmit = function(link)
//	{
//		var value = $("#knownDialog").find("#knownColumn").val();
//		$(self.name).attr('value', value);
//		$(self.value).attr('value', '');
//		$(self.label).html(value);
//		$(self.cell).addClass("green").removeClass("dgreen").removeClass("red");
//		$(self.comment).html("Known column");
//		$(self.unitlink).addClass("invisible");
//		
//		if (value == "ARTICLEID")
//		{
//			$(self.articlelink).removeClass("invisible");
//			$(self.articlelink).html('[article]');			
//		}
//		else
//		{
//			$(self.articlelink).addClass("invisible");
//			$(self.articlelink).html('');			
//		}
//		
//		this.hide();
//	}
//
//	this.handleDialogCancel = function(link)
//	{
//		this.cancel();
//	}

}

function BatchUploadBrowser()
{
	this.controller = "batchuploadnew";
	Browser.call(this);
	this.url = "batchuploadnew/browser_list.do";
	this.actionURL = "batchuploadnew/browser_action.do";
	this.itemElement = "exp-property-wrapper";	
	this.itemTemplate = "js/templates/batchuploadnew-record.ejs";
	var self = this;
	var hide = 0;
	var savedParameters = "";
	this.parentActionQuery = this.getActionQuery;
	this.pager.selectors.pager = ".pgr";
	
	this.onItemDrawn = function()
	{
		this.currentBlock.find('a').tooltip({showURL: false});
	}
	
	this.doDuplicatetoggle = function(link)
	{
		var table = $(link).parents("table").find("tr").get(1);
		$(table).toggleClass("invisible");
		if ($(table).hasClass("invisible"))
			$(link).html("[show&gt;&gt;]");
		else
			$(link).html("[&lt;&lt;hide]");
	}
	
	this.getActionQuery = function(action)
	{
		query = self.parentActionQuery();
		if (action == "state")
		{
			var radioButton = this.currentBlock.find("[action='state']").each(function()
			{
				if ($(this).attr("checked") == true) //doesn't work as a normal selector
				{
					var value =  ($(this).attr('value') == undefined) ? "" : URLEncode($(this).attr('value'));
					query += "&" + $(this).attr('name') + "=" + value;
				}
			});
		} else if (action == "batchstate")
		{
			query += "&"+savedParameters;
		}
		return query;
	}
	
	this.onStateSuccess = function()
	{
	}
	
	this.onBatchstateSuccess = this.onBatchstateFailure = function()
	{
		self.request();
	}
	
	this.doBatchmenu = function(link)
	{
		$("#batchoperations").dialog("open");
		$($("#batchoperations").dialog("widget")).position({ my: "right top", at: "right bottom", of: link } );
		hide = 2;
	}
	
	this.doBatch = function(link)
	{
		savedParameters="btype="+$(link).attr("type")+"&boperation="+$(link).attr("operation")
		self.callAction("batchstate",link,null);
	}
	
	this.dialogs = function()
	{
		$("#batchoperations").dialog({
			autoOpen: false,
			height: 170,
			width: 350,
			modal: false
		});
		
		$(document).click(function() {
			 hide--;
			 if (hide <= 0)
				 $("#batchoperations").dialog("close");
		 });
	}
	
}

function BatchUploadBuBrowser()
{
	this.controller = "batchuploadnew";
	Browser.call(this);
	this.url = "batchuploadnew/bubrowser_list.do";
	this.actionURL = "batchuploadnew/bubrowser_action.do";
	this.itemElement = "batchupload";	
	this.itemTemplate = "js/templates/batchuploadnew-burecord.ejs";
	var self = this;

	this.parentActionQuery = this.getActionQuery;
	this.pager.selectors.pager = ".pgr";
	
	this.onItemDrawn = function()
	{
		this.currentBlock.find('a').tooltip({showURL: false});
	}
	
	this.getActionQuery = function(action)
	{
		query = self.parentActionQuery();
		if (action == "state")
		{
		} 
		return query;
	}
	
	this.onStateSuccess = function()
	{
	}
	
	this.onBatchstate_uploadnewSuccess = this.onBatchstate_skipSuccess = function()
	{
		self.request();
	}
	
}


function ColumnRemappingActionable()
{
	Actionable.call(this);
	EventDispatcher.call(this);
	
	var self = this;
	
	var hide = 0;
	
	this.doMapkc = function(link) //Map known column
	{
		$("#knowncolumns").find("option").each(function(index,element){
			element.handler = link.handler;
		});
		$("#knowncolumns").dialog("open");
	}
	
	this.doMapkcSubmit = function(link) //Map known column
	{
		var handler = link.handler;
		$(handler).html(link.value);
		var td = $(handler).parent();
		td.find("input[type='checkbox']").attr("checked","checked");
		td.find("input[name='name']").attr("value",link.value);
		td.removeClass("dgreen").removeClass("red").addClass("green");
		//$("#knowncolumns").dialog("open");
	}
	
	this.doMappr = function(link) //Map property
	{
		var propWin = openTab("Select property", webRoot + "properties/show.do?render-mode=popup");
		var handler = link.handler;
		propWin.callback = function(property)
		{
			$(handler).html(property["name"]);
			var td = $(handler).parent();
			td.find("input[type='checkbox']").attr("checked","checked");
			td.find("input[name='name']").attr("value",property["name"]);
			td.removeClass("green").removeClass("red").addClass("dgreen");
			propWin.closeTab();
		}		
	}

	this.doMapco = function(link) //Map condition
	{
		var condWin = openTab("Select condition", webRoot + "properties/show.do?render-mode=popup&condition=true");
		var handler = link.handler;
		condWin.callback = function(condition)
		{
			$(handler).html(condition["name"]);
			var td = $(handler).parent();
			td.find("input[type='checkbox']").attr("checked","checked");
			td.find("input[name='name']").attr("value",condition["name"]);
			td.removeClass("green").removeClass("red").addClass("dgreen");
			condWin.closeTab();
		}		
	}
	
	this.doMenu = function(link)
	{
		$("#menu").find("a").each(function(index,element){
			element.handler = link;
		});
		$("#menu").menu().show().position({ my: "left top", at: "right top", of: link } );
		hide = 2;
	}
	
	this.init = function()
	{
		 $(document).click(function() {
			 hide--;
			 if (hide <= 0)
				 $("#menu").hide();
		 });
		 
		$("#knowncolumns").dialog({
			autoOpen: false,
			height: 170,
			width: 280,
			modal: true,
			buttons: {
				"OK": function() {
					self.doMapkcSubmit($(this).find("option:selected").get(0));
					$(this).dialog("close");
				},
				Cancel: function() {
					$(this).dialog("close");
					}
				}
		});
	}

}

function setValidationHook()
{
	$("input[type=file]").change(function(){
		var ext = $(this).val().split('.').pop().toLowerCase();
		if($.inArray(ext, ['csv','xls','sdf','xlsx']) == -1) 
		{
		    $("b[class=error]").html("Filename extension is not one of the following: CSV, XLS, XLSX, SDF");
		    $("input[type=submit]").attr("disabled","disabled");
		} else
		{
			$("b[class=error]").html("");
			$("input[type=submit]").removeAttr("disabled");
		}
	});
}
