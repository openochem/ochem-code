function XlspreviewActionable()
{
	Actionable.call(this); // superclass(es)
	EventDispatcher.call(this);

	var self = this;
	
	this.doKnown = function(link)
	{
		this.knownDialog.show();
	}
	
	this.doProperty = function(link)
	{
		
		var propWin = openTab("Select property", webRoot + "properties/show.do?render-mode=popup");

		propWin.callback = function(newProperty)
		{
			$(self.name).attr('value', newProperty["name"]);
			
			
			var unitName = newProperty.defaultUnit.name;
			
			$(self.value).attr('value', unitName);
			
			if (unitName == '')
				unitName = 'dimensionless';
							
			$(self.unitlink).html('['+unitName+']');
			$(self.articlelink).html('');
			
			$(self.label).html(newProperty["name"]);
			$(self.cell).addClass("dgreen").removeClass("green").removeClass("red");
			$(self.comment).html("Property");
			
			if (newProperty.qualitive == "false")
				$(self.unitlink).removeClass("invisible");
			else
				$(self.unitlink).addClass("invisible");
			
			$(self.articlelink).addClass("invisible");
			propWin.closeTab();
		}
	}

	this.doCondition = function(link)
	{
		var propWin = openTab("Select condition", webRoot + "properties/show.do?render-mode=popup&condition=true");

		propWin.callback = function(newProperty)
		{
			$(self.name).attr('value', newProperty["name"]);
			

			var unitName = newProperty.defaultUnit.name;
			
			$(self.value).attr('value', unitName);
			
			if (unitName == '')
				unitName = 'dimensionless';
							
			$(self.unitlink).html('['+unitName+']');
			$(self.articlelink).html('');
			
			$(self.label).html(newProperty["name"]);			
			$(self.cell).addClass("dgreen").removeClass("green").removeClass("red");
			$(self.comment).html("Condition");

			if (newProperty.qualitive == "false")
				$(self.unitlink).removeClass("invisible");
			else
				$(self.unitlink).addClass("invisible");
			
			$(self.articlelink).addClass("invisible");
			propWin.closeTab();
		}
	}
	
	this.prepareLink = function(link)
	{
		self.link = link;
		self.cell = $(self.link).parents("th");
		self.label = $(self.cell).contents().find("label").get(0);
		
		self.name = $(self.cell).contents().find("input[type='hidden']").get(0);
		self.value = $(self.cell).contents().find("input[type='hidden']").get(1);
		
		self.comment = $(self.cell).contents().find("nobr[name='comment']").get(0);
		
		self.unitlink = $(self.cell).contents().find("a[action='unitmap']");
		self.articlelink = $(self.cell).contents().find("a[action='articlemap']");
	}
	
	this.doRecordmenu = function(link)
	{
		self.prepareLink(link);
		
		self.recordMenu.cfg.setProperty("context", [$(link).attr('id'), "tl", "tr"]); 
		self.recordMenu.show();
	}
	
	this.doUnitmap = function(link)
	{
		self.prepareLink(link);
		
		var unitWin = openTab("Select unit", webRoot + "unit/show.do?render-mode=popup&pname="+URLEncode($(self.name).attr('value')));

		unitWin.callback = function(newUnit)
		{
			$(self.value).attr('value', newUnit["name"]);
			
			var unitName = newUnit.name;
			if (unitName == '')
				unitName = 'dimensionless';
							
			$(self.unitlink).html('['+unitName+']');
			$(self.articlelink).html('');

			unitWin.closeTab();
		}

	}
	
	this.doArticlemap = function(link)
	{
		self.prepareLink(link);
		
		var articleWin = openTab("Select article", webRoot + "article/show.do?render-mode=popup");

		articleWin.callback = function(newArticle)
		{
			var articleName = "A"+newArticle["id"];
			$(self.value).attr('value', articleName);
			
			$(self.unitlink).html();
			$(self.articlelink).html('['+articleName+']');			
			
			articleWin.closeTab();
		}

	}	
	this.handleDialogSubmit = function(link)
	{
		var value = $("#knownDialog").find("#knownColumn").val();
		$(self.name).attr('value', value);
		$(self.value).attr('value', '');
		$(self.label).html(value);
		$(self.cell).addClass("green").removeClass("dgreen").removeClass("red");
		$(self.comment).html("Known column");
		$(self.unitlink).addClass("invisible");
		
		if (value == "ARTICLEID")
		{
			$(self.articlelink).removeClass("invisible");
			$(self.articlelink).html('[article]');			
		}
		else
		{
			$(self.articlelink).addClass("invisible");
			$(self.articlelink).html('');			
		}
		
		this.hide();
	}

	this.handleDialogCancel = function(link)
	{
		this.cancel();
	}

}

var xlsAct = new XlspreviewActionable();

$(document).ready(function() {
	xlsAct.recordMenu = new YAHOO.widget.Menu("basicmenu");
	xlsAct.recordMenu.render();
	
	xlsAct.knownDialog = new YAHOO.widget.Dialog("knownDialog", { 
    	width:"325px", 
		fixedcenter:true, 
		modal:true, 
		visible:false 
    });
	var myButtons = [{ text:"Submit", handler:xlsAct.handleDialogSubmit}, 
	                 { text:"Cancel", handler:xlsAct.handleDialogCancel } ]; 
	xlsAct.knownDialog.cfg.queueProperty("buttons", myButtons);
	xlsAct.knownDialog.render();	
});