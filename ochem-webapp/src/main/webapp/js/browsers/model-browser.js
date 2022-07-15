function ModelBrowser(controller, item)
{
	this.controller = controller;
	Browser.call(this);
	this.itemElement = item;
	this.itemTemplate = 'js/templates/models/model.ejs';
	this.pager.selectors.pager = ".pgr";
	var self = this;
	this.selectedEntities = [];
	this.singleSelect = false;
	
	
	this.filters.setFromUrl();
	
	this.onTogglebasketSuccess = function(xml)
	{
		var newValue = (this.currentEntity.selected == "true") ? "false" : "true";
		this.currentEntity.selected = newValue;
		this.currentBlock.find('img[name="checked"]').attr('src', 
			(newValue == "true") ? "img/icons/checked.gif" : "img/icons/unchecked.gif");
	}
	
	this.onItemDrawn = function()
	{
		var predictionScenario = this.currentEntity.predictionScenario;
		if(predictionScenario >= 1)
			this.currentBlock.addClass("highlightestimated");
		
		if (window.callback) // We draw the javascript selection checkboxes 
		{
			var index = this.findEntitySelected();
			if (index == -1)
			{
				$(this.currentBlock).find("[action='selectmodel']").html("<img name=\"checked\" src=\"img/icons/unchecked.gif\"/>");			
			}
			else
			{
				$(this.currentBlock).find("[action='selectmodel']").html("<img name=\"checked\" src=\"img/icons/checked.gif\"/>");
			}
			
			$(this.currentBlock).find(".apply-model").addClass("invisible");
		}
	}
	
	this.beforeDelete = function()
	{
		return confirm("Are your sure you want to delete this model?\nThis change is permanent and cannot be undone")
	}
	
	this.onDeleteSuccess = function()
	{
		this.deleteRecord();
	}
	
	this.findEntitySelected = function()
	{
		var index = -1;
		for (var i = 0; i < this.selectedEntities.length; i++)
			if (this.currentEntity.id == this.selectedEntities[i].id)
				index = i;
		return index;
	}
	
	this.doSelectmodel = function()
	{
		var index = this.findEntitySelected();
		if (index == -1)
		{
			this.selectedEntities.push(this.currentEntity);
			$(this.currentBlock).find("[action='selectmodel']").html("<img name=\"checked\" src=\"img/icons/checked.gif\"/>");
		}
		else
		{
			this.selectedEntities.splice(index, 1);
			$(this.currentBlock).find("[action='selectmodel']").html("<img name=\"checked\" src=\"img/icons/unchecked.gif\"/>");			
		}
	}
	
	this.doSelectsinglemodel = function()
	{
		window.callback(this.currentEntity);
	}
	
	this.doSelectsubmit = function()
	{
		if (this.selectedEntities.length > 0)
			window.callback(this.selectedEntities);
		else
			window.alert("Please select at least one model!");
	}
	
	this.doShowxml = function()
	{
		$("#xmlDialog pre").html($("div#xml-" + this.currentRecordId).html().replace(/</g, "&lt;").replace(/>/g, "&gt;"));
		xmlDialog.show();
	}
	
	this.listenEvent('items_loaded', function()
	{
		$('span').tooltip({showURL: false});
		if (window.callback) // we have a model selection dialog
			$("#Browser").find("a[action=togglebasket]").remove();
	});
	
	this.originalInitialize = this.initialize;
	
	this.initialize = function(dontRunQuery)
	{
		this.originalInitialize(dontRunQuery);
		if (window.callback) // we have a model selection dialog
		{
			$(".formsubmit").addClass("invisible");
			this.singleSelect = (getParams["single"] != undefined);
			if (!this.singleSelect)
				$("a[action='selectsubmit']").removeClass("invisible");
		}
	}
}