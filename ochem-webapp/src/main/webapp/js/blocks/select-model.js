var ModelsBrowser = function()
{
	var self = this;
	UnifiedBrowser.call(this);
	this.container = "ModelsBrowser";
	this.itemTemplate = "js/templates/models/model-inline.ejs";
	this.minusId = 0;

	this.onItemDrawn = function()
	{
		//this.currentBlock.css({float: "left", width: "auto", padding: "0px 5px 0px 0px"});
		this.mainDiv.find('span').remove();
	}

	this.doDelete = function()
	{
		this.deleteRecord();
	}
	
	this.doAdd = function()
	{
		var selectWin = openTab("Select model", webRoot + "/model/select.do");
		selectWin.callback = function(models)
		{
			var lmodels = array(models);
			for (var i = 0; i < lmodels.length; i++)
			{
				var selectedItem = lmodels[i];
				self.drawFromJSON({id: selectedItem.id, name: selectedItem.name});
			}
			selectWin.closeTab();
			$("#scenario").removeClass("invisible");
		}
	}
	
	this.listenEvent("items_loaded", function() {
		if (self.pager.totalNum == 0)
			self.mainDiv.html("<span>No tags selected</span>");
	})
}
include.plugins("view");
var modelBrowser = new ModelsBrowser();

$(document).ready(function(){
	modelBrowser.initialize();
});
