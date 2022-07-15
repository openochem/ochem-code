var TagsBrowser = function()
{
	var self = this;
	this.controller = "tags";
	UnifiedBrowser.call(this);
	this.container = "TagsBrowser";
	this.itemElement = "tag";
	this.itemTemplate = "js/templates/tag-inline.ejs";
	this.filters.values.filterby = "property";
	this.minusId = 0;
	this.changed = false;
	
	this.onItemDrawn = function()
	{
		this.currentBlock.css({float: "left", width: "auto", padding: "0px 5px 0px 0px"});
		this.mainDiv.find('span').remove();
	}
	
	this.doDelete = function()
	{
		this.changed = true;
		this.deleteRecord();
	}
	
	this.doAdd = function()
	{
		var selectWin = openTab("Select a tag", webRoot + "tags/show.do?render-mode=popup&type=property");
		selectWin.callback = function(selectedItem)
		{  
			tagsBrowser.drawFromJSON({id: selectedItem.id, name: selectedItem.name});
			self.changed = true;
			selectWin.closeTab();
		}
	}
	
	this.listenEvent("items_loaded", function() {
		if (self.pager.totalNum == 0)
			self.mainDiv.html("<span>No tags selected</span>");
	})
	
	this.doTagclick = function()
	{
		if (this.getValue('tag-id') < 0)
		{
			var selectWin = openTab("Edit tag", webRoot + "tags/show.do?render-mode=popup");
			selectWin.moveBy(50, 50);
			selectWin.callback = function(selectedItem)
			{
				self.changed = true;
				self.setValue('tag-id', selectedItem.id);
				self.currentBlock.find('a[name="tag-link"]').html(selectedItem.name);
				selectWin.closeTab();
			}
		}
//		else
//			window.open(webRoot + 'wikipage/action.do?name='+this.currentEntity.name, 'wiki', "location=2,status=0,scrollbars=yes,resizable=yes,width=670,height=565");
	}
}

var tagsBrowser = new TagsBrowser();
$(document).ready(function() {
	tagsBrowser.initialize();
});