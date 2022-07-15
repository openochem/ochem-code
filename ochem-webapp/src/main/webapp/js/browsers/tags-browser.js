var TagsBrowser = function(type)
{
	var self = this;
	this.controller = "tags";
	UnifiedBrowser.call(this);
	this.container = "TagsBrowser-" + type;
	this.scope = ".scope-" + type;
	this.itemElement = "tag";
	this.itemTemplate = "js/templates/tag-inline.ejs";
	this.filters.values.filterby = "filters";
	this.filters.values.type = type;
	this.minusId = 0;
	
	this.onItemDrawn = function()
	{
		this.currentBlock.css({float: "left", width: "auto", padding: "0px 5px 0px 0px"});
		this.mainDiv.find('span').remove();
	}
	this.doDelete = function()
	{
		this.deleteRecord();
		this.callAction("saveall");
	}
	
	this.doAdd = function()
	{
		
		var selectWin = openTab("Edit tag", webRoot + "tags/show.do?render-mode=popup&type=" + this.filters.values.type);
		selectWin.callback = function(selectedItem)
		{  
			appropriateBrowser = tagsBrowser[selectedItem.type];
			appropriateBrowser.drawFromJSON({id: selectedItem.id, name: selectedItem.name});
			//self.setValue('tag-id', selectedItem.id);
			//self.currentBlock.find('a[name="tag-link"]').html(selectedItem.name);
			appropriateBrowser.callAction("saveall");
			selectWin.closeTab();
		}
	}
	
	this.listenEvent("items_loaded", function() {
		if (self.pager.totalNum == 0)
			self.mainDiv.html("<span class='notags'>No tags for "+type+" selected</span>");
	})
	
//	this.doTagclick = function()
//	{
//		openTab("Wiki", webRoot + 'wikipage/action.do?name='+this.currentEntity.name);
//	}
}

var tagsBrowser = {molecule: new TagsBrowser("molecule"), property: new TagsBrowser("property")}
$(document).ready(function() {
	tagsBrowser.molecule.initialize();
	tagsBrowser.property.initialize();
});