
function PropertyBrowser(controller, item)
{
	var self = this;
	this.controller = controller;
	Browser.call(this);
	this.itemElement = item;
	this.itemTemplate = 'js/templates/property.ejs';
	this.pager.selectors.pager = ".pgr";
	
	if (getParams['condition'] == "true")
	{
		this.filters.setValue('condition', 'true');
		this.printedItemName = "condition";
	}
	
	this.beforeDelete = function()
	{
		return window.confirm("You are going to delete the property: \"" + this.currentEntity.name + "\"\nAre you sure you want to continue?");
	}
	
	this.onDeleteSuccess = function()
	{
		this.deleteRecord();
	}	
	
	this.onItemDrawn = function()
	{
		this.currentBlock.find('i[title]').tooltip({showURL: false});
	}
	
	
	this.beforeRemovechild = function()
	{
		return window.confirm("This will remove property \"" + this.currentEntity.name + "\" from group \"" + this.currentEntity.parent.name + "\"\nContinue?");
	}
	
	this.beforeApprove = function()
	{
		return window.confirm("Are you sure you want to approve this property ("+self.currentEntity.name+")?");
	}
	
	this.beforeUnapprove = function()
	{
		return window.confirm("Are you sure you want to un-approve this property?\n This will also unapprove all the records for this property.");
	}
	
	this.onRemovechildSuccess = this.onPublishSuccess = this.onApproveSuccess = this.onUnapproveSuccess = function()
	{
		this.request(false);
	}
	
	this.doAddchildselect = function(link)
	{
		var selectWin = openTab("Add a property in the group &lt;" + this.currentEntity.name + "&gt;", webRoot + "properties/show.do?directories=false");
		selectWin.callback = function(selectedItem)
		{  
			self.callAction("addchild", link, 
					{
						"data": "parent=" + self.currentEntity.id + "&child=" + selectedItem.id,
						"success": function(){
							self.request(false);
							selectWin.closeTab();
						}});
		}
	}
}


