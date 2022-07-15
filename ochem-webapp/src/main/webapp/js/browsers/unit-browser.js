function UnitBrowser(controller, item)
{
	this.controller = controller;
	Browser.call(this);
	this.itemElement = item;
	this.itemTemplate = 'js/templates/unit.ejs';
	this.pager.selectors.pager = ".pgr";
	
	this.onDeleteSuccess = function()
	{
		this.deleteRecord();
	}	
	
	this.onItemDrawn = function()
	{
		this.currentBlock.find('i[title]').tooltip({showURL: false});
	}
}