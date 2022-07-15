//PropertyBrowser.Inherits(Browser);
function AuthorBrowser(controller, item)
{
	this.controller = controller;
	Browser.call(this);
	this.itemElement = item;
	
	this.draw = function(entity)
	{
		return new View({url: 'js/templates/author.ejs'}).render(entity);
	}
	
	this.onDeleteSuccess = function()
	{
		this.deleteRecord();
	}
}
