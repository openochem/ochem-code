function BatchEditBrowser()
{
	this.controller = "batchedit";
	var self = this;
	Browser.call(this);
	this.itemElement = "compressedEP";	
	this.itemTemplate = "js/templates/batchedit-record.ejs";
	
	this.doEdit = function()
	{
		var id = this.currentEntity.id;
		var artid = this.currentEntity.ep.article.id;
		var propid = this.currentEntity.ep.property.id;
		var unitid = (this.currentEntity.ep.unit) ? this.currentEntity.ep.unit.id : "";
		
		var key = this.currentEntity.compressionKey;
		
		var win = openTab("Batch edit", webRoot+this.controller+"/edit.do?render-mode=popup&key="+key);
		win.callback = function(entity) 
		{
			if (self.onItemSaved)
			{
				self.onItemSaved(entity);
			}
				
			self.request(false, function()
			{
				self.setPositionById(entity["id"])
				win.closeTab();
			}); 
		};
			
	}
	
	this.onItemDrawn = function()
	{
		var recStatus = this.currentEntity.epEvidence;
		
		if (recStatus == "3" || recStatus == "4") 	// error or invalid
			this.currentBlock.addClass("highlighterror");
		if (recStatus == "2") 						// to be verified?
			this.currentBlock.addClass("highlightnon-v");
		
		this.currentBlock.removeClass("highlighted");
	}
	
	this.onBatcheditError = function(msg, msgEntity)
	{
		this.ajax.showError(msgEntity);
	}
	
}
