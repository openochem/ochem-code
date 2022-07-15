function TooltipMenu()
{
	this.menuSelector = "#ttmenu";
	this.parentItemDrawn = this.onItemDrawn;
	this.parentInitialize = this.initialize;

	var self = this;
	
	this.initializeMenu = function()
	{
		$(self.menuSelector).hover(
			function(){
				$(self.menuSelector).removeClass("invisible").addClass("softhighlight");
			}, 
			function(){
				$(self.menuSelector).addClass("invisible").removeClass("softhighlight");
			}	
		);
	}
	
	this.menuHoverOn = function(e)
	{
		var t = $(self.menuSelector);
		var p = $(this).position();
		t.css({top: p.top, left: p.left})
		t.removeClass("invisible");		
		self.activeBlock = e.currentTarget;	
	}
	
	this.menuHoverOff = function(e)
	{
		var t = $(self.menuSelector);
		t.addClass("invisible");	
	}
	
	this.onItemDrawn = function()
	{
		//console.log("Item "+self.currentBlock.attr("rec-id")+" drawn in "+self.container+" on "+Date.now());
		self.currentBlock.hover(this.menuHoverOn, this.menuHoverOff);
		if (this.parentItemDrawn)
			this.parentItemDrawn();
	}
	
	this.initialize = function(arg)
	{
		this.initializeMenu();
		
		if (this.parentInitialize)
			this.parentInitialize(arg);
	}
}