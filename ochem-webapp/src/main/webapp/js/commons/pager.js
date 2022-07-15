// Pager : Basic code
// Please think it over before modifying anything here
// YAHOO.History library being available, pager supports browser history buttons 
// Midnighter

// FIXME: Problems with history, pager doent work sometimes

var
	Pager = function(_maxPages, doNotUseHistory)
	{
		this.useHistory = false;//(doNotUseHistory == "undefined");
		
		this.selectors = new Object();
		this.selectors.showed = ".showed";
		this.selectors.total = ".total";
		this.selectors.pager = "#pager";
		this.selectors.strip = ".pager-strip";
		this.selectors.scope = document;
		this.maxPages = _maxPages;
		var self = this;
		
		this.load = function(webList)
		{
			this.totalNum = parseInt(webList.size);
			this.itemsPerPage = parseInt(webList.pageSize);
			this.currentPage = parseInt(webList.pageNum);
			this.firstResult = (this.totalNum == 0) ? this.totalNum : parseInt(webList.firstResult);
			this.lastResult = parseInt(webList.lastResult);
			this.totalPages = Math.floor((1*this.totalNum - 1) / this.itemsPerPage + 1);
			this.queryTime = webList.queryTime;
		}
		
		this.setVisible = function(visible) {
			var scope = this.selectors.scope ? this.selectors.scope : document;
			var pagerBlock = this.pagerBlock = $(scope).find(this.selectors.pager);
			$(scope).find(this.selectors.strip).setClass("invisible", !visible);
			$(scope).find(this.selectors.pager).setClass("invisible", !visible);
		}
		
		this.render = function()
		{
			var scope = this.selectors.scope ? this.selectors.scope : document;
			var pagerBlock = this.pagerBlock = $(scope).find(this.selectors.pager);
			pagerBlock.html('');
			this.visible = false;
			pagerBlock.removeClass('pager');
			pagerBlock.addClass('invisible');
			$(scope).find(this.selectors.showed).html(''+this.firstResult+' - '+this.lastResult);
			
			if (this.queryTime)
				$(scope).find(this.selectors.total).html("" + this.totalNum + " (query time "+this.queryTime+" sec.)")
			else
				$(scope).find(this.selectors.total).html(this.totalNum);
			
			if (this.totalNum > this.itemsPerPage)
			{				
				pagerBlock.addClass("pager");
				pagerBlock.removeClass("invisible");
				this.renderInternal();
				pagerBlock.find('a[page]').attr('href', '#');
				pagerBlock.find('[page]').click(function() 
				{
					self.currentPage = $(this).attr('page');
					self.changeState();
					return false;
				});
			}
			
			if (self.itemsPerPage)
				pagerBlock.find('select').val(self.itemsPerPage);
			
			pagerBlock.find('select').change(function(){
				var oldPageSize = self.itemsPerPage;
				self.itemsPerPage = $(this).val();
				self.currentPage = Math.ceil(self.currentPage*oldPageSize / self.itemsPerPage);
				self.changeState();
				});
		};
		
		this.setVisibility = function(visible)
		{
			var scope = this.selectors.scope ? this.selectors.scope : document;
			var pagerBlock = this.pagerBlock = $(scope).find(this.selectors.pager);
			pagerBlock.setClass("invisible", !visible);
		}
		
		this.renderInternal = function()
		{
			var scope = this.selectors.scope;
			var startPage = Math.max(1, 1*this.currentPage - (this.maxPages-1)/2);
			var endPage = Math.min(1*this.totalPages, 1*startPage + this.maxPages - 1);
			//window.alert(" startpage "+startPage+"  total "+totalPages+"  endpage "+endPage);
			if (1*this.currentPage > 1)
				$(scope).find(this.selectors.pager).append('<a href="#" page="'+(1*this.currentPage - 1)+'">prev</a>');
			for (var i = startPage; i <= endPage; i++)
				if (i != this.currentPage)
					$(scope).find(this.selectors.pager).append('<a href="#" page="'+i+'">'+i+'</a>');
				else
					$(scope).find(this.selectors.pager).append('<b>'+i+'</b>');
			if (1*this.currentPage < this.totalPages)
				$(scope).find(this.selectors.pager).append('<a href="#" page="'+(1*this.currentPage + 1)+'">next</a>');
		}
		
		this.changeState = function()
		{
			if (this.useHistory)
				YAHOO.util.History.navigate("pager", self.currentPage+","+self.itemsPerPage);
			else
				this.onStateChanged();
		}

		
		// Overridable function. Actual action on pager state change		
		this.onStateChanged = function()
		{
			// Stub
			window.alert('Current state: Page '+this.currentPage);
		}
		
		if (this.useHistory)
		{
			this.onHistoryEvent = function(state)
			{
				var parts = state.split(",");
				self.currentPage = parts[0];
				if (self.pagerBlock)
					self.pagerBlock.find('select').val(parts[1]);
				else
					self.itemsPerPage = parts[1];
				self.onStateChanged(); 
			}
			YAHOO.util.History.register("pager", "1,", this.onHistoryEvent);
			YAHOO.util.History.onReady(function () { 
				var currentState = 
					YAHOO.util.History.getBookmarkedState("pager") || YAHOO.util.History.getCurrentState("pager");
				if (currentState)
				{
					//window.alert(currentState);
					//YAHOO.util.History.navigate("pager", currentState);
					self.onHistoryEvent(currentState);
				}
			});
		}
	}
	
	function DirectPager(_maxPages, doNotUseHistory)
	{
		var self = this;
		Pager.call(this, _maxPages, doNotUseHistory);
		
		this.renderInternal = function()
		{
			var scope = this.selectors.scope;
			var pagerBlock = $(scope).find(this.selectors.pager);
			
			if (this.currentPage > 1)
				pagerBlock.append("<a page='1'>&lt;&lt;</a><a page='"+(this.currentPage - 1)+"'>&lt;</a>");
			
			var pageOptions = [5, 10, 15, 30, 50, 100];
			var select = $('<select filter="1" name="pagesize"/>');
			for (var i = 0; i < pageOptions.length; i++)
				select.append('<option value="'+pageOptions[i]+'">'+pageOptions[i]+'</option>');
			pagerBlock.append(select);
			pagerBlock.append(' items on ');
			select.val(this.itemsPerPage);
			
			pagerBlock.append("page <input type='text' id='pageInput' value='" + this.currentPage + "'/> of " + this.totalPages);
			if (this.currentPage < this.totalPages)
				pagerBlock.append("<a page='"+(this.currentPage + 1)+"'>&gt;</a><a page='" + this.totalPages + "'>&gt;&gt;</a>");
			
			pagerBlock.find('input').keypress(
				function(e)
				{
					if (e.which == 13)
					{
						var value = $(this).attr('value');
						if (parseInt(value) != value - 0)
							window.alert('Please enter a number');
						else if (parseInt(value) > self.totalPages && parseInt(value) > 0)
							window.alert('Only ' + self.totalPages + ' pages are available');
						else
						{
							self.currentPage = value;
							self.changeState();
						}
					}	
				}
			);
		}
	}
	
	function DirectPagerLite(_maxPages, doNotUseHistory)
	{
		var self = this;
		DirectPager.call(this, _maxPages, doNotUseHistory);
		
		this.renderInternal = function()
		{
			var scope = this.selectors.scope;
			var pagerBlock = $(scope).find(this.selectors.pager);
			
			if (this.currentPage > 1)
				pagerBlock.append("<a page='1'>&lt;&lt;</a><a page='"+(this.currentPage - 1)+"'>&lt;</a>");
			
			pagerBlock.append("page <input type='text' id='pageInput' value='" + this.currentPage + "'/> of " + this.totalPages);
			if (this.currentPage < this.totalPages)
				pagerBlock.append("<a page='"+(this.currentPage + 1)+"'>&gt;</a><a page='" + this.totalPages + "'>&gt;&gt;</a>");
			
			pagerBlock.find('input').keypress(
				function(e)
				{
					if (e.which == 13)
					{
						var value = $(this).attr('value');
						if (parseInt(value) != value - 0)
							window.alert('Please enter a number');
						else if (parseInt(value) > self.totalPages && parseInt(value) > 0)
							window.alert('Only ' + self.totalPages + ' pages are available');
						else
						{
							self.currentPage = value;
							self.changeState();
						}
					}	
				}
			);
		}
	}