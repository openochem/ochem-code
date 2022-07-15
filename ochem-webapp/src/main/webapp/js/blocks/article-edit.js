
function ArticleForm()
{
	this.scope = ".formscope";
	AjaxForm.call(this);
	
	this.actionURL = "article/action.do";
	this.itemElement = "article";
	this.disabled = false;
	var self = this;
	
	this.doEditjournal = function()
	{
		if (!this.disabled)
		{
			var journalWin = openTab("Select journal", webRoot+"journal/show.do?render-mode=popup&id="+articleForm.getValue('n-journal'));
			journalWin.moveBy(50, 50);
			journalWin.callback = function(newJournal)
			{
				articleForm.journalAutocomplete.setKey(newJournal["id"]);
				self.setValue("n-journal-title", newJournal.title);
				journalWin.closeTab();
			}
		}
	}
	
	this.doIntroducepmid = function()
	{
		$("#inblock").toggleClass("invisible");
	}
	
	this.beforeEdit = function()
	{
		if (articleForm.getValue('media-type') == "book")
			window.alert("book comment");
	}
	
	this.beforesubmit = function()
	{
		$(document).find("[send]").not("[n-disable]").removeAttr("disabled");
		authorsBrowser.disabled = false;
		articleForm.disabled = false;
	}
	
	this.beforeLoad_isbn = function()
	{
		if($("#load-isbn").val() == "")
			window.alert("Please enter ISBN number");
		else
			return true;
	}
	
	this.onLoad_isbnSuccess = function(entity)
	{
		var win = openTab("Edit article", webRoot + "article/edit.do?render-mode=popup&id=-1&fromsession=1");
		win.tmp = window.callback;
		win.entity = entity;
		win.callback = function()
		{
			win.tmp.call(this,win.entity);
			win.closeTab();
		}
		window.closeTab();
	}
	
	this.doSelectparentbook = function()
	{
		var curId = articleForm.getValue('parent-book-id');
		var appendix = curId ? "&selected="+curId+"&identifier=Q"+curId : "";
		var booksWin = openTab("Select book", webRoot+"article/show.do?render-mode=popup&hide-chapters=1&media-type=book"+appendix);
		booksWin.moveBy(50, 50);
		booksWin.callback = function(book)
		{
			self.setValue("parent-book-id", book.id, book.title);
			booksWin.closeTab();
		}
	}
	
	this.onItemSaved = function(entity)
	{
		authorsBrowser.filters.setValue('id', entity.id);
		authorsBrowser.callAction("saveall");
	}
	
	this.doCancel = function()
	{
		window.closeTab();
	}
	
	this.doDeletepdf  = function()
	{
		$("#pdfavail").addClass("hidden");
		$("#nopdfavail").removeClass("hidden");
		this.setValue('n-delpdf','1');
	}
	
	this.doDeleteexcel  = function()
	{
		$("#batchavail").addClass("hidden");
		this.setValue('n-batch','1');
	}
	
	this.onReloadpmSuccess = function(xml)
	{
		this.entity = xml;
		
		
		$('input[name="n-title"]').val(xml.article.title);
		$('textarea[name="n-abstract"]').val(xml.article.articleAbstract);
		$('input[name="n-journal-key"]').val(null);
		$('input[name="n-journal-title"]').val(xml.article.journal.title);
		$('input[name="n-date"]').val(xml.article["publication-date"]["printed-name"]);
		$('input[name="n-volume"]').val(xml.article.volume);
		$('input[name="n-issue"]').val(xml.article.issue);
		$('input[name="n-pages"]').val(xml.article.pageNumbers);
		$('input[name="n-pubmed"]').val(xml.article.pmid);
		$('input[name="n-affiliation"]').val(xml.article.affiliation);
		$('input[name="n-url"]').val(xml.article.url);
		$('input[name="n-doi"]').val(xml.article.doi);
		$('textarea[name="n-comment"]').val(xml.article.comment);
		
		var authorsHtml = "";
		var authors = array(xml.article.authors.author);
		for (var i = 0; i < authors.length; i++) // (for (var i in authors)
		{
			authorsHtml += '<div style="" class="browser-item floating" rec-id="' + (-i-1) + '">' +
				'<input class="narrow" name="name" value="' + authors[i]["printed-name"] + '" send="1" type="text" readonly="readonly">' +
//				'<a href="javascript:void(0)" action="editauthor" title="Edit this author name">' +
//				'<img src="img/icons/edit.gif">' +
//				'</a>' + 
//				'<a href="javascript:void(0)" action="delete" class="delete-link">[x]</a>&nbsp;' +
				'</div>';
		}
		
		$("#ArticleAuthorsBrowser").html(authorsHtml);
		
		$(document).find("[send]").not("[n-disable]").removeAttr("disabled").attr("readonly", "readonly");
		
	}

	this.onReloadisbnSuccess = function(xml)
	{
		this.entity = xml;
		this.onItemSaved.call(this, xml);		
	}
	
	this.onIframeLoaded = function()
	{		
		var data = $("#iframe").contents();
		if ( (data.text() != "null") && (data.text() != "") )
		{			
			var json = eval("(" + data.text() + ")");
			if (typeof(json.message) != "undefined" && (json.message.type == "exception" || json.message.type == "error"))
			{
				if (json.message.title.match("article is duplicate"))
				{
					
					  var re = new RegExp("[0-9]{1,}");
					  var m = re.exec(json.message.title);
					  if (m == null) 
					  {
					    alert("No Article id. Please report the bug");
					  } 
					  else 
					  {
					  //	window.alert("articleid"+m[0]);
					  	$("input[name='art-id']").val(m[0]);
						$("#dupArticle").html(json.message.title);
						this.articleDialog.show();
					  }
				}
				else
					window.alert(json.message.title);
			}
			else
			{
				var item = eval('json["'+this.itemElement+'"]');
				this.entity = item;
				this.onItemSaved.call(this, item);
			}
		}		
	}
	
	this.onJournalChanged = function(journalId)
	{
		$("#journal-info").addClass("invisible");
		if (journalId > 0)
		this.ajax.send(
		{
			url: 'journal/edit.do', 
			data: 'out=json&id='+journalId, 
			success: function(response)
			{
				// Maybe put to separate JMVC template later
				$("#journal-info").html(""
							+ (response.journal.abbreviation ?  response.journal.abbreviation + ", " : "")
							+ "ISSN: "+response.journal.issn+", "+response.journal.articlesCount+" entries");
				$("#journal-info").removeClass("invisible");
			}
		});
	}

	//Currently I don't know how to handle this in multi-part request
	
	this.initialize = function()
	{
		
		this.articleDialog = new YAHOO.widget.Dialog("articleDialog", { 
	    	width:"525px", 
			fixedcenter:true, 
			modal:true, 
			visible:false 
	    });
	    
	    var myButtons = [{ text:"Show article", handler:handleSubmit}, 
						{ text:"Cancel", handler:handleCancel }]; 
		this.articleDialog.cfg.queueProperty("buttons", myButtons);
		this.articleDialog.render();
	}
	
	var handleSubmit = function() { 
		var dupArt = $("input[name='art-id']").val();
		var articleWin = openTab("Duplicate article", webRoot + "article/show.do?render-mode=popup&id="+dupArt);
		this.cancel();
	};
	
	var handleCancel = function() { 
		this.cancel();
	};
	
}

function ArticleAuthorsBrowser()
{
	this.controller = "article";
	this.scope = ".authorsscope";
	UnifiedBrowser.call(this);
	this.itemElement = "author";
	this.itemTemplate = "js/templates/articleauthors.ejs";
	this.url = "article/listauthors.do";
	this.actionURL = "article/saveauthors.do";
	this.container = "ArticleAuthorsBrowser";
	this.minusCounter = -1;
	this.disabled = false;
	this.edit = false;
	
	this.onItemDrawn = function()
	{
		this.currentBlock.addClass("floating").attr("style","");
		if(this.edit)
			$(".narrow").removeAttr("disabled");
	}
	
	this.doNew = function()
	{
		var pmid = articleForm.getValue("n-pubmed");
		if(pmid)
		{
			window.alert("connot be edited");
		}else
		{
			var entity = {"id":this.minusCounter--, "Name":" "};
			this.drawFromJSON(entity);
		}
	}
	
	this.doEditauthor =  function()
	{
		var pmid = articleForm.getValue("n-pubmed");
		if(pmid)
		{
			window.alert("PubMed article connot be edited");
		}else
		{
			var id = this.currentRecordId;
			this.currentBlock.find('[au-id='+id+']').removeAttr("disabled");
		}
	}
	
	this.doDelete = function()
	{
		//alert(this);
		var pmid = articleForm.getValue("n-pubmed");
		if(pmid)
		{
			window.alert("PubMed article connot be edited");
		}else
		{
			this.deleteRecord();
		}
	}
	
	this.getRecordIdentifier = function(entity)
	{
		return entity["id"];
	}
	
	this.beforeSwap = function()
	{
		var pmid = articleForm.getValue("n-pubmed");
		if(pmid)
		{
			window.alert("PubMed article connot be swaped");
		}else{
			return true;
		}
	}
	
	this.onSwapSuccess = function(entity)
	{
		this.edit = true;
		this.mainDiv.html('');
		var items;
		if (entity.list){
			items = entity.list[this.itemElement];
		}
		// Cool drawing of items
		if (items)
		{
			if (items instanceof Array)
			{
				for (var i = 0; i < items.length; i++){
					this.drawFromJSON(items[i]);
				}
			}
			else{
				this.drawFromJSON(items);
			}
		}
	}
	
	this.onSaveallSuccess = function(xml)
	{
		// After authors are saved, close the window
		articleForm.entity.oldId = oldId;
		if (window.callback)
			window.callback({"article": articleForm.entity});
		else
			location.reload();
	}
}

include.plugins('view');
var articleForm = new ArticleForm();
var authorsBrowser = new ArticleAuthorsBrowser();
var oldId = getParams["id"];


$(document).ready(
	function()
	{
		if (articleForm.getValue('media-type') == "article")
		{
			// Its an article
			var value = $("[name='journal-title']").val();
			var id = $("[name='journal-id']").val();
		
			var auto = articleForm.journalAutocomplete = new Autocomplete("n-journal", dsJournalsAC);
			$('#ac-source').append(auto.getJQuery());
			auto.oldSetKey = auto.setKey;
			auto.setKey = function(value)
			{
				auto.oldSetKey(value);
				articleForm.onJournalChanged(value);
			}
		
			auto.setKey(id);
			auto.setValue(value);
			auto.attachHandlers();
			$("[name='n-journal-key']").attr("send",1);
			$("[name='n-journal-title']").attr("send",1);
			var pmid = articleForm.getValue("n-pubmed");
			
			if (pmid)
			{
				$(document).find("[send]").not("[n-disable]").attr("disabled", "disabled");
				authorsBrowser.disabled = true;
				articleForm.disabled = true;
			}
		}
		else
		{
			// Its a book
			$("[name='is-chapter']").change(function(){
				$("#selectparentbook").setClass("invisible", !$(this).is(":checked"));
				$(".book-only").setClass("invisible", $(this).is(":checked"));
			});
			$("[name='is-chapter']").change();
			
			var isbnid = articleForm.getValue("id");
			if(isbnid)
			{
				$(document).find("[n-loadisbn]").attr("disabled", "disabled");
				authorsBrowser.disabled = true;
				articleForm.disabled = true;
			}

		}
		
		authorsBrowser.initialize();
		articleForm.initialize();
		tabView = new YAHOO.widget.TabView("demo");
	}
);

function ismaxlength(obj)
{
	var mlength=obj.getAttribute? parseInt(obj.getAttribute("maxlength")) : ""
	if (obj.getAttribute && obj.value.length>mlength)
	obj.value=obj.value.substring(0,mlength)
}
