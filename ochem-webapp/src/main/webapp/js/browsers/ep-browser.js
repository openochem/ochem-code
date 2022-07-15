//CompoundBrowser.Inherits(Browser);
function CompoundBrowser()
{
	var self = this;
	this.controller = "epbrowser";
	this.editController = "eprecord";
	this.scope = ".browserScope";
	Browser.call(this);
	EventDispatcher.call(this);
	this.itemTemplate = "js/templates/record.ejs";
	this.itemElement = "exp-property";
	this.printedItemName = "record";
	this.basketFilter = new DynamicSelect("basket-select","basket/list.do?no.basket.details=1","basket");
	this.selectedBasketIndex = "Show all";
	this.currentAction = "";
	this.basketSelect = new DynamicSelect('n-basket', 'basket/list.do?no.basket.details=1', 'basket');
	this.parentActionQuery = this.getActionQuery;
	this.pager.selectors.pager = ".pgr";
	this.trash = false;

	CompoundBrowser.hintForDelete = getParams["trash"] ? "Delete this record permanently" : "Move this record to trash";

	this.onBasketFilterLoaded = function()
	{
		var elSel = self.basketFilter.select;
		var elOptOld = elSel.options[0];
		var elOptNew1 = new Option("Show all", -1, false, false);
		var elOptNew2 = new Option("Show records in any basket", "any", false, false);
		if (getParams["compare-baskets"])
		{
			var elOptNew3 = new Option("Compare Basket", getParams["basket-select"], false, false);
		}
		if (getParams["join-baskets"])
		{
			var elOptNew4 = new Option("Join Basket", getParams["basket-select"], false, false);
		}
		try 
		{
			elSel.add(elOptNew2, elOptOld); // standards compliant; doesn't work in IE
			elOptOld = elSel.options[0];
			elSel.add(elOptNew1, elOptOld); // standards compliant; doesn't work in IE
			if (getParams["compare-baskets"])
			{
				elOptOld = elSel.options[0];
				elSel.add(elOptNew3, elOptOld); // standards compliant; doesn't work in IE
			}
			if (getParams["join-baskets"])
			{
				elOptOld = elSel.options[0];
				elSel.add(elOptNew4, elOptOld); // standards compliant; doesn't work in IE
			}
		}
		catch(ex) 
		{
			elSel.add(elOptNew2, 0); // IE only
			elSel.add(elOptNew1, 0); // IE only
			if (getParams["compare-baskets"])
			{
				elSel.add(elOptNew3, 0); // IE only
			}
			if (getParams["join-baskets"])
			{
				elSel.add(elOptNew4, 0); // IE only
			}
		}	 

		for (i=0; i<self.basketFilter.select.options.length; i++)
		{	
			if (self.basketFilter.select.options[i].text == self.selectedBasketIndex)
			{
				self.basketFilter.select.selectedIndex = i;
				break;
			}
		}
		self.request();
	}

	this.parentsetValue = this.filters.setValue;
	this.filters.setValue = function(name, value, title)
	{
		if (name == "basket-select")
		{
			var elSel = self.basketFilter.select;
			var myBasket = false;
			for (var i = 0; i < elSel.options.length; i++)
			{
				if (elSel.options[i].value == value)
					myBasket = true;
			}
			if (!myBasket)
			{
				externalBasket = new Option(title, value, false, false);
				elSel.add(externalBasket, elSel.options[0]);
				self.basketFilter.select.selectedIndex = 0;
			}
		}

		if (name == "introducer")
		{
			if ($("[name=introducer]").find("option[value='"+value+"']").length == 0)
			{
				var option = $("<option/>");
				option.html("User ID: " + value);
				option.attr("value", value);
				$("[name=introducer]").append(option);
			}
		}

		self.parentsetValue.call(self.filters, name, value, title);	
	}

	this.updateBasketFilter = function()
	{
		self.selectedBasketIndex = self.basketFilter.select.options[self.basketFilter.select.selectedIndex].text;
		self.basketFilter.update("pagenum=1&pagesize=500&showsystem=1",null,"basketlist_update");
	}

	this.basketFilter.listenEvent("basketlist_update", this.onBasketFilterLoaded);

//	this.listenEvent("items_load", function(){
//	window.alert(getParams["basket-select"]);
//	if (getParams["basket-select"])
//	$("select[name='basket-select']").append('<li value="sdfsdfsd,yu">Combination of baskets</li>');
//	});

	this.beforeToggleselect = function()
	{
		// Changing checkbox status even before sending request to server
		var newValue = (this.currentEntity.selected == "true") ? "false" : "true";
		this.currentEntity.selected = newValue;
		this.currentBlock.find('img[name="checked"]').attr('src', 
				(newValue == "true") ? "img/icons/checked.gif" : "img/icons/unchecked.gif");
		if (newValue == "true")
			self.selectionSize++;
		else
			self.selectionSize--;
		self.updateSelectionSize();
		return true;
	}

	this.updateSelectionSize = function()
	{
		$(".selection").addClass("invisible");
		if (self.selectionSize)
		{
			if (self.selectionSize > 0)
				$(".selection").removeClass("invisible");
			$("#selection-size").html(self.selectionSize);
		}
	}

	this.onToggleselectSuccess = function(xml)
	{
	}

	this.beforeCleanbasket =
		this.beforeAddselect = 
			this.beforeRemoveselect = 
				this.beforeSelectpage = 
					this.beforeAddtag =
						this.beforeAddbasket =
							this.beforeRestoreselected = 
								this.beforeRemovebasket = 
									this.beforeApprove = 
										this.beforeDisapprove =
											this.beforeUnapprove = 
												this.beforeClearselection
												= function()
												{
		this.waitingDialog.show();
		return true;
												}

	this.doClone = function()
	{
		var id = (this.currentRecordId) ? this.currentRecordId : "-1";
		var win = openTab("Edit cloned record", webRoot + "eprecord/edit.do?render-mode=popup&id=-1&clone="+id);
		//var newEntity = clone(this.currentEntity);
		//newEntity.id = -1;
		//this.drawFromJSON(newEntity, {method: 'prepend'});
		win.callback = function(entity) 
		{
			self.request(false, function()
					{
				self.setPositionById(entity["id"])
				win.closeTab();
					}); 
		};
	}

	this.doDiscuss = function()
	{
		var id = (this.currentRecordId) ? this.currentRecordId : "-1";
		if(id != "-1")
			var win = window.open(webRoot + "properties/action.do?exp-id="+id, this.controller+"Discuss","location=0,status=0,scrollbars="+this.options.editwindow_scrollbars+",resizable=yes,width=200,height=200");
	}

	this.doRecordmenu = function(link)
	{
		this.recordMenu.cfg.setProperty("context", [$(link).attr('id'), "tl", "tr"]); 
		this.recordMenu.show();
	}

	this.doPublishmenu = function(link)
	{
		this.publishMenu.cfg.setProperty("context", [$(link).attr('id'), "tl", "tr"]); 
		this.publishMenu.show();
	}

	this.doRecordchanges = function(link)
	{
		var interval = $(link).attr('interval');
		url = webRoot+"epbrowser/show.do?render-mode=popup&interval="+this.currentEntity.id+","+interval;
		window.location = url;
	}

	this.doShowsimilarcompounds = function()
	{
		openTab("Records with similar compounds", webRoot + "epbrowser/show.do?similarmol="+this.currentEntity.molecule.id);
	}

	this.doShowsamecompounds = function()
	{
		openTab("Records with same compounds", webRoot + "epbrowser/show.do?similarmol="+this.currentEntity.molecule.id+"&stereochemistry=1");
	}

	this.doTrackchanges = function()
	{
		var win = openTab("Record history", "history/show.do?render-mode=popup&record="+this.currentEntity.id);
	}

	this.doReset = function()
	{
		window.location.href = webRoot + "epbrowser/show.do?render-mode=popup";
	}

	this.fragSearch = function(fragId)
	{
		if (fragId == undefined)
		{
			$("#moleculeSearch").html(new View({url: 'js/templates/fragsearch.ejs'}).render({id:-2}));
		}	
		else
		{
			var win = openTab("Draw compound fragment", webRoot+"molecule/edit.do?render-mode=popup&id="+fragId);
			win.callback = function(newId)
			{
				if (newId != fragId)
				{
					self.fragBlock = $("#moleculeSearch").html();
					$("#moleculeSearch").html(new View({url: 'js/templates/fragsearch.ejs'}).render({id:newId}));
					self.page = 1;
					self.request();
				}
			}
		}
	}

	this.cadasterSearch = function(fragId)
	{
		if (fragId == undefined)
		{
			$("#cadasterSearch").html(new View({url: 'js/templates/cadastersearch.ejs'}).render({id:-2}));
		}	
		else
		{
			var win = openTab("Draw compound fragment", webRoot+"molecule/edit.do?render-mode=popup&id="+fragId);
			win.callback = function(newId)
			{
				if (newId != fragId)
				{
					self.fragBlock = $("#cadasterSearch").html();
					$("#cadasterSearch").html(new View({url: 'js/templates/cadastersearch.ejs'}).render({id:newId}));
					self.page = 1;
					self.request();
				}
			}
		}
	}


	this.onItemsLoaded = function()
	{
		var checked = $("input[name='duplicates']").attr("checked");
		self.trash = self.filters.getValue("trash") != undefined;
		$(self.scope).find('a[action="addselect"]')
		.setClass("invisible", self.pager.totalNum > currentUserLimit);
		if (self.filters.getValue("interval"))
			$("#modification-time").removeClass("invisible");
		$("#page-icon").attr('src', self.trash ? 'img/icons/trash.png' : 'img/icons/compounds.png');
		$("h1:first").html("Compounds properties browser" + (self.trash ? ' - Trash' : '') + "<a class=\"infolink\" href=\"https://docs.ochem.eu/display/MAN/Compound+property+browser\" target=\"_blank\"></a>");
		$(self.scope).find(".no-trash").setClass("invisible", self.trash);
		$(self.scope).find(".trash").setClass("invisible", !self.trash);

		$("#duplicate-filter").setClass("invisible", !self.filters.isSet("article") && !self.filters.isSet("basket-select"))

		if(checked)
			$("#clean-basket").removeClass("invisible");
		else
			$("#clean-basket").addClass("invisible");

		$("#unmoderated").setClass("invisible", !$("[name='awaiting-approval']").is(":checked"));
		self.updateSelectionSize();

		if (self.pager.totalNum == 0)
			self.mainDiv.html("<div class='content-message'>Oops.. There are no data matching your filters.<br/>Try to remove some filters</br>Or change Data visibility: All data.</div>");

		//self.pager.setVisibility(self.pager.totalNum > 0);

	}

	this.listenEvent("items_loaded", this.onItemsLoaded);

	this.onItemDrawn = function()
	{
		this.currentBlock.find('[title]').tooltip({showURL: false});

		var recStatus = this.currentEntity.ep_status;

		if(recStatus =="0" || recStatus =="3")
			this.currentBlock.addClass("highlighterror");
		if(recStatus =="1")
			this.currentBlock.addClass("highlightnon-v");

		$(self.scope).find(".trash").setClass("invisible", !self.trash);

		if (this.currentEntity.errorComment != null && this.currentEntity.errorComment != undefined)
		{
			var errorComment = this.currentEntity.errorComment;
			if (errorComment.indexOf("Duplicate") == 0)
			{
				//Replace plain display with link.... tricky :)
				var htmlPiece = this.currentBlock.find(".errorcomment").html();
				var re = new RegExp(" R[0-9]+");
				var m = re.exec(errorComment);
				if (m != null)
				{
					htmlPiece = htmlPiece.replace(/(R[0-9]+)/,"<a action='openduplicate'>$1</a>");
					this.currentBlock.find(".errorcomment").html(htmlPiece);
					this.attachActions(this.currentBlock.find(".errorcomment"));
				}
			}
		}
	}

	this.doOpenduplicate = function()
	{
		this.currentBlock.find(".duplicate").html("<img src='img/roller_small.gif'/>");
		this.ajax.send
		(
				{
					url: this.controller+"/getduplicate.do", 
					data: "id="+this.currentEntity.id, 
					success: function(response)
					{
						var item = response[this.itemElement];
						if (item != undefined)
						{
							item.hidepanel = true;
							var rendered = this.draw(item);
							this.currentBlock.find(".duplicate").html(rendered);
						} else
						{
							this.currentBlock.find(".duplicate").html('<div style="text-align: center; margin: 10px 10px 10px 10px; color: #550000;">Record not found. Please notify us about this behaviour. </div>');
						}
					},

					error: function(e)
					{
						this.currentBlock.find(".duplicate").html('<div style="text-align: center; margin: 10px 10px 10px 10px; color: #550000;">Ajax request failed: '+e+'</div>');
					}
				}
		);
	}

	this.doZoom = function(link)
	{
		var size = Math.round(Math.min(window.innerHeight, window.innerWidth) * 0.9);
		var id = this.currentBlock.find('.block-image IMG').attr("id");
		var img = $(this.molViewDialog.body).find("IMG");
		img.attr("src","");
		img.attr("src","depiction.jsp?w="+size+"&h="+size+"&id="+id);
		img.attr("width",size);
		img.attr("height",size);
		this.molViewDialog.render();
		this.molViewDialog.show();
	}

	this.onRestoreSuccess = function()
	{
		this.deleteRecord();
	}

	this.onRestoreError = function(msg, entity){
		if (entity.attachment)
		{
			$("input[name='rec-id']").val(entity.attachment.content);
			$("#dupRecords").html(msg);
			this.restoreDialog.show();
		}
		else
			window.alert(msg);
	}

	this.onPublishSuccess = this.onPublish_freely_availableSuccess = function()
	{
		this.request();
	}

	//this.beforeDelete = function()
	//{
	//	if (this.filters.getValue("trash"))
	//		return window.confirm('Do you really want to delete given record from the trash?\nIts a permanent change and cannot be undone');
	//	else
	//		return window.confirm('Do you really want to move record to trash?\nThe item will go to the trash of user, who introduced this record');
	//}

	this.doSelectbasket = function()
	{
		var propWin = openTab("Select basket", webRoot + "basket/show.do?render-mode=popup");
		propWin.moveBy(50, 50);
		propWin.callback = function(basket)
		{
			self.filters.setValue("basket", basket.name);
			$("a.mBasket").html(basket.name);
			self.request(true, function(){propWin.closeTab();});
		}
	}

	this.doShowdublicates = function()
	{
		window.openTab("Duplicates", webRoot + "epbrowser/show.do?similarto="+this.currentRecordId);
	}

	this.doEditproperty = function()
	{
		var propWin = openTab("Select property", webRoot + "properties/show.do?render-mode=popup");
		propWin.moveBy(50, 50);
		propWin.callback = function(newProperty)
		{
			self.filters.setValue('property', newProperty.id, newProperty.name);
			//Filter.property.evtChange.fire();
			//Filter.unit.update("property="+newProperty.id);
			self.request(true, function(){propWin.closeTab();});
		}
	}

	this.doEditarticle = function()
	{
		var articleWin = openTab("Select article", webRoot + "article/show.do?render-mode=popup");
		articleWin.moveBy(50, 50);
		articleWin.callback = function(newArticle)
		{
			self.filters.setValue('article', newArticle.id, newArticle.name);
			self.request(true, function(){articleWin.closeTab();});
		}
	}

	this.doHidefilters = function()
	{
		$("#EPFilters").addClass("invisible");
		$(".show-button").removeClass("invisible");
	}

	this.doShowfilters = function()
	{
		$("#EPFilters").removeClass("invisible");
		$(".show-button").addClass("invisible");
	}

	this.parentSetPosition = this.setPosition;
	this.setPosition = function(link)
	{
		// Context menu is outside of the record context
		// So dont change position while clicking menu item link
		// Otherwise you lose focus
		if (!$(link).hasClass("yuimenuitemlabel"))
			this.parentSetPosition(link);
	}

	this.parentInitialize = this.initialize;

	this.propertyACchanged = function()
	{
		//Filter.unit.update("property="+Filter.property.getKey());
	}

	this.structureBasketCbChanged = function(eventObject)
	{
		var div = $("#structure-basket-div");
		if ($(eventObject.target).is(':checked'))
		{
			div.removeClass("invisible");
			self.filters.setValue('structure-basket-dofilter', 1, "");
			if (self.filters.getValue('structure-basket') == -1)
				self.doSelectstructurebasket();
			else
				self.request();
		} else
		{
			div.addClass("invisible");
			self.filters.setValue('structure-basket-dofilter', 0, "");
			self.request();
		}
	}

	this.doSelectstructurebasket = function()
	{
		var structurebasketWin = openTab("Select structures basket", webRoot + "basket/show.do?render-mode=popup");
		structurebasketWin.callback = function(newBasket)
		{
			self.filters.setValue('structure-basket', newBasket.id, newBasket.name);
			$("[name=selectstructurebasket-link]").html(newBasket.name);
			structurebasketWin.closeTab();
			self.request();
		}
	}


	//this.onUnitupdate = function()
	//{
	//	var elSel = Filter.unit.select;
	//	var elOptOld = elSel.options[0];
	//	var elOptNew1 = new Option("No unit filter", -1, false, false);
	//   try 
	//  {
	//   elSel.add(elOptNew1, elOptOld); // standards compliant; doesn't work in IE
	//}
	//catch(ex) 
	//{
	//  elSel.add(elOptNew1, 0); // IE only
	//}	 
	//Filter.unit.select.selectedIndex = 0;
	//}

	this.initialize = function()
	{
		//basket limit increased from 50 to 200
		$("#structure-basket-cb").change(this.structureBasketCbChanged);

		this.basketFilter.update("pagenum=1&pagesize=500&showsystem=1",null,"basketlist_update");	
		this.parentInitialize(true);
		Filter = new Object();
		Filter.property = new Autocomplete("property", dsPropertiesAC);
		Filter.property.evtChange.fire = this.propertyACchanged;
		$('#ac-property').append(Filter.property.getJQuery());
		Filter.property.attachHandlers();

		Filter.source = new Autocomplete("article", dsSourcesAC);
		$('#ac-source').append(Filter.source.getJQuery());
		Filter.source.attachHandlers();

		if (getParams["hidefilters"] == "true" || isMobile)
			$("#EPFilters").addClass("invisible");

		if (isMobile)
			$("#show-filters").removeClass("invisible");

		//Filter.unit = new UnitSelect('unit', 'unit/list.do', 'unit');
		//Filter.unit.scopeBlock = this.currentBlock;
		//Filter.unit.listenEvent("update", this.onUnitupdate);

		$("[name=visibility]").change(function(){
			if ($(this).val() == "private")
				$("[name=approval-status]").attr("disabled", true);
			else
				$("[name=approval-status]").removeAttr("disabled");

		});

		var dlgOptions;

		this.basketDialog = new YAHOO.widget.Dialog("basketDialog", { 
			width:"600px", 
			fixedcenter:true, 
			modal:true, 
			visible:false,
			postmethod: "none",
			buttons:[{text:"Submit", handler:this.handleDialogSubmit, isDefault:true},{text:"Cancel", handler:this.handleDialogCancel}]
		});	
		this.basketDialog.manualSubmitEvent.subscribe(this.handleDialogSubmit);

		this.cleanBasketDialog = new YAHOO.widget.Dialog("cleanBasketDialog", { 
			width:"600px", 
			fixedcenter:true, 
			modal:true, 
			visible:false,
			postmethod: "none",
			buttons:[{text:"Submit", handler:this.cleanDialogSubmit, isDefault:true},{text:"Cancel", handler:this.handleDialogCancel}]
		});
		this.cleanBasketDialog.manualSubmitEvent.subscribe(this.cleanDialogSubmit);

		this.waitingDialog = new YAHOO.widget.Dialog("waitingDialog", { 
			width:"325px", 
			fixedcenter:true, 
			modal:true, 
			visible:false, 
			close: false,
			postmethod: "none"
		});

		this.restoreDialog = new YAHOO.widget.Dialog("restoreDialog", { 
			width:"525px", 
			fixedcenter:true, 
			modal:true, 
			visible:false,
			postmethod: "none",
			buttons:[{text:"Show record", handler:handleSubmit, isDefault:true},{text:"Cancel", handler:handleCancel}]
		});

		this.molViewDialog = new YAHOO.widget.Dialog("molViewDialog", { 
			fixedcenter:true, 
			modal:true, 
			visible:false,
			postmethod: "none" 
		});

		this.restoreDialog.render();
		this.waitingDialog.render();
		this.basketDialog.render();	
		this.cleanBasketDialog.render();
		this.molViewDialog.render();

		this.listenEvent("basket_ready", self.updateBasketFilter);

	}

	this.molViewDialogClose = function()
	{
		this.molViewDialog.cancel();
	}

	this.doCleanbasketaction = function(link)
	{
		this.currentAction = $(link).attr("type");
		this.cleanBasketDialog.show();
	}

	this.handleDialogCancel = function() 
	{
		this.cancel();
	} 

	this.handleDialogSubmit = function() 
	{ 
		self.callAction(self.currentAction, null, null);
	}

	this.cleanDialogSubmit = function() {

		self.callAction(self.currentAction, null, null);
	};

	var handleSubmit = function() { 
		var dupRec = $("input[name='rec-id']").val();
		var recordWin = openTab("Duplicate record", webRoot + "epbrowser/show.do?render-mode=popup&id="+dupRec);
		this.cancel();
	};

	var handleCancel = function() { 
		this.cancel();
	};

	this.onSelectupdate = function()
	{	
		self.basketSelect.select.options[self.basketSelect.select.options.length]
		= new Option("Create new set...", -10, false, false);
		self.onSelectchange();
		////
		var selIndex = self.basketFilter.select.selectedIndex;
		var selOption = self.basketFilter.select.options[selIndex];

		for (i=0; i<self.basketSelect.select.options.length; i++)
		{
			if (self.basketSelect.select.options[i].value == selOption.value)
			{
				self.basketSelect.select.selectedIndex = i;
				break;
			}
		}

	}

	this.basketSelect.listenEvent("update", this.onSelectupdate);

	this.onSelectchange = function()
	{	
		if (self.basketSelect.select.value == -10)
		{
			$("#basketDialog").find("#basketName").show();
			$("#basketDialog").find("#basketOption").hide();

		} else
		{
			$("#basketDialog").find("#basketName").hide();
			$("#basketDialog").find("#basketOption").show();
		}
	}

	this.doBasketaction = function(link)
	{
		this.currentAction = $(link).attr("type");

		this.basketSelect.select.options.length = 1;
		this.basketSelect.select.options[0] = new Option("Loading items...", -11, false, false);

		//molecule set limit has been increased to 500
		this.basketSelect.update("pagenum=1&pagesize=500");
		$("#basketDialog").find("#basketName").attr("style", "display:none;");
		$("#basketDialog").find("#basketNameEdit").val("");
		this.basketDialog.show();
	}

	this.getActionQuery = function(name)
	{
		var data = this.parentActionQuery();

		if (name == "addbasket" || name == "removebasket")
		{
			var optionIndex = self.basketSelect.select.selectedIndex;
			var option = self.basketSelect.select[optionIndex];

			if (option.value == -10)
				data += "&basket-name=" + URLEncode(self.basketDialog.getData()["basket-name"]);
			else
			{
				data += "&basket-name=" + URLEncode(option.text);
				var mol_option = self.basketDialog.getData()["mol-option"];
				data+= "&mol-option="+mol_option;
			}
		}
		else if (name == "selectpage")
		{
			data = data.replace(/id=([0-9]\,?)*/g, ""); // Discard the old ID filter, if there was any
			// Page-targeted action. Probably move this to browser.js, since it looks pretty universal / Midnighter
			$("div[rec-id]").each(function(){
				data += "&id=" + $(this).attr('rec-id');
			});
		}
		else if (name == "cleanbasket")
		{
			data += "&withValue="+self.cleanBasketDialog.getData()["withValue"];
			data += "&withStreo="+self.cleanBasketDialog.getData()["withStreo"];
		}
		return data;
	}

	this.beforeDeleteselected = function()
	{
		if (!this.filters.getValue("trash") || window.confirm('Do you really want to delete all selected records?\nIts a permanent change and cannot be undone'))
		{
			this.waitingDialog.show();
			return true;
		}
	}

	this.onCleanbasketSuccess = function()
	{
		this.waitingDialog.cancel();
		this.cleanBasketDialog.hide();
		self.request();
	}


	this.onApproveSuccess = function()
	{
		this.waitingDialog.cancel();
		var previewWin = openTab("Approve records preview", webRoot + "editor/approvePreview.do");
		previewWin.callback = function()
		{
			self.request();
		}
	}

	this.onDisapproveSuccess = function()
	{
		this.waitingDialog.cancel();
		var previewWin = openTab("Reject records preview", webRoot + "editor/rejectPreview.do");
		previewWin.callback = function()
		{
			self.request();
		}
	}

	this.onUnapproveSuccess = function()
	{
		this.waitingDialog.cancel();
		var previewWin = openTab("Unapprove records preview", webRoot + "editor/unapprovePreview.do");
		previewWin.callback = function()
		{
			self.request();
		}
	}

	this.onRestoreselectedSuccess = 
		this.onAddtagSuccess = 
			this.onRemovetagSuccess =
				this.onAddselectSuccess = 
					this.onRemoveselectSuccess = 
						this.onSelectpageSuccess =
							this.onClearselectionSuccess = 
								function()
								{
		// Unfreeze interface and reload page
		this.waitingDialog.cancel();
		self.request();
								}

	this.onDeleteselectedSuccess = this.onDeleteSuccess = function(response)
	{
		this.filters.setValue("deletion-confirmation", "");
		this.waitingDialog.cancel();
		if (response.message)
			// Ask user confirmation bzw. reason for deletion
			this.overviewDeletion(response.message);
		else
			self.request();	
	}

	this.overviewDeletion = function(recordsSummary)
	{
		if (recordsSummary.context == "provide-message")
		{
			var message = window.prompt(recordsSummary.message);

			if (message)
			{
				message = message.trim();
				this.filters.setValue("deletion-confirmation", message);
			}
		}
		else
			if (window.confirm(recordsSummary.message))
				this.filters.setValue("deletion-confirmation", "OK");

		if (this.filters.getValue("deletion-confirmation"))
			this.callAction("deleteselected");
	}

	this.onBatcheditError =
		this.onCleanbasketError =
			this.onAddtagError = 
				this.onRemovetagError = 
					this.onDeleteselectedError = 
						this.onDeleteError =
							this.onRemoveselectError = 
								this.onSelectpageError = 
									this.onAddbasketError =
										this.onRemovebasketError =
											this.onAddselectError =
												this.onClearselectionError = 
													function(msg, msgEntity)
													{
		this.filters.setValue("deletion-confirmation", "");
		// On error, unfreeze interface and alert this error
		this.waitingDialog.cancel();
		this.ajax.showError(msgEntity);
													}

	this.doSelecttag = function(link)
	{
		var suffix = $(link).attr('removetag') ? "remove" : "add"; 
		var win = openTab("Select tag to " + suffix, webRoot + "tags/show.do?render-mode=popup&type=molecule&force-type=true");
		win.callback = function(tag)
		{
			self.filters.setValue("tagid", tag.id);
			self.callAction(suffix + "tag");
			win.closeTab();			
		}
	}

	this.doBatchedit = function()
	{
		var win = openTab("Batch editing", webRoot + "batchedit/show.do?apply-basket=true&" + this.getActionQuery());
	}

	this.doCompounddetails = function()
	{
		var win = openTab("Molecule profile", webRoot + "molecule/profile.do?depiction="+this.currentEntity.molecule.id);
	}

	this.onAddbasketSuccess =
		this.onRemovebasketSuccess = function()
		{
		this.waitingDialog.cancel();
		this.basketDialog.hide();
		this.fireEvent("basket_ready");
		}

	this.filters.onFilterChange = function(filter)
	{
		if ($(filter).attr("name") == "duplicates")
		{
			var checked = $(filter).attr("checked");
			if (checked)
			{
				$("input[name='nostereo']").attr("disabled", false);
				$("select[name='sortby']").get(0).selectedIndex = 3;
			}
			else
			{
				$("input[name='nostereo']").attr("disabled", true);			
			}
		}
		if ($(filter).attr("name") == "nameerror")
		{
			var checked = $(filter).attr("checked");
			if (checked)
			{
				$("input[name='yesstereo']").attr("disabled", false);
			}
			else
			{
				$("input[name='yesstereo']").attr("disabled", true);			
			}
		}


		self.request(true);
	}

	this.doExpanaloguedetails = function()
	{
		openTab("Experimental measurement", "epbrowser/show.do?id=" + this.currentEntity["exp-analogue"].id);
	}


}

function ConditionsFilter()
{
	this.controller = "properties";
	this.scope = ".conditionsscope";

	UnifiedBrowser.call(this);
	this.itemTemplate = "js/templates/conditionvalue-filter.ejs";
	this.container = "ConditionsFilter";
	this.restriction = "conditionsfilter";
	this.minusCounter = -1;

	var self = this;

	this.doNewfilter = function()
	{
		var propWin = openTab("Select condition", webRoot + "properties/show.do?render-mode=popup&condition=true");
		propWin.callback = function(newProperty)
		{
			self.drawFromJSON(newProperty);
			propWin.closeTab();

			if (newProperty.type == 1)
				setTimeout('conditionsFilter.currentBlock.get(0).optionsSelect.update("id='+newProperty.id+'")', 1);
			else if (newProperty.type == 0)
				setTimeout('conditionsFilter.currentBlock.get(0).unitsSelect.update("category='+newProperty.unitCategory.id+'")', 1);
			self.updateVisibility();
		}	
	}

	this.onItemDrawn = function()
	{
		var dom = this.currentBlock.get(0);
		var ent = this.currentEntity;
		// Create dynamic select for units of condition
		var dynSelect = new UnitSelect('cond-unit'+ent.id, 'unit/list.do?lightweight=1', 'unit');
		dynSelect.scopeBlock = this.currentBlock;
		dom.unitsSelect = dynSelect;

		// Create dynamic select for options for this condition ("qualitive conditions")
		var optSelect = new DynamicSelect('cond-option'+ent.id, 'properties/listoptions.do', 'option');
		optSelect.scopeBlock = this.currentBlock;
		dom.optionsSelect = optSelect;

		// Show/hide appropriate sections
		//dynSelect.fireEvent("update");
		this.updateVisibility();

	}

	this.updateVisibility = function()
	{
		this.currentBlock.find('span').addClass("invisible");
		this.currentBlock.find('span[name="type-'+this.currentEntity.type+'"]').removeClass("invisible");
	}

	this.doDelete = function()
	{
		this.deleteRecord();
	}

}

