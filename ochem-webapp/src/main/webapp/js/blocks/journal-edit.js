
function JournalForm()
{
	this.scope = ".formscope";	
	EditForm.call(this);
	this.itemElement="journal";
	this.actionURL="journal/action.do";
	
	this.beforeEdit = function()
	{
		if(this.getValue('n-title').length < 6)
		{
			window.alert("Journal name is too short. please add correct journal name");
			return false;
		}
		var count = this.getValue('n-issn').split(/-/g).length - 1;
		if(count > 1)
		{
			window.alert("ISSN number cannot have more than one '-'. please provide correct ISSN number.")
			return false;
		}
		if(this.getValue('n-issn').length != 9 || this.getValue('n-issn').match('-') == null)
			return window.confirm("ISSN is missing or wrong format. Please, try to find it in order to avoid duplication of this journal in the database");
		return true;
	}
	
	this.beforeReload = function()
	{
		var count = this.getValue('n-issn').split(/-/g).length - 1;
		if (count > 1)
		{
			window.alert("ISSN number cannot have more than one '-'. please provide correct ISSN number.")
			return false;
		}
		
		if (this.getValue('n-issn').length != 9 || this.getValue('n-issn').match('-') == null)
		{
			window.alert("ISSN is missing or wrong format. Please provide correct ISSN number")
			return false;
		}
		return true;
	}
	
	this.onReloadSuccess = function(xml)
	{
		this.entity = xml;
		$('input[name="n-title"]').val(xml.journal.title);
		$('input[name="n-abbreviation"]').val(xml.journal.abbreviation);
		$('input[name="n-publisher"]').val(xml.journal.publisher);
		$('input[name="n-link"]').val(xml.journal.link);
	}
	
	this.onReloadError = function(xml)
	{
		if(window.confirm("Journal is not found in PubMed journal database. Please, verify the ISSN number and try again."))
		{
			$("input[n-disable]").removeAttr("disabled");
			$("div[n-hidden]").removeClass("hidden");
		}
	}
	
}

include.plugins('view');
var journalForm = new JournalForm();
$(document).ready(
	function()
	{
		if(journalForm.getValue('n-issn').length == 9)
		{
			$(document).find("[n-disable]").attr("disabled", "disabled");
			$("div[n-hidden]").addClass("hidden");
		}
	}
);