var ajaxForm = new EditForm();

ajaxForm.scope = ".EditForm";
ajaxForm.actionURL = "tags/action.do";
ajaxForm.itemElement = 'tag';
ajaxForm.restriction = "form";
ajaxForm.beforeEdit = function()
{	
	var propname = $("input[name='name']").attr("value");
	var desc = $("textarea[name='description']").attr("value");
	desc = desc.trim();

	if(propname.length < 2 || propname.length > 40)
	{
		window.alert("Tag name is too short or too long and can not be added");
		return false;
	}
	if(desc.length < 50 && ajaxForm.getValue("isPublic"))
	{
		window.alert("Tag description is too short and can not be added");
		return false;
	}
	return true;
}

ajaxForm.parentItemSaved = ajaxForm.onItemSaved;
ajaxForm.onItemSaved = function(entity)
{
	ajaxForm.entity = entity;
	if($("#file").val())
	{
		ajaxForm.waitingDialog.show();
		$("#uploadform").submit();
		$("#progress").removeClass("invisible");
	}
	else
		ajaxForm.parentItemSaved(entity);
	
}

ajaxForm.iframeLoaded = function()
{
	$("#progress").addClass("invisible");
	var data = $("#uploadframe").contents();
	if ((data.text() != "null")&&(data.text() != ""))
	{
		ajaxForm.waitingDialog.cancel();
		ajaxForm.parentItemSaved(ajaxForm.entity);
	}
}

$(document).ready(function(){
	if (getParams["type"])
		$("select[name=type]").val(getParams["type"]);
	
	$("select[name=type]").change(function(){ 
		$("#upload").setClass("invisible", $(this).val() != "molecule");
	});	
	$("select[name=type]").change();
	
	ajaxForm.waitingDialog = new YAHOO.widget.Dialog("waitingDialog", { 
	    	width:"325px", 
			fixedcenter:true, 
			modal:true, 
			visible:false, 
			close: false
	    });
	ajaxForm.waitingDialog.render();
	
});	

function checkTagName(obj)
{
	if(obj.value.length <= 1)
	{
		$("#tagName").css("color","red").html("error: tag name length is "+obj.value.length+" and can not be saved (min. 2 characters)");
	}else
	if (obj.value.length > 1 && obj.value.length < 5)
	{
		$("#tagName").css("color","#FA58F4").html("warning: tag name length is "+obj.value.length+" and may not be informative");
	}
	else if(obj.value.length >= 5 && obj.value.length < 40)
	{
		$("#tagName").css("color","green").html("tag name is good and length is "+obj.value.length);
	}
	else if(obj.value.length >= 40)
	{
		$("#tagName").css("color","red").html("error: tag name length is "+obj.value.length+", its should be less than 40 characters and can not be saved (max. 40 characters)");
	}
}

function checkTagDesc(obj)
{
	if(obj.value.length <= 50)
	{
		$("#tagDesc").css("color","red").html("error: tag description length is "+obj.value.length+" and can not be saved (min. 50 characters)");
	}else
	if(obj.value.length > 50 && obj.value.length < 200)
	{
		$("#tagDesc").css("color","#FA58F4").html("warning: tag description length is "+obj.value.length+" and may not be informative");
	}
	else if(obj.value.length > 200)
	{
		$("#tagDesc").css("color","green").html("tag description is good");
	}
}