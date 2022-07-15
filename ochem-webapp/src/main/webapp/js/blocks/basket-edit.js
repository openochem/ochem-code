var ajaxForm = new EditForm();

ajaxForm.scope = ".EditForm";
ajaxForm.actionURL = "basket/action.do";
ajaxForm.itemElement = 'basket';
ajaxForm.restriction = "form";
ajaxForm.beforeEdit = function()
{	
	var basketname = $("input[name='name']").attr("value");

	if(basketname.length < 2)
	{
		window.alert("basket name is too short can not be added");
		return false;
	}
	return true;
}

function checkBasketName(obj)
{
	if(obj.value.length <= 1)
	{
		$("#basketName").css("color","red").html("error: basket name length is "+obj.value.length+" and can not be saved (min. 2 characters)");
	}else
	if (obj.value.length > 1 && obj.value.length < 5)
	{
		$("#basketName").css("color","#FA58F4").html("warning: basket name length is "+obj.value.length+" and may not be informative");
	}
	else if(obj.value.length >= 5 && obj.value.length < 40)
	{
		$("#basketName").css("color","green").html("basket name is good and length is "+obj.value.length);
	}
	else if(obj.value.length >= 40)
	{
		$("#basketName").css("color","red").html("error: basket name length is "+obj.value.length+", its should be less than 40 characters and can not be saved (max. 40 characters)");
	}
}
