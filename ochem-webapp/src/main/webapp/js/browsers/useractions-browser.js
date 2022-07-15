

function loadUserEvents()
{
	var ajax = new QSPR.Ajax();
	ajax.send({
		url: "useractions/list.do",
		data: getParams["login"] ? "login=" + getParams["login"] : "",
		success: function(response){
			var items = response.list["user-event"];
			var block = $("#Browser");
			var user = "";
			var day = "";
			for (var i = 0; i < items.length; i++)
			{
				var item = items[i];
				
				if (item.date != day)
				{
					block.append("<div class='day'>" + item.date + "</div>");
					day = item.date;
					user = "";
				}
				
				if (item.username != user)
				{
					if (item.username.indexOf("Guest") == 0)
						block.append("<div class='user'><img src='img/icons/user-16.png'/>" + item.username + "</div>");
					else
						block.append("<div class='user'><img src='img/icons/user-16.png'/><a href='user/profile.do?login=" + item.username + "' tab='User profile'>" + item.username + "</a></div>");
					user = item.username;
				}
				
				var comment = $("<div class='comment'><span class='time'>" + item.time + "</span><span class='comment-content'></span></div>");
				comment.find(".comment-content").html(item.comment);
				block.append(comment);
			}
			
			if (getParams["login"])
				$("#info").html("Activity of user " + getParams["login"] + ".");
			
			$(document).trigger('DOM_updated');
		}
	});
}