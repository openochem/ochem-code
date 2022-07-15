function User()
{
	AjaxForm.call(this);
	this.controller = "user";
	this.actionURL = "user/action.do";
	var self = this;

	this.doChange = function()
	{
		$("#pwd").parent().html("" +
				"<tr><td width=\"150px\">Password*</td><td><input type=\"password\" name=\"passwd\" send=\"1\" class=\"required\"/></td></tr>" +
				"<tr><td width=\"150px\">Confirm password*</td><td><input type=\"password\" name=\"c_passwd\" send=\"1\" class=\"required\"/></td></tr>"
		);
	}

	this.beforeCheck = function()
	{
		var login = this.getValue('login').trim();

		if (login.length < 3 || login.length > 20)
		{
			$(".message").css("color", "red");
			$("[name=warn_login]").text("Please, enter a valid login!");
			return false;
		}
		return true;
	}

	this.onUpdateSuccess = function(xml)
	{
		$(".redirect").addClass("invisible");
		self.info("Your profile has been successfully updated!");
	}

	this.onCheckSuccess = function(xml)
	{
		$(".message").css("color","red").html(this.getValue('login')+" is not available");
	}

	this.onCheckError = function(xml)
	{
		$(".message").css("color","green").html(this.getValue('login')+" is available");
		$("[name=warn_login]").text("");
	}

	this.warning = function(str)
	{
		$(".warning").removeClass("invisible").css("color", "red").html(str);
	}

	this.info = function(str)
	{
		$(".warning").removeClass("invisible").css("color", "green").html(str);
		$('body,html').animate({ 
			scrollTop: 0
		}, 800);
	}

	this.beforeSubmit =  this.beforeUpdate = function()
	{	
		var login = this.getValue('login').trim();
		var email = this.getValue('emailId').trim();		
		var pwd   = this.getValue('passwd');
		var c_pwd = this.getValue('c_passwd');
		var title   = this.getValue('title');
		var fname = this.getValue('firstname');
		var lname = this.getValue('lastname');
		var orga   = this.getValue('organisation');

		var ok = true;

		// login
		if (login.length < 3 || login.length > 20)
		{
			$(".message").css("color", "red");
			$("[name=warn_login]").text("Please, enter a valid login (between 4 and 20 characters)!");
			ok = false;
		} else
		{
			$("[name=warn_login]").text();
		}

		// email
		if (email == "" || ! echeck(email)) 
		{
			$("[name=warn_email]").text("Please, enter a valid e-mail!");
			ok = false;
		}
		else
		{
			$("[name=warn_email]").text("");
		}

		// passwd
		if (this.getValue('userId') == "") 
		{
			if (pwd == "") {
				$("[name=warn_passwd]").text("Password cannot be empty!");
				ok = false;
			}
			if (pwd.replace(/[A-Z]|[a-z]|[0-9]/g, "") != "") {
				$("[name=warn_passwd]").text("Password can contain only letters and numbers!");
				ok = false;
			}
			if (pwd != c_pwd) {
				$("[name=warn_c_passwd]").text("Passwords do not match");
				ok = false;
			}
		}

		// first name
		if (fname == "") 
		{
			$("[name=warn_firstname]").text("Please, enter your first name");
			ok = false;
		}
		else
		{
			$("[name=warn_firstname]").text("");
		}

		// last name
		if (lname == "") 
		{
			$("[name=warn_lastname]").text("Please, enter your last name");
			ok = false;
		}		
		else
		{
			$("[name=warn_lastname]").text("");
		}

		// organisation
		if (orga == "-- please select --")
		{
			$("[name=warn_orga]").text("Please, select your organization type.");
			ok = false;
		}
		else
		{
			$("[name=warn_orga]").text("");
		}

		// title
		if (orga == "-- please select --")
		{
			$("[name=warn_title]").text("Please, select your title.");
			ok = false;
		}
		else
		{
			$("[name=warn_title]").text("");
		}
		
		return ok;
	}	

	function echeck(str) 
	{
		var at="@"
			var dot="."
				var lat=str.indexOf(at)
				var lstr=str.length
				var ldot=str.indexOf(dot)

				if (
						(str.indexOf(at)==-1) || 
						(str.indexOf(at)==-1 || str.indexOf(at)==0 || str.indexOf(at)==lstr) || 
						(str.indexOf(dot)==-1 || str.indexOf(dot)==0 || str.indexOf(dot)==lstr) || 
						(str.indexOf(at,(lat+1))!=-1) || 
						(str.substring(lat-1,lat)==dot || str.substring(lat+1,lat+2)==dot) || 
						(str.indexOf(dot,(lat+2))==-1) || 
						(str.indexOf(" ")!=-1)
				)
				{
					return false;
				}
		return true;
	}

	this.onSubmitSuccess = function(xml)
	{
		$(".warn").text("");
		window.parent.location.replace(webRoot + "home/show.do?render-mode=full");
	}

}

var user = new User();
$(document).ready(
		function()
		{
			if(user.getValue('userId') != "")
			{
				$("input[n-disable]").attr("disabled", "disabled");
				$("#check").hide();
				$(".message").hide();
			}

			if(user.getValue('organisation') == "-- please select --")
			{
				$("[name=warn_orga]").text("Please, select your organization type.");
			}

			if(user.getValue('title') == "-- please select --")
			{
				$("[name=warn_title]").text("Please, select your title.");
			}
		}
)
