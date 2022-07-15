
function AuthorForm()
{
	this.scope = ".formscope";	
	EditForm.call(this);
	this.itemElement="author";
	this.actionURL="author/action.do";
}

include.plugins('view');
var authorForm = new AuthorForm();
