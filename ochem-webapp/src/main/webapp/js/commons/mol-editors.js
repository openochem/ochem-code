/**
 * MarvinJS abstraction
 * @author midnighter/itetko
 * @returns
 */


function JSME()
{
	var self = this;
	var dfd = new jQuery.Deferred();
	var timer = 0;

	this.setMolecule = function(mol){

		//window.alert(mol)

		if (!self.sketcher){
			self.sketcher = document.getElementById("sketch").contentWindow.jsmeApplet;
		}

		try
		{
			if (self.sketcher)
			{
				self.sketcher.readGenericMolecularInput(mol);
			}
			else
			{
				QSPR._this = this;
				QSPR._mol = mol;
				timer = setTimeout('QSPR._this.setMolecule(QSPR._mol)', 500);
			}
		}
		catch (e)
		{
			QSPR._this = this;
			QSPR._mol = mol;
			timer = setTimeout('QSPR._this.setMolecule(QSPR._mol)', 500);	
		}
	}

	this.getMolecule = function()
	{
		if (!self.sketcher){
			self.sketcher = document.getElementById("sketch").contentWindow.jsmeApplet;
		}

		return self.sketcher.molFile();
	}

	this.promise = function()
	{
		return dfd.promise();
	}

}
