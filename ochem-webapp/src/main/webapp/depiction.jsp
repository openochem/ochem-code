<%@page import="qspr.depiction.*,java.util.*,java.io.*,com.eadmet.mmpa.domain.MMPFragment,qspr.*,org.hibernate.*,qspr.util.*,qspr.workflow.utils.QSPRConstants,qspr.dao.Various,com.eadmet.mmpa.MMPDepiction,com.eadmet.mmpa.MMPFragDepiction"%><%
	try 
	{
		Globals.startAllTransactions();
		
		// initialize the main depiction instance
		MoleculeDepiction depiction = MoleculeDepiction.get(OCHEMConfiguration.getCheminfEngine());
		depiction.configure(request);
	    
	    byte[] imgBinary = null;
	    if (request.getParameter("smiles") != null) 
	    { 
	    	String req = request.getParameter("smiles");
	    	
	    	if (req.equals(""))
	    	{
	    		response.sendRedirect("img/empty-mol.gif");
	    	}

			depiction.appendToStringConfig(",setcolors:a1:#990099");
			depiction.setMolecule(req);
			System.out.println("smiles:\t" + req);
	    }
	    else if ((request.getParameter("id") != null) || (request.getParameter("mp2") != null) || (request.getParameter("mp1") != null))
	    {
  	    	Long depictionId = null;
	    	Integer mp2id = null;
	    	
	    	if (request.getParameter("id") != null)
	    	{
	    		String req = request.getParameter("id");
		    	if (req.equals(""))
		    		response.sendRedirect("img/empty-mol.gif");
		    	
		    	if (req.startsWith("M"))
		    		mp2id = Integer.valueOf(req.substring(1));
		    	else if (req.startsWith("D"))
		    		depictionId = Long.valueOf(req.substring(1));
		    	else
		    		depictionId = Long.valueOf(req);
		    	
		    	if ((depictionId != null && depictionId < 0) || (mp2id != null && mp2id < 0))
		    		response.sendRedirect("img/empty-mol.gif");
		    		
	    	} else if (request.getParameter("mp2") != null)
	    	{
	    		String req = request.getParameter("mp2");
	    		if (req.equals("") || (new Long(req) <= 0))
		    		response.sendRedirect("img/empty-mol.gif");
	    		mp2id = Integer.valueOf(req);
	    	}
	    	
	    	if (depictionId != null) {
	    		qspr.entities.Molecule mol_entity = (qspr.entities.Molecule) Globals.session().get(qspr.entities.Molecule.class, depictionId);
	    		depiction.setMolecule(mol_entity.getData());
	    		if(mol_entity.mapping1.inchi1.length() != 14) {
		    		depiction.setError(true);
		    	}	
	    	}
			else if (mp2id != null)
	    	{
		    	qspr.entities.Mapping2 mp2 = (qspr.entities.Mapping2) Globals.session().get(qspr.entities.Mapping2.class, mp2id);
		    	depiction.setMolecule(mp2.getMolecule().getData());
		    	if(mp2.mapping1.inchi1.length() != 14) {
		    		depiction.setError(true);
		    	}
		    	
		    	//Maximal common subgraph depiction piece
		    	if (request.getParameter("mp2templ") != null || request.getParameter("templ") != null)
		    	{
		    		qspr.entities.Molecule  molTempl;
		    		if (request.getParameter("mp2templ") != null)
		    		{
			    		qspr.entities.Mapping2 mp2teml = (qspr.entities.Mapping2) Globals.session().get(qspr.entities.Mapping2.class,  Integer.valueOf(request.getParameter("mp2templ")));
			    		molTempl = mp2teml.getMolecule();
		    		} else 
		    		{
		    			molTempl = MoleculePeer.fetchFromString(request.getParameter("templ"), null);
		    		}

					MMPDepiction depictionMMP = MMPDepiction.get(
							OCHEMConfiguration.getCheminfEngine(),
							mp2.getMolecule().getData(),
							molTempl.getData()
					);
					depictionMMP.configure(depiction);
					depictionMMP.search();
					imgBinary = depictionMMP.getImage();
		    	}
	    	} else if (request.getParameter("mp1") != null)
	    	{
	    		String req = request.getParameter("mp1");
	    		if (req.equals("") || (new Long(req) <= 0))
		    		response.sendRedirect("img/empty-mol.gif");
	    		
		    	qspr.entities.Mapping1 mp1 = (qspr.entities.Mapping1) Globals.session().get(qspr.entities.Mapping1.class, new Long(req));
		    	qspr.entities.Molecule mp2_mol = mp1.mapping2.get(0).getMolecule();
		    	depiction.setMolecule(mp2_mol.getData());
		    	if(mp2_mol.mapping1.inchi1.length() != 14) {
		    		depiction.setError(true);
		    	}
	    	}
	    } 
// 	    else if (request.getParameter("ecid") != null)
// 	    {
// 	    	String req = request.getParameter("ecid");
// 	    	qspr.ecentities.ECMolecule ecMol = (qspr.ecentities.ECMolecule) Globals.session().get(qspr.ecentities.ECMolecule.class, new Integer(req));	    	
// 	    	if (!customFormat)
//     			imgBinary = ecMol.image;
    	
// 			format += ",#" + validColor;
// 			data = ecMol.smiles;
// 	    } 
	    else if (request.getParameter("mmp_frag") != null) 
	    {
	    	MMPFragment frag = MMPFragment.getFragment(Long.valueOf(request.getParameter("mmp_frag")));
	    	MMPFragDepiction depictionMMPFrag = MMPFragDepiction.get(OCHEMConfiguration.getCheminfEngine());
	    	depictionMMPFrag.configure(depiction);
	    	
	    	imgBinary = depictionMMPFrag.getImage(frag.smiles);
	    }
	    else if (request.getParameter("frag_id") != null)
	    {
	    	String req = request.getParameter("frag_id");
	    	if (req.equals("") || (new Long(req) <= 0))
	    		response.sendRedirect("img/empty-mol.gif");
	    	
	    	qspr.fragmententities.Fragment fragment = (qspr.fragmententities.Fragment) Globals.alternateSession().get(qspr.fragmententities.Fragment.class, new Long(req));
			depiction.setMolecule(fragment.fragment_data);
	    }
// 	    else if (request.getParameter("rgroup") != null)
// 	    {
// 	    	String req = request.getParameter("rgroup");
// 	    	if (req.equals("") || (new Long(req) <= 0))
// 	    		response.sendRedirect("img/empty-mol.gif");
	    	
// 	    	qspr.entities.RGroup rgroup = (qspr.entities.RGroup)Globals.session().get(qspr.entities.RGroup.class, Long.valueOf(req));
// 	    	imgBinary = rgroup.getImage("png:h150,w150");
// 	    }
	    else if (request.getParameter("mol") != null)
	    {
	    	
	    	String req = request.getParameter("mol");
	    	if (req.equals(""))
	    		response.sendRedirect("img/empty-mol.gif"); 
	    	
			String[] inchiKeys = qspr.entities.Molecule.getInChiKeys(Various.molecule.convertToFormat(req, QSPRConstants.SDF));
			if (inchiKeys[0].length() != 14)
				depiction.setError(true);
			depiction.setMolecule(req);
	    }
	    
	    if (imgBinary == null)
		{
			if (depiction.getMolecule() == null)
	    		response.sendRedirect("img/empty-mol.gif"); 
			
			imgBinary = depiction.getImage();
		}
	    
	 	// Set content type and return
	    response.setContentType("image/" + depiction.getFormat());
	 	response.setHeader("Cache-Control","max-age=3600, must-revalidate"); //Experiment with some caching options here
	    ServletOutputStream outs = response.getOutputStream();
	    outs.write(imgBinary, 0, imgBinary.length);
	    outs.flush();	    
	    outs.close();
	} 
	catch(Exception t) 
	{
	   response.sendRedirect("img/invalid-molecule.png"); 
	   t.printStackTrace();
	}
	finally
	{
		Globals.rollbackAllTransactions();
	}
%>