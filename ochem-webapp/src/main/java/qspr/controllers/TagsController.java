/* Copyright (C) 2022 BIGCHEM GmbH <info@bigchem.de>
 *
 * Contact: info@bigchem.de
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License (AGPL)
 * as published by the Free Software Foundation; either version 3.0
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the Affero GNU General Public License for more details.
 *
 * You should have received a copy of the Affero GNU Lesser General Public License
 * along with this program; If not, see <https://www.gnu.org/licenses/>. 
 */

package qspr.controllers;

import java.io.File;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Criteria;
import org.hibernate.HibernateException;
import org.hibernate.criterion.Disjunction;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;
import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.business.WebFilters;
import qspr.dao.Various;
import qspr.entities.Mapping1;
import qspr.entities.Molecule;
import qspr.entities.MoleculeName;
import qspr.entities.Property;
import qspr.entities.Tag;
import qspr.entities.User;
import qspr.export.ExportWriter;
import qspr.export.ExportableColumn;
import qspr.export.ExportableMolecule;
import qspr.export.ExportableSet;
import qspr.export.ExportableSetConfiguration;
import qspr.export.HttpResponseExportWriter;
import qspr.frontend.BrowserModel;
import qspr.frontend.WebList;
import qspr.frontend.WebModel;
import qspr.util.MoleculePeer;

import com.eadmet.exceptions.UserFriendlyException;

@Controller
public class TagsController extends BrowserWrapper
{
	private static transient final Logger logger = LogManager.getLogger(TagsController.class);

	public TagsController()
	{
		sessionRequired = true;
	}

	@Override
	public ModelAndView list(HttpServletRequest req, HttpServletResponse res)
			throws Exception 
	{
		WebList webList = new WebList();
		WebFilters filters = formFilters(req);
		Criteria criteria = Globals.session().createCriteria(Tag.class);


		if (filters.has("query"))
			criteria.add(Restrictions.like("name", "%"+filters.get("query")+"%"));

		if ("property".equals(filters.get("filterby")))
			criteria.createCriteria("properties")
			.add(Restrictions.eq("id", getLongParam("property-id")));
		else if (filters.has("id") && ! "".equals(filters.get("id")))
			criteria.add(Restrictions.eq("id", filters.getLong("id")));

		if (filters.has("type"))
			criteria.add(Restrictions.eq("type", filters.get("type")));

		// Rights restrictions
		Disjunction rights = Restrictions.disjunction();
		rights.add(Restrictions.eq("isPublic", Boolean.TRUE));
		if (Globals.userSession().user != null)
		{
			rights.add(Restrictions.eq("introducer", Globals.userSession().user));
			if (Globals.userSession().user.group != null)
			{
				criteria.createAlias("introducer", "intr");
				rights.add(Restrictions.eq("intr.group", Globals.userSession().user.group));
			}
		}
		criteria.add(rights);

		if ("filters".equals(filters.get("filterby")))
		{
			webList.loadFromSet(Globals.getTaginationFilters(null));
			Iterator<Object> iterator = webList.list.iterator();
			while (iterator.hasNext())
				if (!((Tag)iterator.next()).type.equals(filters.get("type")))
					iterator.remove();
			webList.loadFromList(webList.list); // Update counts
		}
		else
			webList.useEntity(Tag.class).loadDistinctFromCriteria(criteria, getPageNum(), getPageSize(15));

		return new BrowserModel().setFilters(formFilters(req)).setObject(webList).getModelAndView();
	}

	public ModelAndView action(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		if ("edit".equals(getParam("action")))
		{

			User user = Globals.userSession().user;
			if(user == null)
				throw new UserFriendlyException("Guest can not introduce/edit tags. Please register in QSPR site in order to introduce Tags");

			Tag tag;
			if (!assertParam("id") || getLongParam("id") < 0)
				tag = new Tag();
			else
			{
				tag = (Tag) Globals.session().get(Tag.class, Long.valueOf(req.getParameter("id")));
				if (!tag.canEdit())
					throw new UserFriendlyException("You are not allowed to edit this tag");
			}

			if(tag.id != null)
				tag.doCheckRights();

			tag.isPublic = assertParam("isPublic");
			tag.showInBrowser = assertParam("showInBrowser");
			tag.name = req.getParameter("name");
			tag.type = getParam("type");
			tag.description =  req.getParameter("description");
			tag.owner = Globals.userSession().user;
			if(tag.introducer == null)
				tag.introducer = Globals.userSession().user;
			Globals.session().save(tag);
			return new WebModel(tag).getModelAndView();
		}
		else if ("delete".equals(getParam("action")))
		{
			Tag tag = (Tag) Globals.session().get(Tag.class, getLongParam("id"));
			if (tag.id != null)
				tag.doCheckRights();
			Globals.session().createSQLQuery("delete from MoleculeTag where tag_id=:tagId").setParameter("tagId", tag.id).executeUpdate();
			Globals.session().delete(tag);
		}
		else if ("saveall".equals(getParam("action")))
		{
			String[] tagIds = req.getParameterValues("tag-id");
			if ("filters".equals(getParam("filterby")))
			{
				// Set the "area of interest" to the selected tags
				Set<Tag> filterTags = Globals.getTaginationFilters(null);
				Iterator<Tag> iterator = filterTags.iterator();
				while (iterator.hasNext())
					if (iterator.next().type.equals(getParam("type")))
						iterator.remove();
				if (tagIds != null)
				{
					if (getParam("type").equals("molecule"))
						// Only one tag for molecules is allowed
						filterTags.add((Tag) Globals.session().get(Tag.class, Long.valueOf(tagIds[tagIds.length - 1])));
					else
						for (String tagId : tagIds)
							filterTags.add((Tag) Globals.session().get(Tag.class, Long.valueOf(tagId)));
				}
			}
			else
			{
				// Update the tags of a property
				Property property = (Property) Globals.session().get(Property.class, getLongParam("id"));
				property.tags.clear();

				if (tagIds != null)
					for (String tagId : tagIds)
						property.tags.add((Tag) Globals.session().get(Tag.class, Long.valueOf(tagId))); 
				Globals.session().saveOrUpdate(property);
			}
		}
		else if ("uploadfile".equals(getParam("action")))
		{
			Tag tag = (Tag) Globals.session().get(Tag.class, Long.valueOf(getLongParam("id")));
			if (!tag.canEdit())
				throw new UserFriendlyException("You are not allowed to edit this tag");

			File f = Globals.getUploadedFile();

			List<Long> alreadyThere = new ArrayList<Long>();

			int i = 0;
			Set<Long> ids = new HashSet<Long>();
			long time = Calendar.getInstance().getTimeInMillis();
			long skipCount = assertParam("skipCount") ? getLongParam("skipCount") : 0;

			// FIXME: if it is a large data set, this can take A LOT of memory -> implement an iterator in ChemDAO
			List<String> molecules = Various.molecule.readSDFMolsFromFile(f.getAbsolutePath());
			for (String sdfMol : molecules) {
				try
				{
					if (skipCount-- > 0)
						continue;
					Molecule mol = MoleculePeer.getMolecule(sdfMol);
					ids.add(mol.mapping1.id);
					if (i++ % 100 == 0)
					{
						addTag(ids, alreadyThere, tag);
						logger.info("" + i + " compounds processed. Last 100 took "+Math.round(1.0f * (Calendar.getInstance().getTimeInMillis() - time) / 1000)+" sec. Restarting transaction...");
						time = Calendar.getInstance().getTimeInMillis();
						ids.clear();
						Globals.restartAllTransactions(true);
					}
				} catch (HibernateException e)
				{
					Thread.sleep(2000);
					try
					{
						Globals.restartAllTransactions(true);
					} catch (Exception e2){

					}
				}
			}
			addTag(ids, alreadyThere, tag);
			Globals.restartAllTransactions(true);
		}

		return new WebModel().getModelAndView();
	}

	public ModelAndView copy(HttpServletRequest req, HttpServletResponse res)
			throws Exception
	{
		Tag tag = (Tag) Globals.session().get(Tag.class, getLongParam("id"));
		Tag copy = new Tag();
		copy.name = tag.name + " copy";
		copy.introducer = copy.owner = Globals.userSession().user;
		copy.isPublic = false;
		copy.showInBrowser = false;
		copy.type = tag.type;
		Globals.session().save(copy);

		Globals.session().createSQLQuery("insert into MoleculeTag(tag_id, mapping1_id) select "+copy.id+", mapping1_id from MoleculeTag where tag_id=" + tag.id).executeUpdate();

		return redirect("tags/show.do?id=" + copy.id);

	}

	private void addTag(Set<Long> ids, List<Long> alreadyThere, Tag tag)
	{
		if (alreadyThere.size() > 0)
			ids.removeAll(alreadyThere);
		if (ids.size() > 0)
		{
			logger.info("Adding " + ids.size() + " molecules");
			alreadyThere.addAll(ids);
			Globals.session()
			.createSQLQuery("insert ignore into MoleculeTag(tag_id, mapping1_id) select " + tag.id + ", mapping1_id from Mapping1 where mapping1_id in (:ids)")
			.setParameterList("ids", ids).executeUpdate();
		}
	}

	public ModelAndView show(HttpServletRequest req, HttpServletResponse res)
	{
		return new WebModel().setTemplate("tags-browser").getModelAndView();
	}

	public ModelAndView profile(HttpServletRequest req, HttpServletResponse res)
	{
		Tag tag = (Tag) Globals.session().get(Tag.class, getLongParam("id"));
		return new WebModel(tag).setTemplate("tag-profile").getModelAndView();
	} 

	public ModelAndView edit(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		Tag tag;

		if (!assertParam("id") || getLongParam("id") < 0)
			tag = new Tag();
		else
			tag = (Tag) Globals.session().get(Tag.class, Long.valueOf(request.getParameter("id")));

		return new WebModel(tag).setRenderMode("popup").setTemplate("tag-edit").getModelAndView();
	}

	@SuppressWarnings({ "unchecked", "rawtypes" })
	public ModelAndView export(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		ExportableSet eData = new ExportableSet();
		eData.clearColumns();
		eData.addColumn(ExportableColumn.SMILES);
		eData.addColumn(ExportableColumn.CASRN);
		eData.addColumn(ExportableColumn.MOLECULEID);
		eData.addColumn(ExportableColumn.NAMES);
		eData.uncheckColumn(ExportableColumn.NAMES);

		if (assertParam("submit"))
		{
			eData.configure(ExportableSetConfiguration.configureFromDialog(request));

			Tag tag = (Tag) Globals.session().get(Tag.class, getLongParam("id"));
			Criteria criteria = Globals.session().createCriteria(Mapping1.class);
			if(tag.type.equals("property"))
			{
				List<Property> pls = Globals.session().createCriteria(Property.class)
						.createAlias("tags", "t")
						.add(Restrictions.eq("t.id", tag.id)).list();
				criteria.createAlias("molecules", "m").createAlias("m.experimentalProperties", "ep")
				.add(Restrictions.in("ep.property", pls.toArray()));
			} else
			{
				criteria.createAlias("tags", "t")
				.add(Restrictions.eq("t.id", tag.id));
			}
			List ids = criteria
					.setProjection(Projections.distinct(Projections.id()))
					.list();

			for (int fromIndex = 0; fromIndex < ids.size(); fromIndex += 100)
			{
				if (eData.exportableMolecules.size() % 1000 == 0 && eData.exportableMolecules.size() > 0)
				{
					Globals.restartAllTransactions(true);
				}

				int toIndex = fromIndex + 100;
				if (toIndex > ids.size())
					toIndex = ids.size();

				List<Mapping1> ls = Globals.session().createCriteria(Mapping1.class)
						.add(Restrictions.in("id", ids.subList(fromIndex, toIndex))).list();

				for (Mapping1 mp1 : ls) 
				{
					ExportableMolecule eMol = new ExportableMolecule();
					eMol.setMolecule(mp1.molecules.get(0));
					eData.addMolecule(eMol);
					if (eData.selectedColumns.contains(ExportableColumn.NAMES))
					{
						List<MoleculeName> molNames = Globals.session().createCriteria(MoleculeName.class)
								.setResultTransformer(Criteria.DISTINCT_ROOT_ENTITY)
								.createAlias("records", "rec")
								.add(Restrictions.in("rec.molecule", mp1.molecules))
								.setMaxResults(10).list();
						eMol.setMoleculeNames(molNames);
					}

				}
			}

			HttpResponseExportWriter eWriter = new HttpResponseExportWriter(response, ExportWriter.createWriter(getParam("format"), eData, null));
			eWriter.setFileName(tag.name);
			eWriter.write();
			return null;
		}
		else
			return new WebModel(eData).setTemplate("export").getModelAndView();
	}
}
