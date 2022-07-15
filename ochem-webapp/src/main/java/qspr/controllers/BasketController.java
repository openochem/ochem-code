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
import java.util.List;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.hibernate.Criteria;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;
import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.MarshallingOption;
import qspr.business.BasketFilter;
import qspr.business.BasketPeer;
import qspr.business.RandomBasketSplitter;
import qspr.entities.Basket;
import qspr.entities.BasketEntry;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Model;
import qspr.entities.Property;
import qspr.export.ExportableColumn;
import qspr.export.ExportableSet;
import qspr.export.ExportableSetConfiguration;
import qspr.frontend.WebList;
import qspr.frontend.WebModel;
import qspr.util.ExportThread;

import com.eadmet.business.BasketService;
import com.eadmet.business.DiscretizeBasketOptions;
import com.eadmet.business.DiscretizeBasketRunner;
import com.eadmet.business.SolventConvertBasketRunner;
import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.useractions.BasketAction;
import com.eadmet.useractions.BasketAction.BasketActionType;
import com.eadmet.utils.OCHEMUtils;
import com.eadmet.useractions.EventFactory;

@Controller
public class BasketController extends BrowserWrapper 
{
	//	private static transient final Logger logger = Logger.getLogger(BasketController.class);

	public BasketController()
	{
		sessionRequired = true;
	}

	BasketService service = new BasketService();

	public ModelAndView show(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		return new WebModel().setTemplate("basket-browser").getModelAndView();
	}

	public ModelAndView action(HttpServletRequest request, HttpServletResponse response)
	{
		// 		not used Rob		
		//		List<String> baskets = (List<String>) request.getSession().getAttribute("baskets");
		Basket basket;
		boolean newbasket = !assertParam("id") || (getParam("id").equals("-1"));
		if(newbasket)
			basket = new Basket();
		else
			basket = (Basket) Globals.session().get(Basket.class, getLongParam("id"));
		if ("delete".equals(request.getParameter("action")))
		{
			if (!basket.isEditableBy(Globals.userSession()))
				throw new UserFriendlyException("You are not authorized to modify this basket");
			Long count = (Long) Globals.session().createCriteria(Model.class)
					.add(Restrictions.or(Restrictions.eq("trainingSet", basket), Restrictions.eq("validationSet", basket)))
					.setProjection(Projections.countDistinct("id")).list().get(0);

			if (count > 0)
				throw new UserFriendlyException("" + count + " model(s), connected to this basket, are still in pending tasks. Please delete these models first");

			basket.delete();
			EventFactory.document("Basket delete", new BasketAction(basket, 0, BasketActionType.DELETE_BASKET));
		}
		else if ("edit".equals(request.getParameter("action")))
		{
			if (!basket.isEditableBy(Globals.userSession()))
				throw new UserFriendlyException("You are not authorized to modify this basket");
			String name = request.getParameter("name");
			if(basket.id == null)
				basket = Basket.getBasket(Globals.userSession(), name);
			else
			{
				if (name != null && !name.equals(basket.name))
				{
					basket.name = OCHEMUtils.getFilteredBasketName(name);
					basket.basketType = 0L;
					EventFactory.document("Basket rename", new BasketAction(basket, 0, BasketActionType.RENAME));
				}
				if (assertParam("description"))
					basket.description = getParam("description");
				if (assertParam("excludedrecords"))
					basket.description = getParam("excludedrecords");
			}
			Globals.session().saveOrUpdate(basket);
		}
		else if ("fillrecords".equals(request.getParameter("action")))
		{
			BasketEntry be = new BasketEntry();
			be.ep = (ExperimentalProperty) Globals.session().get(ExperimentalProperty.class, getLongParam("record ID"));
		}
		return new WebModel(basket).getModelAndView();
	}

	private BasketFilter getBasketFilterFromRequest()
	{
		BasketFilter filter = new BasketFilter();
		if (assertParam("id"))
			filter.id = getLongParam("id");

		filter.showGroupBaskets = assertParam("group");
		filter.showPublicBaskets = assertParam("public");
		if (assertParam("name"))
			filter.name = getParam("name");

		filter.showSystemBaskets = assertParam("showsystem");
		return filter;

	}

	public ModelAndView list(HttpServletRequest request, HttpServletResponse response)
	{
		if (assertParam("no.basket.details"))
			Globals.setMarshallingOption(MarshallingOption.NO_BASKET_DETAILS);

		Globals.setMarshallingOption(MarshallingOption.BASKET_LIST_MODE);
		int pageSize = 15;
		Criteria basketCriteria = BasketPeer.getListCriteria(getBasketFilterFromRequest());
		WebList webList = new WebList();
		webList.loadFromCriteria(basketCriteria, getPageNum(), getPageSize(pageSize));
		return new WebModel(webList).getModelAndView();
	}

	public ModelAndView listConditions(HttpServletRequest request, HttpServletResponse response)
	{
		Basket basket = Basket.getBasket(Globals.userSession(), getLongParam("id"));
		List<Property> conditionsUsed = service.listUsedConditions(basket);
		WebList wl = new WebList();
		wl.loadFromList(conditionsUsed);
		return new WebModel(wl).getModelAndView();
	}

	public ModelAndView clone(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		Basket b = Basket.getBasket(Globals.userSession(), getLongParam("basket"));
		Basket clone = service.clone(b);		
		return redirect("basket/edit.do?id=" + clone.id + "&render-mode=popup");
	}

	public ModelAndView getbasket(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		Basket bs = null;

		if(request.getParameter("id") != null)
		{
			bs = (Basket) Globals.session().get(Basket.class, getLongParam("id"));
			return new WebModel(bs).setTemplate("basket-record").setRenderMode("popup").getModelAndView();
		}
		return null;
	}

	//add records in basket
	public ModelAndView addRecords(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		Basket targetBasket = (Basket) Globals.session().get(Basket.class, getLongParam("id"));
		targetBasket.getPrivileges().requestModification();
		String action = getParam("action");

		if (!assertParam("type"))
			return new WebModel(targetBasket).setTemplate("basket-record").setRenderMode("popup").getModelAndView();

		if ("includeall".equals(action))
		{
			targetBasket.includeAllEntries();
			return new WebModel(targetBasket).setTemplate("basket-recordadded").setRenderMode("popup").getModelAndView();	
		}

		if (getParam("type").equals("basket"))
		{
			Basket sourceBasket = Basket.getBasket(Globals.userSession(), getLongParam("another-basket-id"));
			service.addRecordsFromAnotherBasket(targetBasket, sourceBasket, action);
		}
		else
		{
			// Work with records listed in an XLS file
			File f = null;
			try 
			{
				f = Globals.getUploadedFile();
			} catch (Exception e) 
			{
				e.printStackTrace();
				return new WebModel(targetBasket).setTemplate("basket-record").setRenderMode("popup").getModelAndView();
			}

			service.addRecordsFromXls(targetBasket, f, action);
		}

		return new WebModel(targetBasket).setTemplate("basket-recordadded").setRenderMode("popup").getModelAndView();		
	}

	//Export basket in excel sheet 
	@SuppressWarnings("unchecked")
	public ModelAndView exportBasket(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		ExportableSet eData = new ExportableSet();
		eData.removeColumn(ExportableColumn.APPLICABILITY_DOMAIN); 
		eData.removeColumn(ExportableColumn.PREDICTED_VALUE); 
		eData.removeColumn(ExportableColumn.DM_VALUE); 
		eData.removeColumn(ExportableColumn.ACCURACY);
		eData.removeColumn(ExportableColumn.DESCRIPTORS);
		eData.removeColumn(ExportableColumn.DESCRIPTORSNAMES);
		eData.removeColumn(ExportableColumn.ERROR);

		if (assertParam("submit"))
		{
			ExportThread eThread = service.getExportThread(getLongParam("id"), ExportableSetConfiguration.configureFromDialog(request), getParam("format"));
			eThread.start();
			return redirect("longoperations/operationWaitingScreen.do?operation-id=" + eThread.operationID);
		}
		else
		{
			Globals.setMarshallingOption(MarshallingOption.PROPERTY_UNITCATEGORY);
			Globals.setMarshallingOption(MarshallingOption.UNITCATEGORY_UNITS);
			List<Property> props = new ArrayList<Property>();
			List<Long> propIds = Globals.session().createCriteria(Property.class)
					.createAlias("experimentalProperties","eps")
					.createAlias("eps.basketEntries","bes")
					.createAlias("bes.basket","b")
					.add(Restrictions.eq("b.id", getLongParam("id")))
					.setProjection(Projections.distinct(Projections.id()))
					.list();
			if (propIds.size() > 0)
				props = Globals.session().createCriteria(Property.class).add(Restrictions.in("id", propIds)).list();

			eData.properties.addAll(props);
			return new WebModel(eData).setTemplate("export").getModelAndView();
		}
	}

	public ModelAndView split(HttpServletRequest request, HttpServletResponse response)
	{
		Basket basket = Basket.getBasket(Globals.userSession(), getLongParam("id"));

		if (!assertParam("splitting-method"))
			return new WebModel(basket).setTemplate("split-basket").getModelAndView();

		RandomBasketSplitter splitter = new RandomBasketSplitter();
		splitter.validationSetPercentage = getIntParam("validation-set-percentage");
		splitter.basket = basket;
		splitter.basketName1 = getParam("basket-name-1");
		splitter.basketName2 = getParam("basket-name-2");

		splitter.start();

		return new WebModel().getModelAndView();
	}

	public ModelAndView findPrimaryRecords(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		Basket basket = Basket.getBasket(Globals.userSession(), getLongParam("basket"));
		Basket primaryBasket = service.getPrimaryBasket(basket);
		return redirect("epbrowser/show.do?basket-select=" + primaryBasket.id + "&render-mode=popup");
	}

	/*
	public ModelAndView convertToSolvent(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		Basket basket = Basket.getBasket(Globals.userSession(), getLongParam("id"));
		basket = service.fillDiscretizeMetadata(basket);
		return new WebModel(basket).setTemplate("discretize-basket").getModelAndView();
	}
	*/

	public ModelAndView solventExractSubmit(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		new SolventConvertBasketRunner(Long.valueOf(request.getParameter("id"))).start();
		return null;
	}

	public ModelAndView discretize(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		Basket basket = Basket.getBasket(Globals.userSession(), getLongParam("id"));
		basket = service.fillDiscretizeMetadata(basket);
		return new WebModel(basket).setTemplate("discretize-basket").getModelAndView();
	}


	public ModelAndView discretizeSubmit(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		new DiscretizeBasketRunner(getOptions(request)).start();
		return null;
	}


	private DiscretizeBasketOptions getOptions(HttpServletRequest request)
	{
		DiscretizeBasketOptions opts = new DiscretizeBasketOptions();
		opts.options = request.getParameterValues("option");
		opts.strThresholds = request.getParameterValues("threshold");
		opts.thresholds = new double[opts.strThresholds.length];
		opts.basketId = Long.valueOf(request.getParameter("id"));
		opts.newBasketName = request.getParameter("basket-name");
		opts.newPropertyName = request.getParameter("property-name");
		opts.propertyId = Long.valueOf(request.getParameter("property"));
		return opts;
	}

	public ModelAndView edit(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		Basket basket = null;

		if (assertParam("id") && getLongParam("id") >= 0)
			basket = (Basket) Globals.session().get(Basket.class, Long.valueOf(request.getParameter("id")));
		else if (assertParam("name"))
			basket = Basket.getBasket(Globals.userSession(), getParam("name"), false);

		if (basket != null && basket.id != null)
			service.fillEditMetadata(basket);

		Globals.setMarshallingOption(MarshallingOption.PROPERTY_OPTIONS_FULL);
		return new WebModel(basket).setRenderMode("popup").setTemplate("basket-edit").getModelAndView();
	}
}