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

package qspr.business;

/**
 * A filter for models. Under construction.
 * @author midnighter
 *
 */
public class ModelFilter
{
	public String query;
	
	public Long publicId;

	public boolean awaitingApproval;

	public String proQuery;

	public long propertyId;

	public long basketId;

	public long trainingSetId;

	public String articleId;

	public boolean privateOnly;

	public boolean publicOnly;

	public String templateName;

	public String order;

	public boolean groupModels = false;
	
	public boolean toBeDeleted;

	public String getQuery()
	{
		return query;
	}

	public void setQuery(String query)
	{
		this.query = query;
	}

	public Long getPublicId()
	{
		return publicId;
	}

	public void setPublicId(Long publicId)
	{
		this.publicId = publicId;
	}

	public boolean isAwaitingApproval()
	{
		return awaitingApproval;
	}

	public void setAwaitingApproval(boolean awaitingApproval)
	{
		this.awaitingApproval = awaitingApproval;
	}

	public String getProQuery()
	{
		return proQuery;
	}

	public void setProQuery(String proQuery)
	{
		this.proQuery = proQuery;
	}

	public long getPropertyId()
	{
		return propertyId;
	}

	public void setPropertyId(long propertyId)
	{
		this.propertyId = propertyId;
	}

	public long getBasketId()
	{
		return basketId;
	}

	public void setBasketId(long basketId)
	{
		this.basketId = basketId;
	}

	public long getTrainingSetId()
	{
		return trainingSetId;
	}

	public void setTrainingSetId(long trainingSetId)
	{
		this.trainingSetId = trainingSetId;
	}

	public String getArticleId()
	{
		return articleId;
	}

	public void setArticleId(String articleId)
	{
		this.articleId = articleId;
	}

	public boolean isPrivateOnly()
	{
		return privateOnly;
	}

	public void setPrivateOnly(boolean privateOnly)
	{
		this.privateOnly = privateOnly;
	}

	public boolean isPublicOnly()
	{
		return publicOnly;
	}

	public void setPublicOnly(boolean publicOnly)
	{
		this.publicOnly = publicOnly;
	}

	public String getTemplateName()
	{
		return templateName;
	}

	public void setTemplateName(String templateName)
	{
		this.templateName = templateName;
	}

	public String getOrder()
	{
		return order;
	}

	public void setOrder(String order)
	{
		this.order = order;
	}

	public boolean isGroupModels()
	{
		return groupModels;
	}

	public void setGroupModels(boolean groupModels)
	{
		this.groupModels = groupModels;
	}
}
