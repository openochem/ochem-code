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
 * A filtering request for experimental records on OCHEM.
 * Under construction!
 * @author midnighter
 *
 */
public class ExperimentalRecordsFilter
{
	public Long propertyId;
	public Long basketId;
	public Long introducerId;
	public Long modifierId;
	public Long articleId;
	public String page;
	public String line;
	public Integer mwMin;
	public Integer mwMax;
	public Long propertyOptionId;
	public String propertyOptionName;
	public SubstructureSearchFilter substructureSearch;
	public String combinedNameFilter;
	
	public boolean publicRecords = true;
	public boolean privateRecords = true;
	public boolean originalRecordsOnly = false;
	public boolean primaryRecordsOnly = false;
	
	public ApprovalStatus dataFromOtherUsers = ApprovalStatus.ALL;
			
	public Long getPropertyId()
	{
		return propertyId;
	}
	public void setPropertyId(Long propertyId)
	{
		this.propertyId = propertyId;
	}
	public Long getBasketId()
	{
		return basketId;
	}
	public void setBasketId(Long basketId)
	{
		this.basketId = basketId;
	}
	public Long getIntroducerId()
	{
		return introducerId;
	}
	public void setIntroducerId(Long introducerId)
	{
		this.introducerId = introducerId;
	}
	public Long getModifierId()
	{
		return modifierId;
	}
	public void setModifierId(Long modifierId)
	{
		this.modifierId = modifierId;
	}
	public Long getArticleId()
	{
		return articleId;
	}
	public void setArticleId(Long articleId)
	{
		this.articleId = articleId;
	}
	public String getPage()
	{
		return page;
	}
	public void setPage(String page)
	{
		this.page = page;
	}
	public String getLine()
	{
		return line;
	}
	public void setLine(String line)
	{
		this.line = line;
	}
	public Integer getMwMin()
	{
		return mwMin;
	}
	public void setMwMin(Integer mwMin)
	{
		this.mwMin = mwMin;
	}
	public Integer getMwMax()
	{
		return mwMax;
	}
	public void setMwMax(Integer mwMax)
	{
		this.mwMax = mwMax;
	}
	public Long getPropertyOptionId()
	{
		return propertyOptionId;
	}
	public void setPropertyOptionId(Long propertyOptionId)
	{
		this.propertyOptionId = propertyOptionId;
	}
	public String getPropertyOptionName()
	{
		return propertyOptionName;
	}
	public void setPropertyOptionName(String propertyOptionName)
	{
		this.propertyOptionName = propertyOptionName;
	}
	public SubstructureSearchFilter getSubstructureSearch()
	{
		return substructureSearch;
	}
	public void setSubstructureSearch(SubstructureSearchFilter substructureSearch)
	{
		this.substructureSearch = substructureSearch;
	}
	public String getCombinedNameFilter()
	{
		return combinedNameFilter;
	}
	public void setCombinedNameFilter(String combinedNameFilter)
	{
		this.combinedNameFilter = combinedNameFilter;
	}
	public boolean isPublicRecords()
	{
		return publicRecords;
	}
	public void setPublicRecords(boolean publicRecords)
	{
		this.publicRecords = publicRecords;
	}
	public boolean isPrivateRecords()
	{
		return privateRecords;
	}
	public void setPrivateRecords(boolean privateRecords)
	{
		this.privateRecords = privateRecords;
	}
	public boolean isOriginalRecordsOnly()
	{
		return originalRecordsOnly;
	}
	public void setOriginalRecordsOnly(boolean originalRecordsOnly)
	{
		this.originalRecordsOnly = originalRecordsOnly;
	}
	public boolean isPrimaryRecordsOnly()
	{
		return primaryRecordsOnly;
	}
	public void setPrimaryRecordsOnly(boolean primaryRecordsOnly)
	{
		this.primaryRecordsOnly = primaryRecordsOnly;
	}
	
	public boolean isErrorRecords()
	{
		return errorRecords;
	}
	public void setErrorRecords(boolean errorRecords)
	{
		this.errorRecords = errorRecords;
	}
	public boolean isErrorInchies()
	{
		return errorInchies;
	}
	public void setErrorInchies(boolean errorInchies)
	{
		this.errorInchies = errorInchies;
	}
	public boolean isEmptyMolecules()
	{
		return emptyMolecules;
	}
	public void setEmptyMolecules(boolean emptyMolecules)
	{
		this.emptyMolecules = emptyMolecules;
	}
	public boolean isMismatchingNames()
	{
		return mismatchingNames;
	}
	public void setMismatchingNames(boolean mismatchingNames)
	{
		this.mismatchingNames = mismatchingNames;
	}
	public static enum ApprovalStatus {ALL, APPROVED_ONLY, AWAITING_ONLY};
	
	public boolean errorRecords;
	public boolean errorInchies;
	public boolean emptyMolecules;
	public boolean mismatchingNames;
}


