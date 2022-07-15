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

package qspr.dao;

/**
 * The repository layer.
 * Holds the instances of the data access objects (DAOs)
 * In future, we might move this logic to a dependency injection container (e.g., Spring)
 * 
 * Only part of the DA operations have been moved here. We will move the rest gradually.
 * This refactoring is required to be able to "mock" data access with stubs, for example for implementing standalone DB-free predictor.
 * 
 * @author midnighter
 *
 */
public class Repository {
	public static ModelTemplateDAO modelTemplate = new ModelTemplateDAOImpl();
	public static PropertyDAO property = new PropertyDAOImpl();
	public static PropertyOptionDAO option = new PropertyOptionDAOImpl();
	public static UnitDAO unit = new UnitDAOImpl();
	public static MoleculeDAO molecule = new MoleculeDAOImpl();
	public static ModelDAO model = new ModelDAOImpl();
	public static UserDAO user = new UserDAOImpl();
	public static BasketDAO basket = new BasketDAOImp();
	public static OthersDAO various = new OthersDAOImp();
	public static RecordDAO record = new RecordDAOImp();
	public static ArticleDAO article = new ArticleDAOImp();
}
