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

package qspr.util;

import java.util.List;

import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.entities.Molecule;

/**
 * Utility to correct errors with incompatible compression of sdfs in Molecules and in RGroups
 * @author itetko
 *
 */
public class MoleculeCleaningUtility
{
	//	private static final org.apache.logging.log4j.Logger logger = org.apache.logging.log4j.Logger.getLogger(MoleculeCleaningUtility.class);
	@SuppressWarnings("unchecked")
	public void execute()
	{
		List<Long> mols = Globals.session().createCriteria(Molecule.class).setProjection(Projections.id()).list();

		while (mols.size() > 0)
		{
			List<Long> batchIDs = mols.subList(0, Math.min(mols.size(), 1000));
			System.out.println("Processing " + batchIDs.size() + " out of " + mols.size());
			List<Molecule> oMols = Globals.session().createCriteria(Molecule.class).add(Restrictions.in("id", batchIDs)).list();
			for (Molecule molecule : oMols)
			{
				String dat = molecule.getData();
				molecule.setData(dat);
				Globals.session().saveOrUpdate(molecule);                       
			}
			batchIDs.clear();
			Globals.restartAllTransactions(true);
		}
	}

/*	
	@SuppressWarnings("unchecked")
	public void executeR() throws IOException
	{
		List<Long> mols = Globals.session().createCriteria(RGroup.class).setProjection(Projections.id()).list();

		while (mols.size() > 0)
		{
			List<Long> batchIDs = mols.subList(0, Math.min(mols.size(), 1000));
			System.out.println("Processing " + batchIDs.size() + " out of " + mols.size());
			List<RGroup> oMols = Globals.session().createCriteria(RGroup.class).add(Restrictions.in("id", batchIDs)).list();
			for (RGroup molecule : oMols)
			{
				String dat = molecule.getData();
				molecule.setData(dat);
				Globals.session().saveOrUpdate(molecule);                       
			}
			batchIDs.clear();
			Globals.restartAllTransactions(true);
		}
	}
*/
	
	public static void main(final String[] args) throws InterruptedException
	{
		final MoleculeCleaningUtility utility = new MoleculeCleaningUtility();
		WrapperThread wt = new WrapperThread()
		{

			@Override
			public void wrapped() throws Exception
			{
				utility.execute();	
//				utility.executeR();	
			}
		};

		wt.start();
		wt.join();
	}
}
