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

package qspr.metaserver.configurations;

import java.io.IOException;

import javax.servlet.http.HttpServletRequest;

/**
 * New configuration in order to allow smooth migration from Old to New Dragon configuration
 * @author itetko
 *
 */

abstract public class DescriptorsAbstractDragonConfiguration extends DescriptorsAbstractConfiguration
{
	private static final long serialVersionUID = 1L;

	// each bit represents one of the descriptor blocks
	public  long  dragonBlocks=Long.MAX_VALUE; // by default all descriptor blocks are selected

	/*
	 * Stores mask indicating blocks required for 3D decsriptor calculations
	 */

	final static public String DRAGONS[] = {DescriptorsConfiguration.RDKIT};

	public DescriptorsAbstractDragonConfiguration()
	{
		//TODO Remove once all these configurations are removed from the database
		try {
			setAllOn();
		}catch(Exception e) {
			//throw new UserFriendlyException(e);
			//Mailer.notifyAdmins("Outdated configuration", "again!");
			System.out.println("Outdated configuration");
		}
	}

	public DescriptorsAbstractDragonConfiguration(long dragonBlocks)
	{
		this.dragonBlocks = dragonBlocks;
	}

	/**
	 *  Checks whether block i is selected
	 * @param i
	 * @return
	 */

	public boolean bit(long i){
		return ((1l<<i)&dragonBlocks)>0;	
	}

	public void setBit(long i) {
		dragonBlocks |= 1l << i; 
	}

	public String toString()
	{
		String s=requires3D()?"3D blocks: ":"2D blocks: ";
		if(requires3D()) {
			long d2 = 0;
			for(long i = 0, block = 1;i<getBlocks(); i++, block*=2)
				if( (block & get3DMask()) == 0) // not 3D bit
					d2 += block;
			if(d2 == 0) s += "(only) ";
		}
		
		boolean first=true;
		for(long i=0;(1l<<i)<=dragonBlocks;i++)if(bit(i)){
			if(first||(1l<<(i+1))>=dragonBlocks||bit(i+1)!=bit(i)){s+=""+(i+1);first=true;}
			if(first&&(1l<<(i+1))<dragonBlocks)s+=bit(i+1)==bit(i)?"-":" ";
			first=false;
		}
		else
			first=true;

		return s;
	}


	abstract public int getBlocks();

	abstract public String getExtension();

	abstract protected long get3DMask();

	@Override
	public DescriptorsAbstractConfiguration setAllOn() throws IOException{
		dragonBlocks = 0;
		for(int i=0;i<getBlocks(); i++)setBit(i);
		return this;
	}

	@Override
	public boolean requires3D() {
		return (dragonBlocks & get3DMask()) != 0;
	}

	@Override
	protected void setAllOn2D() throws IOException{
		dragonBlocks = 0;
		for(long i = 0, block = 1;i<getBlocks(); i++, block*=2)
			if( (block & get3DMask()) == 0) // not 3D bit
				dragonBlocks += block;
	}

	@Override
	public boolean isCachable() {
		return true;
	}

	@Override
	DescriptorsAbstractConfiguration setConfiguration(HttpServletRequest request) {

		dragonBlocks = 0;
		for (long i = 1, block = 1; i <= getBlocks(); i++, block *=2)
		{
			if (request.getParameter(getExtension() + i) != null)
				dragonBlocks |= block;
		}

		if(dragonBlocks == 0) return null;
		
		return this;

	}

}



// Note:
// Descriptor blocks of Dragon 6.0.0
// 1) Constitutional descriptors 	
// 2) Ring descriptors 	
// 3) Topological indices 	
// 4) Walk and path counts 	
// 5) Connectivity indices 	
// 6) Information indices 	
// 7) 2D matrix-based descriptors 	
// 8) 2D autocorrelations 	
// 9) Burden eigenvalues 	
// 10) P_VSA-like descriptors 	
// 11) ETA indices 	
// 12) Edge adjacency indices 	
// 13) Geometrical descriptors 	
// 14) 3D matrix-based descriptors 	
// 15) 3D autocorrelations 	
// 16) RDF descriptors 	
// 17) 3D-MoRSE descriptors 	
// 18) WHIM descriptors 	
// 19) GETAWAY descriptors 	
// 20) Randic molecular profiles 	
// 21) Functional group counts 	
// 22) Atom-centred fragments 	
// 23) Atom-type E-state indices 	
// 24) CATS 2D 	
// 25) 2D Atom Pairs 	
// 26) 3D Atom Pairs 	
// 27) Charge descriptors 	
// 28) Molecular properties 	
// 29) Drug-like indices 	

// Descriptor blocks of Dragon 1.2.4
// 1	Constitutional descriptors	48
// 2	Topological descriptors	119
// 3	Walk and path counts	47
// 4	Connectivity indices	33
// 5	Information indices	47
// 6	2D autocorrelations	96
// 7	Edge adjacency indices	107
// 8	BCUT descriptors	64
// 9	Topological charge indices	21
// 10	Eigenvalue-based indices	44
// 11	Randic molecular profiles	41
// 12	Geometrical descriptors	74
// 13	RDF descriptors	150
// 14	3D-MoRSE descriptors	160
// 15	WHIM descriptors	99
// 16	GETAWAY descriptors	197
// 17	Functional group counts	121
// 18	Atom-centred fragments	120
// 19	Charge descriptors	14
// 20	Molecular properties

