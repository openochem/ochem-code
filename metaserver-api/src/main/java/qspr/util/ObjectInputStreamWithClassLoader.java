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

// Ok guys. This class I have taken from JBoss
// A hardcore solution for unserializing objects using given classloader
// Midnighter

import java.io.IOException;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.io.ObjectStreamClass;
import java.lang.reflect.Proxy;

/**
 * This replaces the EjbossInputStream in the storage package.
 * The input stream will take a class loader in its constructor and look
 * into it to retrieve the class definitions.
 * It is used throughout the server to deserialize parameters and objects
 * whose definition are in a jar and not the global classpath
 * It also has better comments than the previous version.
 *
 * @author  <a href="rickard@dreambean.com">Rickard Oberg</a>
 * @since   Ejboss 0.9
 */


public class ObjectInputStreamWithClassLoader
    extends ObjectInputStream {

    /**
    * The classloader to use when the default classloader cannot find
    * the classes in the stream.
    */
    ClassLoader cl;

    
/******************************************************************************/
/******************************************************************************/
/*
/*   CONSTRUCTORS
/*
/******************************************************************************/
/******************************************************************************/

    /**
    * Construct a new instance with the given classloader and input stream.
    *
    * @param  ClassLoader      classloader to use
    * @param  InputStream      stream to read objects from
    */
    public ObjectInputStreamWithClassLoader(InputStream in, ClassLoader cl)
         throws IOException {

         super(in);

         this.cl = cl;
    }


/******************************************************************************/
/******************************************************************************/
/*
/*   OVERWRITING  <ObjectInputStream>
/*
/******************************************************************************/
/******************************************************************************/

    /**
    * Resolve the class described in the osc parameter. First, try the
    * default classloader (implemented by the super class). If it cannot
    * load the class, try the classloader given to this instance.
    *
    * @param  ObjectStreamClass     class description object
    * @return      the Class corresponding to class description
    * @exception   IOException     if an I/O error occurs
    * @exception   ClassNotFoundException  if the class cannot be found by the classloader
    */
    @SuppressWarnings({ "unchecked", "rawtypes" })
	protected Class resolveClass(ObjectStreamClass osc)
        throws IOException, ClassNotFoundException {

        return
        	Class.forName(osc.getName(), false, cl);
        //cl.loadClass(osc.getName());
    }
    
    @SuppressWarnings({ "unchecked", "rawtypes" })
	protected Class resolveProxyClass( String[] interfaces )
       	throws IOException, ClassNotFoundException {
		
		Class[] interfacesClass = new Class[interfaces.length];
		for( int i=0; i< interfaces.length; i++ )
		{
			interfacesClass[i] = Class.forName(interfaces[i], false, cl);
		}
		
    	return Proxy.getProxyClass(cl, interfacesClass);
    }
}




