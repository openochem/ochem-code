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

package qspr.tests.depiction;

import java.awt.FlowLayout;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.image.BufferedImage;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;

import javax.imageio.ImageIO;
import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JLabel;

import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;

import com.eadmet.exceptions.CriticalException;
import com.eadmet.exceptions.UserFriendlyException;

import junit.framework.TestCase;
import qspr.OCHEMConfiguration;
import qspr.dao.ChemInfEngine;
import qspr.dao.Various;
import qspr.depiction.MoleculeDepiction;

@RunWith(Parameterized.class)
public abstract class DepictionTestCaseBase extends TestCase {

	private static Object lock = new Object();
	protected static final boolean show = false; // set to true to see images
	
	protected static MoleculeDepiction depiction = null;
	protected static ChemInfEngine depictionChoice = null;

	public DepictionTestCaseBase(ChemInfEngine engine) {
        super();
        depiction = MoleculeDepiction.get(engine);
        depictionChoice = engine;
    }
    
    @Parameterized.Parameters
    public static Collection<Object[]> input() {
    	ArrayList<Object[]> engines = new ArrayList<>(Arrays.asList(new Object[][] {
			new Object[] {ChemInfEngine.CDK}
		}));
		
		if (OCHEMConfiguration.chemaxonLicenceAvailable()) {
			try {
				Various.getCheminfImpl(ChemInfEngine.CHEMAXON);
    		} catch (UserFriendlyException e) {
    			System.err.println(ChemInfEngine.CHEMAXON + " is enabled, but no license was found.");
    			throw new CriticalException(e);
    		}
			engines.add(new Object[] {ChemInfEngine.CHEMAXON});
		}
		
		return engines;
    }
	
	protected void showImage(byte[] imgBinary) throws IOException, InterruptedException {
		InputStream is = new ByteArrayInputStream(imgBinary);
        BufferedImage img = ImageIO.read(is);
        ImageIcon icon=new ImageIcon(img);
        JFrame frame=new JFrame();
        frame.setLayout(new FlowLayout());
        frame.setSize(200,300);
        JLabel lbl=new JLabel();
        lbl.setIcon(icon);
        frame.add(lbl);
        frame.setVisible(true);
        frame.setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
        
        Thread t = new Thread() {
            public void run() {
                synchronized(lock) {
                    while (frame.isVisible())
                        try {
                            lock.wait();
                        } catch (InterruptedException e) {
                            e.printStackTrace();
                        }
                    System.out.println("Closed.");
                }
            }
        };
        t.start();

        frame.addWindowListener(new WindowAdapter() {

            @Override
            public void windowClosing(WindowEvent arg0) {
                synchronized (lock) {
                    frame.setVisible(false);
                    lock.notify();
                }
            }

        });

        t.join();
	}
}
