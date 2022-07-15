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

package qspr.entities;

import java.io.File;
import java.sql.Timestamp;
import java.text.SimpleDateFormat;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.persistence.Transient;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import qspr.Globals;
import qspr.OCHEMConfiguration;
import qspr.util.ExportThread;

import com.eadmet.utils.OCHEMUtils;

// A log of the fact that a user has exported some data from  DB
// Intended for controlling and limiting the data export
// Midnighter

@Entity
@XmlRootElement(name = "export-action")
public class ExportAction
{
	@Id
	@GeneratedValue
	@Column(name = "ea_id")
	@XmlAttribute
	public Long id;
	
	@ManyToOne
	@JoinColumn(name = "session_id")
	public Session session;
	
	@Column(name = "exported_columns")
	public String exportedColumns;
	public String format;
	
	@Column(name = "filename")
	private String fileName;
	
	public String getFileName() {
		return fileName;
	}

	public void setFileName(String fileName) 
	{
		fileName = fileName.replaceAll("[^a-zA-Z0-9]", " ").replaceAll("\\s+", "_");
		if (fileName.length() > 98)
			this.fileName = fileName.substring(0, 98);
		else
			this.fileName = fileName;
	}

	public int count; // the total number of exported records
	
	@Column(name = "count_restricted")
	public int countRestricted; // the number of the restricted access records
	
	@XmlTransient
	public Timestamp time;
		
	public String status;
	
	@Column(name = "dataset_summary")
	public String datasetSummary;
	
	public boolean failed;
	
	@Column(name = "file_downloaded")
	public boolean fileDownloaded;
	
	@Column(name = "free_bonus_points_used")
	public double freeBonusPointsUsed;
	
	@Column(name = "total_bonus_points")
	public double totalBonusPoints;
	
	@Transient
	public boolean running = false;
	
	@XmlElement
	public String getFormattedTime()
	{
		return new SimpleDateFormat("HH:mm, dd MMM yy").format(time);
	}
	
	@XmlElement
	public double getFreeBonusPointsAvailable()
	{
			return 0D;
	}
	
	@XmlElement
	public double getEarnedBonusPointsAvailable()
	{
			return 0D;
	}
	
	@XmlElement
	public double getTotalBonusPointsAvailable()
	{
			return 0D;
	}
	
	@XmlElement
	public Boolean getEnoughBonuses()
	{
		return true;
	}
	
	
	@XmlElement
	public String getPublicURL()
	{
		if (fileDownloaded)
			return OCHEMConfiguration.getRootHost()+OCHEMConfiguration.rootDir+"/"+getDownloadUrl();
		else
			return null;
	}
	
	public String getFullFilePath()
	{
		return getExportFolder()+"/"+getFileName()+"."+format;
	}
	
	public String getExportFolder()
	{
		String exportPath = Globals.commonDownloadDirectory;
		
		if (exportPath == null)
			exportPath = new File(ExportThread.class.getClassLoader().getResource("views.properties").getPath()).getParentFile().getAbsolutePath() + "/../../src/main/webapp/exports";
		
		if (Globals.userSession() != null)
			exportPath += (Globals.userSession().user != null) ? "/" + Globals.userSession().user.login : "/anonymous";
		else
			exportPath += "/nosession";
		
		exportPath += "/"+OCHEMUtils.getMD5(getFormattedTime()+"_"+getFileName()+"_");
		File f = new File(exportPath);
		if (!f.exists())
			f.mkdirs();
		return exportPath;
	}
	
	public String getDownloadUrl()
	{
		String path = getFullFilePath();
		int pos = path.indexOf("export");
		return path.substring(pos);
	}
}
