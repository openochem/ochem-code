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

import java.io.Serializable;
import java.sql.Timestamp;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.List;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import qspr.Globals;
import qspr.workflow.datatypes.AbstractDataRow;
import qspr.workflow.datatypes.DataTable;

@Entity
@XmlRootElement(name = "shuffle-key")
public class ShuffleKey 
{
	@Id
	@GeneratedValue
	@Column(name = "shuffle_key_id")
	@XmlAttribute
	public Long id;

	@Column(name="time_created")
	@XmlTransient
	public Timestamp timeCreated;

	@Column
	public String name;

	@ManyToOne
	@JoinColumn(name = "attachment_id")
	public Attachment<ShuffleKeyData> attachment;

	@ManyToOne
	@JoinColumn(name = "session_id")
	public Session session;

	@Column(name = "num_descriptors")
	public int numDescriptors;

	private transient int[] targetDescriptorsPermutation;

	public ShuffleKey()
	{
		timeCreated = new Timestamp(Calendar.getInstance().getTimeInMillis());
		session = Globals.userSession();
	}

	public void createKeyFrom(DataTable dtDescriptors)
	{
		int[] permutation = getRandomPermutation(dtDescriptors.getColumnsSize());
		numDescriptors = permutation.length;
		ShuffleKeyData shuffleData = new ShuffleKeyData();
		for (int i = 0; i < numDescriptors; i++)
			shuffleData.addDescriptor(dtDescriptors.getColumn(permutation[i]), 0.0, 0.0); // TODO: Calculate mean and std
		attachment = new Attachment<ShuffleKey.ShuffleKeyData>(shuffleData, AttachmentSource.ShuffleKey);
	}

	@XmlElement(name = "title")
	public String getTitle()
	{
		SimpleDateFormat sdf = new SimpleDateFormat("HH:mm 'on' dd MMMMM");
		return "KeyID " + id + " (" + numDescriptors + " descriptors, created at " + sdf.format(timeCreated) + ")";
	}

	public void save()
	{
		Globals.session().saveOrUpdate(attachment);
		Globals.session().saveOrUpdate(this);
	}

	public void setTargetDescriptors(DataTable dtTarget)
	{
		targetDescriptorsPermutation = getPermutationFor(dtTarget);
	}

	public Object getDescriptor(AbstractDataRow row, int num)
	{
		if (targetDescriptorsPermutation[num] == -1)
			return null; // a missing descriptor
		else
			return row.getValue(targetDescriptorsPermutation[num]);
	}

	// Get a shuffle-key permutation for a given datatable
	private int[] getPermutationFor(DataTable dtDescriptors)
	{
		List<String> keyDescriptorNames = attachment.getObject().descriptors;
		int[] permutation = new int[keyDescriptorNames.size()];
		for (int i = 0; i < keyDescriptorNames.size(); i++)
		{
			int index = dtDescriptors.indexOfColumn(keyDescriptorNames.get(i));
			// Missing descriptors will be set to NULLs later, no need to fail / Midnighter on Jun 16, 2011
			//if (index == -1)
			//	throw new UserFriendlyException("Key <" + getTitle() + "> is incompatible with given descriptors! Descriptor is " + keyDescriptorNames.get(i) + " is missing in the dataset");
			permutation[i] = index;
		}

		return permutation;
	}

	private int[] getRandomPermutation(int N)
	{
		int[] permutation = new int[N];

		for (int i = 0; i < N; i++)
			permutation[i] = i;

		// Shuffle
		for (int i = 0; i < N; i++)
		{
			int r = (int) (Math.random() * (i + 1));
			int swap = permutation[r];
			permutation[r] = permutation[i];
			permutation[i] = swap;
		}

		return permutation;
	}

	public static class ShuffleKeyData implements Serializable
	{
		private static final long serialVersionUID = 1L;

		public List<String> descriptors = new ArrayList<String>();
		public List<Double> means = new ArrayList<Double>(), stds = new ArrayList<Double>();

		public void addDescriptor(String desc, double mean, double std)
		{
			descriptors.add(desc);
			means.add(mean);
			stds.add(std);
		}
	}
}
