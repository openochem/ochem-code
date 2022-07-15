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

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import qspr.Globals;
import qspr.ThreadScope;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Property;
import qspr.entities.User;
import qspr.metaserver.EventListener;

public class DataApprovalEmailNotifier 
{
	public Map<Long, AffectedRecordsInfo> affectedRecordsInfoMap = new HashMap<Long, AffectedRecordsInfo>();
	public User executingUser;

	public DataApprovalEmailNotifier()
	{
		executingUser = Globals.userSession().user;
		attachListeners();
	}

	public void sendEmail()
	{
		//		for (AffectedRecordsInfo info : affectedRecordsInfoMap.values()) {
		//			info.recordIDs.addAll(info.pendingRecordIDs);
		//			info.pendingRecordIDs.clear();
		//			Mailer.postMailAsync(info.user.email, "Your data has been approved", "Dear " + info.user.login + ",\n\n Your " + info.recordIDs.size() + " records for properties " + info.properties + " have been approved by " + executingUser.login + "\n" +
		//			"You can view the data introduced by you at the following link: "+OCHEMConfiguration.getRootHost()+OCHEMConfiguration.rootDir+"/epbrowser/show.do?introducer="+info.user.id+"\n\n" + 
		//			"Sincerely yours,\nOCHEM Team");
		//		}

	}

	@SuppressWarnings({ "rawtypes", "unchecked" })
	public void attachListeners()
	{
		ThreadScope.get().transactionCommit.addListener(new EventListener() 
		{
			@Override
			public void onEvent(qspr.metaserver.Event event, Object arg) {
				for (AffectedRecordsInfo info : affectedRecordsInfoMap.values()) {
					info.recordIDs.addAll(info.pendingRecordIDs);
					info.pendingRecordIDs.clear();
				}

			}
		});

		ThreadScope.get().recordApproved.addListener(new EventListener<ExperimentalProperty>() 
				{
			@Override
			public void onEvent(qspr.metaserver.Event event, ExperimentalProperty arg) {
				if (arg.introducer != executingUser)
					getInfo(arg.introducer).addRecord(arg);
			}
				});

		ThreadScope.get().transactionRollback.addListener(new EventListener() 
		{
			@Override
			public void onEvent(qspr.metaserver.Event event, Object arg) {
				for (AffectedRecordsInfo info : affectedRecordsInfoMap.values()) {
					info.pendingRecordIDs.clear();
				}
			}
		});

		ThreadScope.get().threadFinished.addListener(new EventListener() 
		{
			@Override
			public void onEvent(qspr.metaserver.Event event, Object arg) {
				sendEmail();
			}
		});
	}

	private AffectedRecordsInfo getInfo(User user)
	{
		AffectedRecordsInfo info = affectedRecordsInfoMap.get(user.id);
		if (info == null)
			affectedRecordsInfoMap.put(user.id, info = new AffectedRecordsInfo(user));

		return info;

	}
}

class AffectedRecordsInfo
{
	public final Set<Long> recordIDs = new HashSet<Long>();
	public final Set<Long> pendingRecordIDs = new HashSet<Long>();
	public final Set<Property> properties = new HashSet<Property>();
	public User user;

	public AffectedRecordsInfo(User user)
	{
		this.user = user;
	}

	public void addRecord(ExperimentalProperty ep)
	{
		if (ep.property != null)
			properties.add(ep.property);
		pendingRecordIDs.add(ep.id);
	}
}
