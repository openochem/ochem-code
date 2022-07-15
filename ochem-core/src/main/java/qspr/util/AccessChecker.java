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

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Criteria;
import org.hibernate.criterion.Conjunction;
import org.hibernate.criterion.Disjunction;
import org.hibernate.criterion.Expression;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.OCHEMConfiguration;
import qspr.entities.Property;
import qspr.entities.User;
import qspr.metaserver.protocol.Task.TaskPriority;

import com.eadmet.exceptions.UserFriendlyException;

// The class defines various global access rules: 
// The logic who can export Dragon descriptors, export Corina structures, etc.
// Later, make something universal out of it / Midnighter on Apr 18, 2011
public class AccessChecker
{

	static String DESCRIPTORS = "</descriptors>";

	private static transient final Logger logger = LogManager.getLogger(AccessChecker.class);

	public static boolean canExportDescriptors(User user, String description) {
		if (user == null) return false;
		if (description == null) return true;
		if (user.isOCHEMDeveloper()) return true;
		if (OCHEMConfiguration.inhouseInstallation) return true;

		description = description.toLowerCase();
		if(description.contains(DESCRIPTORS))description=description.substring(0, description.indexOf(DESCRIPTORS));
		return true;
	}

	public static boolean canExportCorina(User user)
	{
		return user != null && user.isOCHEMDeveloper();
	}

	public static boolean isFromGroup(User owner, User accessingUser) {
		return owner != null && accessingUser != null && accessingUser.group != null && accessingUser.group.equals(owner.group);
	}

	public static boolean isModerator(User user)
	{
		if(user == null) return false;
		if(user.isSuperUser()) return true;
		if(user.isOCHEMDeveloper()) return true;
		return ((Long) Globals.session().createCriteria(Property.class).add(Restrictions.eq("moderator", user)).setProjection(Projections.count("id")).uniqueResult() > 0); 
	}

	public static int getMaximumTaskPriority(User user)
	{
		if (user != null && user.isOCHEMDeveloper()) 
			return TaskPriority.EXTRA_HIGH;
		return TaskPriority.HIGH;  // the current idea is to limit number of high tasks per users
	}

	/**
	 * Checks whether entity can be added to basket for viewing
	 * @param entity
	 * @param requestedRights
	 * @return
	 */
	public static boolean requestAddToBasket(UserContributedEntity entity)
	{
		if (entity.getRights() == Globals.RIGHTS_FREELY_AVAILABLE) return true;
		User accessingUser = Globals.userSession().user;
		if (accessingUser == null) return false;

		if (entity.getOwner().equals(accessingUser) 
				|| entity.getIntroducer().equals(accessingUser)
				|| entity.getOwner().group != null && entity.getOwner().group.equals(accessingUser.group)
				|| entity.getIntroducer().group != null && entity.getIntroducer().group.equals(accessingUser.group)
				)
			return true;

		return false;
	}


	public static void requestPermission(UserContributedEntity entity, Integer requestedRights)
	{
		User accessingUser = Globals.userSession().user;
		Integer rights = entity.getRights();
		if (rights == null)
			rights = Globals.RIGHTS_FREELY_AVAILABLE;

		if (entity.getOwner() == null)
			return;

		if (accessingUser != null)
		{
			if (entity.getOwner().equals(accessingUser))
				return;

			if (entity.getIntroducer().equals(accessingUser) && !entity.getOwner().isPublisher()) // For published data we do not allow to change them
				return;

			if (accessingUser.rank >= entity.getOwner().rank)
				return;

			if (entity.getOwner().group != null && entity.getOwner().group.equals(accessingUser.group))
				return;

			if (accessingUser.rank > entity.getOwner().rank && rights >= requestedRights)
				return;
		}

		// If all possibilities for record edit failed - that it's not permitted after all
		logger.info("You do not have sufficient privileges to modify this entity (introduced by "+entity.getOwner().login+")");
		throw new UserFriendlyException("You do not have sufficient privileges to modify this entity (introduced by "+entity.getOwner().login+")");
	}

	/**
	 * Check if the current user is allowed to modify the entity, throw a runtime exception otherwise
	 * @param entity
	 */
	public static void requestModificationPermission(UserContributedEntity entity)
	{
		requestPermission(entity, Globals.RIGHTS_FREELY_AVAILABLE);
	}

	/**
	 * Check if the current user is allowed to view the entity, throw a runtime exception otherwise
	 * @param entity
	 */
	public static void requestViewingPermission(UserContributedEntity entity)
	{
		requestPermission(entity, Globals.RIGHTS_FREELY_AVAILABLE);
	}

	public static void requestRegisteredUserPrivileges()
	{
		if (Globals.userSession().user == null)
			throw new UserFriendlyException("You need to login to perform this action: sorry, but because of an abuse use we had to restrict an access to this option by Guest users.");
	}

	public static void requestValidatedPrivileges()
	{
		if (!Globals.isValidatedUser())
			throw new UserFriendlyException("You need a validated account to perform this action!");
	}

	public static void requestSuperuserPrivileges()
	{
		if (Globals.userSession().user == null || !Globals.userSession().user.isSuperUser())
			throw new UserFriendlyException("You need superuser privileges to perform this action!");
	}

	public static void requestModeratorPrivileges()
	{
		if (Globals.isSuperUser())
			return;

		// So far, its not clear who is a "moderator". Let it be a super-user or any container user.
		User user = Globals.userSession().user;
		if (user == null || !AccessChecker.isModerator(user))
			throw new UserFriendlyException("You need moderator privileges to perform this action!");
	}

	public static void addAccessRestrictions(Criteria criteria, int requiredPrivileges, User accessingUser, Disjunction disjunctionRestriction, boolean showUnapproved)
	{
		// Apply user rights restriction to a given criteria
		// Midnighter
		User user = Globals.userSession() == null ? null : Globals.userSession().user;
		if (accessingUser != null)
			user = accessingUser;
		criteria.createAlias("introducer", "intr", Criteria.LEFT_JOIN);
		criteria.createAlias("owner", "own", Criteria.LEFT_JOIN);

		if (disjunctionRestriction == null)
			disjunctionRestriction = Restrictions.disjunction();

		if (showUnapproved)
			disjunctionRestriction.add(Restrictions.isNull("intr.id"));

		Conjunction conjunction = Restrictions.conjunction();
		conjunction.add(Expression.ge("rights", requiredPrivileges));
		if (!showUnapproved)
			conjunction.add(Restrictions.eq("approved", Boolean.TRUE)); // the record is approved by administrator
		disjunctionRestriction.add(conjunction);


		if (user != null)
		{

			// Show items that I own or that are owned by a significantly lower-privileged user (difference in ranks >= 2)
			if (showUnapproved)
				disjunctionRestriction.add(Expression.lt("own.rank", user.rank - 1));

			// Show records, owned by me
			disjunctionRestriction.add(Restrictions.eq("own.id", user.id));
			disjunctionRestriction.add(Restrictions.eq("intr.id", user.id));
			if (user.group != null)
			{
				// Show records from my group
				criteria.createAlias("intr.group", "grp", Criteria.LEFT_JOIN);
				disjunctionRestriction.add(Restrictions.eq("grp.id", user.group.id));
			}
		}

		criteria.add(disjunctionRestriction);
	}

	public static String explainUser(User checkUser) {
		User user = Globals.userSession().user;
		if(user.isSuperUser() || user.isOCHEMDeveloper()) 
			return " (user: " + (checkUser == null ?"null" : checkUser.login)+ ")";
		else
			return "";
	}
}
