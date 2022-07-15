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

package com.eadmet.utils;

import java.sql.Timestamp;
import java.util.Calendar;

/**
 * A class to nicely display dates and intervals
 * @author midnighter
 */
public class TimeUtils {

    public final static long ONE_SECOND = 1000;
    public final static long SECONDS = 60;

    public final static long ONE_MINUTE = ONE_SECOND * 60;
    public final static long MINUTES = 60;

    public final static long ONE_HOUR = ONE_MINUTE * 60;
    public final static long HOURS = 24;

    public final static long ONE_DAY = ONE_HOUR * 24;
    public final static long ONE_WEEK = ONE_DAY * 7;
    public final static long ONE_MONTH = ONE_DAY * 30;

    private TimeUtils() {
    }

    // Nice time display // Midnighter on Oct 11, 2011
    public static String niceTimeDifference(long duration) 
    {
      if (duration <= 31 * ONE_SECOND)
    	  return "seconds";
      else if (duration <= 1.51 * ONE_MINUTE)
    	  return "about a minute";
      else if (duration <= 50 * ONE_MINUTE)
    	  return ""  + Math.round(duration / ONE_MINUTE) + " minutes";
      else if (duration <= 1.5 * ONE_HOUR)
    	  return "about an hour";
      else if (duration <= 20 * ONE_HOUR)
    	  return "about " + Math.round(duration / ONE_HOUR) + " hours";
      else if (duration <= 1.5 * ONE_DAY)
    	  return "about a day";
      else if (duration <= 6 * ONE_DAY)
    	  return "" +  Math.round(duration / ONE_DAY) + " hours";
      else if (duration <= 9 * ONE_DAY)
    	  return "about a week";
      else if (duration <= 3.5 * ONE_WEEK)
    	  return "several weeks";
      else if (duration <= 6 * ONE_WEEK)
    	  return "about a month";
      else if (duration <= 10 * ONE_MONTH)
    	  return "" +  Math.round(duration / ONE_MONTH) + " months";
      else if (duration <= 14 * ONE_MONTH)
    	  return "about a year";
      else
    	  return "more than a year";
    }
    
    public static String ago(Timestamp when)
    {
    	return niceTimeDifference(Calendar.getInstance().getTimeInMillis() - when.getTime()) + " ago";
    }
}
