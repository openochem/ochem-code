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

import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.List;

import com.eadmet.exceptions.UserFriendlyException;

public class DynaWrap {

    private Object obj;

    public DynaWrap(Object obj) {
        this.obj = obj;
    }
    
    public void setField(String name, Object obj) {
    	try {
			Field field = this.obj.getClass().getField(name);
			field.set(this.obj, obj);
		} catch (NoSuchFieldException | SecurityException | IllegalArgumentException | IllegalAccessException e) {
			throw new UserFriendlyException(e);
		}
    }

    public Object get(String path) {
        String[] steps = path.split("\\.");
        Object result = this.obj;
        try {
            for (String step : steps) {
                try {
                	Method method = result.getClass().getMethod(createGetterName(step));
                    result = method.invoke(result);
                } catch (NoSuchMethodException exp) {
                	Field field = result.getClass().getField(step);
                	result = field.get(result);
                }
            }
            return result;
        }
        catch (Exception e) {
            throw new RuntimeException(e);
        }
    }
    
    public String discoverClassName() {
        return this.obj.getClass().getName();
    }
    
    public List discoverProperties() {
        return discoverProperties(this.obj);
    }
    
    private List<String> discoverProperties(Object obj) {
        Method[] methods = obj.getClass().getDeclaredMethods();
        List<String> result = new ArrayList<>();
        for (Method method : methods) {
            if (method.getName().startsWith("get") && method.getParameters().length == 0) {
                String propName = createPropertyName(method.getName());
                result.add(propName);
                try {
                    if (method.invoke(obj) instanceof Object) {
                        Object subObj = method.invoke(obj);
                        List<String> subProps = discoverProperties(subObj);
                        for (String subProp : subProps) {
                            result.add(propName + "." + subProp);
                        }
                    }
                }
                catch (Exception e) {
                    // ignore
                }
            }
        }
        return result;
    }

    private String createGetterName(String name) {
        StringBuilder sb = new StringBuilder("get");
        sb.append(name.substring(0, 1).toUpperCase());
        sb.append(name.substring(1));
        return sb.toString();
    }
    
    private String createPropertyName(String name) {
        return name.substring(3, 4).toLowerCase() + name.substring(4);
    }

    
    // convenience methods to avoid casting
    public boolean getBoolean(String path) {
        return (Boolean)get(path);
    }
    public byte getByte(String path) {
        return (Byte)get(path);
    }
    public short getShort(String path) {
        return (Short)get(path);
    }
    public int getInt(String path) {
        return (Integer)get(path);
    }
    public long getLong(String path) {
        return (Long)get(path);
    }
    public float getFloat(String path) {
        return (Float)get(path);
    }
    public double getDouble(String path) {
        return (Double)get(path);
    }
    public String getString(String path) {
        return (String)get(path);
    }

}
