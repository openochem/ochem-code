import numpy as np

class hive ():
    

    def func (self, vec):
        key = tuple(vec)
        if self.vals.get(key) is None:
            self.vals[key] = self.fn(vec)
            self.fn_counter += 1
        return self.vals[key]
        

    def iter (self):
        self.prev_points = self.points.copy()
        self.points = self.best_points(np.vstack((self.random_search(), self.informed_search(self.points))))
        
    def random_search (self):
        return np.random.rand(self.g_density, len(self.points[0]))*(self.g_range[:,1]-self.g_range[:,0])+self.g_range[:,0]

    def informed_search (self, points):
        cur_density = self.density
        inc = cur_density//len(points)
        vals = np.zeros(len(points))
        for i in range(len(points)):
            cur_density += inc*i
            candidates = self.explore_point(points[i], cur_density)
            points[i] = candidates[np.argmin(map(self.func, candidates))]
        # remove duplicates
        new_points = []
        vals = []
        for p in points:
            v = self.func(p)
            if v not in vals:
                vals.append(self.func(p))
                new_points.append(p)
        return np.array(new_points)[np.argsort(vals)]

    def explore_point (self, point, cur_density):
        points = np.zeros((cur_density+1, len(point)), np.float)
        for i in range(cur_density):
            points[i] = point + np.random.rand(len(point))*self.range*2 - self.range
        points[-1] = point
        return points
        
    def best_points (self, points):
        vals = map(self.func, points)
        sorted_points = np.array(points)[np.argsort(vals)]
        return sorted_points[:self.ngood]
    
    def mod_params (self, j):
        for i in range(len(self.range)):
            if self.range[i] > abs(self.prev_points[0][i]-self.points[0][i]):
                self.range[i] /= 1.5
        self.g_density -= int(round(self.g_density*((j+1)/float(self.maxiter))**2))

    def print_stats (self, j):
        print ("Iter %d - best: %f, range: %s, global_density: %f" %\
          (j, self.func(self.points[0]), np.array2string(self.range), self.g_density))
        
    def run (self):
        for i in range(self.maxiter):
            self.iter()
            if i % int(self.maxiter/10) == 0:
                self.print_stats(i)
                self.mod_params(i)
        return self.points[0]

    def __init__ (self, fn, x0, g_range, g_density=10, range=None, density=3, ngood=None, maxiter=10):
        self.g_range = np.array(g_range, np.float)
        self.g_density = g_density
        self.range = np.array(range, np.float) if range is not None else (self.g_range[:,1]-self.g_range[:,0])*0.1
        self.density = density
        self.ngood = ngood if ngood is not None else self.density+2
        self.maxiter = maxiter if maxiter > 10 else 10
        self.points = np.array([x0], np.float)
        self.vals = {}
        self.prev_points = None
        self.fn = fn
        self.fn_counter = 0
        
def griewank (x):
    sum=0
    prod=1
    for i in range(len(x)):
        sum += x[i]**2
        prod *= np.cos(x[i]/np.sqrt(i+1))
    return -(1/(sum/4000. - prod + 1.1))

#H = hive(griewank, [4,4], g_range=np.tile([-10,10],(2,1)), maxiter=100)
#print (H.run(), H.fn_counter)
