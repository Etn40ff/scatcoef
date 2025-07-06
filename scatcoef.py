from sage.plot.colors import ColorsDict
from collections import defaultdict

class vars(defaultdict):

    def __init__(self, name):
        def index(i):
            idx = hash(i)
            if idx < 0:
                return -2*idx+1
            else:
                return 2*idx
        super().__init__(lambda i: SR.var(name+'_{}'.format(index(i)),latex_name=name+'_{'+'{},{}'.format(i[0],i[1])+'}'))

    def __missing__(self, key):
        if self.default_factory is not None:
            self[key] = self.default_factory(key)
            return self[key]
        else:
            raise KeyError(key)

def _gcd_eval(self, a, b):
    if a == 0:
        return b
    elif b == 0:
        return a
    elif a in ZZ and b in ZZ:
        return gcd(a,b)
    a_factors = list(factor(a)) if a in ZZ else a.factor_list()
    b_factors = list(factor(b)) if b in ZZ else b.factor_list()
    symbolic_gcd = 1
    for (x,px) in a_factors:
        for (y,py) in b_factors:
            if x==y and x not in ZZ:
                symbolic_gcd *= x**min(px,py)
            elif x in ZZ and y in ZZ:
                symbolic_gcd *= gcd(x**px, y**py)
    if symbolic_gcd == 1:
        pass
    else:
        aa = a/symbolic_gcd
        bb = b/symbolic_gcd
        return symbolic_gcd * self(aa,bb)

my_gcd = sage.symbolic.function_factory.function('gcd',nargs=2, eval_func=_gcd_eval)

class ScatteringCoefficients(SageObject):
    
    def __init__(self, a, b):
        self.a = a
        self.b = b
        self.x0 = SR('x0')
        self.x1 = SR('x1')
        self.y0 = SR('y0')
        self.y1 = SR('y1')
        self.c = vars('c')
        self.ck = {}
        self.c[0,0] = 1
        self.c[1,0] = 1
        self.c[0,1] = 1
        self._explored_depth=0
        g = my_gcd(a,b)
        self.hyperbola = 2*a/g*self.x0**2 -2*a*b/g*self.x0*self.x1 + 2*b/g*self.x1**2

    def load_from_savefile(self, savefile=None):
        a = SR('a')
        b = SR('b')
        # HACK: avoid segfaults due to #30018
        gcd = sage.symbolic.function_factory.function('gcd',nargs=2)(a,b)
        
        if not savefile:
            savefile = "precomputed_c.sobj"
        depth, precomputed_c = load(savefile)

        self.c.update([(key,value.subs({a:self.a,b:self.b}) if value not in QQ else value) for key,value in precomputed_c])
        self._explored_depth = max(depth,self._explored_depth)

    def save_to_savefile(self, savefile=None, force=False):
        if not (self.a in SR and self.b in SR) and not force:
            raise ValueError("You are trying to save a numerical computation, If sure force it with force=True")
        if not savefile:
            savefile = "precomputed_c.sobj"
        save((self._explored_depth, list(self.c.items())),savefile)

    def wall_hom(self, normal, f, degree):
        normal = vector(normal)
        if normal == vector((1,0)):
            scat_term = 1+self.c[(1,0)]*self.y0
        elif normal == vector((0,1)):
            scat_term = 1+self.c[(0,1)]*self.y1
        else:
            scat_term = 1 + sum( self.c[r,s]*self.y0**r*self.y1**s for (r,s) in [ i*normal for i in range(1,ceil(degree/(sum(normal)))+1) ] )
       
        e = vector((self.a*normal[0], self.b*normal[1]))
        e = e/my_gcd(*e)

        return taylor(f.subs({
            self.x0 : self.x0*scat_term**(-e[0]),
            self.x1 : self.x1*scat_term**(-e[1]),
            self.y0 : self.y0*scat_term**(self.a*e[1]),
            self.y1 : self.y1*scat_term**(-self.b*e[0])
            }), (self.y0,0), (self.y1,0), degree)

    def gamma0(self, f, degree):
        return self.wall_hom( (0,1), self.wall_hom((1,0), f, degree), degree)

    def wall_list(self, degree):
        walls = [ (1,0), (0,1) ] + [ (i,j) for i in range(1,degree+1) for j in range(1,degree-i+1) if gcd(i,j)==1 ]
        return sorted(walls, key=lambda i: Infinity if i[1] == 0 else i[0]/i[1])

    def gamma1(self, f, degree):
        for wall in self.wall_list(degree):
            f = self.wall_hom(wall, f, degree)
        return f

    def compute_scattering_coefficients(self, degree):
        if self._explored_depth >= degree:
            return

        if degree > 1:
            self.compute_scattering_coefficients(degree-1)

        lhs0 = self.gamma0(self.x0**(-1), degree)
        lhs1 = self.gamma0(self.x1**(-1), degree)
        rhs0 = self.gamma1(self.x0**(-1), degree)
        rhs1 = self.gamma1(self.x1**(-1), degree)
        
        variables = []
        equations = []
        for i in range(1, degree):
            variables.append( self.c[i,degree-i] )
            equations.append( lhs0.coefficient(self.y0**i).coefficient(self.y1**(degree-i)) - rhs0.coefficient(self.y0**i).coefficient(self.y1**(degree-i)) )
            equations.append( lhs1.coefficient(self.y0**i).coefficient(self.y1**(degree-i)) - rhs1.coefficient(self.y0**i).coefficient(self.y1**(degree-i)) )

        solution = solve(equations,variables,solution_dict=True)[0]

        for i in range(1, degree):
            self.c[i,degree-i] = self.c[i,degree-i].subs(solution)

        if degree > 1:
            self.c[degree,0] = self.c[0,degree] = 0

        self._explored_depth = degree

    def polynomials(self, i, j):
        gij = gcd(i,j)
        g = my_gcd(i/gij*a, j/gij*b)
        return [ f.coefficient(g**k).factor() for f in [self.c[(i,j)].collect(g)] for k in range(gij+1)]

    def compute_polynomials(self, degree):
        self.compute_scattering_coefficients(degree)
        for i in range(1,degree):
            for j in range(1,degree-i):
                self.ck[(i,j)] = self.polynomials(i,j)

    def show_k(self, degree, k, a, b,  color_numbers=False, fontsize='medium'):
        self.compute_polynomials(degree)
        G = Graphics()
        
        for i in range(1, degree):
            for j in range(1, degree-i):
                G += text(self.ck[(i,j)][k] if len(self.ck[(i,j)]) > k else 0, (i,j+0.2*(i%3) ), color=colors_dict[lengths_dict[i,j]] if color_numbers else 'black', fontsize=fontsize)
        G.set_aspect_ratio(1)
        G.SHOW_OPTIONS['gridlines'] = [list(range(degree))]*2
        return G


    def show(self, degree, color_numbers=False, fontsize='medium'):
        self.compute_scattering_coefficients(degree)
        G = Graphics()
        
        if self.a in ZZ and self.b in ZZ:
            # plot the region where representatives of roots are located
            G += Polyhedron(ieqs=[(0,-2,self.b),(0,self.a,-2),(degree+0.5,-1,0), (degree+0.5,0,1)]).plot(wireframe=False,fill='lightgreen')

            lengths_dict = {}
            for i in range(degree):
                for j in range(degree-i):
                    if self.c[i,j] in ZZ and self.c[i,j] != 0:
                        lengths_dict[i,j] = self.hyperbola.subs({self.x0:i,self.x1:j})

            lengths = set(lengths_dict.values())
            colors = [ sage.plot.colors.Color(rgb) for rgb in rainbow(len(lengths), format='rgbtuple') ]
            colors_dict = { j:colors[i] for i,j in enumerate(lengths) } 
        
            for l in lengths:
                G += implicit_plot( self.hyperbola-l, (self.x0, -0.5, degree+0.5), (self.x1, -0.5, degree+0.5), color=colors_dict[l].lighter(0.5) )

            for (i,j) in lengths_dict:
                G += text(self.c[(i,j)], (i,j), color=colors_dict[lengths_dict[i,j]] if color_numbers else 'black', fontsize=fontsize)
            
            G.set_axes_range(-0.5, max( i for i,_ in lengths_dict )+0.5, -0.5, max( j for _,j in lengths_dict )+0.5 )

        else:
            for i in range(degree):
                for j in range(degree-i):
                    if self.c[(i,j)] != 0:
                        G += text('$'+latex(factor(self.c[(i,j)]))+'$', (i,j), color='blue', fontsize=fontsize)

        G.set_aspect_ratio(1)
        G.SHOW_OPTIONS['gridlines'] = [list(range(degree))]*2
        return G

    def squared_root_length(self, alpha):
        return self.hyperbola.subs({self.x0:alpha[0], self.x1:alpha[1]})

    def primitive_imaginary(self, degree):
        return [ (i,j) for i in range(1,degree) for j in range(1,degree-i+1) if gcd(i,j) == 1 and self.squared_root_length((i,j)) < 0 ]

    def oeis(self, degree):
        out = {}
        self.compute_scattering_coefficients(degree)
        primitive = [ (i,j) for (i,j) in self.c if self.c[i,j] in ZZ and self.c[i,j] != 0 and gcd(i,j) == 1 and 0 not in (i,j)]
        first_points = {}
        for (i,j) in primitive:
            srl = self.squared_root_length((i,j))
            if sum(first_points.get(srl, [Infinity,Infinity])) > i+j:
                first_points[srl] = (i,j)
        for (i,j) in first_points.values():
            k = 0
            seq = []
            while self.c[k*i,k*j] in ZZ:
                seq.append(self.c[k*i,k*j])
                k+=1
            if len([s for s in seq if s!=0 ]) > 2:
                out[(i,j)] = (seq, len(oeis(seq)))
        return out
