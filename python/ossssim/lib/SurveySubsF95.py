from __future__ import print_function, absolute_import, division
from . import _SurveySubsF95
import f90wrap.runtime
import logging

class Poly_Dec(f90wrap.runtime.FortranModule):
    """
    Module poly_dec
    
    
    Defined at poly_dec.fpp lines 5-11
    
    """
    @f90wrap.runtime.register_class("SurveySubsF95.t_polygon")
    class t_polygon(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=t_polygon)
        
        
        Defined at poly_dec.fpp lines 9-11
        
        """
        def __init__(self, handle=None):
            """
            self = T_Polygon()
            
            
            Defined at poly_dec.fpp lines 9-11
            
            
            Returns
            -------
            this : T_Polygon
            	Object to be constructed
            
            
            Automatically generated constructor for t_polygon
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _SurveySubsF95.f90wrap_t_polygon_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class T_Polygon
            
            
            Defined at poly_dec.fpp lines 9-11
            
            Parameters
            ----------
            this : T_Polygon
            	Object to be destructed
            
            
            Automatically generated destructor for t_polygon
            """
            if self._alloc:
                _SurveySubsF95.f90wrap_t_polygon_finalise(this=self._handle)
        
        @property
        def n(self):
            """
            Element n ftype=integer  pytype=int
            
            
            Defined at poly_dec.fpp line 10
            
            """
            return _SurveySubsF95.f90wrap_t_polygon__get__n(self._handle)
        
        @n.setter
        def n(self, n):
            _SurveySubsF95.f90wrap_t_polygon__set__n(self._handle, n)
        
        @property
        def x(self):
            """
            Element x ftype=real(kind=8) pytype=float
            
            
            Defined at poly_dec.fpp line 11
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _SurveySubsF95.f90wrap_t_polygon__array__x(self._handle)
            if array_handle in self._arrays:
                x = self._arrays[array_handle]
            else:
                x = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _SurveySubsF95.f90wrap_t_polygon__array__x)
                self._arrays[array_handle] = x
            return x
        
        @x.setter
        def x(self, x):
            self.x[...] = x
        
        @property
        def y(self):
            """
            Element y ftype=real(kind=8) pytype=float
            
            
            Defined at poly_dec.fpp line 11
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _SurveySubsF95.f90wrap_t_polygon__array__y(self._handle)
            if array_handle in self._arrays:
                y = self._arrays[array_handle]
            else:
                y = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _SurveySubsF95.f90wrap_t_polygon__array__y)
                self._arrays[array_handle] = y
            return y
        
        @y.setter
        def y(self, y):
            self.y[...] = y
        
        def __str__(self):
            ret = ['<t_polygon>{\n']
            ret.append('    n : ')
            ret.append(repr(self.n))
            ret.append(',\n    x : ')
            ret.append(repr(self.x))
            ret.append(',\n    y : ')
            ret.append(repr(self.y))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @property
    def n_e_max(self):
        """
        Element n_e_max ftype=integer pytype=int
        
        
        Defined at poly_dec.fpp line 7
        
        """
        return _SurveySubsF95.f90wrap_poly_dec__get__n_e_max()
    
    def __str__(self):
        ret = ['<poly_dec>{\n']
        ret.append('    n_e_max : ')
        ret.append(repr(self.n_e_max))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

poly_dec = Poly_Dec()

class Poly_Lib(f90wrap.runtime.FortranModule):
    """
    Module poly_lib
    
    
    Defined at poly_lib.fpp lines 5-208
    
    """
    @staticmethod
    def point_in_polygon(p, poly):
        """
        point_in_polygon = point_in_polygon(p, poly)
        
        
        Defined at poly_lib.fpp lines 22-86
        
        Parameters
        ----------
        p : float array
        poly : T_Polygon
        
        Returns
        -------
        point_in_polygon : int
        
        """
        point_in_polygon = _SurveySubsF95.f90wrap_point_in_polygon(p=p, \
            poly=poly._handle)
        return point_in_polygon
    
    @staticmethod
    def calc_walk_summand(p1, p2):
        """
        calc_walk_summand = calc_walk_summand(p1, p2)
        
        
        Defined at poly_lib.fpp lines 88-167
        
        Parameters
        ----------
        p1 : float array
        p2 : float array
        
        Returns
        -------
        calc_walk_summand : int
        
        """
        calc_walk_summand = _SurveySubsF95.f90wrap_calc_walk_summand(p1=p1, p2=p2)
        return calc_walk_summand
    
    @staticmethod
    def check_polygon(self):
        """
        check_polygon(self)
        
        
        Defined at poly_lib.fpp lines 169-208
        
        Parameters
        ----------
        poly : T_Polygon
        
        """
        _SurveySubsF95.f90wrap_check_polygon(poly=self._handle)
    
    _dt_array_initialisers = []
    

poly_lib = Poly_Lib()

class Datadec(f90wrap.runtime.FortranModule):
    """
    Module datadec
    
    
    Defined at datadec.fpp lines 5-53
    
    """
    @f90wrap.runtime.register_class("SurveySubsF95.t_ratecut")
    class t_ratecut(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=t_ratecut)
        
        
        Defined at datadec.fpp lines 20-21
        
        """
        def __init__(self, handle=None):
            """
            self = T_Ratecut()
            
            
            Defined at datadec.fpp lines 20-21
            
            
            Returns
            -------
            this : T_Ratecut
            	Object to be constructed
            
            
            Automatically generated constructor for t_ratecut
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _SurveySubsF95.f90wrap_t_ratecut_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class T_Ratecut
            
            
            Defined at datadec.fpp lines 20-21
            
            Parameters
            ----------
            this : T_Ratecut
            	Object to be destructed
            
            
            Automatically generated destructor for t_ratecut
            """
            if self._alloc:
                _SurveySubsF95.f90wrap_t_ratecut_finalise(this=self._handle)
        
        @property
        def min_bn(self):
            """
            Element min_bn ftype=real(kind=8) pytype=float
            
            
            Defined at datadec.fpp line 21
            
            """
            return _SurveySubsF95.f90wrap_t_ratecut__get__min_bn(self._handle)
        
        @min_bn.setter
        def min_bn(self, min_bn):
            _SurveySubsF95.f90wrap_t_ratecut__set__min_bn(self._handle, min_bn)
        
        @property
        def max_bn(self):
            """
            Element max_bn ftype=real(kind=8) pytype=float
            
            
            Defined at datadec.fpp line 21
            
            """
            return _SurveySubsF95.f90wrap_t_ratecut__get__max_bn(self._handle)
        
        @max_bn.setter
        def max_bn(self, max_bn):
            _SurveySubsF95.f90wrap_t_ratecut__set__max_bn(self._handle, max_bn)
        
        @property
        def angle(self):
            """
            Element angle ftype=real(kind=8) pytype=float
            
            
            Defined at datadec.fpp line 21
            
            """
            return _SurveySubsF95.f90wrap_t_ratecut__get__angle(self._handle)
        
        @angle.setter
        def angle(self, angle):
            _SurveySubsF95.f90wrap_t_ratecut__set__angle(self._handle, angle)
        
        @property
        def hwidth(self):
            """
            Element hwidth ftype=real(kind=8) pytype=float
            
            
            Defined at datadec.fpp line 21
            
            """
            return _SurveySubsF95.f90wrap_t_ratecut__get__hwidth(self._handle)
        
        @hwidth.setter
        def hwidth(self, hwidth):
            _SurveySubsF95.f90wrap_t_ratecut__set__hwidth(self._handle, hwidth)
        
        def __str__(self):
            ret = ['<t_ratecut>{\n']
            ret.append('    min_bn : ')
            ret.append(repr(self.min_bn))
            ret.append(',\n    max_bn : ')
            ret.append(repr(self.max_bn))
            ret.append(',\n    angle : ')
            ret.append(repr(self.angle))
            ret.append(',\n    hwidth : ')
            ret.append(repr(self.hwidth))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @f90wrap.runtime.register_class("SurveySubsF95.t_orb_m")
    class t_orb_m(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=t_orb_m)
        
        
        Defined at datadec.fpp lines 23-24
        
        """
        def __init__(self, handle=None):
            """
            self = T_Orb_M()
            
            
            Defined at datadec.fpp lines 23-24
            
            
            Returns
            -------
            this : T_Orb_M
            	Object to be constructed
            
            
            Automatically generated constructor for t_orb_m
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _SurveySubsF95.f90wrap_t_orb_m_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class T_Orb_M
            
            
            Defined at datadec.fpp lines 23-24
            
            Parameters
            ----------
            this : T_Orb_M
            	Object to be destructed
            
            
            Automatically generated destructor for t_orb_m
            """
            if self._alloc:
                _SurveySubsF95.f90wrap_t_orb_m_finalise(this=self._handle)
        
        @property
        def a(self):
            """
            Element a ftype=real(kind=8) pytype=float
            
            
            Defined at datadec.fpp line 24
            
            """
            return _SurveySubsF95.f90wrap_t_orb_m__get__a(self._handle)
        
        @a.setter
        def a(self, a):
            _SurveySubsF95.f90wrap_t_orb_m__set__a(self._handle, a)
        
        @property
        def e(self):
            """
            Element e ftype=real(kind=8) pytype=float
            
            
            Defined at datadec.fpp line 24
            
            """
            return _SurveySubsF95.f90wrap_t_orb_m__get__e(self._handle)
        
        @e.setter
        def e(self, e):
            _SurveySubsF95.f90wrap_t_orb_m__set__e(self._handle, e)
        
        @property
        def inc(self):
            """
            Element inc ftype=real(kind=8) pytype=float
            
            
            Defined at datadec.fpp line 24
            
            """
            return _SurveySubsF95.f90wrap_t_orb_m__get__inc(self._handle)
        
        @inc.setter
        def inc(self, inc):
            _SurveySubsF95.f90wrap_t_orb_m__set__inc(self._handle, inc)
        
        @property
        def node(self):
            """
            Element node ftype=real(kind=8) pytype=float
            
            
            Defined at datadec.fpp line 24
            
            """
            return _SurveySubsF95.f90wrap_t_orb_m__get__node(self._handle)
        
        @node.setter
        def node(self, node):
            _SurveySubsF95.f90wrap_t_orb_m__set__node(self._handle, node)
        
        @property
        def peri(self):
            """
            Element peri ftype=real(kind=8) pytype=float
            
            
            Defined at datadec.fpp line 24
            
            """
            return _SurveySubsF95.f90wrap_t_orb_m__get__peri(self._handle)
        
        @peri.setter
        def peri(self, peri):
            _SurveySubsF95.f90wrap_t_orb_m__set__peri(self._handle, peri)
        
        @property
        def m(self):
            """
            Element m ftype=real(kind=8) pytype=float
            
            
            Defined at datadec.fpp line 24
            
            """
            return _SurveySubsF95.f90wrap_t_orb_m__get__m(self._handle)
        
        @m.setter
        def m(self, m):
            _SurveySubsF95.f90wrap_t_orb_m__set__m(self._handle, m)
        
        def __str__(self):
            ret = ['<t_orb_m>{\n']
            ret.append('    a : ')
            ret.append(repr(self.a))
            ret.append(',\n    e : ')
            ret.append(repr(self.e))
            ret.append(',\n    inc : ')
            ret.append(repr(self.inc))
            ret.append(',\n    node : ')
            ret.append(repr(self.node))
            ret.append(',\n    peri : ')
            ret.append(repr(self.peri))
            ret.append(',\n    m : ')
            ret.append(repr(self.m))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @f90wrap.runtime.register_class("SurveySubsF95.t_orb_p")
    class t_orb_p(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=t_orb_p)
        
        
        Defined at datadec.fpp lines 26-27
        
        """
        def __init__(self, handle=None):
            """
            self = T_Orb_P()
            
            
            Defined at datadec.fpp lines 26-27
            
            
            Returns
            -------
            this : T_Orb_P
            	Object to be constructed
            
            
            Automatically generated constructor for t_orb_p
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _SurveySubsF95.f90wrap_t_orb_p_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class T_Orb_P
            
            
            Defined at datadec.fpp lines 26-27
            
            Parameters
            ----------
            this : T_Orb_P
            	Object to be destructed
            
            
            Automatically generated destructor for t_orb_p
            """
            if self._alloc:
                _SurveySubsF95.f90wrap_t_orb_p_finalise(this=self._handle)
        
        @property
        def a(self):
            """
            Element a ftype=real(kind=8) pytype=float
            
            
            Defined at datadec.fpp line 27
            
            """
            return _SurveySubsF95.f90wrap_t_orb_p__get__a(self._handle)
        
        @a.setter
        def a(self, a):
            _SurveySubsF95.f90wrap_t_orb_p__set__a(self._handle, a)
        
        @property
        def e(self):
            """
            Element e ftype=real(kind=8) pytype=float
            
            
            Defined at datadec.fpp line 27
            
            """
            return _SurveySubsF95.f90wrap_t_orb_p__get__e(self._handle)
        
        @e.setter
        def e(self, e):
            _SurveySubsF95.f90wrap_t_orb_p__set__e(self._handle, e)
        
        @property
        def inc(self):
            """
            Element inc ftype=real(kind=8) pytype=float
            
            
            Defined at datadec.fpp line 27
            
            """
            return _SurveySubsF95.f90wrap_t_orb_p__get__inc(self._handle)
        
        @inc.setter
        def inc(self, inc):
            _SurveySubsF95.f90wrap_t_orb_p__set__inc(self._handle, inc)
        
        @property
        def node(self):
            """
            Element node ftype=real(kind=8) pytype=float
            
            
            Defined at datadec.fpp line 27
            
            """
            return _SurveySubsF95.f90wrap_t_orb_p__get__node(self._handle)
        
        @node.setter
        def node(self, node):
            _SurveySubsF95.f90wrap_t_orb_p__set__node(self._handle, node)
        
        @property
        def peri(self):
            """
            Element peri ftype=real(kind=8) pytype=float
            
            
            Defined at datadec.fpp line 27
            
            """
            return _SurveySubsF95.f90wrap_t_orb_p__get__peri(self._handle)
        
        @peri.setter
        def peri(self, peri):
            _SurveySubsF95.f90wrap_t_orb_p__set__peri(self._handle, peri)
        
        @property
        def tperi(self):
            """
            Element tperi ftype=real(kind=8) pytype=float
            
            
            Defined at datadec.fpp line 27
            
            """
            return _SurveySubsF95.f90wrap_t_orb_p__get__tperi(self._handle)
        
        @tperi.setter
        def tperi(self, tperi):
            _SurveySubsF95.f90wrap_t_orb_p__set__tperi(self._handle, tperi)
        
        def __str__(self):
            ret = ['<t_orb_p>{\n']
            ret.append('    a : ')
            ret.append(repr(self.a))
            ret.append(',\n    e : ')
            ret.append(repr(self.e))
            ret.append(',\n    inc : ')
            ret.append(repr(self.inc))
            ret.append(',\n    node : ')
            ret.append(repr(self.node))
            ret.append(',\n    peri : ')
            ret.append(repr(self.peri))
            ret.append(',\n    tperi : ')
            ret.append(repr(self.tperi))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @f90wrap.runtime.register_class("SurveySubsF95.t_v3d")
    class t_v3d(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=t_v3d)
        
        
        Defined at datadec.fpp lines 29-30
        
        """
        def __init__(self, handle=None):
            """
            self = T_V3D()
            
            
            Defined at datadec.fpp lines 29-30
            
            
            Returns
            -------
            this : T_V3D
            	Object to be constructed
            
            
            Automatically generated constructor for t_v3d
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _SurveySubsF95.f90wrap_t_v3d_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class T_V3D
            
            
            Defined at datadec.fpp lines 29-30
            
            Parameters
            ----------
            this : T_V3D
            	Object to be destructed
            
            
            Automatically generated destructor for t_v3d
            """
            if self._alloc:
                _SurveySubsF95.f90wrap_t_v3d_finalise(this=self._handle)
        
        @property
        def x(self):
            """
            Element x ftype=real(kind=8) pytype=float
            
            
            Defined at datadec.fpp line 30
            
            """
            return _SurveySubsF95.f90wrap_t_v3d__get__x(self._handle)
        
        @x.setter
        def x(self, x):
            _SurveySubsF95.f90wrap_t_v3d__set__x(self._handle, x)
        
        @property
        def y(self):
            """
            Element y ftype=real(kind=8) pytype=float
            
            
            Defined at datadec.fpp line 30
            
            """
            return _SurveySubsF95.f90wrap_t_v3d__get__y(self._handle)
        
        @y.setter
        def y(self, y):
            _SurveySubsF95.f90wrap_t_v3d__set__y(self._handle, y)
        
        @property
        def z(self):
            """
            Element z ftype=real(kind=8) pytype=float
            
            
            Defined at datadec.fpp line 30
            
            """
            return _SurveySubsF95.f90wrap_t_v3d__get__z(self._handle)
        
        @z.setter
        def z(self, z):
            _SurveySubsF95.f90wrap_t_v3d__set__z(self._handle, z)
        
        def __str__(self):
            ret = ['<t_v3d>{\n']
            ret.append('    x : ')
            ret.append(repr(self.x))
            ret.append(',\n    y : ')
            ret.append(repr(self.y))
            ret.append(',\n    z : ')
            ret.append(repr(self.z))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @f90wrap.runtime.register_class("SurveySubsF95.t_obspos")
    class t_obspos(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=t_obspos)
        
        
        Defined at datadec.fpp lines 32-34
        
        """
        def __init__(self, handle=None):
            """
            self = T_Obspos()
            
            
            Defined at datadec.fpp lines 32-34
            
            
            Returns
            -------
            this : T_Obspos
            	Object to be constructed
            
            
            Automatically generated constructor for t_obspos
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _SurveySubsF95.f90wrap_t_obspos_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class T_Obspos
            
            
            Defined at datadec.fpp lines 32-34
            
            Parameters
            ----------
            this : T_Obspos
            	Object to be destructed
            
            
            Automatically generated destructor for t_obspos
            """
            if self._alloc:
                _SurveySubsF95.f90wrap_t_obspos_finalise(this=self._handle)
        
        @property
        def jday(self):
            """
            Element jday ftype=real(kind=8) pytype=float
            
            
            Defined at datadec.fpp line 33
            
            """
            return _SurveySubsF95.f90wrap_t_obspos__get__jday(self._handle)
        
        @jday.setter
        def jday(self, jday):
            _SurveySubsF95.f90wrap_t_obspos__set__jday(self._handle, jday)
        
        @property
        def r(self):
            """
            Element r ftype=real(kind=8) pytype=float
            
            
            Defined at datadec.fpp line 33
            
            """
            return _SurveySubsF95.f90wrap_t_obspos__get__r(self._handle)
        
        @r.setter
        def r(self, r):
            _SurveySubsF95.f90wrap_t_obspos__set__r(self._handle, r)
        
        @property
        def pos(self):
            """
            Element pos ftype=type(t_v3d) pytype=T_V3D
            
            
            Defined at datadec.fpp line 34
            
            """
            pos_handle = _SurveySubsF95.f90wrap_t_obspos__get__pos(self._handle)
            if tuple(pos_handle) in self._objs:
                pos = self._objs[tuple(pos_handle)]
            else:
                pos = datadec.t_v3d.from_handle(pos_handle)
                self._objs[tuple(pos_handle)] = pos
            return pos
        
        @pos.setter
        def pos(self, pos):
            pos = pos._handle
            _SurveySubsF95.f90wrap_t_obspos__set__pos(self._handle, pos)
        
        def __str__(self):
            ret = ['<t_obspos>{\n']
            ret.append('    jday : ')
            ret.append(repr(self.jday))
            ret.append(',\n    r : ')
            ret.append(repr(self.r))
            ret.append(',\n    pos : ')
            ret.append(repr(self.pos))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @f90wrap.runtime.register_class("SurveySubsF95.t_eff_r")
    class t_eff_r(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=t_eff_r)
        
        
        Defined at datadec.fpp lines 36-39
        
        """
        def __init__(self, handle=None):
            """
            self = T_Eff_R()
            
            
            Defined at datadec.fpp lines 36-39
            
            
            Returns
            -------
            this : T_Eff_R
            	Object to be constructed
            
            
            Automatically generated constructor for t_eff_r
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _SurveySubsF95.f90wrap_t_eff_r_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class T_Eff_R
            
            
            Defined at datadec.fpp lines 36-39
            
            Parameters
            ----------
            this : T_Eff_R
            	Object to be destructed
            
            
            Automatically generated destructor for t_eff_r
            """
            if self._alloc:
                _SurveySubsF95.f90wrap_t_eff_r_finalise(this=self._handle)
        
        @property
        def min_bn(self):
            """
            Element min_bn ftype=real(kind=8) pytype=float
            
            
            Defined at datadec.fpp line 37
            
            """
            return _SurveySubsF95.f90wrap_t_eff_r__get__min_bn(self._handle)
        
        @min_bn.setter
        def min_bn(self, min_bn):
            _SurveySubsF95.f90wrap_t_eff_r__set__min_bn(self._handle, min_bn)
        
        @property
        def max_bn(self):
            """
            Element max_bn ftype=real(kind=8) pytype=float
            
            
            Defined at datadec.fpp line 37
            
            """
            return _SurveySubsF95.f90wrap_t_eff_r__get__max_bn(self._handle)
        
        @max_bn.setter
        def max_bn(self, max_bn):
            _SurveySubsF95.f90wrap_t_eff_r__set__max_bn(self._handle, max_bn)
        
        @property
        def mag_lim(self):
            """
            Element mag_lim ftype=real(kind=8) pytype=float
            
            
            Defined at datadec.fpp line 37
            
            """
            return _SurveySubsF95.f90wrap_t_eff_r__get__mag_lim(self._handle)
        
        @mag_lim.setter
        def mag_lim(self, mag_lim):
            _SurveySubsF95.f90wrap_t_eff_r__set__mag_lim(self._handle, mag_lim)
        
        @property
        def n(self):
            """
            Element n ftype=integer  pytype=int
            
            
            Defined at datadec.fpp line 38
            
            """
            return _SurveySubsF95.f90wrap_t_eff_r__get__n(self._handle)
        
        @n.setter
        def n(self, n):
            _SurveySubsF95.f90wrap_t_eff_r__set__n(self._handle, n)
        
        @property
        def b(self):
            """
            Element b ftype=real(kind=8) pytype=float
            
            
            Defined at datadec.fpp line 39
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _SurveySubsF95.f90wrap_t_eff_r__array__b(self._handle)
            if array_handle in self._arrays:
                b = self._arrays[array_handle]
            else:
                b = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _SurveySubsF95.f90wrap_t_eff_r__array__b)
                self._arrays[array_handle] = b
            return b
        
        @b.setter
        def b(self, b):
            self.b[...] = b
        
        @property
        def e(self):
            """
            Element e ftype=real(kind=8) pytype=float
            
            
            Defined at datadec.fpp line 39
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _SurveySubsF95.f90wrap_t_eff_r__array__e(self._handle)
            if array_handle in self._arrays:
                e = self._arrays[array_handle]
            else:
                e = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _SurveySubsF95.f90wrap_t_eff_r__array__e)
                self._arrays[array_handle] = e
            return e
        
        @e.setter
        def e(self, e):
            self.e[...] = e
        
        def __str__(self):
            ret = ['<t_eff_r>{\n']
            ret.append('    min_bn : ')
            ret.append(repr(self.min_bn))
            ret.append(',\n    max_bn : ')
            ret.append(repr(self.max_bn))
            ret.append(',\n    mag_lim : ')
            ret.append(repr(self.mag_lim))
            ret.append(',\n    n : ')
            ret.append(repr(self.n))
            ret.append(',\n    b : ')
            ret.append(repr(self.b))
            ret.append(',\n    e : ')
            ret.append(repr(self.e))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @f90wrap.runtime.register_class("SurveySubsF95.t_charact")
    class t_charact(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=t_charact)
        
        
        Defined at datadec.fpp lines 41-45
        
        """
        def __init__(self, handle=None):
            """
            self = T_Charact()
            
            
            Defined at datadec.fpp lines 41-45
            
            
            Returns
            -------
            this : T_Charact
            	Object to be constructed
            
            
            Automatically generated constructor for t_charact
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _SurveySubsF95.f90wrap_t_charact_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class T_Charact
            
            
            Defined at datadec.fpp lines 41-45
            
            Parameters
            ----------
            this : T_Charact
            	Object to be destructed
            
            
            Automatically generated destructor for t_charact
            """
            if self._alloc:
                _SurveySubsF95.f90wrap_t_charact_finalise(this=self._handle)
        
        @property
        def r_cut(self):
            """
            Element r_cut ftype=type(t_ratecut) pytype=T_Ratecut
            
            
            Defined at datadec.fpp line 42
            
            """
            r_cut_handle = _SurveySubsF95.f90wrap_t_charact__get__r_cut(self._handle)
            if tuple(r_cut_handle) in self._objs:
                r_cut = self._objs[tuple(r_cut_handle)]
            else:
                r_cut = datadec.t_ratecut.from_handle(r_cut_handle)
                self._objs[tuple(r_cut_handle)] = r_cut
            return r_cut
        
        @r_cut.setter
        def r_cut(self, r_cut):
            r_cut = r_cut._handle
            _SurveySubsF95.f90wrap_t_charact__set__r_cut(self._handle, r_cut)
        
        @property
        def mag_er(self):
            """
            Element mag_er ftype=real(kind=8) pytype=float
            
            
            Defined at datadec.fpp line 43
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _SurveySubsF95.f90wrap_t_charact__array__mag_er(self._handle)
            if array_handle in self._arrays:
                mag_er = self._arrays[array_handle]
            else:
                mag_er = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _SurveySubsF95.f90wrap_t_charact__array__mag_er)
                self._arrays[array_handle] = mag_er
            return mag_er
        
        @mag_er.setter
        def mag_er(self, mag_er):
            self.mag_er[...] = mag_er
        
        @property
        def photf(self):
            """
            Element photf ftype=real(kind=8) pytype=float
            
            
            Defined at datadec.fpp line 43
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _SurveySubsF95.f90wrap_t_charact__array__photf(self._handle)
            if array_handle in self._arrays:
                photf = self._arrays[array_handle]
            else:
                photf = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _SurveySubsF95.f90wrap_t_charact__array__photf)
                self._arrays[array_handle] = photf
            return photf
        
        @photf.setter
        def photf(self, photf):
            self.photf[...] = photf
        
        @property
        def track(self):
            """
            Element track ftype=real(kind=8) pytype=float
            
            
            Defined at datadec.fpp line 43
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _SurveySubsF95.f90wrap_t_charact__array__track(self._handle)
            if array_handle in self._arrays:
                track = self._arrays[array_handle]
            else:
                track = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _SurveySubsF95.f90wrap_t_charact__array__track)
                self._arrays[array_handle] = track
            return track
        
        @track.setter
        def track(self, track):
            self.track[...] = track
        
        @property
        def f(self):
            """
            Element f ftype=integer  pytype=int
            
            
            Defined at datadec.fpp line 44
            
            """
            return _SurveySubsF95.f90wrap_t_charact__get__f(self._handle)
        
        @f.setter
        def f(self, f):
            _SurveySubsF95.f90wrap_t_charact__set__f(self._handle, f)
        
        @property
        def nr(self):
            """
            Element nr ftype=integer  pytype=int
            
            
            Defined at datadec.fpp line 44
            
            """
            return _SurveySubsF95.f90wrap_t_charact__get__nr(self._handle)
        
        @nr.setter
        def nr(self, nr):
            _SurveySubsF95.f90wrap_t_charact__set__nr(self._handle, nr)
        
        def init_array_eff_p(self):
            self.eff_p = f90wrap.runtime.FortranDerivedTypeArray(self,
                                            _SurveySubsF95.f90wrap_t_charact__array_getitem__eff_p,
                                            _SurveySubsF95.f90wrap_t_charact__array_setitem__eff_p,
                                            _SurveySubsF95.f90wrap_t_charact__array_len__eff_p,
                                            """
            Element eff_p ftype=type(t_eff_r) pytype=T_Eff_R
            
            
            Defined at datadec.fpp line 45
            
            """, Datadec.t_eff_r)
            return self.eff_p
        
        def __str__(self):
            ret = ['<t_charact>{\n']
            ret.append('    r_cut : ')
            ret.append(repr(self.r_cut))
            ret.append(',\n    mag_er : ')
            ret.append(repr(self.mag_er))
            ret.append(',\n    photf : ')
            ret.append(repr(self.photf))
            ret.append(',\n    track : ')
            ret.append(repr(self.track))
            ret.append(',\n    f : ')
            ret.append(repr(self.f))
            ret.append(',\n    nr : ')
            ret.append(repr(self.nr))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = [init_array_eff_p]
        
    
    @f90wrap.runtime.register_class("SurveySubsF95.t_pointing")
    class t_pointing(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=t_pointing)
        
        
        Defined at datadec.fpp lines 47-53
        
        """
        def __init__(self, handle=None):
            """
            self = T_Pointing()
            
            
            Defined at datadec.fpp lines 47-53
            
            
            Returns
            -------
            this : T_Pointing
            	Object to be constructed
            
            
            Automatically generated constructor for t_pointing
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _SurveySubsF95.f90wrap_t_pointing_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class T_Pointing
            
            
            Defined at datadec.fpp lines 47-53
            
            Parameters
            ----------
            this : T_Pointing
            	Object to be destructed
            
            
            Automatically generated destructor for t_pointing
            """
            if self._alloc:
                _SurveySubsF95.f90wrap_t_pointing_finalise(this=self._handle)
        
        @property
        def ff(self):
            """
            Element ff ftype=real(kind=8) pytype=float
            
            
            Defined at datadec.fpp line 48
            
            """
            return _SurveySubsF95.f90wrap_t_pointing__get__ff(self._handle)
        
        @ff.setter
        def ff(self, ff):
            _SurveySubsF95.f90wrap_t_pointing__set__ff(self._handle, ff)
        
        @property
        def code(self):
            """
            Element code ftype=integer  pytype=int
            
            
            Defined at datadec.fpp line 49
            
            """
            return _SurveySubsF95.f90wrap_t_pointing__get__code(self._handle)
        
        @code.setter
        def code(self, code):
            _SurveySubsF95.f90wrap_t_pointing__set__code(self._handle, code)
        
        @property
        def efnam(self):
            """
            Element efnam ftype=character(80) pytype=str
            
            
            Defined at datadec.fpp line 50
            
            """
            return _SurveySubsF95.f90wrap_t_pointing__get__efnam(self._handle)
        
        @efnam.setter
        def efnam(self, efnam):
            _SurveySubsF95.f90wrap_t_pointing__set__efnam(self._handle, efnam)
        
        def init_array_o_pos(self):
            self.o_pos = f90wrap.runtime.FortranDerivedTypeArray(self,
                                            _SurveySubsF95.f90wrap_t_pointing__array_getitem__o_pos,
                                            _SurveySubsF95.f90wrap_t_pointing__array_setitem__o_pos,
                                            _SurveySubsF95.f90wrap_t_pointing__array_len__o_pos,
                                            """
            Element o_pos ftype=type(t_obspos) pytype=T_Obspos
            
            
            Defined at datadec.fpp line 51
            
            """, Datadec.t_obspos)
            return self.o_pos
        
        @property
        def poly(self):
            """
            Element poly ftype=type(t_polygon) pytype=T_Polygon
            
            
            Defined at datadec.fpp line 52
            
            """
            poly_handle = _SurveySubsF95.f90wrap_t_pointing__get__poly(self._handle)
            if tuple(poly_handle) in self._objs:
                poly = self._objs[tuple(poly_handle)]
            else:
                poly = poly_dec.t_polygon.from_handle(poly_handle)
                self._objs[tuple(poly_handle)] = poly
            return poly
        
        @poly.setter
        def poly(self, poly):
            poly = poly._handle
            _SurveySubsF95.f90wrap_t_pointing__set__poly(self._handle, poly)
        
        @property
        def c(self):
            """
            Element c ftype=type(t_charact) pytype=T_Charact
            
            
            Defined at datadec.fpp line 53
            
            """
            c_handle = _SurveySubsF95.f90wrap_t_pointing__get__c(self._handle)
            if tuple(c_handle) in self._objs:
                c = self._objs[tuple(c_handle)]
            else:
                c = datadec.t_charact.from_handle(c_handle)
                self._objs[tuple(c_handle)] = c
            return c
        
        @c.setter
        def c(self, c):
            c = c._handle
            _SurveySubsF95.f90wrap_t_pointing__set__c(self._handle, c)
        
        def __str__(self):
            ret = ['<t_pointing>{\n']
            ret.append('    ff : ')
            ret.append(repr(self.ff))
            ret.append(',\n    code : ')
            ret.append(repr(self.code))
            ret.append(',\n    efnam : ')
            ret.append(repr(self.efnam))
            ret.append(',\n    poly : ')
            ret.append(repr(self.poly))
            ret.append(',\n    c : ')
            ret.append(repr(self.c))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = [init_array_o_pos]
        
    
    @f90wrap.runtime.register_class("SurveySubsF95.T_Pointing_Xn_Sur_Max_Array")
    class T_Pointing_Xn_Sur_Max_Array(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=t_pointing_xn_sur_max_array)
        
        
        Defined at datadec.fpp lines 47-53
        
        super-type
        Automatically generated to handle derived type arrays as a new derived type
        """
        def __init__(self, handle=None):
            """
            self = T_Pointing_Xn_Sur_Max_Array()
            
            
            Defined at datadec.fpp lines 47-53
            
            
            Returns
            -------
            this : T_Pointing_Xn_Sur_Max_Array
            	Object to be constructed
            
            
            Automatically generated constructor for t_pointing_xn_sur_max_array
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _SurveySubsF95.f90wrap_t_pointing_xn_sur_max_array_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class T_Pointing_Xn_Sur_Max_Array
            
            
            Defined at datadec.fpp lines 47-53
            
            Parameters
            ----------
            this : T_Pointing_Xn_Sur_Max_Array
            	Object to be destructed
            
            
            Automatically generated destructor for t_pointing_xn_sur_max_array
            """
            if self._alloc:
                _SurveySubsF95.f90wrap_t_pointing_xn_sur_max_array_finalise(this=self._handle)
        
        def init_array_items(self):
            self.items = f90wrap.runtime.FortranDerivedTypeArray(self,
                                            _SurveySubsF95.f90wrap_t_pointing_xn_sur_max_array__array_getitem__items,
                                            _SurveySubsF95.f90wrap_t_pointing_xn_sur_max_array__array_setitem__items,
                                            _SurveySubsF95.f90wrap_t_pointing_xn_sur_max_array__array_len__items,
                                            """
            Element items ftype=type(t_pointing) pytype=T_Pointing
            
            
            Defined at  line 0
            
            """, Datadec.t_pointing)
            return self.items
        
        _dt_array_initialisers = [init_array_items]
        
    
    @property
    def n_sur_max(self):
        """
        Element n_sur_max ftype=integer pytype=int
        
        
        Defined at datadec.fpp line 9
        
        """
        return _SurveySubsF95.f90wrap_datadec__get__n_sur_max()
    
    @property
    def n_bin_max(self):
        """
        Element n_bin_max ftype=integer pytype=int
        
        
        Defined at datadec.fpp line 9
        
        """
        return _SurveySubsF95.f90wrap_datadec__get__n_bin_max()
    
    @property
    def n_r_max(self):
        """
        Element n_r_max ftype=integer pytype=int
        
        
        Defined at datadec.fpp line 9
        
        """
        return _SurveySubsF95.f90wrap_datadec__get__n_r_max()
    
    @property
    def nw_max(self):
        """
        Element nw_max ftype=integer pytype=int
        
        
        Defined at datadec.fpp line 9
        
        """
        return _SurveySubsF95.f90wrap_datadec__get__nw_max()
    
    @property
    def pi(self):
        """
        Element pi ftype=real(kind=8) pytype=float
        
        
        Defined at datadec.fpp line 12
        
        """
        return _SurveySubsF95.f90wrap_datadec__get__pi()
    
    @property
    def drad(self):
        """
        Element drad ftype=real(kind=8) pytype=float
        
        
        Defined at datadec.fpp line 12
        
        """
        return _SurveySubsF95.f90wrap_datadec__get__drad()
    
    @property
    def twohours(self):
        """
        Element twohours ftype=real(kind=8) pytype=float
        
        
        Defined at datadec.fpp line 12
        
        """
        return _SurveySubsF95.f90wrap_datadec__get__twohours()
    
    @property
    def twopi(self):
        """
        Element twopi ftype=real(kind=8) pytype=float
        
        
        Defined at datadec.fpp line 12
        
        """
        return _SurveySubsF95.f90wrap_datadec__get__twopi()
    
    @property
    def eps(self):
        """
        Element eps ftype=real(kind=8) pytype=float
        
        
        Defined at datadec.fpp line 12
        
        """
        return _SurveySubsF95.f90wrap_datadec__get__eps()
    
    @property
    def gmb(self):
        """
        Element gmb ftype=real(kind=8) pytype=float
        
        
        Defined at datadec.fpp line 15
        
        """
        return _SurveySubsF95.f90wrap_datadec__get__gmb()
    
    @property
    def om_lim_low(self):
        """
        Element om_lim_low ftype=real(kind=8) pytype=float
        
        
        Defined at datadec.fpp line 17
        
        """
        return _SurveySubsF95.f90wrap_datadec__get__om_lim_low()
    
    @om_lim_low.setter
    def om_lim_low(self, om_lim_low):
        _SurveySubsF95.f90wrap_datadec__set__om_lim_low(om_lim_low)
    
    @property
    def om_lim_high(self):
        """
        Element om_lim_high ftype=real(kind=8) pytype=float
        
        
        Defined at datadec.fpp line 17
        
        """
        return _SurveySubsF95.f90wrap_datadec__get__om_lim_high()
    
    @om_lim_high.setter
    def om_lim_high(self, om_lim_high):
        _SurveySubsF95.f90wrap_datadec__set__om_lim_high(om_lim_high)
    
    def __str__(self):
        ret = ['<datadec>{\n']
        ret.append('    n_sur_max : ')
        ret.append(repr(self.n_sur_max))
        ret.append(',\n    n_bin_max : ')
        ret.append(repr(self.n_bin_max))
        ret.append(',\n    n_r_max : ')
        ret.append(repr(self.n_r_max))
        ret.append(',\n    nw_max : ')
        ret.append(repr(self.nw_max))
        ret.append(',\n    pi : ')
        ret.append(repr(self.pi))
        ret.append(',\n    drad : ')
        ret.append(repr(self.drad))
        ret.append(',\n    twohours : ')
        ret.append(repr(self.twohours))
        ret.append(',\n    twopi : ')
        ret.append(repr(self.twopi))
        ret.append(',\n    eps : ')
        ret.append(repr(self.eps))
        ret.append(',\n    gmb : ')
        ret.append(repr(self.gmb))
        ret.append(',\n    om_lim_low : ')
        ret.append(repr(self.om_lim_low))
        ret.append(',\n    om_lim_high : ')
        ret.append(repr(self.om_lim_high))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

datadec = Datadec()

class Ioutils(f90wrap.runtime.FortranModule):
    """
    Module ioutils
    
    
    Defined at ioutils.fpp lines 5-285
    
    """
    @staticmethod
    def trim(base_name, len_bn):
        """
        i1, i2, finished = trim(base_name, len_bn)
        
        
        Defined at ioutils.fpp lines 8-60
        
        Parameters
        ----------
        base_name : str
        len_bn : int
        
        Returns
        -------
        i1 : int
        i2 : int
        finished : bool
        
        """
        i1, i2, finished = _SurveySubsF95.f90wrap_trim(base_name=base_name, \
            len_bn=len_bn)
        return i1, i2, finished
    
    @staticmethod
    def read_file_name(base_name, len_bn):
        """
        i1, i2, finished = read_file_name(base_name, len_bn)
        
        
        Defined at ioutils.fpp lines 62-115
        
        Parameters
        ----------
        base_name : str
        len_bn : int
        
        Returns
        -------
        i1 : int
        i2 : int
        finished : bool
        
        """
        i1, i2, finished = _SurveySubsF95.f90wrap_read_file_name(base_name=base_name, \
            len_bn=len_bn)
        return i1, i2, finished
    
    @staticmethod
    def parse(comd, nwmax, word, lw):
        """
        nw = parse(comd, nwmax, word, lw)
        
        
        Defined at ioutils.fpp lines 119-196
        
        Parameters
        ----------
        comd : str
        nwmax : int
        word : str array
        lw : int array
        
        Returns
        -------
        nw : int
        
        """
        nw = _SurveySubsF95.f90wrap_parse(comd=comd, nwmax=nwmax, word=word, lw=lw)
        return nw
    
    @staticmethod
    def format(angle, incode, outcod):
        """
        string_bn, ierr = format(angle, incode, outcod)
        
        
        Defined at ioutils.fpp lines 198-273
        
        Parameters
        ----------
        angle : float
        incode : int
        outcod : int
        
        Returns
        -------
        string_bn : str
        ierr : int
        
        """
        string_bn, ierr = _SurveySubsF95.f90wrap_format(angle=angle, incode=incode, \
            outcod=outcod)
        return string_bn, ierr
    
    @staticmethod
    def strip_comment(str):
        """
        strip_comment = strip_comment(str)
        
        
        Defined at ioutils.fpp lines 275-285
        
        Parameters
        ----------
        str : str
        
        Returns
        -------
        strip_comment : str
        
        """
        strip_comment = _SurveySubsF95.f90wrap_strip_comment(str=str)
        return strip_comment
    
    _dt_array_initialisers = []
    

ioutils = Ioutils()

class Effut(f90wrap.runtime.FortranModule):
    """
    Module effut
    
    
    Defined at effut.fpp lines 5-192
    
    """
    pass
    _dt_array_initialisers = []
    

effut = Effut()

class Elemutils(f90wrap.runtime.FortranModule):
    """
    Module elemutils
    
    
    Defined at elemutils.fpp lines 5-425
    
    """
    @staticmethod
    def coord_cart(mu, o_m):
        """
        p, v = coord_cart(mu, o_m)
        
        
        Defined at elemutils.fpp lines 9-129
        
        Parameters
        ----------
        mu : float
        o_m : T_Orb_M
        
        Returns
        -------
        p : T_V3D
        v : T_V3D
        
        """
        p, v = _SurveySubsF95.f90wrap_coord_cart(mu=mu, o_m=o_m._handle)
        p = f90wrap.runtime.lookup_class("SurveySubsF95.t_v3d").from_handle(p, \
            alloc=True)
        v = f90wrap.runtime.lookup_class("SurveySubsF95.t_v3d").from_handle(v, \
            alloc=True)
        return p, v
    
    @staticmethod
    def pos_cart(self):
        """
        p = pos_cart(self)
        
        
        Defined at elemutils.fpp lines 132-240
        
        Parameters
        ----------
        o_m : T_Orb_M
        
        Returns
        -------
        p : T_V3D
        
        """
        p = _SurveySubsF95.f90wrap_pos_cart(o_m=self._handle)
        p = f90wrap.runtime.lookup_class("SurveySubsF95.t_v3d").from_handle(p, \
            alloc=True)
        return p
    
    @staticmethod
    def pq_cart(inc, node, peri, p, q, r):
        """
        pq_cart(inc, node, peri, p, q, r)
        
        
        Defined at elemutils.fpp lines 243-307
        
        Parameters
        ----------
        inc : float
        node : float
        peri : float
        p : float array
        q : float array
        r : float array
        
        """
        _SurveySubsF95.f90wrap_pq_cart(inc=inc, node=node, peri=peri, p=p, q=q, r=r)
    
    @staticmethod
    def osc_el(mu, p, v):
        """
        o_m = osc_el(mu, p, v)
        
        
        Defined at elemutils.fpp lines 310-425
        
        Parameters
        ----------
        mu : float
        p : T_V3D
        v : T_V3D
        
        Returns
        -------
        o_m : T_Orb_M
        
        """
        o_m = _SurveySubsF95.f90wrap_osc_el(mu=mu, p=p._handle, v=v._handle)
        o_m = f90wrap.runtime.lookup_class("SurveySubsF95.t_orb_m").from_handle(o_m, \
            alloc=True)
        return o_m
    
    _dt_array_initialisers = []
    

elemutils = Elemutils()

class Rot(f90wrap.runtime.FortranModule):
    """
    Module rot
    
    
    Defined at rot.fpp lines 5-763
    
    """
    @staticmethod
    def equ_ecl(epsilon, poseq, posecl):
        """
        equ_ecl(epsilon, poseq, posecl)
        
        
        Defined at rot.fpp lines 9-40
        
        Parameters
        ----------
        epsilon : float
        poseq : float array
        posecl : float array
        
        """
        _SurveySubsF95.f90wrap_equ_ecl(epsilon=epsilon, poseq=poseq, posecl=posecl)
    
    @staticmethod
    def ecl_equ(epsilon, posecl, poseq):
        """
        ecl_equ(epsilon, posecl, poseq)
        
        
        Defined at rot.fpp lines 42-73
        
        Parameters
        ----------
        epsilon : float
        posecl : float array
        poseq : float array
        
        """
        _SurveySubsF95.f90wrap_ecl_equ(epsilon=epsilon, posecl=posecl, poseq=poseq)
    
    @staticmethod
    def rotx(alpha, posin, posout):
        """
        rotx(alpha, posin, posout)
        
        
        Defined at rot.fpp lines 75-107
        
        Parameters
        ----------
        alpha : float
        posin : float array
        posout : float array
        
        """
        _SurveySubsF95.f90wrap_rotx(alpha=alpha, posin=posin, posout=posout)
    
    @staticmethod
    def roty(alpha, posin, posout):
        """
        roty(alpha, posin, posout)
        
        
        Defined at rot.fpp lines 109-141
        
        Parameters
        ----------
        alpha : float
        posin : float array
        posout : float array
        
        """
        _SurveySubsF95.f90wrap_roty(alpha=alpha, posin=posin, posout=posout)
    
    @staticmethod
    def rotz(alpha, posin, posout):
        """
        rotz(alpha, posin, posout)
        
        
        Defined at rot.fpp lines 143-175
        
        Parameters
        ----------
        alpha : float
        posin : float array
        posout : float array
        
        """
        _SurveySubsF95.f90wrap_rotz(alpha=alpha, posin=posin, posout=posout)
    
    @staticmethod
    def equat_ecl(ieqec, v_in):
        """
        v_out, ierr = equat_ecl(ieqec, v_in)
        
        
        Defined at rot.fpp lines 177-238
        
        Parameters
        ----------
        ieqec : int
        v_in : T_V3D
        
        Returns
        -------
        v_out : T_V3D
        ierr : int
        
        """
        v_out, ierr = _SurveySubsF95.f90wrap_equat_ecl(ieqec=ieqec, v_in=v_in._handle)
        v_out = f90wrap.runtime.lookup_class("SurveySubsF95.t_v3d").from_handle(v_out, \
            alloc=True)
        return v_out, ierr
    
    @staticmethod
    def invar_ecl(ieqec, v_in):
        """
        v_out, ierr = invar_ecl(ieqec, v_in)
        
        
        Defined at rot.fpp lines 240-311
        
        Parameters
        ----------
        ieqec : int
        v_in : T_V3D
        
        Returns
        -------
        v_out : T_V3D
        ierr : int
        
        """
        v_out, ierr = _SurveySubsF95.f90wrap_invar_ecl(ieqec=ieqec, v_in=v_in._handle)
        v_out = f90wrap.runtime.lookup_class("SurveySubsF95.t_v3d").from_handle(v_out, \
            alloc=True)
        return v_out, ierr
    
    @staticmethod
    def invar_ecl_osc(ieqec, o_mi):
        """
        o_mo, ierr = invar_ecl_osc(ieqec, o_mi)
        
        
        Defined at rot.fpp lines 313-351
        
        Parameters
        ----------
        ieqec : int
        o_mi : T_Orb_M
        
        Returns
        -------
        o_mo : T_Orb_M
        ierr : int
        
        """
        o_mo, ierr = _SurveySubsF95.f90wrap_invar_ecl_osc(ieqec=ieqec, \
            o_mi=o_mi._handle)
        o_mo = f90wrap.runtime.lookup_class("SurveySubsF95.t_orb_m").from_handle(o_mo, \
            alloc=True)
        return o_mo, ierr
    
    @staticmethod
    def invar_ecl_inc_node(ieqec, ii, noi):
        """
        io, noo, ierr = invar_ecl_inc_node(ieqec, ii, noi)
        
        
        Defined at rot.fpp lines 353-383
        
        Parameters
        ----------
        ieqec : int
        ii : float
        noi : float
        
        Returns
        -------
        io : float
        noo : float
        ierr : int
        
        """
        io, noo, ierr = _SurveySubsF95.f90wrap_invar_ecl_inc_node(ieqec=ieqec, ii=ii, \
            noi=noi)
        return io, noo, ierr
    
    @staticmethod
    def ref_ecl(ieqec, v_in, eps, om):
        """
        v_out, ierr = ref_ecl(ieqec, v_in, eps, om)
        
        
        Defined at rot.fpp lines 385-468
        
        Parameters
        ----------
        ieqec : int
        v_in : T_V3D
        eps : float
        om : float
        
        Returns
        -------
        v_out : T_V3D
        ierr : int
        
        """
        v_out, ierr = _SurveySubsF95.f90wrap_ref_ecl(ieqec=ieqec, v_in=v_in._handle, \
            eps=eps, om=om)
        v_out = f90wrap.runtime.lookup_class("SurveySubsF95.t_v3d").from_handle(v_out, \
            alloc=True)
        return v_out, ierr
    
    @staticmethod
    def ref_ecl_osc(ieqec, o_mi, eps, om):
        """
        o_mo, ierr = ref_ecl_osc(ieqec, o_mi, eps, om)
        
        
        Defined at rot.fpp lines 470-513
        
        Parameters
        ----------
        ieqec : int
        o_mi : T_Orb_M
        eps : float
        om : float
        
        Returns
        -------
        o_mo : T_Orb_M
        ierr : int
        
        """
        o_mo, ierr = _SurveySubsF95.f90wrap_ref_ecl_osc(ieqec=ieqec, o_mi=o_mi._handle, \
            eps=eps, om=om)
        o_mo = f90wrap.runtime.lookup_class("SurveySubsF95.t_orb_m").from_handle(o_mo, \
            alloc=True)
        return o_mo, ierr
    
    @staticmethod
    def ref_ecl_inc_node(ieqec, ii, noi, eps, om):
        """
        io, noo, ierr = ref_ecl_inc_node(ieqec, ii, noi, eps, om)
        
        
        Defined at rot.fpp lines 515-550
        
        Parameters
        ----------
        ieqec : int
        ii : float
        noi : float
        eps : float
        om : float
        
        Returns
        -------
        io : float
        noo : float
        ierr : int
        
        """
        io, noo, ierr = _SurveySubsF95.f90wrap_ref_ecl_inc_node(ieqec=ieqec, ii=ii, \
            noi=noi, eps=eps, om=om)
        return io, noo, ierr
    
    @staticmethod
    def forced_plane(a):
        """
        ifd, omfd = forced_plane(a)
        
        
        Defined at rot.fpp lines 552-647
        
        Parameters
        ----------
        a : float
        
        Returns
        -------
        ifd : float
        omfd : float
        
        """
        ifd, omfd = _SurveySubsF95.f90wrap_forced_plane(a=a)
        return ifd, omfd
    
    @staticmethod
    def forced_plane_damp(a, inc):
        """
        ifd, omfd = forced_plane_damp(a, inc)
        
        
        Defined at rot.fpp lines 649-743
        
        Parameters
        ----------
        a : float
        inc : float
        
        Returns
        -------
        ifd : float
        omfd : float
        
        """
        ifd, omfd = _SurveySubsF95.f90wrap_forced_plane_damp(a=a, inc=inc)
        return ifd, omfd
    
    @staticmethod
    def ztopi(var):
        """
        ztopi(var)
        
        
        Defined at rot.fpp lines 745-763
        
        Parameters
        ----------
        var : float
        
        """
        _SurveySubsF95.f90wrap_ztopi(var=var)
    
    _dt_array_initialisers = []
    

rot = Rot()

class Xvutils(f90wrap.runtime.FortranModule):
    """
    Module xvutils
    
    
    Defined at xvutils.fpp lines 5-813
    
    """
    @staticmethod
    def baryxv(jday):
        """
        pos, vel, ierr = baryxv(jday)
        
        
        Defined at xvutils.fpp lines 9-88
        
        Parameters
        ----------
        jday : float
        
        Returns
        -------
        pos : T_V3D
        vel : T_V3D
        ierr : int
        
        """
        pos, vel, ierr = _SurveySubsF95.f90wrap_baryxv(jday=jday)
        pos = f90wrap.runtime.lookup_class("SurveySubsF95.t_v3d").from_handle(pos, \
            alloc=True)
        vel = f90wrap.runtime.lookup_class("SurveySubsF95.t_v3d").from_handle(vel, \
            alloc=True)
        return pos, vel, ierr
    
    @staticmethod
    def distsunecl(jday, pos):
        """
        r = distsunecl(jday, pos)
        
        
        Defined at xvutils.fpp lines 90-129
        
        Parameters
        ----------
        jday : float
        pos : T_V3D
        
        Returns
        -------
        r : float
        
        """
        r = _SurveySubsF95.f90wrap_distsunecl(jday=jday, pos=pos._handle)
        return r
    
    @staticmethod
    def obspos(code, t):
        """
        pos, vel, r, ierr = obspos(code, t)
        
        
        Defined at xvutils.fpp lines 131-217
        
        Parameters
        ----------
        code : int
        t : float
        
        Returns
        -------
        pos : T_V3D
        vel : T_V3D
        r : float
        ierr : int
        
        """
        pos, vel, r, ierr = _SurveySubsF95.f90wrap_obspos(code=code, t=t)
        pos = f90wrap.runtime.lookup_class("SurveySubsF95.t_v3d").from_handle(pos, \
            alloc=True)
        vel = f90wrap.runtime.lookup_class("SurveySubsF95.t_v3d").from_handle(vel, \
            alloc=True)
        return pos, vel, r, ierr
    
    @staticmethod
    def planetelem(ind, jday):
        """
        o_m, ierr = planetelem(ind, jday)
        
        
        Defined at xvutils.fpp lines 219-344
        
        Parameters
        ----------
        ind : int
        jday : float
        
        Returns
        -------
        o_m : T_Orb_M
        ierr : int
        
        """
        o_m, ierr = _SurveySubsF95.f90wrap_planetelem(ind=ind, jday=jday)
        o_m = f90wrap.runtime.lookup_class("SurveySubsF95.t_orb_m").from_handle(o_m, \
            alloc=True)
        return o_m, ierr
    
    @staticmethod
    def planetxv(ind, jday):
        """
        pos, vel, ierr = planetxv(ind, jday)
        
        
        Defined at xvutils.fpp lines 346-422
        
        Parameters
        ----------
        ind : int
        jday : float
        
        Returns
        -------
        pos : T_V3D
        vel : T_V3D
        ierr : int
        
        """
        pos, vel, ierr = _SurveySubsF95.f90wrap_planetxv(ind=ind, jday=jday)
        pos = f90wrap.runtime.lookup_class("SurveySubsF95.t_v3d").from_handle(pos, \
            alloc=True)
        vel = f90wrap.runtime.lookup_class("SurveySubsF95.t_v3d").from_handle(vel, \
            alloc=True)
        return pos, vel, ierr
    
    @staticmethod
    def newcomb(dj, rout):
        """
        newcomb(dj, rout)
        
        
        Defined at xvutils.fpp lines 432-813
        
        Parameters
        ----------
        dj : float
        rout : float array
        
        """
        _SurveySubsF95.f90wrap_newcomb(dj=dj, rout=rout)
    
    _dt_array_initialisers = []
    

xvutils = Xvutils()

class Getsur(f90wrap.runtime.FortranModule):
    """
    Module getsur
    
    
    Defined at getsur.fpp lines 5-875
    
    """
    @staticmethod
    def create_ears(ra, dec):
        """
        poly = create_ears(ra, dec)
        
        
        Defined at getsur.fpp lines 12-139
        
        Parameters
        ----------
        ra : float
        dec : float
        
        Returns
        -------
        poly : T_Polygon
        
        """
        poly = _SurveySubsF95.f90wrap_create_ears(ra=ra, dec=dec)
        poly = f90wrap.runtime.lookup_class("SurveySubsF95.t_polygon").from_handle(poly, \
            alloc=True)
        return poly
    
    @staticmethod
    def create_rectangle(w, h, ra, dec):
        """
        poly = create_rectangle(w, h, ra, dec)
        
        
        Defined at getsur.fpp lines 141-180
        
        Parameters
        ----------
        w : float
        h : float
        ra : float
        dec : float
        
        Returns
        -------
        poly : T_Polygon
        
        """
        poly = _SurveySubsF95.f90wrap_create_rectangle(w=w, h=h, ra=ra, dec=dec)
        poly = f90wrap.runtime.lookup_class("SurveySubsF95.t_polygon").from_handle(poly, \
            alloc=True)
        return poly
    
    @staticmethod
    def create_poly(ra, dec):
        """
        poly = create_poly(ra, dec)
        
        
        Defined at getsur.fpp lines 182-211
        
        Parameters
        ----------
        ra : float
        dec : float
        
        Returns
        -------
        poly : T_Polygon
        
        """
        poly = _SurveySubsF95.f90wrap_create_poly(ra=ra, dec=dec)
        poly = f90wrap.runtime.lookup_class("SurveySubsF95.t_polygon").from_handle(poly, \
            alloc=True)
        return poly
    
    @staticmethod
    def hms(str):
        """
        val = hms(str)
        
        
        Defined at getsur.fpp lines 213-271
        
        Parameters
        ----------
        str : str
        
        Returns
        -------
        val : float
        
        """
        val = _SurveySubsF95.f90wrap_hms(str=str)
        return val
    
    @staticmethod
    def read_eff(filen, lun_in):
        """
        c, ierr = read_eff(filen, lun_in)
        
        
        Defined at getsur.fpp lines 273-548
        
        Parameters
        ----------
        filen : str
        lun_in : int
        
        Returns
        -------
        c : T_Charact
        ierr : int
        
        """
        c, ierr = _SurveySubsF95.f90wrap_read_eff(filen=filen, lun_in=lun_in)
        c = f90wrap.runtime.lookup_class("SurveySubsF95.t_charact").from_handle(c, \
            alloc=True)
        return c, ierr
    
    @staticmethod
    def read_sur(dirn, lun_in):
        """
        point, ierr = read_sur(dirn, lun_in)
        
        
        Defined at getsur.fpp lines 550-737
        
        Parameters
        ----------
        dirn : str
        lun_in : int
        
        Returns
        -------
        point : T_Pointing
        ierr : int
        
        """
        point, ierr = _SurveySubsF95.f90wrap_read_sur(dirn=dirn, lun_in=lun_in)
        point = \
            f90wrap.runtime.lookup_class("SurveySubsF95.t_pointing").from_handle(point, \
            alloc=True)
        return point, ierr
    
    @staticmethod
    def getsurvey(survey, lun_s, sur_mm):
        """
        n_sur, points, ierr = getsurvey(survey, lun_s, sur_mm)
        
        
        Defined at getsur.fpp lines 739-875
        
        Parameters
        ----------
        survey : str
        lun_s : int
        sur_mm : float array
        
        Returns
        -------
        n_sur : int
        points : T_Pointing_Xn_Sur_Max_Array
        	super-type
        
        ierr : int
        
        """
        n_sur, points, ierr = _SurveySubsF95.f90wrap_getsurvey(survey=survey, \
            lun_s=lun_s, sur_mm=sur_mm)
        points = \
            f90wrap.runtime.lookup_class("SurveySubsF95.T_Pointing_Xn_Sur_Max_Array").from_handle(points, \
            alloc=True)
        return n_sur, points, ierr
    
    _dt_array_initialisers = []
    

getsur = Getsur()

class Numutils(f90wrap.runtime.FortranModule):
    """
    Module numutils
    
    
    Defined at numutils.fpp lines 5-507
    
    """
    @staticmethod
    def appmag(r, delta, robs, h, g):
        """
        alpha, mag, ierr = appmag(r, delta, robs, h, g)
        
        
        Defined at numutils.fpp lines 10-59
        
        Parameters
        ----------
        r : float
        delta : float
        robs : float
        h : float
        g : float
        
        Returns
        -------
        alpha : float
        mag : float
        ierr : int
        
        """
        alpha, mag, ierr = _SurveySubsF95.f90wrap_appmag(r=r, delta=delta, robs=robs, \
            h=h, g=g)
        return alpha, mag, ierr
    
    @staticmethod
    def absmag(r, delta, robs, mag, g):
        """
        alpha, h, ierr = absmag(r, delta, robs, mag, g)
        
        
        Defined at numutils.fpp lines 61-105
        
        Parameters
        ----------
        r : float
        delta : float
        robs : float
        mag : float
        g : float
        
        Returns
        -------
        alpha : float
        h : float
        ierr : int
        
        """
        alpha, h, ierr = _SurveySubsF95.f90wrap_absmag(r=r, delta=delta, robs=robs, \
            mag=mag, g=g)
        return alpha, h, ierr
    
    @staticmethod
    def dgauss(i):
        """
        y = dgauss(i)
        
        
        Defined at numutils.fpp lines 107-150
        
        Parameters
        ----------
        i : int
        
        Returns
        -------
        y : float
        
        """
        y = _SurveySubsF95.f90wrap_dgauss(i=i)
        return y
    
    @staticmethod
    def magran(mag_t, mag_er):
        """
        seed, mag, magerr = magran(mag_t, mag_er)
        
        
        Defined at numutils.fpp lines 152-215
        
        Parameters
        ----------
        mag_t : float
        mag_er : float array
        
        Returns
        -------
        seed : int
        mag : float
        magerr : float
        
        """
        seed, mag, magerr = _SurveySubsF95.f90wrap_magran(mag_t=mag_t, mag_er=mag_er)
        return seed, mag, magerr
    
    @staticmethod
    def latlong(self):
        """
        long_bn, lat, r = latlong(self)
        
        
        Defined at numutils.fpp lines 217-247
        
        Parameters
        ----------
        pos : T_V3D
        
        Returns
        -------
        long_bn : float
        lat : float
        r : float
        
        """
        long_bn, lat, r = _SurveySubsF95.f90wrap_latlong(pos=self._handle)
        return long_bn, lat, r
    
    @staticmethod
    def radececlxv(self, obspos):
        """
        delta, ra, dec = radececlxv(self, obspos)
        
        
        Defined at numutils.fpp lines 249-289
        
        Parameters
        ----------
        pos : T_V3D
        obspos : T_V3D
        
        Returns
        -------
        delta : float
        ra : float
        dec : float
        
        """
        delta, ra, dec = _SurveySubsF95.f90wrap_radececlxv(pos=self._handle, \
            obspos=obspos._handle)
        return delta, ra, dec
    
    @staticmethod
    def ran3(idum):
        """
        ran3 = ran3(idum)
        
        
        Defined at numutils.fpp lines 291-331
        
        Parameters
        ----------
        idum : int
        
        Returns
        -------
        ran3 : float
        
        """
        ran3 = _SurveySubsF95.f90wrap_ran3(idum=idum)
        return ran3
    
    @staticmethod
    def zero2pi(var):
        """
        zero2pi(var)
        
        
        Defined at numutils.fpp lines 333-361
        
        Parameters
        ----------
        var : float
        
        """
        _SurveySubsF95.f90wrap_zero2pi(var=var)
    
    @staticmethod
    def cal2jul(iyyy, mm, dd):
        """
        jul = cal2jul(iyyy, mm, dd)
        
        
        Defined at numutils.fpp lines 363-404
        
        Parameters
        ----------
        iyyy : int
        mm : int
        dd : float
        
        Returns
        -------
        jul : float
        
        """
        jul = _SurveySubsF95.f90wrap_cal2jul(iyyy=iyyy, mm=mm, dd=dd)
        return jul
    
    @staticmethod
    def julday(mm, id, iyyd):
        """
        julday = julday(mm, id, iyyd)
        
        
        Defined at numutils.fpp lines 406-426
        
        Parameters
        ----------
        mm : int
        id : int
        iyyd : int
        
        Returns
        -------
        julday : int
        
        """
        julday = _SurveySubsF95.f90wrap_julday(mm=mm, id=id, iyyd=iyyd)
        return julday
    
    @staticmethod
    def objabs(self, jday, mag, code, gb):
        """
        alpha, h, ra, dec, ierr = objabs(self, jday, mag, code, gb)
        
        
        Defined at numutils.fpp lines 428-507
        
        Parameters
        ----------
        o_p : T_Orb_P
        jday : float
        mag : float
        code : int
        gb : float
        
        Returns
        -------
        alpha : float
        h : float
        ra : float
        dec : float
        ierr : int
        
        """
        alpha, h, ra, dec, ierr = _SurveySubsF95.f90wrap_objabs(o_p=self._handle, \
            jday=jday, mag=mag, code=code, gb=gb)
        return alpha, h, ra, dec, ierr
    
    _dt_array_initialisers = []
    

numutils = Numutils()

class Surveysub(f90wrap.runtime.FortranModule):
    """
    Module surveysub
    
    
    Defined at surveysub.fpp lines 5-454
    
    """
    @staticmethod
    def detos1(self, jday, hx, color, gb, ph, period, amp, surnam, seed):
        """
        flag, ra, dec, d_ra, d_dec, r, delta, m_int, m_rand, eff, isur, Mt, jdayp, ic, \
            surna, h_rand, ierr = detos1(self, jday, hx, color, gb, ph, period, amp, \
            surnam, seed)
        
        
        Defined at surveysub.fpp lines 15-454
        
        Parameters
        ----------
        o_m : T_Orb_M
        jday : float
        hx : float
        color : float array
        gb : float
        ph : float
        period : float
        amp : float
        surnam : str
        seed : int
        
        Returns
        -------
        flag : int
        ra : float
        dec : float
        d_ra : float
        d_dec : float
        r : float
        delta : float
        m_int : float
        m_rand : float
        eff : float
        isur : int
        Mt : float
        jdayp : float
        ic : int
        surna : str
        h_rand : float
        ierr : int
        
        """
        flag, ra, dec, d_ra, d_dec, r, delta, m_int, m_rand, eff, isur, Mt, jdayp, ic, \
            surna, h_rand, ierr = _SurveySubsF95.f90wrap_detos1(o_m=self._handle, \
            jday=jday, hx=hx, color=color, gb=gb, ph=ph, period=period, amp=amp, \
            surnam=surnam, seed=seed)
        return flag, ra, dec, d_ra, d_dec, r, delta, m_int, m_rand, eff, isur, Mt, \
            jdayp, ic, surna, h_rand, ierr
    
    _dt_array_initialisers = []
    

surveysub = Surveysub()

