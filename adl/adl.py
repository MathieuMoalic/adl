import inspect

from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np

from . import antidot
from . import lattice

antidots = {
    "square": antidot.square,
    "circle": antidot.circle,
    "triangle": antidot.triangle,
    "diamond": antidot.diamond,
}
lattices = {
    "square": lattice.square_lattice,
    "hexagonal": lattice.hexagonal_lattice,
    "rectangular": lattice.rectangular_lattice,
    "honeycomb": lattice.honeycomb_lattice,
    "octagonal": lattice.octagonal_lattice,
}


class adl:
    def __init__(
        self,
        lattice_param,
        ad_size,
        ring,
        lattice,
        antidot,
        lattice_param2=0,
        ad_size2=0,
        dx=1,
        dy=1,
        dz=13.2,
        pbc=32,
        angle=10,
        B0=0.223,
        amps="B0/100",
        f_cut="25e9",
        delay="2.5 / f_cut",
        t_sampl="0.5 / (f_cut * 1.4)",
        maxerr="0.25e-7",
        minimizerstop="5e-6",
        maxdt="t_sampl/100",
        msat="810e3",
        aex="13e-12",
        ku1="453195",
        anisu="vector(0,0,1)",
        alpha="0.000000001",
        gammall="187e9",
        m="uniform(0, 0, 1)",
        autosave="autosave(m,t_sampl)",
        bextx="B0*sin(angle)",
        bexty="0",
        bextz="amps*sinc(2*pi*f_cut*(t-delay-tp)) + B0*cos(angle)",
    ):
        # input parameters, can be changed
        self.lattice_param = lattice_param
        self.ad_size = ad_size
        self.ring_width = ring
        self.lattice_name = lattice
        self.antidot_name = antidot
        if lattice == "rectangular":
            if lattice_param2 == 0:
                raise ValueError(
                    "'lattice_param2' is needed when using a rectangular lattice"
                )
            else:
                self.lattice_param2 = lattice_param2
        self.ad_size2 = ad_size2
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.PBC = pbc
        self.angle = angle
        self.B0 = B0
        self.amps = amps
        self.f_cut = f_cut
        self.delay = delay
        self.t_sampl = t_sampl
        self.maxerr = maxerr
        self.minimizerstop = minimizerstop
        self.maxdt = maxdt
        self.msat = msat
        self.aex = aex
        self.ku1 = ku1
        self.anisu = anisu
        self.alpha = alpha
        self.gammall = gammall
        self.m = m
        self.autosave = autosave
        self.bextx = bextx
        self.bexty = bexty
        self.bextz = bextz

        self.make()

    def make(self):
        self._s = ""
        self._lattice = lattices[self.lattice_name](self)
        self._Nx = int(self._lattice.xsize / self.dx)
        self._Ny = int(self._lattice.ysize / self.dy)
        self._Nz = 1
        new_nx = self._Nx - self._Nx % 5
        self.dx = self._Nx / new_nx
        self._Nx = new_nx
        new_ny = self._Ny - self._Ny % 5
        self.dy = self._Ny / new_ny
        self._Ny = new_ny
        self._antidot = antidots[self.antidot_name](self)
        self.pre()
        self.geom()
        self.post()

    def save(self, path):
        self._s = inspect.cleandoc(self._s)  # removes padding
        if path[-4:] == ".mx3":
            with open(path, "w") as f:
                f.writelines(self._s)
        else:
            name = f"{self.lattice_name}_lat_{self.antidot_name}_ad_ring_{self.ring_width}nm_ad_{self.ad_size}nm_a_{self.lattice_param}nm"
            with open(f"{path}/{name}.mx3", "w") as f:
                f.writelines(self._s)

    def pre(self):
        self._s += f"""
        // lattice:             {self.lattice_name}
        // antidot:             {self.antidot_name}
        // antidot size:        {self.ad_size} nm
        // ring width:          {self.ring_width} nm
        // lattice parameter:   {self.lattice_param} nm 
        
        setgridsize({self._Nx},{self._Ny}, {self._Nz})
        setcellsize({self.dx:.5f}e-9, {self.dy:.5f}e-9, {self.dz}e-9)
        setPBC({self.PBC}, {self.PBC}, 0)
        edgesmooth=0

        // CoPd film
        adl := Universe()
        m = {self.m}
        msat = {self.msat}
        aex = {self.aex}
        ku1 = {self.ku1}
        anisu = {self.anisu}
        alpha = {self.alpha}
        gammall = {self.gammall}

        """

    def geom(self):
        self._s += self._antidot.s
        new_origin_x = self._Nx * self.dx / 2
        new_origin_y = self._Ny * self.dy / 2

        self._s += """
        // Lattice
            """
        for i, (x, y) in enumerate(self._lattice.coordinates):
            if i == 0:
                q = ":"
            else:
                q = ""
            new_x = f"{(x-new_origin_x):.0f}e-9"
            new_y = f"{(y-new_origin_y):.0f}e-9"
            self._s += f"""
        // {self.antidot_name.capitalize()} at ({x:.0f}nm,{y:.0f}nm)
        inner_dot {q}= inner_geom.transl({new_x},{new_y},0)
        outer_dot {q}= outer_geom.transl({new_x},{new_y},0)
        adl = adl.add(outer_dot).sub(inner_dot)
        m.setInShape(outer_dot, vortex(1, -1).transl({new_x},{new_y},0))
        defregion(1, outer_dot)
        Ku1.SetRegion(1, 0)
                    """.replace(
                ".transl(0e-9,0e-9,0)", ""
            )

    def post(self):
        self._s += f"""
        setgeom(adl)
        
        // Static Field
        angle := {self.angle} * pi / 180
        B0 := {self.B0}
        B_ext = vector(B0*sin(angle), 0, B0*cos(angle))

        // Dynamics Params
        amps := {self.amps}
        f_cut := {self.f_cut}
        delay := {self.delay}
        t_sampl := {self.t_sampl}

        // Solver Params
        maxerr = {self.maxerr}
        minimizerstop = {self.minimizerstop}
        maxdt = {self.maxdt}

        // Relaxation
        relax()
        run(1e-10)
        relax()
        minimize()
        relax()
        tp:=t
        saveas(m,"stable")
        snapshotas(m,"stable.png")
        
        // Dynamics
        B_ext = vector( B0*sin(angle), amps*sinc(2*pi*f_cut*(t-delay-tp)), B0*cos(angle))
        tableadd(B_ext)
        tableautosave(t_sampl)
        {self.autosave}
        run(1500 * t_sampl)
        """
