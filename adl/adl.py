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
        dx=0.75,
        dy=0.75,
        dz=13.2,
        pbc=32,
        angle="0.0001",
        B0="0.223",
        amps="9e-2",
        f_cut="25e9",
        t0="t + 10/f_cut",
        t_sampl="0.5 / (f_cut * 1.5)",
        maxerr="0.25e-7",
        minimizerstop="5e-6",
        maxdt="t_sampl/100",
        msat="810e3",
        aex="13e-12",
        ku1="453195",
        anisu="vector(0,0,1)",
        alpha="0.02",
        gammall="187e9",
        m="uniform(1e-5, 1e-5, 1)",
        autosave="autosave(m,t_sampl)",
        bextx="B0*sin(angle)",
        bexty="0",
        bextz="amps*sinc(2*pi*f_cut*(t-delay-tp)) + B0*cos(angle)",
        edgesmooth=3,
        trun="1500",
        static="",
        dyn="",
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
        self.edgesmooth = edgesmooth
        self.trun = trun
        self.t0 = t0
        self.static = static
        self.dyn = dyn

        self.make()

    def make(self):
        self._s = ""
        self._lattice = lattices[self.lattice_name](self)
        self._Nx = int(self._lattice.xsize / self.dx)
        self._Ny = int(self._lattice.ysize / self.dy)
        self._Nz = 1
        new_nx = self._Nx - self._Nx % 5
        self.dx = self._Nx / new_nx * self.dx
        self._Nx = new_nx
        new_ny = self._Ny - self._Ny % 5
        self.dy = self._Ny / new_ny * self.dy
        self._Ny = new_ny
        self._antidot = antidots[self.antidot_name](self)
        self.pre()
        self.geom()
        self.add_static()
        self.add_dyn()

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
        Nx := {self._Nx}
        Ny := {self._Ny}
        Nz := {self._Nz}
        SetMesh(Nx,Ny,Nz,{self.dx:.5f}e-9, {self.dy:.5f}e-9, {self.dz}e-9,{self.PBC}, {self.PBC}, 0)
        edgesmooth={self.edgesmooth}

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

        self._s += """
        // Lattice
            """
        for i, (x, y) in enumerate(self._lattice.coordinates):
            if i == 0:
                q = ":"
            else:
                q = ""
            self._s += f"""
        // {self.antidot_name.capitalize()} at ({x:.0f}nm,{y:.0f}nm)
        inner_dot {q}= inner_geom.transl({x}e-9,{y}e-9,0)
        outer_dot {q}= outer_geom.transl({x}e-9,{y}e-9,0)
        adl = adl.add(outer_dot).sub(inner_dot)
        m.setInShape(outer_dot, vortex(1, 1).transl({x}e-9,{y}e-9,0))
        defregion(1, outer_dot)
        Ku1.SetRegion(1, 0)
                    """.replace(
                ".transl(0e-9,0e-9,0)", ""
            )

    def add_static(self):
        if self.static != "":
            self._s += self.static
        else:
            self._s += f"""
        // Static
        angle := 0 * pi / 180
        B0 := 0.223
        B_ext = vector(B0*sin(angle), 0, B0*cos(angle))

        // Relaxation
        MaxErr = 4e-6
        minimizerstop = 1e-8
        RelaxTorqueThreshold = 1e-6
        minimize()
        saveas(m,"stable")
        snapshotas(m,"stable.png")
        """

    def add_dyn(self):
        if self.dyn != "":
            self._s += self.dyn
        else:
            self._s += f"""
        // Dynamics
        f_cut := 20e9
        t_sampl := 0.5 / (f_cut * 1.5)
        amps:= B0/2

        setsolver(5)
        t0 := t + 10/f_cut
        maxdt = 1.42e-12
        mindt = 1e-14
        maxerr = 1e-7

        Bphi := 10
        Bvalue := 1.0
        maskB := newSlice(3, Nx, Ny, Nz)
        theta := 0.0
        for x:=0; x<Nx; x++{{
            for y:=0; y<Ny; y++ {{
                for z:=0; z<Nz; z++ {{
                    r := index2coord(x, y, 0)
                    x2 := r.X()
                    y2 := r.Y()
                    theta = atan(y2/x2)*sign(x2)
                    maskB.set(0, x, y, z, 0)
                    maskB.set(1, x, y, z, Bvalue*sin(theta))
                    maskB.set(2, x, y, z, 0)
                }}
            }}
        }}


        B_ext.add(mask_MF, amps*sinc(2*pi*f_cut*(t-t0)))
        tableadd(B_ext)
        tableautosave(t_sampl)
        autosave(m,t_sampl)
        run(400 * t_sampl)
        """
