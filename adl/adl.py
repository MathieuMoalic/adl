import inspect
import os

from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np

from . import antidot, lattice, parms

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
    def __init__(self, parms: parms):
        self._lattice = lattices[parms.lattice](parms)
        self._s = ""
        parms = self.add_dimensions(parms)
        self._antidot = antidots[parms.antidot](parms)
        self.add_mesh(parms)
        self.add_material(parms)
        self.add_geom(parms)
        self.add_static(parms)
        self.add_dynamics(parms)

    def add_dimensions(self, parms):
        parms.Nx = int(self._lattice.xsize / parms.dx)
        parms.Ny = int(self._lattice.ysize / parms.dy)
        parms.Nz = 1
        new_nx = parms.Nx - parms.Nx % 5
        parms.dx = parms.Nx / new_nx * parms.dx
        parms.Nx = new_nx
        new_ny = parms.Ny - parms.Ny % 5
        parms.dy = parms.Ny / new_ny * parms.dy
        parms.Ny = new_ny
        self._s += f"""
        Nx := {parms.Nx}
        Ny := {parms.Ny}
        Nz := {parms.Nz}"""
        return parms

    def add_mesh(self, parms):
        if parms.mesh == "":
            self._s += f"""
        SetMesh(Nx,Ny,Nz,{parms.dx:.5f}e-9, {parms.dy:.5f}e-9, {parms.dz}e-9,{parms.PBC}, {parms.PBC}, 0)
        edgesmooth={parms.edgesmooth}
        """
        else:
            self._s += parms.mesh

    def add_material(self, parms):
        if parms.material == "":
            self._s += f"""
        // CoPd film
        msat = {parms.msat}
        aex = {parms.aex}
        ku1 = {parms.ku1}
        anisu = {parms.anisu}
        alpha = {parms.alpha}
        gammall = {parms.gammall}

        // Geom
        adl := Universe()
        m = {parms.m}
        """
        else:
            self._s += parms.material

    def add_geom(self, parms):
        self._s += self._antidot.s
        for i, (x, y) in enumerate(self._lattice.coordinates):
            if i == 0:
                q = ":"
            else:
                q = ""
            self._s += f"""
        inner_dot {q}= inner_geom.transl({x}e-9,{y}e-9,0)
        outer_dot {q}= outer_geom.transl({x}e-9,{y}e-9,0)
        adl = adl.add(outer_dot).sub(inner_dot)
        m.setInShape(outer_dot, vortex(1, 1).transl({x}e-9,{y}e-9,0))
        defregion(1, outer_dot)
        Ku1.SetRegion(1, 0)
                    """.replace(
                ".transl(0e-9,0e-9,0)", ""
            )
        self._s += """
        setgeom(adl)
        """

    def add_static(self, parms):
        if parms.static == "":
            self._s += f"""
        // Static
        angle := {parms.angle} * pi / 180
        B0 := {parms.B0}
        B_ext = vector(B0*sin(angle), 0, B0*cos(angle))

        // Relaxation
        maxerr = {parms.maxerr_s}
        minimizerstop = {parms.minimizerstop}
        relaxtorquethreshold = {parms.relaxtorquethreshold}
        minimize()
        saveas(m,"stable")
        snapshotas(m,"stable.png")
        """
        else:
            self._s += parms.static

    def add_dynamics(self, parms):
        if parms.dynamics == "":
            self._s += f"""
        // Dynamics
        setsolver({parms.solver})
        maxdt = {parms.maxdt}
        mindt = {parms.mindt}
        maxerr = {parms.maxerr_d}
        amps:= {parms.amps}
        f_cut := {parms.f_cut}
        t_sampl := {parms.t_sampl}
        t0 := {parms.t0}
        """
            if parms.Bmask == "":
                self._s += """
        // Bmask
        B_mask:=newSlice(3, Nx, Ny, Nz)
        Bxyz:= 0.0
        for x:=0; x<Nx; x++{{
            for y:=0; y<Ny; y++{{
                for z:=0; z<Nz; z++{{
                    Bxyz= randNorm()
                    B_mask.set(0, x, y, z, Bxyz)
                    B_mask.set(1, x, y, z, Bxyz)
                    B_mask.set(2, x, y, z, Bxyz)
                }}
            }}
        }}
        B_ext.add(B_mask, amps*sinc(2*pi*f_cut*(t-t0)))
        """
            else:
                self._s += parms.Bmask
            self._s += f"""
        // Saving
        tableadd(B_ext)
        tableautosave(t_sampl)
        {parms.autosave}
        run({parms.trun} * t_sampl)
        """
        else:
            self._s += parms.dynamics

    def save(self, path):
        import time

        self._s = inspect.cleandoc(self._s)  # removes padding
        if path[-4:] == ".mx3":
            with open(path, "w") as f:
                f.writelines(self._s)
        else:
            i = 0
            while True:
                mx3_path = f"{path}/adl_{i}.mx3"
                if not os.path.exists(mx3_path):
                    print(f"Saved as '{mx3_path}'")
                    with open(mx3_path, "w") as f:
                        f.writelines(self._s)
                    break
                i += 1
