import inspect
import os

from . import antidot, lattice, parms

antidots = {
    "square": antidot.square,
    "circle": antidot.circle,
    "triangle": antidot.triangle,
    "diamond": antidot.diamond,
    "squ": antidot.square,
    "cir": antidot.circle,
    "tri": antidot.triangle,
    "dia": antidot.diamond,
}
lattices = {
    "square": lattice.square_lattice,
    "hexagonal": lattice.hexagonal_lattice,
    "rectangular": lattice.rectangular_lattice,
    "honeycomb": lattice.honeycomb_lattice,
    "octagonal": lattice.octagonal_lattice,
    "squ": lattice.square_lattice,
    "hex": lattice.hexagonal_lattice,
    "rec": lattice.rectangular_lattice,
    "hon": lattice.honeycomb_lattice,
    "oct": lattice.octagonal_lattice,
}


class adl:
    def __init__(self, p: parms):
        self._lattice = lattices[p.lattice](p)
        self._s = ""
        p = self.add_dimensions(p)
        self._antidot = antidots[p.antidot]()
        self.add_mesh(p)
        self.add_material(p)
        self.add_geom()
        self.add_static(p)
        self.add_dynamics(p)

    def add_dimensions(self, p: parms):
        def closest_gridsize(size, step):
            out = [(2**i, size / (2**i)) for i in range(13)]
            out = [i for i in out if i[1] > 0.7]
            out = [(i[0], abs(i[1] - step)) for i in out]
            out = sorted(out, key=lambda x: x[1])
            return out[0][0]

        p.Nx = closest_gridsize(self._lattice.xsize, p.dx)
        p.Ny = closest_gridsize(self._lattice.ysize, p.dy)
        p.Nz = 1
        return p

    def add_mesh(self, p: parms):
        if p.mesh == "":
            self._s += f"""
        lattice_param := {p.lattice_param}e-9
        ring_size := {p.ring}e-9
        ad_size := {p.ad_size}e-9"""
            if p.lattice == "rectangular":
                self._s += f"""
        lattice_param2 := {p.lattice_param2}e-9"""
            if p.antidot == "diamond":
                self._s += f"""
        ad_size2 := {p.ad_size2}e-9"""

            self._s += f"""
        PBC := {p.PBC}
        Tx := {self._lattice.xsize:.5f}e-9
        Ty := {self._lattice.ysize:.5f}e-9
        Tz := {p.Nz*p.dz:.5f}e-9
        Nx := {p.Nx}
        Ny := {p.Ny}
        Nz := {p.Nz}
        setgridsize(Nx,Ny,Nz)
        setcellsize(Tx/Nx,Ty/Ny,Tz/Nz)
        setpbc(PBC,PBC,0)
        edgesmooth={p.edgesmooth}
            """
        else:
            self._s += p.mesh

    def add_material(self, p: parms):
        if p.material == "":
            self._s += f"""
        // CoPd film
        msat = {p.msat}
        aex = {p.aex}
        ku1 = {p.ku1}
        anisu = {p.anisu}
        alpha = {p.alpha}
        gammall = {p.gammall}

        """
        else:
            self._s += p.material

    def add_geom(self):
        self._s += self._antidot.s
        self._s += self._lattice.s
        self._s += """
        bulk := universe().sub(rings).sub(ads)
        m.setInShape(bulk,uniform(0,0,1))
        defregion(201,rings)
        defregion(202,ads)
        defregion(203,bulk)
        ku1.setregion(201,0)
        setgeom(bulk.add(rings).sub(ads))
        """

    def add_static(self, p: parms):
        if p.static == "":
            self._s += f"""
        // Static
        phi := {p.phi} * pi / 180
        theta := {p.theta} * pi / 180
        B0 := {p.B0}
        B_ext = vector(B0*sin(theta)*cos(phi), B0*sin(theta)*sin(phi), B0*cos(theta))

        // Relaxation
        maxerr = {p.maxerr_s}
        minimizerstop = {p.minimizerstop}
        relax()
        minimize()
        """
        else:
            self._s += p.static

    def add_dynamics(self, p: parms):
        if p.dynamics == "":
            self._s += f"""
        // Dynamics
        setsolver({p.solver})
        maxdt = {p.maxdt}
        mindt = {p.mindt}
        maxerr = {p.maxerr_d}
        amps:= {p.amps}
        f_cut := {p.f_cut}
        t_sampl := {p.t_sampl}
        t0 := {p.t0}
        """
            if p.Bmask == "":
                self._s += """
        // Bmask
        grainSize  := 20e-9 
        randomSeed := 1234567
        maxRegion  := 60
        ext_makegrains(grainSize, maxRegion, randomSeed)
        for i:=0; i<maxRegion; i++{
            b:=0.1*rand()*1/f_cut
            B_ext.setregion(i, vector(
                B0*sin(theta)*cos(phi)+amps*sinc(2*pi*f_cut*(t-t0+b)),
                B0*sin(theta)*sin(phi)+amps*sinc(2*pi*f_cut*(t-t0+b)),
                B0*cos(theta),
            ))
            defregion(i+maxRegion, ShapeFromRegion(i).intersect(ring))
            Ku1.SetRegion(i+maxRegion, 0)
        }
        saveas(m,"stable")
        saveas(regions,"regions")
        saveas(ku1,"ku1")
        saveas(torque,"torque")
        saveas(b_ext,"b_ext")
        """
            else:
                self._s += """
        // Bmask"""
                self._s += p.Bmask
            self._s += f"""
        // Saving
        run(20/f_cut)
        B_ext.RemoveExtraTerms( )
        B_ext = vector(B0*sin(theta)*cos(phi), B0*sin(theta)*sin(phi), B0*cos(theta))
        run(5/f_cut)
        t = 0
        tableadd(m)
        tableadd(B_ext)
        tableautosave(t_sampl)
        {p.autosave}
        run({p.trun} * t_sampl)
        """
        else:
            self._s += p.dynamics

    def save(self, path: str):
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
