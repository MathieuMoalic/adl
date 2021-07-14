import inspect

from matplotlib import pyplot as plt
import matplotlib as mpl

import antidot
import lattice

antidots = {
    "square": antidot.square,
    "circle": antidot.circle,
    "triangle": antidot.triangle,
    "diamond": antidot.diamond,
}
lattices = {
    "square": lattice.square_lattice,
    "honeycomb": lattice.honeycomb_lattice,
    "hexagonal": lattice.hexagon_lattice,
}


class adl:
    def __init__(
        self, a=100, ad=40, ring=10, reps=4, lattice="square", antidot="square"
    ):
        self.s = ""
        self.a = a
        self.ad_size = ad
        self.ring_width = ring
        self.reps = reps
        self.lattice_name = lattice
        self.antidot_name = antidot
        self.lattice = lattices[lattice](self.a, self.reps)
        self.antidot = antidots[antidot](self.ad_size, self.ring_width)
        self._Nx = self.reps * self.a
        self._Ny = self.reps * self.a
        self._Nz = 1
        self.dx = 1
        self.dy = 1
        self.dz = 13.4
        self.PBC = 128
        self.make()

    def make(self):
        self.pre()
        self.geom()
        self.post()

    def save(self, dir):
        name = f"{self.lattice_name}_lat_{self.antidot_name}_ad_{self.reps}x{self.reps}_ring_{self.ring_width}nm_ad_{self.ad_size}nm_a_{self.a}nm"
        self.s = inspect.cleandoc(self.s)
        with open(f"{dir}/{name}.mx3", "w") as f:
            f.writelines(self.s)

    def preview(self):
        fig, ax = plt.subplots(figsize=(3, 3))
        ax.set_xlim(0, self.a * self.reps)
        ax.set_ylim(0, self.a * self.reps)
        ax.patches.append(
            mpl.patches.Rectangle(
                (0, 0),
                self.a * self.reps,
                self.a * self.reps,
                color="k",
                transform=ax.transData,
            )
        )
        for x, y in self.lattice.coordinates:
            patch = self.antidot.get_patch(x, y, self.ad_size + self.ring_width)
            patch.set(transform=ax.transData, lw=0, facecolor="gray")
            ax.patches.append(patch)
        for x, y in self.lattice.coordinates:
            patch = self.antidot.get_patch(x, y, self.ad_size)
            patch.set(transform=ax.transData, lw=0, facecolor="w")
            ax.patches.append(patch)
        ax.set(xlabel="x (nm)", ylabel="y (nm)")
        ax.tick_params(which="both", direction="out")
        fig.tight_layout()

    def pre(self):
        self.s += f"""
        dx := {self.dx}e-9
        dy := {self.dy}e-9
        dz := {self.dz}e-9
        Nx := {self._Nx}
        Ny := {self._Ny}
        Nz := {self._Nz}
        setgridsize(Nx,Ny, Nz)
        setcellsize(dx, dy, dz)
        setPBC({self.PBC}, {self.PBC}, 0)

        // CoPd stripe
        adl := rect(Nx*dx,Ny*dy)
        m = uniform(0, 0, 1)
        Msat = 810e3
        aex = 13e-12
        alpha = 0.01
        Ku1 = 453195
        anisU = vector(0, 0, 1)
        alpha = 0.01

        """

    def geom(self):
        self.s += self.antidot.s

        self.s += """
        inner_dot := inner_geom // just initializing the variable
        outer_dot := outer_geom // just initializing the variable
            """
        for (x, y) in self.lattice.coordinates:
            new_origin = self._Nx * self.dx / 2
            x = f"{(x-new_origin):.0f}e-9"
            y = f"{(new_origin-y):.0f}e-9"
            self.s += f"""
        inner_dot = inner_geom.transl({x},{y},0)
        outer_dot = outer_geom.transl({x},{y},0)
        adl = adl.add(outer_dot).sub(inner_dot)
        m.setInShape(outer_dot, vortex(1, -1).transl({x},{y},0))
        defregion(1, outer_dot)
        Ku1.SetRegion(1, 0)

                    """

    def post(self):
        self.s += """
        Amps := 0.5e-4
        f_cut := 20e09
        t_sampl := 0.5 / (f_cut * 1.4)
        t0 := 50 / f_cut
        Maxdt = t_sampl / 100

        angle := 0 * pi / 180
        B0 := 0.223
        B_ext = vector(B0*sin(angle), 0, B0*cos(angle))


        delay := 2.5 / f_cut     

        setgeom(adl)

        maxerr = 0.0000010976
        MinimizerStop = 5e-7

        minimize()
        relax()
        tableadd(B_ext)
        tableadd(e_total)
        TableAutoSave(t_sampl)
        alpha=0.5
        setsolver(4)
        run(2.5e-9)

        setsolver(3)
        relax()
        alpha = 0.01
        saveAs(m, "stable")
        setsolver(4)
        tp:=t


        B_ext = vector( B0*sin(angle) + Amps * sinc(2*pi*f_cut*(t-delay-tp)),  Amps * sinc(2*pi*f_cut*(t-delay-tp)), B0*cos(angle))
        AutoSave(m, t_sampl)
        run(1500 * t_sampl)
        """
