import inspect

from matplotlib import pyplot as plt
import matplotlib as mpl

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
        a=100,
        ad=40,
        ring=10,
        lattice="square",
        antidot="square",
        dx=1,
        dy=1,
        dz=13.2,
        pbc=32,
        b=None,
        **kwargs,
    ):
        # input parameters, can be changed
        self.a = a
        self.b = b
        self.ad_size = ad
        self.ring_width = ring
        self.lattice_name = lattice
        self.antidot_name = antidot
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.PBC = pbc
        self.kwargs = kwargs

        self.make()

    def make(self):
        self._s = ""
        self._lattice = lattices[self.lattice_name](self.a, self.b)
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
            name = f"{self.lattice_name}_lat_{self.antidot_name}_ad_ring_{self.ring_width}nm_ad_{self.ad_size}nm_a_{self.a}nm"
            with open(f"{path}/{name}.mx3", "w") as f:
                f.writelines(self._s)

    def preview(self):
        figshape = (3, self._lattice.ysize / self._lattice.xsize * 3)
        fig, ax = plt.subplots(figsize=figshape)
        ax.set_xlim(0, self._lattice.xsize)
        ax.set_ylim(0, self._lattice.ysize)
        ax.patches.append(
            mpl.patches.Rectangle(
                (0, 0),
                self._lattice.xsize,
                self._lattice.ysize,
                color="k",
                transform=ax.transData,
            )
        )
        for x, y in self._lattice.coordinates:
            patch = self._antidot.get_patch(x, y, self.ad_size + self.ring_width)
            patch.set(transform=ax.transData, lw=0, facecolor="gray")
            ax.add_patch(patch)
        for x, y in self._lattice.coordinates:
            patch = self._antidot.get_patch(x, y, self.ad_size)
            patch.set(transform=ax.transData, lw=0, facecolor="w")
            ax.add_patch(patch)
        ax.set(xlabel="x (nm)", ylabel="y (nm)")
        ax.tick_params(which="both", direction="out")
        fig.tight_layout()

    def pre(self):
        self._s += f"""
        // lattice:             {self.lattice_name}
        // antidot:             {self.antidot_name}
        // antidot size:        {self.ad_size} nm
        // ring width:          {self.ring_width} nm
        // lattice parameter:   {self.a} nm 
        
        setgridsize({self._Nx},{self._Ny}, {self._Nz})
        setcellsize({self.dx:.5f}e-9, {self.dy:.5f}e-9, {self.dz}e-9)
        setPBC({self.PBC}, {self.PBC}, 0)
        edgesmooth=0

        // CoPd stripe
        adl := Universe()
        m = uniform(0, 0, 1)
        Msat = 810e3
        aex = 13e-12
        Ku1 = 453195
        anisU = vector(0, 0, 1)
        alpha = 0.000000001
        gammaLL = 187e9

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
        angle := {self.kwargs.get('angle',"0")} * pi / 180
        B0 := {self.kwargs.get('B0',"0.223")}
        B_ext = vector(B0*sin(angle), 0, B0*cos(angle))

        // Dynamics Params
        amps := {self.kwargs.get('amps',"B0 / 100")}
        f_cut := {self.kwargs.get('f_cut',"25e9")}
        delay := {self.kwargs.get('delay',"2.5 / f_cut     ")}
        t_sampl := {self.kwargs.get('t_sampl',"0.5 / (f_cut * 1.4)")}

        // Solver Params
        maxerr = {self.kwargs.get('maxerr',"0.25e-7")}
        minimizerStop = {self.kwargs.get('minimizerStop',"5e-6")}
        maxdt = {self.kwargs.get('maxdt',"t_sampl / 100")}

        // Relaxation
        relax()
        run(1e-10)
        relax()
        minimize()
        relax()
        tp:=t
        saveas(m,"stable")
        
        // Dynamics
        B_ext = vector( B0*sin(angle), 0, amps*sinc(2*pi*f_cut*(t-delay-tp)) + B0*cos(angle))
        tableadd(B_ext)
        tableautosave(t_sampl)
        autosave(m,t_sampl)
        run(1500 * t_sampl)
        """
