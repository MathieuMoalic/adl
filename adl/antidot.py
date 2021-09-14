from math import sqrt
import matplotlib as mpl


class square:
    def __init__(self, adl):

        self.s = f"""// Square Antidot
        inner_geom := rect({adl.ad_size:.0f}e-9,{adl.ad_size:.0f}e-9)
        outer_geom := rect({adl.ad_size+adl.ring_width*2:.0f}e-9,{adl.ad_size+adl.ring_width*2:.0f}e-9)
        """


class diamond:
    def __init__(self, adl):
        self.s = f"""// Diamond Antidot
        inner_geom := rect({adl.ad_size:.0f}e-9,{adl.ad_size:.0f}e-9).RotZ(pi/4)
        outer_geom := rect({adl.ad_size+adl.ring_width*2:.0f}e-9,{adl.ad_size+adl.ring_width*2:.0f}e-9).RotZ(pi/4)
        """


class circle:
    def __init__(self, adl):
        self.s = f"""// Circle Antidot
        inner_geom := cylinder({adl.ad_size:.0f}e-9,{adl.dz*adl._Nz:.0f}e-9)
        outer_geom := cylinder({adl.ad_size+adl.ring_width*2:.0f}e-9,{adl.dz*adl._Nz:.0f}e-9)
        """


class triangle:
    def __init__(self, adl):
        self.s = f"""// Equilateral Triangle Antidot
        inner_geom := triangle({adl.ad_size:.0f}e-9)
        outer_geom := triangle({adl.ad_size+adl.ring_width*2:.0f}e-9)
        """
