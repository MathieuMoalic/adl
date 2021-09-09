from math import sqrt
import matplotlib as mpl


class square:
    def __init__(self, adl):
        self.s = f"""// Square Antidot
        inner_geom := rect({adl.ad_size:.0f}e-9,{adl.ad_size:.0f}e-9)
        outer_geom := rect({adl.ad_size+adl.ring_width*2:.0f}e-9,{adl.ad_size+adl.ring_width*2:.0f}e-9)
        """

    def get_patch(self, center_x, center_y, ad_size):
        return mpl.patches.Rectangle(
            (center_x - ad_size / 2, center_y - ad_size / 2), ad_size, ad_size
        )


class diamond:
    def __init__(self, adl):
        self.s = f"""// Diamond Antidot
        inner_geom := rect({adl.ad_size:.0f}e-9,{adl.ad_size:.0f}e-9).RotZ(pi/4)
        outer_geom := rect({adl.ad_size+adl.ring_width*2:.0f}e-9,{adl.ad_size+adl.ring_width*2:.0f}e-9).RotZ(pi/4)
        """

    def get_patch(self, center_x, center_y, ad_size):
        r = sqrt(2) * ad_size
        verts = [
            [center_x, center_y + r / 2],
            [center_x - r / 2, center_y],
            [center_x, center_y - r / 2],
            [center_x + r / 2, center_y],
            [center_x, center_y + r / 2],
        ]
        return mpl.patches.PathPatch(mpl.path.Path(verts, closed=True))


class circle:
    def __init__(self, adl):
        self.s = f"""// Circle Antidot
        inner_geom := cylinder({adl.ad_size:.0f}e-9,{adl.dz*adl._Nz:.0f}e-9)
        outer_geom := cylinder({adl.ad_size+adl.ring_width*2:.0f}e-9,{adl.dz*adl._Nz:.0f}e-9)
        """

    def get_patch(self, center_x, center_y, ad_size):
        return mpl.patches.Circle((center_x, center_y), radius=ad_size / 2)


class triangle:
    def __init__(self, adl):
        self.s = f"""// Triangle Antidot
        ad_size := {adl.ad_size:.0f}e-9
        ring_width := {adl.ring_width:.0f}e-9
        ad_size_outer := ad_size + ring_width * 2
        rec_right := rect(ad_size, ad_size * 2).RotZ(pi/6).transl(3/4*ad_size,(2-sqrt(3))*ad_size/4,0)
        rec_left := rect(ad_size, ad_size * 2).RotZ(-pi/6).transl(-3/4*ad_size,(2-sqrt(3))*ad_size/4,0)
        inner_geom := rect(ad_size, ad_size).sub(rec_right).sub(rec_left)
        rec_right = rect(ad_size_outer, ad_size_outer * 2).RotZ(pi/6).transl(3/4*ad_size_outer,(2-sqrt(3))*ad_size_outer/4,0)
        rec_left = rect(ad_size_outer, ad_size_outer * 2).RotZ(-pi/6).transl(-3/4*ad_size_outer,(2-sqrt(3))*ad_size_outer/4,0)
        outer_geom := rect(ad_size_outer, ad_size_outer).sub(rec_right).sub(rec_left).transl(0,4e-9,0)
        """

    def get_patch(self, center_x, center_y, ad_size):
        r = sqrt(3) * ad_size / 3
        verts = [
            [center_x, center_y + r],
            [center_x - ad_size / 2, center_y - r / 2],
            [center_x + ad_size / 2, center_y - r / 2],
            [center_x, center_y + r],
        ]
        return mpl.patches.PathPatch(mpl.path.Path(verts, closed=True))
