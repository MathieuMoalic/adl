from math import sqrt
import matplotlib as mpl


class square:
    def __init__(self, ad_size, ring_width):
        self.s = f"""// Square
        ad_size := {ad_size:.0f}e-9
        ring_width := {ring_width:.0f}e-9
        inner_geom := rect(ad_size,ad_size)
        outer_geom := rect(ad_size + ring_width * 2,ad_size + ring_width * 2)
        """

    def get_patch(self, center_x, center_y, ad_size):
        return mpl.patches.Rectangle(
            (center_x - ad_size / 2, center_y - ad_size / 2), ad_size, ad_size
        )


class diamond:
    def __init__(self, ad_size, ring_width):
        self.s = f"""// Diamond
        ad_size := {ad_size:.0f}e-9
        ring_width := {ring_width:.0f}e-9
        inner_geom := rect(ad_size,ad_size).RotZ(pi/4)
        outer_geom := rect(ad_size + ring_width * 2,ad_size + ring_width * 2).RotZ(pi/4)
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
    def __init__(self, ad_size, ring_width):
        self.s = f"""// Circle
        ad_size := {ad_size:.0f}e-9
        ring_width := {ring_width:.0f}e-9
        inner_geom := cylinder(ad_size,Nz*dz)
        outer_geom := cylinder(ad_size + ring_width * 2,Nz*dz)
        """

    def get_patch(self, center_x, center_y, ad_size):
        return mpl.patches.Circle((center_x, center_y), radius=ad_size / 2)


class triangle:
    def __init__(self, ad_size, ring_width):
        self.s = f"""// Triangle
        ad_size := {ad_size:.0f}e-9
        ring_width := {ring_width:.0f}e-9
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
