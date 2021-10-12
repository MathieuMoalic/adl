import math


class square_lattice:
    def __init__(self, parms):
        self.xsize = parms.lattice_param
        self.ysize = parms.lattice_param
        self.coordinates = [(0, 0)]


class rectangular_lattice:
    def __init__(self, parms):
        self.xsize = parms.lattice_param
        self.ysize = parms.lattice_param2
        self.coordinates = [(0, 0)]


class honeycomb_lattice:
    def __init__(self, parms):
        h = parms.lattice_param * math.sin(60 * math.pi / 180)
        self.coordinates = [
            (parms.lattice_param, 0),
            (2 * parms.lattice_param, 0),
            (parms.lattice_param / 2, h),
            (2.5 * parms.lattice_param, h),
            (parms.lattice_param, 2 * h),
            (2 * parms.lattice_param, 2 * h),
        ]
        self.xsize = 3 * parms.lattice_param
        self.ysize = 2 * h


class hexagonal_lattice:
    def __init__(self, parms):
        h = parms.lattice_param * math.sin(60 * math.pi / 180)
        self.coordinates = [
            (0.25 * parms.lattice_param, 1.5 * h),
            (0.75 * parms.lattice_param, 0.5 * h),
            (
                1.25 * parms.lattice_param,
                1.5 * h,
            ),  # these extra 2 nodes prevent bad clipping
            (-0.25 * parms.lattice_param, 0.5 * h),  # needs top and bottom ones too
        ]
        self.xsize = parms.lattice_param
        self.ysize = 2 * h


class octagonal_lattice:
    def __init__(self, parms):
        self.coordinates = [
            (0.5 * parms.lattice_param, 0),
            (1.5 * parms.lattice_param, 0),
            (0.5 * parms.lattice_param, 2 * parms.lattice_param),
            (1.5 * parms.lattice_param, 2 * parms.lattice_param),
            (0, 0.5 * parms.lattice_param),
            (0, 1.5 * parms.lattice_param),
            (2 * parms.lattice_param, 0.5 * parms.lattice_param),
            (2 * parms.lattice_param, 1.5 * parms.lattice_param),
        ]
        self.xsize = 2 * parms.lattice_param
        self.ysize = 2 * parms.lattice_param
