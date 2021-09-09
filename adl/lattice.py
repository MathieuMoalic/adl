import math


class square_lattice:
    def __init__(self, a, b):
        self.coordinates = [(a / 2, a / 2)]
        self.xsize = a
        self.ysize = a


class rectangular_lattice:
    def __init__(self, a, b):
        self.coordinates = [(a / 2, b / 2)]
        self.xsize = a
        self.ysize = b


class honeycomb_lattice:
    def __init__(self, a, b):
        h = a * math.sin(60 * math.pi / 180)
        self.coordinates = [
            (a, 0),
            (2 * a, 0),
            (a / 2, h),
            (2.5 * a, h),
            (a, 2 * h),
            (2 * a, 2 * h),
        ]
        self.xsize = 3 * a
        self.ysize = 2 * h


class hexagonal_lattice:
    def __init__(self, a, b):
        h = a * math.sin(60 * math.pi / 180)
        self.coordinates = [
            (0.25 * a, 1.5 * h),
            (0.75 * a, 0.5 * h),
            (1.25 * a, 1.5 * h),  # these extra 2 nodes prevent bad clipping
            (-0.25 * a, 0.5 * h),  # needs top and bottom ones too
        ]
        self.xsize = a
        self.ysize = 2 * h


class octagonal_lattice:
    def __init__(self, a, b):
        self.coordinates = [
            (0.5 * a, 0),
            (1.5 * a, 0),
            (0.5 * a, 2 * b),
            (1.5 * a, 2 * b),
            (0, 0.5 * b),
            (0, 1.5 * b),
            (2 * a, 0.5 * b),
            (2 * a, 1.5 * b),
        ]
        self.xsize = 2 * a
        self.ysize = 2 * b
