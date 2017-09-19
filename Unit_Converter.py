# Converts num from original unit to millimeters


def convert(num, unit):
    unit_conversions = {
        "em": 4.233333,
        "ex": 2.116666,
        "px": 0.282222,
        "pt": 0.352777,
        "pc": 4.233333,
        "cm": 10,
        "mm": 1,
        "in": 25.4,
        "m": 1000,
        "ft": 304.8,
        "yd": 914.4,
    }
    return unit_conversions[unit] * num
