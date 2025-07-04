# Copyright (c) 2010 Andrew Brown <brownan@cs.duke.edu, brownan@gmail.com>
# See LICENSE.txt for license terms

class GF256int(int):
    """Instances of this object are elements of the field GF(2^8)
    Instances are integers in the range 0 to 255
    This field is defined using the irreducable polynomial
    x^8 + x^4 + x^3 + x + 1
    and using 3 as the generator for the exponent table and log table.
    """
    # Maps integers to GF256int instances
    cache = {}
    # Exponent table for 3, a generator for GF(256)
    exptable = (1, 3, 5, 15, 17, 51, 85, 255, 26, 46, 114, 150, 161, 248, 19,
            53, 95, 225, 56, 72, 216, 115, 149, 164, 247, 2, 6, 10, 30, 34,
            102, 170, 229, 52, 92, 228, 55, 89, 235, 38, 106, 190, 217, 112,
            144, 171, 230, 49, 83, 245, 4, 12, 20, 60, 68, 204, 79, 209, 104,
            184, 211, 110, 178, 205, 76, 212, 103, 169, 224, 59, 77, 215, 98,
            166, 241, 8, 24, 40, 120, 136, 131, 158, 185, 208, 107, 189, 220,
            127, 129, 152, 179, 206, 73, 219, 118, 154, 181, 196, 87, 249, 16,
            48, 80, 240, 11, 29, 39, 105, 187, 214, 97, 163, 254, 25, 43, 125,
            135, 146, 173, 236, 47, 113, 147, 174, 233, 32, 96, 160, 251, 22,
            58, 78, 210, 109, 183, 194, 93, 231, 50, 86, 250, 21, 63, 65, 195,
            94, 226, 61, 71, 201, 64, 192, 91, 237, 44, 116, 156, 191, 218,
            117, 159, 186, 213, 100, 172, 239, 42, 126, 130, 157, 188, 223,
            122, 142, 137, 128, 155, 182, 193, 88, 232, 35, 101, 175, 234, 37,
            111, 177, 200, 67, 197, 84, 252, 31, 33, 99, 165, 244, 7, 9, 27,
            45, 119, 153, 176, 203, 70, 202, 69, 207, 74, 222, 121, 139, 134,
            145, 168, 227, 62, 66, 198, 81, 243, 14, 18, 54, 90, 238, 41, 123,
            141, 140, 143, 138, 133, 148, 167, 242, 13, 23, 57, 75, 221, 124,
            132, 151, 162, 253, 28, 36, 108, 180, 199, 82, 246, 1)

    # Logarithm table, base 3
    logtable = (None, 0, 25, 1, 50, 2, 26, 198, 75, 199, 27, 104, 51, 238, 223,
            3, 100, 4, 224, 14, 52, 141, 129, 239, 76, 113, 8, 200, 248, 105,
            28, 193, 125, 194, 29, 181, 249, 185, 39, 106, 77, 228, 166, 114,
            154, 201, 9, 120, 101, 47, 138, 5, 33, 15, 225, 36, 18, 240, 130,
            69, 53, 147, 218, 142, 150, 143, 219, 189, 54, 208, 206, 148, 19,
            92, 210, 241, 64, 70, 131, 56, 102, 221, 253, 48, 191, 6, 139, 98,
            179, 37, 226, 152, 34, 136, 145, 16, 126, 110, 72, 195, 163, 182,
            30, 66, 58, 107, 40, 84, 250, 133, 61, 186, 43, 121, 10, 21, 155,
            159, 94, 202, 78, 212, 172, 229, 243, 115, 167, 87, 175, 88, 168,
            80, 244, 234, 214, 116, 79, 174, 233, 213, 231, 230, 173, 232, 44,
            215, 117, 122, 235, 22, 11, 245, 89, 203, 95, 176, 156, 169, 81,
            160, 127, 12, 246, 111, 23, 196, 73, 236, 216, 67, 31, 45, 164,
            118, 123, 183, 204, 187, 62, 90, 251, 96, 177, 134, 59, 82, 161,
            108, 170, 85, 41, 157, 151, 178, 135, 144, 97, 190, 220, 252, 188,
            149, 207, 205, 55, 63, 91, 209, 83, 57, 132, 60, 65, 162, 109, 71,
            20, 42, 158, 93, 86, 242, 211, 171, 68, 17, 146, 217, 35, 32, 46,
            137, 180, 124, 184, 38, 119, 153, 227, 165, 103, 74, 237, 222, 197,
            49, 254, 24, 13, 99, 140, 128, 192, 247, 112, 7)

    def __new__(cls, value):
        # Check cache
        # Caching sacrifices a bit of speed for less memory usage. This way,
        # there are only a max of 256 instances of this class at any time.
        try:
            return GF256int.cache[value]
        except KeyError:
            if value > 255 or value < 0:
                raise ValueError("Field elements of GF(2^8) are between 0 and 255. Cannot be %s" % value)

            newval = int.__new__(cls, value)
            GF256int.cache[int(value)] = newval
            return newval

    def __add__(self, b):
        "Addition in GF(2^8) is the xor of the two"
        return GF256int(self ^ b)
    __sub__ = __add__
    __radd__ = __add__
    __rsub__ = __add__
    def __neg__(self):
        return self

    def __mul__(self, b):
        "Multiplication in GF(2^8)"
        if self == 0 or b == 0:
            return GF256int(0)
        x = GF256int.logtable[self]
        y = GF256int.logtable[b]
        z = (x + y) % 255
        return GF256int(GF256int.exptable[z])
    __rmul__ = __mul__

    def __pow__(self, power): #type: ignore
        if isinstance(power, GF256int):
            raise TypeError("Raising a Field element to another Field element is not defined. power must be a regular integer")
        x = GF256int.logtable[self]
        z = (x * power) % 255
        return GF256int(GF256int.exptable[z])

    def inverse(self):
        """Returns the multiplicative inverse of this field element"""
        if self == 0:
            raise ZeroDivisionError("Cannot find the inverse of 0")

        # Check if element exists in logtable
        # Patch for GF256int(255) - it should map to log value 7
        if int(self) == 255:
            return GF256int(GF256int.exptable[255 - 7])

        if self not in GF256int.logtable:
            # Fallback implementation if logtable lookup fails
            # Find inverse by brute force (slow but reliable)
            for i in range(1, 256):
                if (self * GF256int(i)) % 256 == 1:
                    return GF256int(i)
            raise ValueError(f"Element {self} has no inverse in GF(256)")

        e = GF256int.logtable[self]
        if e is None:
            raise ValueError(f"Element {self} has no log value in table")

        return GF256int(GF256int.exptable[255 - e])

    def __div__(self, other):
        return self * GF256int(other).inverse()
    def __rdiv__(self, other):
        return self.inverse() * other

    def __repr__(self):
        n = self.__class__.__name__
        return "%s(%r)" % (n, int(self))

    def multiply(self, other):
        """A slow multiply method. This method gives the same results as the
        other multiply method, but is implemented to illustrate how it works
        and how the above tables were generated.

        This procedure is called Peasant's Algorithm (I believe)
        """
        a = int(self)
        b = int(other)

        p = a
        r = 0
        while b:
            if b & 1: r = r ^ p
            b = b >> 1
            p = p << 1
            if p & 0x100: p = p ^ 0x11b

        return GF256int(r)
