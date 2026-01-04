import math

EPS = 1e-9

def is_close(a, b, eps=EPS):
    return abs(a - b) <= eps

class Point:
    def __init__(self, x, y):
        self.x = float(x)
        self.y = float(y)

    def __repr__(self):
        return f"Point({self.x:.6f}, {self.y:.6f})"

    def distance_to(self, other):
        dx = self.x - other.x
        dy = self.y - other.y
        return math.hypot(dx, dy)

    def translate(self, dx, dy):
        return Point(self.x + dx, self.y + dy)

    def scale(self, k, center=None):
        if center is None:
            center = Point(0, 0)
        return Point(center.x + k * (self.x - center.x),
                     center.y + k * (self.y - center.y))

    def rotate(self, theta_rad, center=None):
        if center is None:
            center = Point(0, 0)
        x = self.x - center.x
        y = self.y - center.y
        c = math.cos(theta_rad)
        s = math.sin(theta_rad)
        xr = c * x - s * y
        yr = s * x + c * y
        return Point(xr + center.x, yr + center.y)

class Line:
    # Ax + By + C = 0
    def __init__(self, A, B, C):
        self.A = float(A)
        self.B = float(B)
        self.C = float(C)
        if abs(self.A) < EPS and abs(self.B) < EPS:
            raise ValueError("Invalid line: A and B cannot both be 0.")

    def __repr__(self):
        return f"Line({self.A:.6f}x + {self.B:.6f}y + {self.C:.6f} = 0)"

    @staticmethod
    def through(p1, p2):
        x1, y1 = p1.x, p1.y
        x2, y2 = p2.x, p2.y
        if is_close(x1, x2) and is_close(y1, y2):
            raise ValueError("Cannot make a line through two identical points.")
        A = y1 - y2
        B = x2 - x1
        C = x1 * y2 - x2 * y1
        return Line(A, B, C)

    def direction_vector(self):
        # For Ax + By + C = 0, a direction vector is (B, -A)
        return (self.B, -self.A)

    def normal_vector(self):
        return (self.A, self.B)

    def contains_point(self, p):
        return abs(self.A * p.x + self.B * p.y + self.C) <= 1e-7

    def intersection_with_line(self, other):
        A1, B1, C1 = self.A, self.B, self.C
        A2, B2, C2 = other.A, other.B, other.C
        det = A1 * B2 - A2 * B1
        if abs(det) < EPS:
            return None
        x = (B1 * C2 - B2 * C1) / det
        y = (C1 * A2 - C2 * A1) / det
        return Point(x, y)

    def perpendicular_through(self, p):
        # line has normal (A,B). A perpendicular line has direction same as normal => normal' = (B, -A)
        A2 = self.B
        B2 = -self.A
        C2 = -(A2 * p.x + B2 * p.y)
        return Line(A2, B2, C2)

    def projection_of_point(self, p):
        # Foot of perpendicular from p to this line
        A, B, C = self.A, self.B, self.C
        denom = A * A + B * B
        t = (A * p.x + B * p.y + C) / denom
        x = p.x - A * t
        y = p.y - B * t
        return Point(x, y)

class Circle:
    def __init__(self, center, radius):
        self.center = center
        self.radius = float(radius)
        if self.radius < 0:
            raise ValueError("Radius must be non-negative.")

    def __repr__(self):
        return f"Circle(center={self.center}, radius={self.radius:.6f})"

    def intersection_with_circle(self, other):
        c1, r1 = self.center, self.radius
        c2, r2 = other.center, other.radius
        d = c1.distance_to(c2)

        if d < EPS and is_close(r1, r2):
            return []  # infinite intersections (same circle) -> return empty for assignment
        if d > r1 + r2 + EPS:
            return []
        if d < abs(r1 - r2) - EPS:
            return []
        if d < EPS:
            return []

        a = (r1 * r1 - r2 * r2 + d * d) / (2 * d)
        h2 = r1 * r1 - a * a
        if h2 < -EPS:
            return []
        h = math.sqrt(max(0.0, h2))

        x0 = c1.x + a * (c2.x - c1.x) / d
        y0 = c1.y + a * (c2.y - c1.y) / d

        rx = -(c2.y - c1.y) * (h / d)
        ry = (c2.x - c1.x) * (h / d)

        p1 = Point(x0 + rx, y0 + ry)
        p2 = Point(x0 - rx, y0 - ry)

        if p1.distance_to(p2) < 1e-7:
            return [p1]
        return [p1, p2]

    def intersection_with_line(self, line):
        # Use projection then check distance from center to line
        foot = line.projection_of_point(self.center)
        dist = foot.distance_to(self.center)

        if dist > self.radius + EPS:
            return []
        if is_close(dist, self.radius):
            return [foot]

        dx, dy = line.direction_vector()
        norm = math.hypot(dx, dy)
        dx /= norm
        dy /= norm

        t = math.sqrt(max(0.0, self.radius * self.radius - dist * dist))
        p1 = Point(foot.x + dx * t, foot.y + dy * t)
        p2 = Point(foot.x - dx * t, foot.y - dy * t)
        return [p1, p2]

class Triangle:
    def __init__(self, p1, p2, p3):
        self.p1, self.p2, self.p3 = p1, p2, p3

    def __repr__(self):
        return f"Triangle({self.p1}, {self.p2}, {self.p3})"

    def side_lengths(self):
        a = self.p2.distance_to(self.p3)
        b = self.p1.distance_to(self.p3)
        c = self.p1.distance_to(self.p2)
        return (a, b, c)

    def translate(self, dx, dy):
        return Triangle(self.p1.translate(dx, dy),
                        self.p2.translate(dx, dy),
                        self.p3.translate(dx, dy))

    def scale(self, k, center=None):
        return Triangle(self.p1.scale(k, center),
                        self.p2.scale(k, center),
                        self.p3.scale(k, center))

    def rotate(self, theta_rad, center=None):
        return Triangle(self.p1.rotate(theta_rad, center),
                        self.p2.rotate(theta_rad, center),
                        self.p3.rotate(theta_rad, center))

def verify_pythagorean(tri, right_vertex_index=None, tol=1e-6):
    pts = [tri.p1, tri.p2, tri.p3]

    def sq(x): return x * x

    if right_vertex_index is None:
        # Try all vertices as the right angle point
        for i in range(3):
            A = pts[i]
            B = pts[(i + 1) % 3]
            C = pts[(i + 2) % 3]
            AB2 = sq(A.distance_to(B))
            AC2 = sq(A.distance_to(C))
            BC2 = sq(B.distance_to(C))
            if abs((AB2 + AC2) - BC2) < tol:
                return True, i
        return False, None

    i = right_vertex_index
    A = pts[i]
    B = pts[(i + 1) % 3]
    C = pts[(i + 2) % 3]
    AB2 = sq(A.distance_to(B))
    AC2 = sq(A.distance_to(C))
    BC2 = sq(B.distance_to(C))
    return (abs((AB2 + AC2) - BC2) < tol), i

if __name__ == "__main__":
    # ----- Example demo to match assignment -----

    # 1) Define objects
    P = Point(0, 0)
    Q = Point(4, 0)
    L = Line.through(P, Q)  # x-axis

    # 2) Intersections
    L2 = Line.through(Point(2, -2), Point(2, 3))  # vertical line x=2
    inter_LL = L.intersection_with_line(L2)
    print("Line-Line intersection:", inter_LL)

    c1 = Circle(Point(0, 0), 5)
    c2 = Circle(Point(6, 0), 5)
    print("Circle-Circle intersections:", c1.intersection_with_circle(c2))

    print("Line-Circle intersections:", c1.intersection_with_line(L2))

    # 3) Perpendicular from an external point
    X = Point(1, 3)  # point not on x-axis
    perp = L.perpendicular_through(X)
    foot = L.projection_of_point(X)
    print("Perpendicular line:", perp)
    print("Foot of perpendicular:", foot)

    # 4) Verify Pythagorean using triangle (point on line, outside point, foot)
    A = Point(3, 0)     # on line L
    B = X               # outside point
    C = L.projection_of_point(B)  # foot on line
    tri = Triangle(A, B, C)
    ok, idx = verify_pythagorean(tri)
    print("Triangle:", tri)
    print("Side lengths (a,b,c):", tri.side_lengths())
    print("Pythagorean holds?", ok, "Right angle at vertex index:", idx)

    # 6) Transformations
    tri2 = tri.translate(2, -1)
    tri3 = tri.scale(1.5, center=Point(0, 0))
    tri4 = tri.rotate(math.radians(30), center=Point(0, 0))
    print("Translated:", tri2)
    print("Scaled:", tri3)
    print("Rotated:", tri4)
