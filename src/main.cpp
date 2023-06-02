////////////////////////////////////////////////////////////////////////////////
#include <algorithm>
#include <complex>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>
////////////////////////////////////////////////////////////////////////////////

const std::string root_path = DATA_DIR;

typedef std::complex<double> Point;
typedef std::vector<Point> Polygon;


double inline det(const Point& u, const Point& v)
{
    double determinant = u.real() * v.imag() - u.imag() * v.real();
    return determinant;
}

// Return true iff [a,b] intersects [c,d], and store the intersection in ans
bool intersect_segment(const Point &a, const Point &b, const Point &c, const Point &d, Point &ans)
{
    Point direction_ab = b - a;
    Point direction_cd = d - c;
    double detvalue = det(direction_ab, direction_cd);
    if (detvalue == 0)
    {
        return false;
    }

    double t = det(a - c, direction_cd) / detvalue;
    double u = det(a - c, direction_ab) / detvalue;

    if (t >= 0 && t <= 1 && u >= 0 && u <= 1)
    {
      
        ans = a + direction_ab * t;
        return true; 
    }
    return false;
    
}

////////////////////////////////////////////////////////////////////////////////

bool is_inside(const Polygon& poly, const Point& query)
{
    double minX = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::lowest();
    double minY = std::numeric_limits<double>::max();
    double maxY = std::numeric_limits<double>::lowest();

    for (const auto& point : poly)
    {
        minX = std::min(minX, point.real());
        maxX = std::max(maxX, point.real());
        minY = std::min(minY, point.real());
        maxY = std::max(maxY, point.real());
    }

    Point outside(maxX + 1, maxY + 1);
    Point ans(0,0);

    int count = 0;
    size_t n = poly.size();
    for (size_t i = 0; i < n; ++i)
    {
        const Point& current = poly[i];
        const Point& next = poly[(i + 1) % n];

        if (intersect_segment(query, outside, current, next, ans))
        {
            count++;
        }
    }
    return count % 2 == 1;
}
////////////////////////////////////////////////////////////////////////////////
struct Compare
{
    Point p0; // Leftmost point of the poly
    bool operator()(const Point &p1, const Point &p2)
    {
        return det(p1, p2) - det(p0, p2) + det(p0, p1) > 0 || (det(p1, p2) - det(p0, p2) + det(p0, p1) == 0 && abs(p2 - p0) > abs(p1 - p0));
    }
};

bool inline salientAngle(Point &a, Point &b, Point &c)
{
    return (b.real() - a.real()) * (c.imag() - a.imag()) - (b.imag() - a.imag()) * (c.real() - a.real()) < 0;
}

Polygon convex_hull(std::vector<Point> &points)
{
    Compare order;
    Point p0 = points[0];

    for (int i = 1; i < points.size(); i++) {
        if (points[i].imag() < p0.imag() || (points[i].imag() == p0.imag() && points[i].real() < p0.real())) {
            p0 = points[i];
        }
    }
    
    order.p0 = p0;

    std::sort(points.begin(), points.end(), order);
    Polygon hull;
    
    for (int i = 0; i < points.size(); i++) {
        hull.push_back(points[i]);
        while (hull.size() > 2 && salientAngle(hull[hull.size() - 3], hull[hull.size() - 2], hull[hull.size() - 1])) {
            hull.pop_back();
            hull.pop_back();
            hull.push_back(points[i]);
        }
    }
    return hull;
}
std::vector<Point> load_xyz(const std::string& filename)
{
    std::ifstream in(filename);
    if (!in.is_open())
    {
        throw std::runtime_error("Unable to open file: " + filename);
    }

    std::vector<Point> points;
    std::string line;
    while (std::getline(in, line))
    {
        std::istringstream iss(line);
        double x, y;
        if (iss >> x >> y)
        {
            points.emplace_back(x, y);
        }
    }

    return points;
}

void save_xyz(const std::string& filename, const std::vector<Point>& points)
{
    std::ofstream out(filename);
    if (!out.is_open())
    {
        throw std::runtime_error("Unable to open file: " + filename);
    }

    out << std::fixed << std::setprecision(6);
    for (const Point& point : points)
    {
        out << point.real() << ' ' << point.imag() << '\n';
    }
}

Polygon load_obj(const std::string& filename)
{
    std::ifstream in(filename);
    if (!in.is_open())
    {
        throw std::runtime_error("Unable to open file: " + filename);
    }

    Polygon poly;
    std::string line;
    while (std::getline(in, line))
    {
        std::istringstream iss(line);
        std::string type;
        double x, y, z;
        if (iss >> type >> x >> y >> z && type == "v")
        {
            poly.emplace_back(x, y);
        }
    }

    return poly;
}
void save_obj(const std::string &filename, Polygon &poly)
{
    std::ofstream out(filename);
    if (!out.is_open())
    {
        throw std::runtime_error("failed to open file " + filename);
    }
    out << std::fixed;
    for (const auto &v : poly)
    {
        out << "v " << v.real() << ' ' << v.imag() << " 0\n";
    }
    for (size_t i = 0; i < poly.size(); ++i)
    {
        out << "l " << i + 1 << ' ' << 1 + (i + 1) % poly.size() << "\n";
    }
    out << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
    const std::string points_path = root_path + "/points.xyz";
    const std::string poly_path = root_path + "/polygon.obj";

    std::vector<Point> points = load_xyz(points_path);

    ////////////////////////////////////////////////////////////////////////////////
    //Point in polygon
    Polygon poly = load_obj(poly_path);
    std::vector<Point> result;
    for (size_t i = 0; i < points.size(); ++i)
    {
        if (is_inside(poly, points[i]))
        {
            result.push_back(points[i]);
        }
    }
    save_xyz("output.xyz", result);

    ////////////////////////////////////////////////////////////////////////////////
    //Convex hull
    Polygon hull = convex_hull(points);
    save_obj("output.obj", hull);
    save_obj("points.obj", result);

    return 0;
}
