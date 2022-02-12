#pragma once

#include <vector>
#include <map>
#include <string>

struct Point2D
{
    double X;
    double Y;
};

struct Line
{
    Point2D Start;
    Point2D End;
};

class Contours
{
public:

    Contours();

    ~Contours() = default;

    void SetNX(int nx) { m_NX = nx; }
    int GetNX() const { return m_NX; }

    void SetNY(int ny) { m_NY = ny; }
    int GetNY() const { return m_NY; }

    void SetScale(int scale) { m_SCALE = scale; }
    int GetScale() const { return m_SCALE; }

    void InitAxis();

    void InitLevels();

    void ReadDataFromFile(const std::string& fileName);

    void CONREC();

    const std::map<int, std::vector<Line>>& contour() const { return m_contours; }

    void SmoothData();

private:
    std::map<int, std::vector<Line>> m_contours;

    int m_NX;
    int m_NY;
    int m_SCALE;
    int m_nContour;

    std::vector<std::vector<double>> m_datas;
    std::vector<double> m_xaxis;
    std::vector<double> m_yaxis;

    std::vector<double> m_levels;

    double m_xmin, m_xmax;
    double m_ymin, m_ymax;
    double m_zmin, m_zmax;
};
