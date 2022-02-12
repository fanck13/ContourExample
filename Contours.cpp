#include "Contours.h"

#include <stdio.h>
#include <stdlib.h>

Contours::Contours() : 
    m_SCALE(5), 
    m_nContour(5),
    m_NX(0), m_NY(0),
    m_xmin(1e32), m_xmax(-1e32),
    m_ymin(1e32), m_ymax(-1e32),
    m_zmin(1e32), m_zmax(-1e32)
{
    
}

void Contours::InitAxis()
{
    m_xaxis.resize(m_NX * m_SCALE);

    for (int i = 0; i < m_SCALE * m_NX; i++)
    {
        m_xaxis[i] = i;
    }

    m_yaxis.resize(m_NY * m_SCALE);

    for (int i = 0; i < m_SCALE * m_NY; i++)
    {
        m_yaxis[i] = i;
    }
}

void Contours::InitLevels()
{
    m_levels.resize(m_nContour);
    double lMin = m_zmin + (m_zmax - m_zmin) / ((m_nContour - 1) * 2);
    double lMax = m_zmax - (m_zmax - m_zmin) / ((m_nContour - 1) * 2);

    for (int i = 0; i < m_nContour; i++)
    {
        m_levels[i] = lMin + (lMax - lMin) / (m_nContour - 1) * i;
    }
}

void Contours::ReadDataFromFile(const std::string& fileName)
{
    FILE* fptr;
    if ((fptr = fopen(fileName.c_str(), "r")) == nullptr)
    {
        fprintf(stderr, "Failed to open data file\n");
        exit(-1);
    }

    if (fscanf(fptr, "%d %d", &m_NX, &m_NY) != 2)
    {
        fprintf(stderr, "Expected to be able to read width and height of grid\n");
        exit(-1);
    }

    if (m_NX < 4 || m_NY < 4 || m_NX > 1000 || m_NY > 1000)
    {
        fprintf(stderr, "Got an unexpected grid size (%d x %d)\n", m_NX, m_NY);
        exit(-1);
    }

    m_datas.resize(m_SCALE * m_NX);

    for (int i = 0; i < m_SCALE * m_NX; i++)
    {
        m_datas[i].resize(m_SCALE * m_NY);
    }

    for (int i = 0; i < m_SCALE * m_NX; i++)
    {
        for (int j = 0; j < m_SCALE * m_NY; j++)
        {
            m_datas[i][j] = 0;
        }
    }

    /*
       Read the data.
       This is VERY specific to the details of this example data.
    */
    double x, y, z;
    while (fscanf(fptr, "%lf %lf %lf", &x, &y, &z) == 3)
    {
        /* Spread the data over a block, for future smoothing */
        for (int ii = 0; ii < m_SCALE; ii++)
        {
            for (int jj = 0; jj < m_SCALE; jj++)
            {
                m_datas[m_SCALE * (int)x + ii][m_SCALE * (int)y + jj] = z;
            }
        }

        /* Check the bounds - debugging purposes */
        m_xmin = std::min(m_xmin, x);
        m_xmax = std::max(m_xmax, x);
        m_ymin = std::min(m_ymin, y);
        m_ymax = std::max(m_ymax, y);
        m_zmin = std::min(m_zmin, z);
        m_zmax = std::max(m_zmax, z);
    }

    fclose(fptr);
}

void Contours::CONREC()
{
#define xsect(p1,p2) (h[p2]*xh[p1]-h[p1]*xh[p2])/(h[p2]-h[p1])
#define ysect(p1,p2) (h[p2]*yh[p1]-h[p1]*yh[p2])/(h[p2]-h[p1])

    int m1, m2, m3, case_value;
    double dmin, dmax, x1, x2, y1, y2;
    double h[5];
    int sh[5];
    double xh[5], yh[5];
    int im[4] = { 0,1,1,0 }, jm[4] = { 0,0,1,1 };

    int castab[3][3][3] =
    {
      { {0,0,8},{0,2,5},{7,6,9} },
      { {0,3,4},{1,3,1},{4,3,0} },
      { {9,6,7},{5,2,0},{8,0,0} }
    };

    double temp1, temp2;

    int ilb = 0;
    int iub = m_SCALE * m_NX - 1;
    int jlb = 0;
    int jub = m_SCALE * m_NY - 1;

    for (int j = (jub - 1); j >= jlb; j--)
    {
        for (int i = ilb; i <= iub - 1; i++)
        {
            temp1 = std::min(m_datas[i][j], m_datas[i][j + 1]);
            temp2 = std::min(m_datas[i + 1][j], m_datas[i + 1][j + 1]);
            dmin = std::min(temp1, temp2);
            temp1 = std::max(m_datas[i][j], m_datas[i][j + 1]);
            temp2 = std::max(m_datas[i + 1][j], m_datas[i + 1][j + 1]);
            dmax = std::max(temp1, temp2);

            if (dmax < m_levels[0] || dmin > m_levels[m_nContour - 1])
                continue;

            for (int k = 0; k < m_nContour; k++)
            {
                if (m_levels[k] < dmin || m_levels[k] > dmax)
                    continue;

                for (int m = 4; m >= 0; m--)
                {
                    if (m > 0)
                    {
                        h[m] = m_datas[i + im[m - 1]][j + jm[m - 1]] - m_levels[k];
                        xh[m] = m_xaxis[i + im[m - 1]];
                        yh[m] = m_yaxis[j + jm[m - 1]];
                    }
                    else
                    {
                        h[0] = 0.25 * (h[1] + h[2] + h[3] + h[4]);
                        xh[0] = 0.50 * (m_xaxis[i] + m_xaxis[i + 1]);
                        yh[0] = 0.50 * (m_yaxis[j] + m_yaxis[j + 1]);
                    }

                    if (h[m] > 0.0)
                    {
                        sh[m] = 1;
                    }
                    else if (h[m] < 0.0)
                    {
                        sh[m] = -1;
                    }
                    else
                    {
                        sh[m] = 0;
                    }
                }

                /*
                   Note: at this stage the relative heights of the corners and the
                   centre are in the h array, and the corresponding coordinates are
                   in the xh and yh arrays. The centre of the box is indexed by 0
                   and the 4 corners by 1 to 4 as shown below.
                   Each triangle is then indexed by the parameter m, and the 3
                   vertices of each triangle are indexed by parameters m1,m2,and m3.
                   It is assumed that the centre of the box is always vertex 2
                   though this isimportant only when all 3 vertices lie exactly on
                   the same contour level, in which case only the side of the box
                   is drawn.
                      vertex 4 +-------------------+ vertex 3
                               | \               / |
                               |   \    m-3    /   |
                               |     \       /     |
                               |       \   /       |
                               |  m=2    X   m=2   |       the centre is vertex 0
                               |       /   \       |
                               |     /       \     |
                               |   /    m=1    \   |
                               | /               \ |
                      vertex 1 +-------------------+ vertex 2
                */
                /* Scan each triangle in the box */

                for (int m = 1; m <= 4; m++)
                {
                    m1 = m;
                    m2 = 0;
                    if (m != 4)
                        m3 = m + 1;
                    else
                        m3 = 1;
                    if ((case_value = castab[sh[m1] + 1][sh[m2] + 1][sh[m3] + 1]) == 0)
                        continue;
                    switch (case_value)
                    {
                    case 1: /* Line between vertices 1 and 2 */
                        x1 = xh[m1];
                        y1 = yh[m1];
                        x2 = xh[m2];
                        y2 = yh[m2];
                        break;
                    case 2: /* Line between vertices 2 and 3 */
                        x1 = xh[m2];
                        y1 = yh[m2];
                        x2 = xh[m3];
                        y2 = yh[m3];
                        break;
                    case 3: /* Line between vertices 3 and 1 */
                        x1 = xh[m3];
                        y1 = yh[m3];
                        x2 = xh[m1];
                        y2 = yh[m1];
                        break;
                    case 4: /* Line between vertex 1 and side 2-3 */
                        x1 = xh[m1];
                        y1 = yh[m1];
                        x2 = xsect(m2, m3);
                        y2 = ysect(m2, m3);
                        break;
                    case 5: /* Line between vertex 2 and side 3-1 */
                        x1 = xh[m2];
                        y1 = yh[m2];
                        x2 = xsect(m3, m1);
                        y2 = ysect(m3, m1);
                        break;
                    case 6: /* Line between vertex 3 and side 1-2 */
                        x1 = xh[m1];
                        y1 = yh[m1];
                        x2 = xsect(m1, m2);
                        y2 = ysect(m1, m2);
                        break;
                    case 7: /* Line between sides 1-2 and 2-3 */
                        x1 = xsect(m1, m2);
                        y1 = ysect(m1, m2);
                        x2 = xsect(m2, m3);
                        y2 = ysect(m2, m3);
                        break;
                    case 8: /* Line between sides 2-3 and 3-1 */
                        x1 = xsect(m2, m3);
                        y1 = ysect(m2, m3);
                        x2 = xsect(m3, m1);
                        y2 = ysect(m3, m1);
                        break;
                    case 9: /* Line between sides 3-1 and 1-2 */
                        x1 = xsect(m3, m1);
                        y1 = ysect(m3, m1);
                        x2 = xsect(m1, m2);
                        y2 = ysect(m1, m2);
                        break;
                    default:
                        break;
                    }

                    if (m_contours.count(k) == 0)
                    {
                        std::vector<Line> lines{ {{x1, y1}, {x2, y2}} };
                        m_contours[k] = lines;
                    }
                    else
                    {
                        m_contours[k].push_back({ {x1, y1}, {x2, y2} });
                    }
                } /* m */
            } /* k - contour */
        } /* i */
    } /* j */
}

void Contours::SmoothData()
{
    int  n = 0;
    int sum = 0;
    for (int i = 0; i < m_SCALE * m_NX; i++)
    {
        for (int j = 0; j < m_NY * m_SCALE; j++)
        {
            n = 0;
            sum = 0;
            for (int ii = -4; ii <= 4; ii++)
            {
                for (int jj = -4; jj <= 4; jj++)
                {
                    if (i + ii < 0 || i + ii >= m_SCALE * m_NX)
                        continue;
                    if (j + jj < 0 || j + jj >= m_SCALE * m_NY)
                        continue;
                    sum += m_datas[i + ii][j + jj];
                    n++;
                }
            }

            if (n <= 0)
            {
                fprintf(stderr, "No cells averaged, this shouldn't happen!\n");
                exit(-1);
            }
            m_datas[i][j] = sum / n;
        }
    }
}