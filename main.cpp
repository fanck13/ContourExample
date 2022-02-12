#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "paulslib.h"
#include "bitmaplib.h"

#include <vector>
#include <map>
#include <iostream>



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

    Contours() :m_SCALE(5), m_nContour(5)
    {
        m_xmin = 1e32;
        m_xmax = -1e32;
        m_ymin = 1e32;
        m_ymax = -1e32;
        m_zmin = 1e32;
        m_zmax = -1e32;
    }
    void SetNX(int nx) { m_NX = nx; }
    int GetNX() const { return m_NX; }

    void SetNY(int ny) { m_NY = ny; }
    int GetNY() const { return m_NY; }

    void SetScale(int scale) { m_SCALE = scale; }
    int GetScale() const { return m_SCALE; }

    void InitAxis()
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

    void InitLevels()
    {
        m_levels.resize(m_nContour);
        m_levels[0] = m_zmin + (m_zmax - m_zmin) / 8;
        m_levels[1] = m_zmin + (m_zmax - m_zmin) / 4;
        m_levels[2] = (m_zmax + m_zmin) / 2;
        m_levels[3] = m_zmax - (m_zmax - m_zmin) / 4;
        m_levels[4] = m_zmax - (m_zmax - m_zmin) / 8;
    }

    void ReadDataFromFile(const std::string& fileName)
    {
        FILE* fptr;
        if ((fptr = fopen("D:\\GitHub\\ContourExample\\example.data", "r")) == NULL)
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

        /*
           Malloc space for the two dimensional data
           Initialise it to 0
        */

        /*if ((data = (double**)malloc(SCALE * NX * sizeof(double*))) == NULL)
        {
            fprintf(stderr, "Failed to malloc space for the data\n");
            exit(-1);
        }*/

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
            //lines++;
        }

        fclose(fptr);
        //fprintf(stderr, "Read %d lines\n", lines);
        fprintf(stderr, "X Range: %lf -> %lf\n", m_xmin, m_xmax);
        fprintf(stderr, "Y Range: %lf -> %lf\n", m_ymin, m_ymax);
        fprintf(stderr, "Z Range: %lf -> %lf\n", m_zmin, m_zmax);
    }

    void CONREC()
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

                        /*if (k != 1)
                            continue;*/

                        /* Finally draw the line */
                        //ConrecLine(x1, y1, x2, y2, z[k]);
                    } /* m */
                } /* k - contour */
            } /* i */
        } /* j */
    }

    const std::map<int, std::vector<Line>>& contour() const { return m_contours; }

    void SmoothData()
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

std::map<int, std::vector<Line>> CONTOURS;

/* The dimension of the grid */
int NX = 0;
int NY = 0;

/*
   The data will be scaled up by this amount to give some
   resolution for the contour lines.
*/
#define SCALE 5

/* Two dimensional array of data */
double** datas;
//std::vector<std::vector<double>> datas;

/* Arrays for the x axis and y axis coordinates */
double* xaxis, * yaxis;
//std::vector<double> xaxis;
//std::vector<double> yaxis;

/* Array for the contour levels, 5 of them */
#define NCONTOUR 5
double contours[NCONTOUR];

/* Image on which the contours will be drawn, see bitmaplib.c */
BITMAP4* image;

/* Prototype for CONREC and the line drawing function */
void CONREC(double**, int, int, int, int, double*, double*, int, double*,
    void (*drawline)(double, double, double, double, double));
void drawline(double, double, double, double, double);

/* Debugging - count the number of line segments drawn */
int vectorsdrawn = 0;

int main(int argc, char** argv)
{
    COLOUR colour;
    BITMAP4 col, grey = { 128,128,128,0 };
    FILE* fptr;

    {
        Contours con;
        con.ReadDataFromFile("D:\\GitHub\\ContourExample\\example.data");
        con.InitAxis();
        con.InitLevels();
        con.SmoothData();
        con.CONREC();
        auto c = con.contour();

        if ((image = Create_Bitmap(con.GetNX() * con.GetScale(), con.GetNY() * con.GetScale())) == NULL)
        {
            fprintf(stderr, "Malloc of bitmap failed\n");
            exit(-1);
        }

        Erase_Bitmap(image, con.GetNX() * con.GetScale(), con.GetNY() * con.GetScale(), grey); /* Not strictly necessary */

        for (int j = 0; j < con.GetScale() * con.GetNY(); j++)
        {
            for (int i = 0; i < con.GetScale() * con.GetNX(); i++)
            {
                //colour = GetColour(data[i][j], zmin, zmax, 1);
                col.r = /* colour.r * */ 255;
                col.g = /* colour.g * */ 255;
                col.b = /* colour.b * */ 255;
                Draw_Pixel(image, con.GetScale() * con.GetNX(), con.GetScale() * con.GetNY(), (double)i, (double)j, col);
            }
        }

        BITMAP4 black = { 0,0,0,0 };

        for (const auto& contour : c)
        {
            //std::cout << contour.first << std::endl;

            for (const auto& line : contour.second)
            {
                Draw_Line(image, con.GetScale() * con.GetNX(), con.GetScale() * con.GetNY(), (int)line.Start.X, (int)line.Start.Y, (int)line.End.X, (int)line.End.Y, black);
                vectorsdrawn++;

                //std::cout << "    " << "{" << "(" << line.Start.X << ", " << line.Start.Y << "), " << "(" << line.Start.X << ", " << line.Start.Y << ")}" << std::endl;
            }
        }

        

        if ((fptr = fopen("image1.BMP", "w")) == NULL) {
            fprintf(stderr, "Failed to open output image\n");
            exit(-1);
        }
        Write_Bitmap(fptr, image, con.GetScale() * con.GetNX(), con.GetScale() * con.GetNY(), 9);
        fclose(fptr);
    }

    int i, j, ii, jj;
    int n = 0, lines = 0;
    double sum;
    double x, y, z;
    double xmin = 1e32, xmax = -1e32;
    double ymin = 1e32, ymax = -1e32;
    double zmin = 1e32, zmax = -1e32;
    //COLOUR colour;
    //BITMAP4 col, grey = { 128,128,128,0 };
    //FILE* fptr;

    /*
       Open the data file and read the header.
       This is specific to the test dataset being used here,
       the header consists of the x and y dimensions of the data grid
       (which need not be fill)
    */
    if ((fptr = fopen("D:\\GitHub\\ContourExample\\example.data", "r")) == NULL)
    {
        fprintf(stderr, "Failed to open data file\n");
        exit(-1);
    }

    if (fscanf(fptr, "%d %d", &NX, &NY) != 2)
    {
        fprintf(stderr, "Expected to be able to read width and height of grid\n");
        exit(-1);
    }

    if (NX < 4 || NY < 4 || NX > 1000 || NY > 1000)
    {
        fprintf(stderr, "Got an unexpected grid size (%d x %d)\n", NX, NY);
        exit(-1);
    }

    /*
       Malloc space for the two dimensional data
       Initialise it to 0
    */
    
    if ((datas = (double**)malloc(SCALE * NX * sizeof(double*))) == NULL)
    {
        fprintf(stderr, "Failed to malloc space for the data\n");
        exit(-1);
    }

    if ((datas = (double**)malloc(SCALE * NX * sizeof(double*))) == NULL) {
        fprintf(stderr, "Failed to malloc space for the data\n");
        exit(-1);
    }
    for (i = 0; i < SCALE * NX; i++) {
        if ((datas[i] = (double*)malloc(SCALE * NY * sizeof(double))) == NULL) {
            fprintf(stderr, "Failed to malloc space for the data\n");
            exit(-1);
        }
    }
    for (i = 0; i < SCALE * NX; i++)
        for (j = 0; j < SCALE * NY; j++)
            datas[i][j] = 0;


    /*
       Read the data.
       This is VERY specific to the details of this example data.
    */
    while (fscanf(fptr, "%lf %lf %lf", &x, &y, &z) == 3)
    {
        /* Spread the data over a block, for future smoothing */
        for (ii = 0; ii < SCALE; ii++)
        {
            for (jj = 0; jj < SCALE; jj++)
            {
                datas[SCALE * (int)x + ii][SCALE * (int)y + jj] = z;
            }
        }

        /* Check the bounds - debugging purposes */
        xmin = std::min(xmin, x);
        xmax = std::max(xmax, x);
        ymin = std::min(ymin, y);
        ymax = std::max(ymax, y);
        zmin = std::min(zmin, z);
        zmax = std::max(zmax, z);
        lines++;
    }

    fclose(fptr);
    fprintf(stderr, "Read %d lines\n", lines);
    fprintf(stderr, "X Range: %lf -> %lf\n", xmin, xmax);
    fprintf(stderr, "Y Range: %lf -> %lf\n", ymin, ymax);
    fprintf(stderr, "Z Range: %lf -> %lf\n", zmin, zmax);

    /*
       Smooth the data
       This is just a simple 4x4 rectangular filter, it isn't
       actually necessary but makes the result sexier.
    */
    for (i = 0; i < SCALE * NX; i++)
    {
        for (j = 0; j < NY * SCALE; j++)
        {
            n = 0;
            sum = 0;
            for (ii = -4; ii <= 4; ii++)
            {
                for (jj = -4; jj <= 4; jj++)
                {
                    if (i + ii < 0 || i + ii >= SCALE * NX)
                        continue;
                    if (j + jj < 0 || j + jj >= SCALE * NY)
                        continue;
                    sum += datas[i + ii][j + jj];
                    n++;
                }
            }

            if (n <= 0)
            {
                fprintf(stderr, "No cells averaged, this shouldn't happen!\n");
                exit(-1);
            }
            datas[i][j] = sum / n;
        }
    }

    /*
       Set up the axis coordinates
       If helps to do this with thought so the line segments drawn
       by CONREC are in the most convenient coordinate system.
    */
    if ((xaxis = (double*)malloc(NX * SCALE * sizeof(double))) == NULL)
    {
        fprintf(stderr, "Failed to malloc space for the xaxis data\n");
        exit(-1);
    }

    for (i = 0; i < SCALE * NX; i++)
    {
        xaxis[i] = i;
    }

    if ((yaxis = (double*)malloc(NY * SCALE * sizeof(double))) == NULL)
    {
        fprintf(stderr, "Failed to malloc space for the yaxis data\n");
        exit(-1);
    }
    for (i = 0; i < SCALE * NY; i++)
    {
        yaxis[i] = i;
    }

    double lMin = zmin + (zmax - zmin) / ((NCONTOUR - 1) * 2);
    double lMax = zmax - (zmax - zmin) / ((NCONTOUR - 1) * 2);
    /* Set up the contour levels */

    for (i = 0; i < NCONTOUR; i++)
    {
        contours[i] = lMin + (lMax - lMin) / (NCONTOUR - 1) * i;
    }

    contours[0] = zmin + (zmax - zmin) / 8;
    contours[1] = zmin + (zmax - zmin) / 4;
    contours[2] = (zmax + zmin) / 2;
    contours[3] = zmax - (zmax - zmin) / 4;
    contours[4] = zmax - (zmax - zmin) / 8;

    /*
       Create the image in memory and draw a colour ramped version
       GetColour() just performs a colour ramp from blue to red.
    */
    if ((image = Create_Bitmap(NX * SCALE, NY * SCALE)) == NULL)
    {
        fprintf(stderr, "Malloc of bitmap failed\n");
        exit(-1);
    }

    Erase_Bitmap(image, NX * SCALE, NY * SCALE, grey); /* Not strictly necessary */

    for (j = 0; j < SCALE * NY; j++)
    {
        for (i = 0; i < SCALE * NX; i++)
        {
            //colour = GetColour(data[i][j], zmin, zmax, 1);
            col.r = /* colour.r * */ 255;
            col.g = /* colour.g * */ 255;
            col.b = /* colour.b * */ 255;
            Draw_Pixel(image, SCALE * NX, SCALE * NY, (double)i, (double)j, col);
        }
    }

    /* Finally do the contouring */
    CONREC(datas, 0, SCALE * NX - 1, 0, SCALE * NY - 1,
        xaxis, yaxis, NCONTOUR, contours, drawline);
    fprintf(stderr, "Drew %d vectors\n", vectorsdrawn);

    /*
       Write the image as a TGA file
       See bitmaplib.c for more details, or write "image"
       in your own prefered format.
    */
    if ((fptr = fopen("image.BMP", "w")) == NULL) {
        fprintf(stderr, "Failed to open output image\n");
        exit(-1);
    }
    Write_Bitmap(fptr, image, SCALE * NX, SCALE * NY, 9);
    fclose(fptr);

    exit(0);
}

void drawline(double x1, double y1, double x2, double y2, double z)
{
    BITMAP4 black = { 0,0,0,0 };

    if (x1 < 0 || x1 >= NX * SCALE || x2 < 0 || x2 > NX * SCALE)
        fprintf(stderr, "Shouldn't get here, x out of bounds: %g %g\n", x1, x2);
    if (y1 < 0 || y1 >= NY * SCALE || y2 < 0 || y2 > NY * SCALE)
        fprintf(stderr, "Shouldn't get here, y out of bounds: %g %g\n", y1, y2);
    Draw_Line(image, SCALE * NX, SCALE * NY, (int)x1, (int)y1, (int)x2, (int)y2, black);
    vectorsdrawn++;
}

/*
   Derivation from CONREC
   d               ! matrix of data to contour
   ilb,iub,jlb,jub ! index bounds of data matrix
   x               ! data matrix column coordinates
   y               ! data matrix row coordinates
   nc              ! number of contour levels
   z               ! contour levels in increasing order
*/
void CONREC(double** data, int ilb, int iub, int jlb, int jub,
    double* x, double* y, int nc, double* z,
    void (*ConrecLine)(double, double, double, double, double))
{
#define xsect(p1,p2) (h[p2]*xh[p1]-h[p1]*xh[p2])/(h[p2]-h[p1])
#define ysect(p1,p2) (h[p2]*yh[p1]-h[p1]*yh[p2])/(h[p2]-h[p1])

    int m1, m2, m3, case_value;
    double dmin, dmax, x1, x2, y1, y2;
    int i, j, k, m;
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

    for (j = (jub - 1); j >= jlb; j--) 
    {
        for (i = ilb; i <= iub - 1; i++) 
        {
            temp1 = std::min(data[i][j], data[i][j + 1]);
            temp2 = std::min(data[i + 1][j], data[i + 1][j + 1]);
            dmin = std::min(temp1, temp2);
            temp1 = std::max(data[i][j], data[i][j + 1]);
            temp2 = std::max(data[i + 1][j], data[i + 1][j + 1]);
            dmax = std::max(temp1, temp2);

            if (dmax < z[0] || dmin > z[nc - 1])
                continue;

            for (k = 0; k < nc; k++) 
            {
                if (z[k] < dmin || z[k] > dmax)
                    continue;
                
                for (m = 4; m >= 0; m--) 
                {
                    if (m > 0) 
                    {
                        h[m] = data[i + im[m - 1]][j + jm[m - 1]] - z[k];
                        xh[m] = x[i + im[m - 1]];
                        yh[m] = y[j + jm[m - 1]];
                    }
                    else 
                    {
                        h[0] = 0.25 * (h[1] + h[2] + h[3] + h[4]);
                        xh[0] = 0.50 * (x[i] + x[i + 1]);
                        yh[0] = 0.50 * (y[j] + y[j + 1]);
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

                for (m = 1; m <= 4; m++) 
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

                    if (CONTOURS.count(k) == 0)
                    {
                        std::vector<Line> lines{ {{x1, y1}, {x2, y2}} };
                        CONTOURS[k] = lines;
                    }
                    else
                    {
                        CONTOURS[k].push_back({ {x1, y1}, {x2, y2} });
                    }

                    if (k != 1)
                        continue;

                    /* Finally draw the line */
                    ConrecLine(x1, y1, x2, y2, z[k]);
                } /* m */
            } /* k - contour */
        } /* i */
    } /* j */

    /*for (const auto& contour : CONTOURS)
    {
        std::cout << contour.first << std::endl;

        for (const auto& line : contour.second)
        {
            std::cout << "    " << "{" << "(" << line.Start.X << ", " << line.Start.Y << "), " << "(" << line.Start.X << ", " << line.Start.Y << ")}" << std::endl;
        }
    }*/
}
