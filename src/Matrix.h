#ifndef __MATRIX_H__
#define __MATRIX_H__
#include "EasyBMP.h"

#include <stdio.h>
#include <stdlib.h>


#ifdef OMP
#pragma omp parallel for
#else
#endif


struct Pixel
{
    static Pixel Zero;

    Pixel(){ r=g=b=0;}
    Pixel(double r, double g, double b)
    {
        this->r = r;
        this->g = g;
        this->b = b;
    }
    
    double r;
    double g;
    double b;

    friend Pixel operator - (const Pixel& a, const Pixel& b)
    {
        return Pixel(a.r - b.r,
        a.g - b.g,
        a.b - b.b);
    }

    Pixel& operator += (const Pixel& other)
    {
        r+=other.r;
        g+=other.g;
        b+=other.b;

        return *this;
    }
    
    Pixel& operator += (double eps)
    {
        r+=eps;
        g+=eps;
        b+=eps;

        return *this;
    }
    Pixel& operator -= (const Pixel& other)
    {
        r-=other.r;
        g-=other.g;
        b-=other.b;

        return *this;
    }
    Pixel& operator *= (const Pixel& other)
    {
        r*=other.r;
        g*=other.g;
        b*=other.b;

        return *this;
    }
    Pixel& operator *= (double v)
    {
        r*=v;
        g*=v;
        b*=v;

        return *this;
    }
    
    Pixel& operator /= (const Pixel& other)
    {
        r/=other.r;
        g/=other.g;
        b/=other.b;

        return *this;
    }
};

Pixel Pixel::Zero = Pixel(0,0,0);

struct Matrix
{
private:
    int m_h;
    int m_w;
    bool m_isGray;
    Pixel *p;
    
public:
    Matrix(BMP& bmp) 
    {
        p = NULL;
        init(bmp.TellHeight(), bmp.TellWidth(), bmp.TellBitDepth()==8);

        for(int h = 0; h < m_h; h++ )
        {
          for( int w = 0; w < m_w ; w++ )
          {
            RGBApixel rgb = bmp.GetPixel(w, h);
            Pixel pixel( double(rgb.Red)/255.0, 
                        double(rgb.Green)/255.0,
                        double(rgb.Blue)/255.0);
            set(h, w, pixel);
          }
        }
    }
    Matrix(const Matrix& other)
    {
        p = NULL;
        init(other.m_h, other.m_w, other.m_isGray);

    //int begin = clock();
       // memcpy(p, other.p, sizeof(Pixel)*m_h*m_w);
       #pragma omp parallel for
        for (int i = 0; i < m_h*m_w; i++)
        {
            p[i] = other.p[i];
        }
     //   printf("memcpy1 %d\n", clock()-begin);
    }
    Matrix& operator = (const Matrix& ref)
    {
        if (this == &ref)
        {
            return *this;
        }
        
        init(ref.m_h, ref.m_w, ref.m_isGray);
     //   int begin = clock();
        //memcpy(p, ref.p, sizeof(Pixel)*m_h*m_w);

        #pragma omp parallel for
        for (int i = 0; i < m_h*m_w; i++)
        {
            p[i] = ref.p[i];
        }

        
      //  printf("memcpy2 %d\n", clock()-begin);
        return *this;
    }

    Matrix& operator += (const Matrix& v)
    {
        #pragma omp parallel for
        for (int pos = 0; pos < m_h*m_w; pos++)
        {
            p[pos] += v.p[pos];
        }

        return *this;
    }
    Matrix& operator += (double eps)
    {
    #pragma omp parallel for
        for (int pos = 0; pos < m_h*m_w; pos++)
        {
        //if (pos ==m_h*m_w/7 ||pos ==m_h*m_w/7*2|| pos ==m_h*m_w/7*3||pos ==m_h*m_w/7*4||pos ==m_h*m_w/5||pos ==m_h*m_w/7*6)
        //{
        //    printf("pos(%d) %d %d %d\n", pos, omp_get_num_threads(),  omp_get_max_threads(),  omp_get_thread_num() );
        //}

        p[pos] += eps;
        }
        return *this;
    }

    Matrix& operator *= (const Matrix& v)
    {
    #pragma omp parallel for
        for (int pos = 0; pos < m_h*m_w; pos++)
        {
            p[pos] *= v.p[pos];
        }


        return *this;
    }
    Matrix& operator *= (double v)
    {
        for (int pos = 0; pos < m_h*m_w; pos++)
        {
            p[pos] *= v;
        }

        return *this;
    }
    
    Matrix& operator /= (const Matrix& v)
    {
    #pragma omp parallel for
        for (int pos = 0; pos < m_h*m_w; pos++)
        {
            p[pos] /= v.p[pos];
        }

        return *this;
    }
    Matrix& operator -= (const Matrix& v)
    {
    
    
    #pragma omp parallel for
        for (int pos = 0; pos < m_h*m_w; pos++)
        {
            p[pos] -= v.p[pos];
        }

        return *this;
    }

    friend Matrix operator + (const Matrix& a, const Matrix& b)
    {
        return Matrix(a) += b;
    }
    friend Matrix operator + (const Matrix& a, double eps)
    {
        Matrix r(a);
       return r += eps;
    }
    friend Matrix operator - (const Matrix& a, const Matrix& b)
    {
    Matrix r(a);
        
        return r -= b;
    }
    friend Matrix operator * (const Matrix& a, const Matrix& b)
    {
        Matrix r(a);
         return r *= b;
    }
    friend Matrix operator * (const Matrix& a, double eps)
    {
       Matrix r(a);
       return r *= eps;
    }
    
    friend Matrix operator / (const Matrix& a, const Matrix& b)
    {
        return Matrix(a) /= b;
    }
    
    Matrix(int h, int w, bool isGray)
    {
        init(h, w, isGray);
    }
    ~Matrix()
    {

        if (p)
            free(p); 

    }

    bool toBMP(BMP& bmp) const
    {
        bmp.SetBitDepth(24);
        bmp.SetSize(m_w, m_h);

        for (int w = 0; w < m_w; w++)
        {
            for (int h = 0; h < m_h; h++)
            {
                Pixel p = get(h, w);
            
                RGBApixel temp;
                temp.Red = round(p.r*255);
                temp.Green = round(p.g*255);
                temp.Blue = round(p.b*255);
                temp.Alpha = 0;
                    
                bmp.SetPixel(w, h, temp);
            }
        }


        return true;
    }

    void show(const char* p=NULL ) const
    {
        return ;
        printf("%s = [\n", p);
    
        for(int i=0; i < m_h; i++ )
        {
          for( int j=0; j < m_w ; j++ )
          {
              Pixel p = get(i, j);
              printf("%f,%f,%f ", p.b, p.g, p.r);
          }
          printf("; -------------\n");
        }
        printf("]\n\n\n");
    }

    void newsize(int h, int w)
    {
        if (p) free(p);
        
        p = NULL;
        init(h, w, m_isGray);
    }

    int Height() const {return m_h;} 
    int Width() const {return m_w;}

    void cum_for_every_h()
    {    
        //每行积分上一行的
        for (int i = 1; i < m_h; i++)
        {
            for (int j = 0; j < m_w; j++)
            {
                get(i,j) += get(i-1, j);
            }
        }
    }
    void cum_for_every_w()
    {    
        //每行积分上一行的
        for (int i = 1; i < m_w; i++)
        {
            for (int j = 0; j < m_h; j++)
            {
                get(j,i) += get(j, i-1);
            }
        }
    }

    void set_all_value_1()
    {
        for (int pos = 0; pos < m_h*m_w; pos++)
        {
            Pixel& t = p[pos];
            t.r = t.g = t.b = 1;
        }
    }

    void cum_diff_r_for_every_h(const Matrix& m, int r)
    {
        for (int i = 0; i < m_h; i++)
        {
            int min = i - r - 1;
            int max = i + r;
            if (max >= m_h)
            {
                max = m_h - 1;
            }
        
            for (int j = 0; j < m_w; j++)
            {
                get(i, j) = m.get(max, j) - ((min>=0)?m.get(min, j):Pixel::Zero);
            }
        }
    }

    void cum_diff_r_for_every_w(const Matrix& m, int r)
    {
        for (int i = 0; i < m_w; i++)
        {
            int min = i - r - 1;
            int max = i + r;
            if (max >= m_w)
            {
                max = m_w - 1;
            }
        
            for (int j = 0; j < m_h; j++)
            {
                get(j, i) = m.get(j, max) - ((min>=0)?m.get(j, min):Pixel::Zero);
            }
        }
    }
    
protected:
    const Pixel& get(int i , int j) const
    {
        return p[i*m_w+j];
    }

    Pixel& get(int i , int j) 
    {
        return p[i*m_w+j];
    }

    inline void set(int h, int w, Pixel v)
    {
        p[ h * m_w + w ] = v;
    }



private:
    Matrix();

    
    void init(int h, int w, bool isGray) {
        m_h = h;
        m_w = w;
        m_isGray = isGray;

        if (p)
        {
            //delete[] p;
            free(p);
        }
        //p = new Pixel[m_h*m_w];

     //   int begin = clock();
        p = (Pixel*)malloc(m_h*m_w*sizeof(Pixel));

      //  printf("malloc: %d  \n", clock()-begin);
    }    





};

#endif /* __MATRIX_H__ */

