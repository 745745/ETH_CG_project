/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Romain Pr¨¦vost

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/emitter.h>
#include <nori/warp.h>
#include <nori/shape.h>
#include <nori/bitmap.h>

NORI_NAMESPACE_BEGIN

class EnvMapper : public Emitter {
public:
    EnvMapper(const PropertyList& props) {
        m_envfile = props.getString("envfile","D:/ETH_CG_project/nori/external_image/test.exr");
        isEnv = true;
        bitmap = Bitmap(m_envfile);
        width = bitmap.cols();
        height = bitmap.rows();
        scalarMap = MatrixXf(height, width);
        marginalPDF = Eigen::VectorXf(height);
        marginalCDF = Eigen::VectorXf(height);
        conditionalPDF = MatrixXf(height, width);
        conditionalCDF = MatrixXf(height, width);
        for(int i=0;i<height;i++)
            for (int j = 0; j < width; j++)
            {
                Color3f pixel = bitmap.coeff(i, j);
                scalarMap(i, j) = pixel.getLuminance();
            }

        //avoid 0
        scalarMap.array() += 1e-8;      
        for (int i = 0; i < height; i++) 
        {
            scalarMap.row(i) *= sin(i*M_PI / height);
                                  
        }
        scalarMap.array() += 1e-8;
        for (int i = 0; i < height; i++)
        {
            marginalPDF(i) = scalarMap.row(i).sum();
        }

        marginalPDF /= marginalPDF.sum();
        marginalCDF(0) = marginalPDF(0);
        for (int i = 1; i < height; i++)
        {
            marginalCDF(i) = marginalCDF(i - 1) + marginalPDF(i);
        }

        for (int i = 0; i < height; i++) 
        {
            float rowSum = scalarMap.row(i).sum();            
            for (int j = 0; j < width; j++) 
                conditionalPDF(i, j) = scalarMap.coeff(i, j) / rowSum;            
        }

        for (int i = 0; i < height; i++)
        {
            for (int j = 0; j < width; j++)
                if (j == 0)
                    conditionalCDF(i, j) = conditionalPDF(i, j);
                else
                    conditionalCDF(i, j) = conditionalCDF(i, j - 1) + conditionalPDF(i, j);
        }

        //fullwhite();
        //test();

    }

    virtual std::string toString() const override {
        return tfm::format(
            "EnvLight[\n"
            "  envFile = %s,\n"
            "]",
            m_envfile);
    }

    int sampleFromMargin(const float sample) const
    {
        float sum = 0;
        for (int i = 0; i < height; i++)
        {
            if(marginalCDF(i) > sample)
                return i;
        }
        return height - 1;
    }
    
    int sampleFromCondition(const float sample,int row) const
    {
        float sum = 0;
        for (int i = 0; i < width; i++)
        {
            if (conditionalCDF(row,i) > sample)
                return i;
        }  
        return width - 1;
    }

    virtual Color3f eval(const EmitterQueryRecord& lRec) const override {
        Vector3f wi = lRec.wi.normalized();
        float theta = std::acos(wi.z()); 
        float phi = std::atan2(wi.y(), wi.x()); 
        if (phi < 0) phi += 2 * M_PI; 


        int i = int(theta / M_PI * height);
        int j = int(phi / (2 * M_PI) * width);

        if (i < 0)
            i = 0;
        if (j < 0)
            j = 0;

        if (i == height)
            i--;
        if (j == width)
            j--;
        return bitmap(i, j);

    }

    virtual Color3f sample(EmitterQueryRecord& lRec, const Point2f& sample) const override {
        int i = sampleFromMargin(sample.x());
        int j = sampleFromCondition(sample.y(), i);
        float theta = M_PI*((i + 0.5) / height);
        float phi = 2 * M_PI * ((j + 0.5) / width);

        Vector3f wi = Vector3f(sin(theta) * cos(phi),sin(theta) * sin(phi),cos(theta));

        lRec.wi = wi;
        lRec.shadowRay= Ray3f(lRec.ref, lRec.wi);
        Color3f color = bitmap(i, j);
        if (pdf(lRec) > 1e-9)
            return eval(lRec) / pdf(lRec);
        else
            return Color3f(0.f);
    }

    virtual float pdf(const EmitterQueryRecord& lRec) const override {
        Vector3f wi = lRec.wi.normalized();
        float theta = std::acos(wi.z());
        float phi = std::atan2(wi.y(), wi.x());
        if (phi < 0) phi += 2 * M_PI;


        int i = int(theta *INV_PI * height);
        int j = int(phi *INV_PI * width /2.f);


        float marginal = marginalPDF(i);
        float conditional = conditionalPDF(i, j);

        return marginal * conditional *INV_PI*INV_PI /(sin(theta));
    }


    void test()
    {       
        for (int i = 0; i < 5000; i++)
        {
            srand(i*i);
            float x = rand() / double(RAND_MAX);
            float y= rand() / double(RAND_MAX);
            int a = sampleFromMargin(x);
            int b = sampleFromCondition(y, a);

            bitmap(a, b) = Color3f(1, 0, 0);
        }
        bitmap.savePNG("D:/test");
    }

    void fullwhite()
    {
        bitmap = Bitmap(m_envfile);
        width = bitmap.cols();
        height = bitmap.rows();
        for (int i = 0; i < height; i++)
            for (int j = 0; j < width; j++)
            {
                bitmap(i, j) = 1;
            }
        bitmap.saveEXR("D:/test");
    }

protected:
    
    std::string m_envfile;
    Eigen::MatrixXf scalarMap;      
    Eigen::VectorXf marginalPDF;    
    Eigen::VectorXf marginalCDF;
    Eigen::MatrixXf conditionalPDF; 
    Eigen::MatrixXf conditionalCDF;
    int height;
    int width;        

    public:
        Bitmap bitmap;
};

NORI_REGISTER_CLASS(EnvMapper, "EnvMapper")
NORI_NAMESPACE_END